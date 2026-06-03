#!/usr/bin/env python3
"""
Re-run label transfer with hierarchical mode on QC-passing cells.

For each sample:
  1. Load the annotated h5ad (which now has QC columns)
  2. Subset to qc_pass == True cells
  3. Run correlation_label_transfer with hierarchical=True
  4. Re-run depth prediction on re-labeled data
  5. Save updated h5ad (preserving QC columns and including QC-fail cells)

Usage:
    python run_hierarchical_relabel.py
"""

import os
import sys
import time
import glob
import numpy as np
import pandas as pd
import anndata as ad
from multiprocessing import Pool, current_process

# Add modules to path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MODULES_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "modules")
sys.path.insert(0, MODULES_DIR)

from loading import load_xenium_sample, discover_samples
from label_transfer import (
    load_reference, correlation_label_transfer,
    get_seaad_colors, build_hierarchy_table, TAXONOMY_LEVELS
)
from depth_model import load_model, predict_depth
from plotting import plot_summary

# Paths
BASE_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))
H5AD_DIR = os.path.join(BASE_DIR, "output", "h5ad")
OUTPUT_DIR = os.path.join(BASE_DIR, "output")
REFERENCE_PATH = os.path.join(BASE_DIR, "data", "reference",
                               "Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad")
DEPTH_MODEL_PATH = os.path.join(BASE_DIR, "output", "depth_model_normalized.pkl")
DATA_DIR = os.path.join(BASE_DIR, "data", "raw")  # Raw h5 files in data/raw/

N_WORKERS = 4

# Globals for multiprocessing
_ref = None
_colors = None
_depth_model = None


def _init_worker():
    """Initialize worker: load reference and depth model once."""
    global _ref, _colors, _depth_model
    pid = current_process().pid

    print(f"  [Worker {pid}] Loading reference...")
    _ref = load_reference(REFERENCE_PATH)
    _colors = get_seaad_colors(_ref)

    if os.path.exists(DEPTH_MODEL_PATH):
        try:
            print(f"  [Worker {pid}] Loading depth model...")
            _depth_model = load_model(DEPTH_MODEL_PATH)
        except Exception as e:
            print(f"  [Worker {pid}] WARNING: Could not load depth model ({e}), "
                  f"will use existing depth predictions from h5ad")
            _depth_model = None
    else:
        _depth_model = None
        print(f"  [Worker {pid}] No depth model found.")

    print(f"  [Worker {pid}] Ready.")


def _process_one_sample(sample_info):
    """
    Re-label one sample:
    1. Load raw h5 (for expression data)
    2. Load existing annotated h5ad (for QC columns)
    3. Run hierarchical label transfer on QC-pass cells
    4. Re-run depth prediction
    5. Merge QC-fail cells back (with 'Unassigned' labels)
    6. Save
    """
    global _ref, _colors, _depth_model
    sid = sample_info["sample_id"]
    pid = current_process().pid
    t0 = time.time()

    try:
        # Load raw Xenium data (for fresh expression matrix)
        adata_raw = load_xenium_sample(
            sample_info["h5_path"], sample_info["boundaries_path"]
        )

        # Load existing annotated h5ad (for QC columns)
        h5ad_path = os.path.join(H5AD_DIR, f"{sid}_annotated.h5ad")
        adata_old = ad.read_h5ad(h5ad_path)

        # Store old labels for comparison
        old_labels = {}
        for col in ["class_label", "subclass_label", "supertype_label"]:
            if col in adata_old.obs.columns:
                old_labels[col] = adata_old.obs[col].values.copy()

        n_total = adata_raw.shape[0]

        # Get QC pass mask from old h5ad
        if "qc_pass" not in adata_old.obs.columns:
            print(f"  [{pid}] {sid}: WARNING - no qc_pass column, using all cells")
            qc_mask = np.ones(n_total, dtype=bool)
        else:
            qc_mask = adata_old.obs["qc_pass"].values.astype(bool)

        n_pass = qc_mask.sum()
        n_fail = n_total - n_pass
        print(f"  [{pid}] {sid}: {n_total:,} total, {n_pass:,} QC pass, "
              f"{n_fail:,} QC fail ({100*n_fail/n_total:.1f}%)")

        # Subset to QC-passing cells for label transfer
        adata_pass = adata_raw[qc_mask].copy()

        # Run hierarchical label transfer
        print(f"  [{pid}] {sid}: Running hierarchical label transfer...")
        annotated_pass = correlation_label_transfer(
            adata_pass, _ref, hierarchical=True
        )

        # Run depth prediction on pass cells
        if _depth_model is not None:
            print(f"  [{pid}] {sid}: Predicting depth...")
            pred_depth = predict_depth(annotated_pass, _depth_model)
            annotated_pass.obs['predicted_norm_depth'] = pred_depth
        elif "predicted_norm_depth" in adata_old.obs.columns:
            # Use existing depth predictions from old h5ad
            print(f"  [{pid}] {sid}: Using existing depth predictions from h5ad...")
            old_depth = adata_old.obs["predicted_norm_depth"].values
            annotated_pass.obs['predicted_norm_depth'] = old_depth[qc_mask]

        # Now build the full output: all cells, with QC-fail cells getting
        # 'Unassigned' labels
        # Start with all cells from raw
        # We'll rebuild the obs DataFrame for ALL cells

        # Create output adata from the pass cells' gene set
        shared_genes = list(annotated_pass.var_names)
        adata_out = adata_raw[:, shared_genes].copy()

        # Initialize all labels as Unassigned
        for col in ["class_label", "subclass_label", "supertype_label"]:
            adata_out.obs[col] = "Unassigned"
            adata_out.obs[f"{col}_confidence"] = 0.0
        adata_out.obs["predicted_norm_depth"] = np.nan

        # Copy in the pass-cell annotations
        pass_indices = np.where(qc_mask)[0]
        for col in ["class_label", "subclass_label", "supertype_label",
                     "class_label_confidence", "subclass_label_confidence",
                     "supertype_label_confidence", "predicted_norm_depth"]:
            if col in annotated_pass.obs.columns:
                adata_out.obs.iloc[pass_indices,
                                   adata_out.obs.columns.get_loc(col)] = \
                    annotated_pass.obs[col].values

        # Copy QC columns from old h5ad
        qc_cols = [c for c in adata_old.obs.columns
                   if c.startswith("fail_") or c.startswith("qc_") or
                   c in ["total_counts", "n_genes", "neg_probe_sum",
                          "neg_codeword_sum", "unassigned_sum"]]
        for col in qc_cols:
            adata_out.obs[col] = adata_old.obs[col].values

        # Add sample_id
        adata_out.obs["sample_id"] = sid

        # Copy lognorm and counts layers from annotated pass cells
        if "lognorm" in annotated_pass.layers:
            from scipy.sparse import issparse, lil_matrix
            lognorm_full = np.zeros((n_total, len(shared_genes)), dtype=np.float32)
            lognorm_pass = annotated_pass.layers["lognorm"]
            if issparse(lognorm_pass):
                lognorm_pass = lognorm_pass.toarray()
            lognorm_full[pass_indices] = lognorm_pass
            adata_out.layers["lognorm"] = lognorm_full

        if "counts" in annotated_pass.layers:
            from scipy.sparse import issparse as is_sparse
            counts_full = np.zeros((n_total, len(shared_genes)), dtype=np.float32)
            counts_pass = annotated_pass.layers["counts"]
            if is_sparse(counts_pass):
                counts_pass = counts_pass.toarray()
            counts_full[pass_indices] = counts_pass
            adata_out.layers["counts"] = counts_full

        # Save
        adata_out.write_h5ad(h5ad_path)

        # Count label changes (only for pass cells that had old labels)
        n_changed = {}
        for col, old_vals in old_labels.items():
            old_pass = old_vals[qc_mask]
            new_pass = annotated_pass.obs[col].values if col in annotated_pass.obs.columns else old_pass
            n_diff = (old_pass != new_pass).sum()
            n_changed[col] = n_diff

        # Summary plot
        plot_dir = os.path.join(OUTPUT_DIR, "plots")
        os.makedirs(plot_dir, exist_ok=True)
        plot_path = os.path.join(plot_dir, f"{sid}_summary.png")
        try:
            plot_summary(adata_out[qc_mask], plot_path, model_bundle=_depth_model)
        except Exception as e:
            print(f"  [{pid}] {sid}: Plot failed: {e}")

        elapsed = time.time() - t0
        depth_range = "N/A"
        if "predicted_norm_depth" in annotated_pass.obs.columns:
            d = annotated_pass.obs["predicted_norm_depth"]
            depth_range = f"[{d.min():.3f}, {d.max():.3f}]"

        print(f"  [{pid}] {sid}: Done in {elapsed:.0f}s, depth {depth_range}")
        print(f"  [{pid}] {sid}: Label changes (pass cells): {n_changed}")

        return {
            "sample_id": sid, "status": "success",
            "n_total": n_total, "n_pass": int(n_pass),
            "n_fail": int(n_fail),
            "n_changed": n_changed, "time": elapsed,
        }

    except Exception as e:
        import traceback
        elapsed = time.time() - t0
        print(f"  [{pid}] {sid}: FAILED after {elapsed:.0f}s - {e}")
        traceback.print_exc()
        return {
            "sample_id": sid, "status": "error",
            "error": str(e), "time": elapsed,
        }


def main():
    t_start = time.time()

    # Discover all samples
    all_samples = discover_samples(DATA_DIR)
    print(f"Found {len(all_samples)} samples, using {N_WORKERS} workers")
    print(f"Mode: HIERARCHICAL label transfer on QC-passing cells")

    # Process in parallel
    with Pool(
        processes=N_WORKERS,
        initializer=_init_worker,
    ) as pool:
        results = pool.map(_process_one_sample, all_samples)

    # Summary
    total_time = time.time() - t_start
    print(f"\n{'='*60}")
    print(f"SUMMARY - {total_time:.0f}s total")
    print(f"{'='*60}")

    n_ok = sum(1 for r in results if r["status"] == "success")
    total_cells = sum(r.get("n_total", 0) for r in results)
    total_pass = sum(r.get("n_pass", 0) for r in results)
    total_fail = sum(r.get("n_fail", 0) for r in results)

    print(f"  Success: {n_ok}/{len(results)}")
    print(f"  Total cells: {total_cells:,}")
    print(f"  QC pass: {total_pass:,} ({100*total_pass/total_cells:.1f}%)")
    print(f"  QC fail: {total_fail:,} ({100*total_fail/total_cells:.1f}%)")

    # Label change summary
    total_changed = {}
    for r in results:
        if r["status"] == "success" and "n_changed" in r:
            for level, n in r["n_changed"].items():
                total_changed[level] = total_changed.get(level, 0) + n

    print(f"\n  Label changes (QC-pass cells, old → new):")
    for level, n in sorted(total_changed.items()):
        pct = 100 * n / total_pass if total_pass > 0 else 0
        print(f"    {level}: {n:,} cells changed ({pct:.1f}%)")

    for r in results:
        status = "✓" if r["status"] == "success" else "✗"
        time_s = r.get("time", 0)
        n = r.get("n_total", 0)
        n_p = r.get("n_pass", 0)
        print(f"  {status} {r['sample_id']:10s}: {n:>6,} total, "
              f"{n_p:>6,} pass, {time_s:>4.0f}s")

    for r in results:
        if r["status"] == "error":
            print(f"  FAILED: {r['sample_id']}: {r.get('error', 'unknown')}")

    # Merge all into combined h5ad
    print("\nMerging all samples into combined h5ad...")
    t_merge = time.time()
    adatas = []
    for r in results:
        if r["status"] == "success":
            path = os.path.join(H5AD_DIR, f"{r['sample_id']}_annotated.h5ad")
            adatas.append(ad.read_h5ad(path))
    combined = ad.concat(adatas, join='outer')
    combined_path = os.path.join(OUTPUT_DIR, "all_samples_annotated.h5ad")
    combined.write_h5ad(combined_path)
    file_size = os.path.getsize(combined_path) / 1e9
    print(f"  Saved: {combined_path} ({file_size:.2f} GB)")
    print(f"  Shape: {combined.shape}")
    print(f"  Merge took {time.time()-t_merge:.0f}s")

    print(f"\nTotal time: {total_time:.0f}s")


if __name__ == "__main__":
    main()
