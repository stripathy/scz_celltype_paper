#!/usr/bin/env python3
"""
Step 5: BANKSY spatial domain classification + depth-based layer assignment.

For each sample:
  1. Load annotated h5ad (with QC, hierarchical labels, depth predictions)
  2. Subset to corr_qc_pass cells (from step 02b; falls back to qc_pass)
  3. Run BANKSY clustering (λ=0.8, res=0.3) for spatially coherent domains
  4. Classify domains: Cortical / Vascular / WM, with L1 border flag
  5. Assign layers: depth-bin layers for all cells, Vascular cells overridden
  6. Spatially smooth layers (within-domain majority vote + vascular trim
     + BANKSY-anchored L1 contiguity)
  7. Save updated h5ad with new columns

BANKSY (Nature Genetics 2024) augments gene expression with spatial
neighbor expression, producing spatially coherent clusters that improve
domain classification vs the older K-NN composition approach:
  - L1 border cells correctly identified as Cortical (not "Extra-cortical")
  - White matter detected via oligo + depth thresholds
  - Vascular threshold lowered to 0.50 (clusters are spatially coherent)

Spatial smoothing pipeline (3 steps):
  - Within-domain majority vote (k=30, 2 rounds): smooths cortical layer
    boundaries without crossing BANKSY domain borders
  - Vascular border trim: reassigns border Vascular cells to cortical layers
    when >33% of neighbors are in L2/3–L6
  - BANKSY-anchored L1 contiguity: promotes banksy_is_l1 cells with shallow
    depth to L1, removes isolated non-BANKSY L1 cells

Requires: Step 02b (corr_qc_pass), Step 04 (depth predictions), pybanksy.

New/updated h5ad columns:
  - banksy_cluster:  int   — raw BANKSY cluster ID
  - banksy_domain:   str   — Cortical / Vascular / WM
  - banksy_is_l1:    bool  — True if cell in L1 border cluster
  - spatial_domain:  str   — Cortical / Vascular / WM
  - layer:           str   — final spatially-smoothed layer assignment

Usage:
    python3 -u 05_run_spatial_domains.py
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import anndata as ad

# Pipeline config
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipeline_config import H5AD_DIR, OUTPUT_DIR, MODULES_DIR, N_WORKERS

# Modules
sys.path.insert(0, MODULES_DIR)
from depth_model import assign_discrete_layers, smooth_layers_spatial, LAYER_BINS
from banksy_domains import (
    preprocess_for_banksy, run_banksy, classify_banksy_domains,
)


def _process_one_sample(h5ad_path):
    """Process one sample: BANKSY → domain classification → layer assignment."""
    t0 = time.time()
    sample_id = os.path.basename(h5ad_path).replace("_annotated.h5ad", "")

    try:
        adata = ad.read_h5ad(h5ad_path)
        n_total = adata.shape[0]

        # Use corr_qc_pass (spatial QC + margin filter + doublet exclusion);
        # fall back to qc_pass if step 02b not run
        if "corr_qc_pass" in adata.obs.columns:
            qc_mask = adata.obs["corr_qc_pass"].values.astype(bool)
        elif "qc_pass" in adata.obs.columns:
            qc_mask = adata.obs["qc_pass"].values.astype(bool)
        else:
            qc_mask = np.ones(n_total, dtype=bool)

        n_pass = qc_mask.sum()
        adata_pass = adata[qc_mask].copy()
        print(f"  [{sample_id}] {n_pass:,} / {n_total:,} QC-pass cells")

        # ── BANKSY clustering ─────────────────────────────────
        t_banksy = time.time()
        adata_b = preprocess_for_banksy(adata_pass)
        banksy_labels = run_banksy(adata_b)
        del adata_b  # free memory
        n_clusters = len(np.unique(banksy_labels))
        print(f"  [{sample_id}] BANKSY: {n_clusters} clusters "
              f"({time.time()-t_banksy:.0f}s)")

        # ── Domain classification ─────────────────────────────
        domains, is_l1, cluster_info = classify_banksy_domains(
            adata_pass, banksy_labels
        )
        domain_counts = {
            d: int((domains == d).sum())
            for d in ["Cortical", "Vascular", "WM"]
        }
        n_l1 = int(is_l1.sum())
        print(f"  [{sample_id}] Domains: "
              + ", ".join(f"{d}={n:,}" for d, n in domain_counts.items())
              + f", L1_border={n_l1:,}")

        # ── Layer assignment ──────────────────────────────────
        pred_depth = adata_pass.obs['predicted_norm_depth'].values
        depth_layers = assign_discrete_layers(pred_depth)

        # Hybrid: depth bins for all, override Vascular only
        combined_layers = depth_layers.copy()
        combined_layers[domains == 'Vascular'] = 'Vascular'

        # ── Spatial smoothing ───────────────────────────────────
        # 3-step pipeline: within-domain vote → vascular trim → L1 contiguity
        t_smooth = time.time()
        smoothed_layers = smooth_layers_spatial(
            coords=adata_pass.obsm['spatial'],
            layers=combined_layers,
            domains=domains,
            is_l1_banksy=is_l1,
            depths=pred_depth,
            verbose=True,
        )
        print(f"  [{sample_id}] Spatial smoothing: {time.time()-t_smooth:.0f}s")

        # ── spatial_domain column ──────────────────────────────
        # Map BANKSY domains to canonical spatial_domain values:
        #   Cortical → Cortical
        #   Vascular → Vascular
        #   WM → WM (white matter cells excluded from cortical analyses)
        # Also: cells with smoothed layer == 'WM' get spatial_domain = 'WM'
        # (handles edge cases where BANKSY says Cortical but depth says WM)
        spatial_domain_compat = domains.copy()
        spatial_domain_compat[smoothed_layers == 'WM'] = 'WM'

        # ── Write back to full adata ──────────────────────────
        # BANKSY columns
        full_banksy_cluster = np.full(n_total, -1, dtype=int)
        full_banksy_cluster[qc_mask] = banksy_labels
        adata.obs["banksy_cluster"] = full_banksy_cluster

        full_banksy_domain = np.full(n_total, "", dtype=object)
        full_banksy_domain[qc_mask] = domains
        adata.obs["banksy_domain"] = full_banksy_domain

        full_banksy_is_l1 = np.zeros(n_total, dtype=bool)
        full_banksy_is_l1[qc_mask] = is_l1
        adata.obs["banksy_is_l1"] = full_banksy_is_l1

        # Backward-compatible columns (QC-pass cells)
        full_spatial_domain = np.full(n_total, 'Unassigned', dtype=object)
        full_spatial_domain[qc_mask] = spatial_domain_compat

        full_layer = np.full(n_total, 'Unassigned', dtype=object)
        full_layer[qc_mask] = smoothed_layers

        # ── Impute layer/spatial_domain for QC-fail cells via KNN ──
        # Use spatial proximity to QC-pass neighbors so every cell
        # gets a layer assignment (needed for viewer rendering).
        n_fail = n_total - n_pass
        if n_fail > 0:
            from scipy.spatial import cKDTree
            coords_all = adata.obsm["spatial"][:, :2]
            pass_idx = np.where(qc_mask)[0]
            fail_idx = np.where(~qc_mask)[0]
            tree = cKDTree(coords_all[pass_idx])
            _, nn_idx = tree.query(coords_all[fail_idx], k=min(15, len(pass_idx)))
            # nn_idx indexes into pass_idx
            for col_full, col_pass in [
                (full_spatial_domain, spatial_domain_compat),
                (full_layer, smoothed_layers),
            ]:
                for j, fi in enumerate(fail_idx):
                    neighbors = nn_idx[j] if nn_idx.ndim > 1 else [nn_idx[j]]
                    # Majority vote among k nearest QC-pass neighbors
                    votes = {}
                    for ni in neighbors:
                        v = col_pass[ni]
                        votes[v] = votes.get(v, 0) + 1
                    col_full[fi] = max(votes, key=votes.get)
            print(f"  [{sample_id}] Imputed layer/domain for {n_fail:,} QC-fail cells via KNN")

        adata.obs['spatial_domain'] = full_spatial_domain
        adata.obs['layer'] = full_layer

        # ── Drop stale columns from archived exploration scripts ──
        stale_cols = [
            'cortical_strip_id', 'cortical_strip_tier', 'in_cortical_strip',
            'curved_strip_id', 'in_curved_strip', 'curved_strip_tier',
            'curved_strip_bank',
            'layer_unsmoothed', 'layer_depth_only',
        ]
        dropped = []
        for col in stale_cols:
            if col in adata.obs.columns:
                del adata.obs[col]
                dropped.append(col)
        if dropped:
            print(f"  [{sample_id}] Dropped {len(dropped)} stale columns: "
                  f"{', '.join(dropped)}")

        # Save
        adata.write_h5ad(h5ad_path)

        elapsed = time.time() - t0
        print(f"  [{sample_id}] Saved h5ad ({elapsed:.0f}s)")

        # Layer distribution for summary (using smoothed layers)
        layer_dist = {}
        for lname in list(LAYER_BINS.keys()) + ['Vascular']:
            layer_dist[lname] = int((smoothed_layers == lname).sum())

        return {
            "sample_id": sample_id, "status": "success",
            "n_total": n_total, "n_pass": int(n_pass),
            "n_cortical": domain_counts.get("Cortical", 0),
            "n_vascular": domain_counts.get("Vascular", 0),
            "n_wm": domain_counts.get("WM", 0),
            "n_l1_border": n_l1,
            "n_clusters": n_clusters,
            "layer_dist": layer_dist,
            "time": elapsed,
        }

    except Exception as e:
        import traceback
        elapsed = time.time() - t0
        print(f"  [{sample_id}] FAILED - {e}")
        traceback.print_exc()
        return {
            "sample_id": sample_id, "status": "error",
            "error": str(e), "time": elapsed,
        }


def main():
    t_start = time.time()

    # Find all h5ad files
    h5ad_files = sorted(
        os.path.join(H5AD_DIR, f)
        for f in os.listdir(H5AD_DIR)
        if f.endswith("_annotated.h5ad")
    )
    print(f"Found {len(h5ad_files)} samples")
    print(f"Pipeline: BANKSY domain classification (λ=0.8, res=0.3) "
          f"+ depth layers + spatial smoothing")
    print(f"{'='*60}")

    # Process sequentially (BANKSY is memory-intensive)
    results = []
    for i, h5ad_path in enumerate(h5ad_files):
        sid = os.path.basename(h5ad_path).replace("_annotated.h5ad", "")
        print(f"\n{'='*60}")
        print(f"[{i+1}/{len(h5ad_files)}] {sid}")
        print(f"{'='*60}")
        result = _process_one_sample(h5ad_path)
        results.append(result)

    # ── Summary ────────────────────────────────────────────────
    total_time = time.time() - t_start
    print(f"\n{'='*60}")
    print(f"SUMMARY - {total_time:.0f}s total")
    print(f"{'='*60}")

    n_ok = sum(1 for r in results if r["status"] == "success")
    print(f"  Success: {n_ok}/{len(results)}")

    # Aggregate stats
    success = [r for r in results if r["status"] == "success"]
    if success:
        total_pass = sum(r["n_pass"] for r in success)
        total_cortical = sum(r["n_cortical"] for r in success)
        total_vascular = sum(r["n_vascular"] for r in success)
        total_wm = sum(r["n_wm"] for r in success)
        total_l1 = sum(r["n_l1_border"] for r in success)

        print(f"\n  Total QC-pass cells: {total_pass:,}")
        print(f"  Cortical: {total_cortical:,} "
              f"({100*total_cortical/total_pass:.1f}%)")
        print(f"    of which L1 border: {total_l1:,} "
              f"({100*total_l1/total_pass:.1f}%)")
        print(f"  Vascular: {total_vascular:,} "
              f"({100*total_vascular/total_pass:.1f}%)")
        print(f"  WM: {total_wm:,} "
              f"({100*total_wm/total_pass:.1f}%)")

        # Aggregate layer distribution
        agg_layers = {}
        for r in success:
            for lname, n in r["layer_dist"].items():
                agg_layers[lname] = agg_layers.get(lname, 0) + n

        print(f"\n  Combined layer distribution:")
        for lname in ['L1', 'L2/3', 'L4', 'L5', 'L6', 'WM', 'Vascular']:
            n = agg_layers.get(lname, 0)
            print(f"    {lname:20s}: {n:>8,} ({100*n/total_pass:5.1f}%)")

        # Per-sample table
        print(f"\n  Per-sample breakdown:")
        print(f"  {'Sample':>10s} {'QC-pass':>8s} {'Cort%':>6s} "
              f"{'Vasc%':>6s} {'WM%':>5s} {'L1%':>5s} {'Time':>5s}")
        for r in sorted(success, key=lambda x: x['sample_id']):
            n = r["n_pass"]
            print(f"  {r['sample_id']:>10s} {n:>8,} "
                  f"{100*r['n_cortical']/n:>5.1f}% "
                  f"{100*r['n_vascular']/n:>5.1f}% "
                  f"{100*r['n_wm']/n:>4.1f}% "
                  f"{100*r['n_l1_border']/n:>4.1f}% "
                  f"{r['time']:>4.0f}s")

    failed = [r for r in results if r["status"] != "success"]
    if failed:
        print(f"\nFailed ({len(failed)}):")
        for r in failed:
            print(f"  {r['sample_id']}: {r.get('error', 'unknown')}")

    # Save summary CSV
    if success:
        summary_rows = []
        for r in success:
            row = {
                'sample_id': r['sample_id'],
                'n_pass': r['n_pass'],
                'n_cortical': r['n_cortical'],
                'n_vascular': r['n_vascular'],
                'n_wm': r['n_wm'],
                'n_l1_border': r['n_l1_border'],
                'n_clusters': r['n_clusters'],
                'pct_cortical': 100 * r['n_cortical'] / r['n_pass'],
                'pct_vascular': 100 * r['n_vascular'] / r['n_pass'],
                'pct_l1_border': 100 * r['n_l1_border'] / r['n_pass'],
            }
            row.update(r['layer_dist'])
            summary_rows.append(row)

        df = pd.DataFrame(summary_rows)
        csv_path = os.path.join(OUTPUT_DIR, "spatial_domain_summary.csv")
        df.to_csv(csv_path, index=False)
        print(f"\n  Saved: {csv_path}")

    # ── Merge into combined h5ad ───────────────────────────────
    print("\nMerging all samples into combined h5ad...")
    t_merge = time.time()
    adatas = []
    for r in success:
        path = os.path.join(H5AD_DIR, f"{r['sample_id']}_annotated.h5ad")
        adatas.append(ad.read_h5ad(path))
    combined = ad.concat(adatas, join='outer')
    combined_path = os.path.join(OUTPUT_DIR, "all_samples_annotated.h5ad")
    combined.write_h5ad(combined_path)
    file_size = os.path.getsize(combined_path) / 1e9
    print(f"  Saved: {combined_path} ({file_size:.2f} GB, {combined.shape})")
    print(f"  Merge took {time.time()-t_merge:.0f}s")

    print(f"\nTotal time: {time.time()-t_start:.0f}s")


if __name__ == "__main__":
    main()
