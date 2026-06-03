#!/usr/bin/env python3
"""
Step 2b: Two-stage hierarchical correlation classifier.

Reclassifies all QC-pass cells using Pearson correlation against centroids
built from high-confidence HANN-labeled exemplar cells (top-100 per type).

Stage 1: Classify against 24 subclass centroids.
Stage 2: Within the assigned subclass, classify against supertype centroids.

QC 1: Flag bottom N% by subclass correlation margin per sample (default 5%, set in pipeline_config.py)
QC 2: Flag suspected spatial doublets via marker co-expression

New columns added to each h5ad:
  - corr_subclass, corr_supertype, corr_class (str/categorical)
  - corr_subclass_corr, corr_subclass_margin, corr_supertype_corr (float32)
  - corr_qc_pass (bool)
  - doublet_suspect (bool)
  - doublet_type (str: '', 'Glut+GABA', 'GABA+GABA')

Requires: Step 02 (MapMyCells HANN labels) must have been run first.
The original HANN columns (subclass_label, supertype_label, class_label)
are preserved for comparison.

Usage:
    python3 -u 02b_run_correlation_classifier.py
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import anndata as ad

# Pipeline config
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipeline_config import (
    H5AD_DIR, MODULES_DIR, CENTROID_PATH,
    CORR_CLASSIFIER_TOP_N, CORR_CLASSIFIER_QC_PERCENTILE,
    L6B_MARGIN_THRESHOLD,
)

# Shared constants
sys.path.insert(0, MODULES_DIR)
from constants import SAMPLE_TO_DX, SUBCLASS_TO_CLASS

# Correlation classifier module
sys.path.insert(0, MODULES_DIR)
from correlation_classifier import (
    build_subclass_centroids,
    build_supertype_centroids,
    run_two_stage_classifier,
    flag_low_margin_cells,
    flag_doublet_cells,
)

UPPER_LAYERS = {"L1", "L2/3", "L4"}


def main():
    t_start = time.time()
    top_n = CORR_CLASSIFIER_TOP_N
    qc_pctl = CORR_CLASSIFIER_QC_PERCENTILE

    print("=" * 70)
    print("Step 2b: Two-Stage Hierarchical Correlation Classifier")
    print(f"  Top-N exemplars: {top_n}")
    print(f"  QC percentile: {qc_pctl}%")
    print("=" * 70)

    # ── Load all samples ──
    sample_ids = sorted(SAMPLE_TO_DX.keys())

    print(f"\nLoading {len(sample_ids)} samples...")
    adatas = []
    for i, sid in enumerate(sample_ids):
        fpath = os.path.join(H5AD_DIR, f"{sid}_annotated.h5ad")
        if not os.path.exists(fpath):
            print(f"  WARNING: {fpath} not found, skipping")
            continue
        adata = ad.read_h5ad(fpath)

        # Filter to QC-pass
        if "qc_pass" in adata.obs.columns:
            adata = adata[adata.obs["qc_pass"].values.astype(bool)].copy()

        print(f"  [{i+1:2d}/{len(sample_ids)}] {sid}: {adata.n_obs:,} cells",
              flush=True)
        adatas.append(adata)

    print("  Concatenating...", flush=True)
    combined = ad.concat(adatas, join='outer')
    print(f"  Total: {combined.n_obs:,} cells x {combined.n_vars} genes")

    # ── Build centroids ──
    print(f"\nBuilding subclass centroids (top-{top_n})...")
    sub_centroids, sub_counts, gene_names = build_subclass_centroids(
        combined, top_n=top_n)

    print(f"\nBuilding supertype centroids (top-{top_n})...")
    sup_centroids, sup_counts = build_supertype_centroids(
        combined, top_n=top_n)

    # ── Save centroids for optional nuclear resolution ──
    import pickle
    centroid_bundle = {
        'sub_centroids': sub_centroids,
        'sup_centroids': sup_centroids,
        'gene_names': gene_names,
        'sub_counts': sub_counts,
        'sup_counts': sup_counts,
    }
    with open(CENTROID_PATH, 'wb') as f:
        pickle.dump(centroid_bundle, f)
    print(f"Saved centroids: {CENTROID_PATH}")

    # ── Run two-stage classifier ──
    print("\nRunning two-stage hierarchical classifier...")
    results = run_two_stage_classifier(
        combined, sub_centroids, sup_centroids, gene_names)

    # ── Derive corr_class from corr_subclass ──
    results['corr_class'] = results['corr_subclass'].map(
        lambda x: SUBCLASS_TO_CLASS.get(x, 'Unknown'))

    # ── QC flagging: low margin ──
    print(f"\nFlagging bottom {qc_pctl}% margin per sample...")
    sample_ids_arr = combined.obs['sample_id'].astype(str).values
    margins = results['corr_subclass_margin'].values
    corr_qc_pass, thresholds = flag_low_margin_cells(
        margins, sample_ids_arr, percentile=qc_pctl,
        subclass_labels=results['corr_subclass'].values,
        l6b_margin_threshold=L6B_MARGIN_THRESHOLD)

    # ── QC flagging: spatial doublets ──
    print(f"\nDetecting spatial doublets...")
    class_labels = results['corr_class'].values
    doublet_suspect, doublet_type, doublet_stats = flag_doublet_cells(
        combined, class_labels, SUBCLASS_TO_CLASS)

    # Merge: doublets also fail QC
    n_doublet_only = (doublet_suspect & corr_qc_pass).sum()
    corr_qc_pass = corr_qc_pass & (~doublet_suspect)
    print(f"  Doublets not already low-margin: {n_doublet_only:,}")

    results['corr_qc_pass'] = corr_qc_pass
    results['doublet_suspect'] = doublet_suspect
    results['doublet_type'] = doublet_type

    # ── Summary statistics ──
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    # Overall subclass agreement
    hann_sub = combined.obs['subclass_label'].astype(str).values
    corr_sub = results['corr_subclass'].values
    agree_sub = (hann_sub == corr_sub).mean()
    print(f"\n  Subclass agreement (HANN vs Corr): {100*agree_sub:.1f}%")

    # L6b evaluation
    hann_l6b = hann_sub == 'L6b'
    corr_l6b = corr_sub == 'L6b'
    layer = combined.obs['layer'].astype(str).values

    if hann_l6b.sum() > 0:
        hann_l6b_upper = hann_l6b & np.isin(layer, list(UPPER_LAYERS))
        print(f"\n  L6b in upper layers:")
        print(f"    HANN: {hann_l6b_upper.sum():,}/{hann_l6b.sum():,} "
              f"({100*hann_l6b_upper.sum()/hann_l6b.sum():.1f}%)")

    if corr_l6b.sum() > 0:
        corr_l6b_upper = corr_l6b & np.isin(layer, list(UPPER_LAYERS))
        print(f"    Corr: {corr_l6b_upper.sum():,}/{corr_l6b.sum():,} "
              f"({100*corr_l6b_upper.sum()/corr_l6b.sum():.1f}%)")

    # QC summary
    n_fail = (~corr_qc_pass).sum()
    print(f"\n  QC flagged (total): {n_fail:,} / {len(corr_qc_pass):,} "
          f"({100*n_fail/len(corr_qc_pass):.1f}%)")
    print(f"    Low margin:   {(~corr_qc_pass & ~doublet_suspect).sum():,}")
    print(f"    Doublet only: {n_doublet_only:,}")
    print(f"    Both:         {(~results['corr_qc_pass'].values & doublet_suspect).sum() - n_doublet_only:,}")

    # Per-sample QC
    print(f"\n  Per-sample QC thresholds:")
    for sid in sorted(thresholds.keys()):
        sid_mask = sample_ids_arr == sid
        sid_fail = sid_mask & (~corr_qc_pass)
        print(f"    {sid}: margin threshold = {thresholds[sid]:.4f}, "
              f"flagged = {sid_fail.sum():,}/{sid_mask.sum():,}")

    # ── Write results back to individual h5ad files ──
    print(f"\n{'='*70}")
    print("Writing results to h5ad files...")
    print(f"{'='*70}")

    # Build index mapping: concatenated position -> (sample_id, cell_barcode)
    # We need to track which cells in the concatenated object belong to which sample
    cell_sample_ids = combined.obs['sample_id'].astype(str).values
    cell_barcodes = combined.obs.index.values

    # Group results by sample
    for sid in sorted(set(cell_sample_ids)):
        fpath = os.path.join(H5AD_DIR, f"{sid}_annotated.h5ad")
        adata = ad.read_h5ad(fpath)

        # Get indices of this sample's cells in the concatenated object
        sid_mask = cell_sample_ids == sid
        sid_barcodes = cell_barcodes[sid_mask]

        # Get QC-pass mask in the original h5ad
        if "qc_pass" in adata.obs.columns:
            qc_mask = adata.obs["qc_pass"].values.astype(bool)
        else:
            qc_mask = np.ones(adata.n_obs, dtype=bool)

        # Initialize columns with defaults for non-QC-pass cells
        adata.obs['corr_subclass'] = 'Unassigned'
        adata.obs['corr_supertype'] = 'Unassigned'
        adata.obs['corr_class'] = 'Unassigned'
        adata.obs['corr_subclass_corr'] = np.float32(0.0)
        adata.obs['corr_subclass_margin'] = np.float32(0.0)
        adata.obs['corr_supertype_corr'] = np.float32(0.0)
        adata.obs['corr_qc_pass'] = False
        adata.obs['doublet_suspect'] = False
        adata.obs['doublet_type'] = ''

        # Fill in QC-pass cells with classifier results
        sid_results = results.loc[sid_mask]

        # Match by cell barcode position
        qc_indices = np.where(qc_mask)[0]
        for col in ['corr_subclass', 'corr_supertype', 'corr_class',
                     'corr_subclass_corr', 'corr_subclass_margin',
                     'corr_supertype_corr', 'corr_qc_pass',
                     'doublet_suspect', 'doublet_type']:
            vals = sid_results[col].values
            if col in ['corr_subclass', 'corr_supertype', 'corr_class',
                        'doublet_type']:
                adata.obs[col] = adata.obs[col].astype(str)
                adata.obs.iloc[qc_indices,
                               adata.obs.columns.get_loc(col)] = vals.astype(str)
            elif col in ['corr_qc_pass', 'doublet_suspect']:
                adata.obs.iloc[qc_indices,
                               adata.obs.columns.get_loc(col)] = vals.astype(bool)
            else:
                adata.obs[col] = adata.obs[col].astype(np.float32)
                adata.obs.iloc[qc_indices,
                               adata.obs.columns.get_loc(col)] = vals.astype(np.float32)

        adata.write_h5ad(fpath)
        n_qc = qc_mask.sum()
        n_fail_s = (~sid_results['corr_qc_pass'].values).sum()
        n_doublet_s = sid_results['doublet_suspect'].values.sum()
        print(f"  {sid}: {n_qc:,} cells written "
              f"({n_fail_s} QC-flagged, {n_doublet_s} doublets)")

    elapsed = time.time() - t_start
    print(f"\nTotal elapsed: {elapsed:.0f}s ({elapsed/60:.1f} min)")


if __name__ == "__main__":
    main()
