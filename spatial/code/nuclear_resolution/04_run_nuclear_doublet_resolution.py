#!/usr/bin/env python3
"""
Optional (unnumbered): Hybrid nuclear doublet resolution.

This step is not part of the core numbered pipeline but can be run after
step 03 (transcript export) for additional doublet filtering.

For each sample with transcript data, builds a nuclear-only count matrix
(transcripts within nucleus polygons), runs doublet detection on nuclear
counts, and classifies each whole-cell doublet as:
  - "resolved"     — WC doublet=True, nuclear doublet=False (cytoplasmic spillover)
  - "persistent"   — WC doublet=True, nuclear doublet=True (likely real doublet)
  - "insufficient" — WC doublet=True, nuclear UMI < min threshold
  - "nuclear_only" — WC doublet=False, nuclear doublet=True (new catch)
  - "clean"        — not doublet in either assay

New columns added to each h5ad:
  - nuclear_total_counts (int)
  - nuclear_n_genes (int)
  - nuclear_fraction (float32)
  - nuclear_doublet_suspect (bool)
  - nuclear_doublet_type (str: '', 'Glut+GABA', 'GABA+GABA')
  - nuclear_doublet_status (str: 'resolved', 'persistent', etc.)
  - hybrid_qc_pass (bool)

Requires: Step 02b (correlation classifier) and Step 03 (transcript export).

Usage:
    python3 -u 06_run_nuclear_doublet_resolution.py              # all available
    python3 -u 06_run_nuclear_doublet_resolution.py Br8667       # specific sample
"""

import os
import sys
import glob
import re
import json
import time
import numpy as np
import pandas as pd
import anndata as ad
from multiprocessing import Pool, current_process

# Pipeline config
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipeline_config import (
    RAW_DIR, H5AD_DIR, MODULES_DIR, TRANSCRIPT_DIR, CENTROID_PATH,
    N_WORKERS, CORR_CLASSIFIER_TOP_N, CORR_CLASSIFIER_QC_PERCENTILE,
    NUCLEAR_CHUNK_SIZE, NUCLEAR_MIN_UMI,
)

# Shared constants + modules
sys.path.insert(0, MODULES_DIR)
from constants import SAMPLE_TO_DX, SUBCLASS_TO_CLASS
from nuclear_counts import (
    load_nucleus_polygons,
    build_nuclear_count_matrix,
    build_nuclear_adata,
)
from hybrid_qc import infer_class_from_markers, compute_hybrid_qc_pass
from correlation_classifier import (
    flag_doublet_cells,
    build_subclass_centroids,
    build_supertype_centroids,
    run_two_stage_classifier,
)


# ──────────────────────────────────────────────────────────────────────
# Discovery helpers
# ──────────────────────────────────────────────────────────────────────

def discover_nucleus_boundaries():
    """Find nucleus boundary CSVs and map sample_id → path."""
    result = {}
    for p in glob.glob(os.path.join(RAW_DIR, "GSM*-nucleus_boundaries.csv.gz")):
        m = re.search(r'(Br\d+)-nucleus_boundaries\.csv\.gz$', p)
        if m:
            result[m.group(1)] = p
    return result


def discover_transcript_dirs():
    """Find exported transcript directories (from step 03)."""
    result = {}
    if not os.path.isdir(TRANSCRIPT_DIR):
        return result
    for d in os.listdir(TRANSCRIPT_DIR):
        gene_idx = os.path.join(TRANSCRIPT_DIR, d, 'gene_index.json')
        if os.path.isfile(gene_idx):
            result[d] = os.path.join(TRANSCRIPT_DIR, d)
    return result


def find_ready_samples(requested=None):
    """Find samples that have all prerequisites for nuclear doublet resolution.

    Returns list of (sample_id, h5ad_path, nuc_boundary_path, transcript_dir).
    """
    nuc_boundaries = discover_nucleus_boundaries()
    transcript_dirs = discover_transcript_dirs()

    sample_ids = sorted(SAMPLE_TO_DX.keys())
    if requested:
        sample_ids = [s for s in requested if s in sample_ids]

    ready = []
    skipped = []
    for sid in sample_ids:
        h5ad_path = os.path.join(H5AD_DIR, f"{sid}_annotated.h5ad")
        if not os.path.exists(h5ad_path):
            skipped.append((sid, "no h5ad"))
            continue
        if sid not in nuc_boundaries:
            skipped.append((sid, "no nucleus boundaries"))
            continue
        if sid not in transcript_dirs:
            skipped.append((sid, "no transcript data (run step 03)"))
            continue
        ready.append((sid, h5ad_path, nuc_boundaries[sid], transcript_dirs[sid]))

    return ready, skipped


# ──────────────────────────────────────────────────────────────────────
# Worker initializer (centroids for high-UMI cell classification)
# ──────────────────────────────────────────────────────────────────────

_sub_centroids = None
_sup_centroids = None
_centroid_gene_names = None


def _init_worker(sub_c, sup_c, gn):
    """Initialize worker process with pre-built centroids."""
    global _sub_centroids, _sup_centroids, _centroid_gene_names
    _sub_centroids = sub_c
    _sup_centroids = sup_c
    _centroid_gene_names = gn


# ──────────────────────────────────────────────────────────────────────
# Per-sample processing
# ──────────────────────────────────────────────────────────────────────

def process_sample(args):
    """Process one sample: build nuclear counts, detect doublets, classify.

    Runs doublet detection on ALL cells (not just QC-pass), using marker-based
    class inference for cells that lack a corr_class assignment. This is more
    rigorous than UMI-based QC filtering — nuclear evidence directly tests
    whether marker co-expression is real (nuclear) or artifactual (cytoplasmic
    spillover).

    Parameters
    ----------
    args : tuple
        (sample_id, h5ad_path, nuc_boundary_path, transcript_dir)

    Returns
    -------
    dict with sample-level summary statistics
    """
    sample_id, h5ad_path, nuc_boundary_path, transcript_dir = args
    pid = current_process().pid
    t0 = time.time()

    print(f"\n[{sample_id}] Starting (PID {pid})...", flush=True)

    # 1. Load annotated h5ad (whole-cell)
    print(f"  [{sample_id}] Loading h5ad...", flush=True)
    adata = ad.read_h5ad(h5ad_path)
    n_cells = adata.n_obs
    print(f"  [{sample_id}] {n_cells:,} cells × {adata.n_vars} genes", flush=True)

    # 1b. Classify high-UMI cells that lack correlation assignments
    #     These cells were excluded by step 01's fail_total_counts_high filter,
    #     so step 02b never ran the correlation classifier on them. We classify
    #     them now so they can be properly evaluated by the hybrid QC.
    n_newly_classified = 0
    n_hu_margin_pass = 0
    if ('fail_total_counts_high' in adata.obs.columns and
            _sub_centroids is not None):
        high_umi = adata.obs['fail_total_counts_high'].values.astype(bool)
        # Cells that failed ONLY for high UMI (no other basic QC failure)
        other_fails = np.zeros(n_cells, dtype=bool)
        for fail_col in ['fail_neg_probe', 'fail_neg_codeword',
                         'fail_unassigned', 'fail_n_genes_low',
                         'fail_total_counts_low']:
            if fail_col in adata.obs.columns:
                other_fails |= adata.obs[fail_col].values.astype(bool)
        high_umi_only = high_umi & ~other_fails

        # Check which need classification (Unassigned = never classified by 02b)
        existing_sub = (adata.obs['corr_subclass'].astype(str).values
                        if 'corr_subclass' in adata.obs.columns
                        else np.full(n_cells, 'Unassigned'))
        needs_class = high_umi_only & np.isin(existing_sub, ['Unassigned', 'nan'])

        if needs_class.sum() > 0:
            print(f"  [{sample_id}] Classifying {needs_class.sum():,} "
                  f"high-UMI cells with correlation classifier...", flush=True)

            # Ensure string columns for assignment
            for col in ['corr_subclass', 'corr_supertype', 'corr_class']:
                if col in adata.obs.columns:
                    adata.obs[col] = adata.obs[col].astype(str)

            adata_hu = adata[needs_class].copy()
            cls_results = run_two_stage_classifier(
                adata_hu, _sub_centroids, _sup_centroids,
                _centroid_gene_names)
            cls_results['corr_class'] = cls_results['corr_subclass'].map(
                lambda x: SUBCLASS_TO_CLASS.get(x, 'Unknown'))

            # Write classifications back to adata
            for col in ['corr_subclass', 'corr_supertype', 'corr_class',
                        'corr_subclass_corr', 'corr_subclass_margin',
                        'corr_supertype_corr']:
                adata.obs.loc[cls_results.index, col] = cls_results[col].values

            # Apply margin filter using threshold from existing classified cells
            orig_classified = (~high_umi &
                               ~np.isin(adata.obs['corr_subclass'].astype(str).values,
                                        ['Unassigned', 'nan']))
            if orig_classified.sum() > 0:
                orig_margins = adata.obs.loc[orig_classified,
                                             'corr_subclass_margin'].values.astype(float)
                margin_threshold = np.percentile(orig_margins,
                                                 CORR_CLASSIFIER_QC_PERCENTILE)
            else:
                margin_threshold = 0.0

            new_margins = cls_results['corr_subclass_margin'].values
            margin_pass = new_margins >= margin_threshold

            # Update corr_qc_pass for newly classified cells
            if 'corr_qc_pass' not in adata.obs.columns:
                adata.obs['corr_qc_pass'] = False
            adata.obs.loc[cls_results.index, 'corr_qc_pass'] = pd.Series(
                margin_pass, index=cls_results.index)

            n_newly_classified = int(needs_class.sum())
            n_hu_margin_pass = int(margin_pass.sum())

            # Print subclass distribution of margin-passing cells
            if n_hu_margin_pass > 0:
                pass_idx = cls_results.index[margin_pass]
                pass_subs = adata.obs.loc[pass_idx, 'corr_subclass']
                sub_counts = pass_subs.value_counts()
                top_types = ', '.join(f"{k}: {v}" for k, v in
                                      sub_counts.head(6).items())
                print(f"    Margin pass: {n_hu_margin_pass:,}/{n_newly_classified:,} "
                      f"({100*n_hu_margin_pass/n_newly_classified:.1f}%) | "
                      f"Threshold: {margin_threshold:.4f}", flush=True)
                print(f"    Top types: {top_types}", flush=True)
            else:
                print(f"    Margin pass: 0/{n_newly_classified:,} | "
                      f"Threshold: {margin_threshold:.4f}", flush=True)

    # 2. Load nucleus polygons
    print(f"  [{sample_id}] Loading nucleus polygons...", flush=True)
    cell_ids, polygons, cell_id_to_poly_idx = load_nucleus_polygons(
        nuc_boundary_path)

    # 3. Load gene index
    gene_index_path = os.path.join(transcript_dir, 'gene_index.json')
    with open(gene_index_path) as f:
        gene_index = json.load(f)
    print(f"  [{sample_id}] Gene index: {gene_index['n_genes']} genes, "
          f"{gene_index['total_transcripts']:,} transcripts", flush=True)

    # 4. Build nuclear count matrix
    print(f"  [{sample_id}] Building nuclear count matrix...", flush=True)
    nuc_counts, gene_stats = build_nuclear_count_matrix(
        polygons, cell_ids, cell_id_to_poly_idx,
        adata, transcript_dir, gene_index,
        chunk_size=NUCLEAR_CHUNK_SIZE)

    # 5. Compute nuclear metadata
    nuc_total = np.array(nuc_counts.sum(axis=1)).flatten()
    nuc_ngenes = np.array((nuc_counts > 0).sum(axis=1)).flatten()
    wc_total = adata.obs['total_counts'].values.astype(float)
    nuc_fraction = np.where(wc_total > 0, nuc_total / wc_total, 0.0).astype(np.float32)

    # 6. Build class labels for ALL cells
    # For cells with corr_class, use it; for others, infer from marker expression
    print(f"  [{sample_id}] Assigning class labels for all cells...", flush=True)

    if 'corr_class' in adata.obs.columns:
        wc_class_labels = adata.obs['corr_class'].astype(str).values.copy()
        unassigned_mask = (wc_class_labels == 'Unassigned') | (wc_class_labels == 'nan')
        n_unassigned = unassigned_mask.sum()

        if n_unassigned > 0:
            print(f"    {n_unassigned:,} cells lack class labels — inferring from markers",
                  flush=True)
            inferred = infer_class_from_markers(adata[unassigned_mask])
            wc_class_labels[unassigned_mask] = inferred
    else:
        print(f"    No corr_class column — inferring all from markers", flush=True)
        wc_class_labels = infer_class_from_markers(adata)

    # 7. Run WHOLE-CELL doublet detection on ALL cells
    print(f"  [{sample_id}] Running WC doublet detection (all cells)...", flush=True)
    wc_doublet_all, wc_doublet_type_all, wc_dbl_stats = flag_doublet_cells(
        adata, wc_class_labels, SUBCLASS_TO_CLASS)

    # 8. Run NUCLEAR doublet detection on ALL cells
    print(f"  [{sample_id}] Running nuclear doublet detection (all cells)...", flush=True)

    # Build temporary AnnData for nuclear counts
    adata_nuc_tmp = ad.AnnData(
        X=nuc_counts,
        obs=adata.obs[['sample_id']].copy(),
        var=adata.var.copy(),
    )

    nuc_doublet_suspect, nuc_doublet_type, nuc_doublet_stats = flag_doublet_cells(
        adata_nuc_tmp, wc_class_labels, SUBCLASS_TO_CLASS)

    # 9. Classify doublets into resolution categories
    status = np.full(n_cells, 'clean', dtype=object)

    wc_dbl_mask = wc_doublet_all
    nuc_dbl_mask = nuc_doublet_suspect

    # Insufficient nuclear evidence
    insufficient_mask = wc_dbl_mask & (nuc_total < NUCLEAR_MIN_UMI)
    status[insufficient_mask] = 'insufficient'

    # Resolved: WC doublet but NOT nuclear doublet (and sufficient UMI)
    resolved_mask = wc_dbl_mask & ~nuc_dbl_mask & (nuc_total >= NUCLEAR_MIN_UMI)
    status[resolved_mask] = 'resolved'

    # Persistent: doublet in both WC and nuclear
    persistent_mask = wc_dbl_mask & nuc_dbl_mask & (nuc_total >= NUCLEAR_MIN_UMI)
    status[persistent_mask] = 'persistent'

    # Nuclear-only: NOT WC doublet but IS nuclear doublet
    nuclear_only_mask = ~wc_dbl_mask & nuc_dbl_mask
    status[nuclear_only_mask] = 'nuclear_only'

    # 10. Build hybrid_qc_pass — nuclear evidence replaces blunt UMI filtering
    hybrid_qc, basic_qc, n_high_umi_rescued = compute_hybrid_qc_pass(adata, status)
    if n_high_umi_rescued > 0:
        print(f"    Re-evaluating {n_high_umi_rescued:,} cells that failed only due to high UMI",
              flush=True)

    # 11. Write results back to h5ad
    print(f"  [{sample_id}] Writing results to h5ad...", flush=True)

    # Update doublet columns to reflect ALL-cell detection
    adata.obs['doublet_suspect'] = wc_doublet_all
    adata.obs['doublet_type'] = wc_doublet_type_all.astype(str)

    adata.obs['nuclear_total_counts'] = nuc_total.astype(int)
    adata.obs['nuclear_n_genes'] = nuc_ngenes.astype(int)
    adata.obs['nuclear_fraction'] = nuc_fraction
    adata.obs['nuclear_doublet_suspect'] = nuc_doublet_suspect
    adata.obs['nuclear_doublet_type'] = nuc_doublet_type.astype(str)
    adata.obs['nuclear_doublet_status'] = status.astype(str)
    adata.obs['hybrid_qc_pass'] = hybrid_qc

    adata.write_h5ad(h5ad_path)

    # 12. Compute summary statistics (report on all cells, not just old qc_pass)
    # Use basic_qc as the denominator (cells that pass non-doublet QC filters)
    n_basic = basic_qc.sum()

    n_wc_doublets_all = wc_dbl_mask.sum()
    n_wc_doublets = wc_dbl_mask[basic_qc].sum()
    n_resolved = (status[basic_qc] == 'resolved').sum()
    n_persistent = (status[basic_qc] == 'persistent').sum()
    n_insufficient = (status[basic_qc] == 'insufficient').sum()
    n_nuclear_only = (status[basic_qc] == 'nuclear_only').sum()

    # Type-specific resolution
    n_gg_wc = ((wc_doublet_type_all == 'Glut+GABA') & basic_qc).sum()
    n_gg_resolved = ((wc_doublet_type_all == 'Glut+GABA') & (status == 'resolved') & basic_qc).sum()
    n_gaga_wc = ((wc_doublet_type_all == 'GABA+GABA') & basic_qc).sum()
    n_gaga_resolved = ((wc_doublet_type_all == 'GABA+GABA') & (status == 'resolved') & basic_qc).sum()

    # Hybrid QC impact vs old corr_qc_pass
    old_qc_mask = adata.obs['qc_pass'].values.astype(bool)
    n_old_corr_pass = 0
    if 'corr_qc_pass' in adata.obs.columns:
        n_old_corr_pass = int((adata.obs['corr_qc_pass'].values.astype(bool) & old_qc_mask).sum())
    n_hybrid_qc_pass = int(hybrid_qc.sum())
    net_rescued = n_hybrid_qc_pass - n_old_corr_pass

    # Nuclear fraction stats
    valid_nf = nuc_fraction[basic_qc & (nuc_total > 0)]
    median_nf = float(np.median(valid_nf)) if len(valid_nf) > 0 else 0.0

    # High-UMI cell rescue stats
    n_high_umi_total = 0
    n_high_umi_rescued = 0
    high_umi_subclass_counts = {}
    if 'fail_total_counts_high' in adata.obs.columns:
        high_umi_mask = adata.obs['fail_total_counts_high'].values.astype(bool)
        n_high_umi_total = int(high_umi_mask.sum())
        high_umi_rescued_mask = high_umi_mask & hybrid_qc
        n_high_umi_rescued = int(high_umi_rescued_mask.sum())

        if n_high_umi_rescued > 0:
            rescued_subs = adata.obs.loc[high_umi_rescued_mask, 'corr_subclass']
            high_umi_subclass_counts = rescued_subs.value_counts().to_dict()

    elapsed = time.time() - t0

    stats = {
        'sample_id': sample_id,
        'status': 'success',
        'n_cells': n_cells,
        'n_basic_qc': int(n_basic),
        'n_wc_doublets': int(n_wc_doublets),
        'n_wc_doublets_all': int(n_wc_doublets_all),
        'n_resolved': int(n_resolved),
        'n_persistent': int(n_persistent),
        'n_insufficient': int(n_insufficient),
        'n_nuclear_only': int(n_nuclear_only),
        'n_glut_gaba_wc': int(n_gg_wc),
        'n_glut_gaba_resolved': int(n_gg_resolved),
        'n_gaba_gaba_wc': int(n_gaga_wc),
        'n_gaba_gaba_resolved': int(n_gaga_resolved),
        'resolution_rate': round(n_resolved / n_wc_doublets * 100, 1) if n_wc_doublets > 0 else 0,
        'n_old_corr_qc_pass': n_old_corr_pass,
        'n_hybrid_qc_pass': n_hybrid_qc_pass,
        'net_rescued': net_rescued,
        'median_nuclear_fraction': round(median_nf, 3),
        'n_high_umi_classified': n_newly_classified,
        'n_high_umi_margin_pass': n_hu_margin_pass,
        'n_high_umi_total': n_high_umi_total,
        'n_high_umi_rescued': n_high_umi_rescued,
        'high_umi_subclass_counts': high_umi_subclass_counts,
        'elapsed_sec': round(elapsed, 1),
    }

    print(f"\n  [{sample_id}] DONE in {elapsed:.0f}s ({elapsed/60:.1f} min)", flush=True)
    print(f"    All cells: {n_cells:,} | Basic QC: {n_basic:,}", flush=True)
    print(f"    WC doublets (all cells): {n_wc_doublets_all:,} | "
          f"(basic QC): {n_wc_doublets:,}", flush=True)
    print(f"    Resolved: {n_resolved:,} | Persistent: {n_persistent:,} | "
          f"Insufficient: {n_insufficient:,} | Nuclear-only: {n_nuclear_only:,}", flush=True)
    print(f"    Resolution rate: {stats['resolution_rate']:.1f}%", flush=True)
    print(f"    Old corr_qc_pass: {n_old_corr_pass:,} → hybrid_qc_pass: "
          f"{n_hybrid_qc_pass:,} ({net_rescued:+,})", flush=True)
    if n_high_umi_rescued > 0:
        top_types = ', '.join(f"{k}: {v}" for k, v in
                              sorted(high_umi_subclass_counts.items(),
                                     key=lambda x: -x[1])[:6])
        print(f"    High-UMI rescued: {n_high_umi_rescued:,}/{n_high_umi_total:,} | "
              f"Types: {top_types}", flush=True)

    return stats


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────

def main():
    t_start = time.time()

    print("=" * 70)
    print("Optional: Hybrid Nuclear Doublet Resolution")
    print(f"  Workers: {N_WORKERS}")
    print(f"  Nuclear min UMI: {NUCLEAR_MIN_UMI}")
    print(f"  Chunk size: {NUCLEAR_CHUNK_SIZE:,}")
    print("=" * 70)

    # Parse optional sample ID arguments
    requested = sys.argv[1:] if len(sys.argv) > 1 else None

    # Discover ready samples
    ready, skipped = find_ready_samples(requested)

    print(f"\nReady samples: {len(ready)}")
    for sid, _, _, _ in ready:
        print(f"  ✓ {sid}")

    if skipped:
        print(f"\nSkipped samples: {len(skipped)}")
        for sid, reason in skipped:
            print(f"  ✗ {sid}: {reason}")

    if not ready:
        print("\nNo samples ready for processing. Ensure step 03 "
              "(transcript export) has been run.")
        return

    # Load pre-built centroids (from step 02b) or build from scratch
    import pickle
    if os.path.exists(CENTROID_PATH):
        print(f"\nLoading pre-built centroids from {CENTROID_PATH}")
        t_cent = time.time()
        with open(CENTROID_PATH, 'rb') as f:
            centroid_bundle = pickle.load(f)
        sub_centroids = centroid_bundle['sub_centroids']
        sup_centroids = centroid_bundle['sup_centroids']
        gene_names = centroid_bundle['gene_names']
        print(f"  Loaded in {time.time()-t_cent:.1f}s "
              f"({len(sub_centroids)} subclass, {len(sup_centroids)} supertype centroids)")
    else:
        print(f"\nNo pre-built centroids found at {CENTROID_PATH}")
        print("Building from scratch (run step 02b first to pre-compute)...")
        t_cent = time.time()
        adatas_for_centroids = []
        for sid, h5ad_path, _, _ in ready:
            adata_tmp = ad.read_h5ad(h5ad_path)
            if 'qc_pass' in adata_tmp.obs.columns:
                adata_tmp = adata_tmp[adata_tmp.obs['qc_pass'].values.astype(bool)].copy()
            adatas_for_centroids.append(adata_tmp)
            print(f"  Loaded {sid}: {adata_tmp.n_obs:,} QC-pass cells", flush=True)

        combined = ad.concat(adatas_for_centroids, join='outer')
        del adatas_for_centroids
        print(f"  Combined: {combined.n_obs:,} cells x {combined.n_vars} genes")

        sub_centroids, _, gene_names = build_subclass_centroids(
            combined, top_n=CORR_CLASSIFIER_TOP_N)
        sup_centroids, _ = build_supertype_centroids(
            combined, top_n=CORR_CLASSIFIER_TOP_N)
        del combined
        print(f"  Centroids built in {time.time()-t_cent:.0f}s")

    # Process samples
    if N_WORKERS > 1 and len(ready) > 1:
        print(f"\nProcessing {len(ready)} samples with {N_WORKERS} workers...")
        with Pool(N_WORKERS, initializer=_init_worker,
                  initargs=(sub_centroids, sup_centroids, gene_names)) as pool:
            all_stats = pool.map(process_sample, ready)
    else:
        _init_worker(sub_centroids, sup_centroids, gene_names)
        print(f"\nProcessing {len(ready)} samples sequentially...")
        all_stats = [process_sample(args) for args in ready]

    # ── Summary ──
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    successful = [s for s in all_stats if s.get('status') == 'success']
    if not successful:
        print("No samples processed successfully.")
        return

    # Print per-sample table
    print(f"\n{'Sample':<10} {'Cells':>8} {'WC Dbl':>8} {'Resolv':>8} {'Persist':>8} "
          f"{'Insuff':>8} {'NucOnly':>8} {'Rate':>7} {'OldQC':>8} {'HybQC':>8} "
          f"{'Net±':>8}")
    print("-" * 108)

    total = {k: 0 for k in ['n_basic_qc', 'n_wc_doublets', 'n_resolved',
                             'n_persistent', 'n_insufficient', 'n_nuclear_only',
                             'n_old_corr_qc_pass', 'n_hybrid_qc_pass', 'net_rescued']}

    for s in sorted(successful, key=lambda x: x['sample_id']):
        print(f"  {s['sample_id']:<8} {s['n_basic_qc']:>8,} {s['n_wc_doublets']:>8,} "
              f"{s['n_resolved']:>8,} {s['n_persistent']:>8,} {s['n_insufficient']:>8,} "
              f"{s['n_nuclear_only']:>8,} {s['resolution_rate']:>6.1f}% "
              f"{s['n_old_corr_qc_pass']:>8,} {s['n_hybrid_qc_pass']:>8,} "
              f"{s['net_rescued']:>+8,}")
        for k in total:
            total[k] += s.get(k, 0)

    print("-" * 108)
    n_dbl = total['n_wc_doublets']
    rate = total['n_resolved'] / n_dbl * 100 if n_dbl > 0 else 0
    print(f"  {'TOTAL':<8} {total['n_basic_qc']:>8,} {n_dbl:>8,} "
          f"{total['n_resolved']:>8,} {total['n_persistent']:>8,} "
          f"{total['n_insufficient']:>8,} {total['n_nuclear_only']:>8,} "
          f"{rate:>6.1f}% {total['n_old_corr_qc_pass']:>8,} "
          f"{total['n_hybrid_qc_pass']:>8,} {total['net_rescued']:>+8,}")

    # High-UMI cell rescue summary
    total_hu_classified = sum(s.get('n_high_umi_classified', 0) for s in successful)
    total_hu_margin = sum(s.get('n_high_umi_margin_pass', 0) for s in successful)
    total_hu_total = sum(s.get('n_high_umi_total', 0) for s in successful)
    total_hu_rescued = sum(s.get('n_high_umi_rescued', 0) for s in successful)
    if total_hu_classified > 0:
        print(f"\n  High-UMI cell rescue:")
        print(f"    Total high-UMI cells:  {total_hu_total:,}")
        print(f"    Classified by corr:    {total_hu_classified:,}")
        print(f"    Margin pass:           {total_hu_margin:,} "
              f"({100*total_hu_margin/total_hu_classified:.1f}%)")
        print(f"    Rescued (hybrid_qc):   {total_hu_rescued:,}")

        # Aggregate subclass counts across all samples
        from collections import Counter
        all_sub_counts = Counter()
        for s in successful:
            for sub, cnt in s.get('high_umi_subclass_counts', {}).items():
                all_sub_counts[sub] += cnt
        if all_sub_counts:
            print(f"    Rescued by cell type:")
            for sub, cnt in all_sub_counts.most_common(10):
                print(f"      {sub:20s}: {cnt:,}")

    # Save summary CSV
    stats_df = pd.DataFrame(successful)
    # Drop dict column for CSV compatibility
    if 'high_umi_subclass_counts' in stats_df.columns:
        stats_df = stats_df.drop(columns=['high_umi_subclass_counts'])
    out_csv = os.path.join(os.path.dirname(H5AD_DIR), "presentation",
                           "nuclear_doublet_resolution_all_samples.csv")
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    stats_df.to_csv(out_csv, index=False)
    print(f"\nSaved: {out_csv}")

    elapsed = time.time() - t_start
    print(f"\nTotal elapsed: {elapsed:.0f}s ({elapsed/60:.1f} min)")


if __name__ == "__main__":
    main()
