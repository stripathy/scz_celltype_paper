#!/usr/bin/env python3
"""
MERSCOPE annotation pipeline using SEA-AD MTG taxonomy.

Adapts the SCZ Xenium annotation pipeline (steps 00-05) for MERSCOPE data
from Fang et al. (Dryad: doi:10.5061/dryad.x3ffbg7mw). Reuses the core
annotation modules (MapMyCells HANN + correlation classifier + depth model
+ BANKSY domains) while handling MERSCOPE-specific data format (triplet
CSVs instead of Xenium HDF5).

Supports both 250-gene and 4000-gene MERSCOPE panels. When processing both
panel sizes, the correlation classifier runs separately per panel size
(since gene sets differ, centroids must be built independently).

Pipeline steps:
  1. Load MERSCOPE samples from triplet CSVs into h5ad
  2. Basic QC (n_genes / total_counts thresholds, no control probes available)
  3. MapMyCells HANN mapping (SEA-AD MTG taxonomy) — reused from 02
  4. Two-stage correlation classifier — reused from 02b (grouped by panel size)
  5. Depth-from-pia prediction using trained GBR model
  6. BANKSY spatial domains + layer assignment + spatial smoothing
  7. Validation against pre-existing MERSCOPE cluster labels (L1/L2/L3)

Usage:
    # Run on a single sample (test mode)
    python3 -u merscope_annotate.py --sample H18.06.006_MTG_4000_rep1

    # Run on all 4K samples
    python3 -u merscope_annotate.py --panel-size 4000

    # Run on both 250 and 4K panels
    python3 -u merscope_annotate.py --panel-size all

    # Skip steps already completed
    python3 -u merscope_annotate.py --skip-mapmycells --skip-load --skip-qc

Requires:
  - cell_type_mapper (MapMyCells): pip install cell_type_mapper
  - SEA-AD precomputed stats: data/reference/precomputed_stats.20231120.sea_ad.MTG.h5
  - Comprehensive gene mapping: data/reference/gene_symbol_to_ensembl_comprehensive.json
  - Trained depth model: output/depth_model_normalized.pkl
"""

import os
import sys
import time
import json
import argparse
import tempfile
import pickle
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import sparse
from collections import defaultdict

# Pipeline config
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipeline_config import (
    BASE_DIR, OUTPUT_DIR, MODULES_DIR,
    PRECOMPUTED_STATS_PATH,
    GENE_MAPPING_COMPREHENSIVE_PATH,
    TAXONOMY_LEVELS,
    MAPMYCELLS_BOOTSTRAP_ITER, MAPMYCELLS_BOOTSTRAP_FACTOR,
    MAPMYCELLS_N_PER_UTILITY,
)

# Modules
sys.path.insert(0, MODULES_DIR)
from merscope_loading import (
    load_merscope_sample,
    discover_merscope_samples,
)
from correlation_classifier import (
    build_subclass_centroids,
    build_supertype_centroids,
    run_two_stage_classifier,
    flag_low_margin_cells,
    flag_doublet_cells,
)
from constants import SUBCLASS_TO_CLASS
from depth_model import (
    predict_depth,
    load_model as load_depth_model,
    assign_discrete_layers,
    smooth_layers_spatial,
    LAYER_BINS,
)
from banksy_domains import (
    preprocess_for_banksy,
    run_banksy,
    classify_banksy_domains,
)

# ──────────────────────────────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────────────────────────────

MERSCOPE_DIR = os.path.join(BASE_DIR, "data", "merscope_4k_probe_testing")
MERSCOPE_H5AD_DIR = os.path.join(OUTPUT_DIR, "merscope_h5ad")
MERSCOPE_OUTPUT_DIR = os.path.join(OUTPUT_DIR, "merscope_annotation")

# QC thresholds (simple, no control probes)
MIN_GENES = 10       # minimum genes detected per cell
MIN_COUNTS = 20      # minimum total counts per cell
MAX_COUNTS_PCTL = 99  # flag cells above this percentile of total_counts

# Depth model (reuse trained model from Xenium pipeline)
DEPTH_MODEL_PATH = os.path.join(OUTPUT_DIR, "depth_model_normalized.pkl")

# HANN class label → standard class name mapping
# (HANN taxonomy uses verbose names like "Neuronal: Glutamatergic")
HANN_CLASS_TO_CLASS = {
    'Neuronal: Glutamatergic': 'Glutamatergic',
    'Neuronal: GABAergic': 'GABAergic',
    'Non-neuronal and Non-neural': 'Non-neuronal',
}

# MERSCOPE L1 → standard class mapping
L1_TO_CLASS = {
    'EXC': 'Glutamatergic', 'INC': 'GABAergic',
    'lASC': 'Non-neuronal', 'lMGC': 'Non-neuronal',
    'lOGC': 'Non-neuronal', 'lOPC': 'Non-neuronal',
    'oENDO': 'Non-neuronal', 'oMURAL': 'Non-neuronal',
}

# MERSCOPE L2 → SEA-AD subclass mapping
L2_TO_SUBCLASS = {
    'eL2/3.IT': 'L2/3 IT', 'eL4/5.IT': 'L4 IT',
    'eL5.ET': 'L5 ET', 'eL5.IT': 'L5 IT',
    'eL5/6.NP': 'L5/6 NP', 'eL6.CT': 'L6 CT',
    'eL6.IT': 'L6 IT', 'eL6.IT.CAR3': 'L6 IT Car3',
    'eL6b': 'L6b',
    'iLAMP5': 'Lamp5', 'iSST': 'Sst', 'iVIP': 'Vip',
    'iPVALB': 'Pvalb',
    'lASC': 'Astrocyte', 'lMGC': 'Microglia-PVM',
    'lOGC': 'Oligodendrocyte', 'lOPC': 'OPC',
    'oENDO': 'Endothelial', 'oMURAL': 'VLMC',
}


# ──────────────────────────────────────────────────────────────────────
# Gene mapping
# ──────────────────────────────────────────────────────────────────────

def load_gene_mapping():
    """Load comprehensive gene symbol → Ensembl ID mapping (36K genes)."""
    path = GENE_MAPPING_COMPREHENSIVE_PATH
    if not os.path.exists(path):
        print(f"ERROR: Comprehensive gene mapping not found at {path}")
        print("Build it by aligning snRNAseq reference var_names with "
              "precomputed stats col_names.")
        sys.exit(1)
    with open(path) as f:
        mapping = json.load(f)
    return mapping


def convert_genes_to_ensembl(adata, gene_mapping):
    """Convert gene symbols to Ensembl IDs, dropping unmapped genes."""
    mappable = [g for g in adata.var_names if g in gene_mapping]
    adata_sub = adata[:, mappable].copy()
    adata_sub.var_names = [gene_mapping[g] for g in mappable]
    adata_sub.var_names_make_unique()
    return adata_sub


# ──────────────────────────────────────────────────────────────────────
# Step 1: Load data and create h5ad
# ──────────────────────────────────────────────────────────────────────

def step1_create_h5ad(samples):
    """Load MERSCOPE samples from CSV triplets and save as h5ad."""
    print("\n" + "=" * 70)
    print("STEP 1: Load MERSCOPE data → h5ad")
    print("=" * 70)

    os.makedirs(MERSCOPE_H5AD_DIR, exist_ok=True)

    for i, info in enumerate(samples):
        sid = info["sample_id"]
        t0 = time.time()
        print(f"\n[{i+1}/{len(samples)}] {sid}")

        try:
            adata = load_merscope_sample(info["prefix"])
            # Store panel_size in obs for later grouping
            adata.obs["panel_size"] = info["panel_size"]

            out_path = os.path.join(MERSCOPE_H5AD_DIR, f"{sid}_annotated.h5ad")
            adata.write_h5ad(out_path)

            elapsed = time.time() - t0
            file_mb = os.path.getsize(out_path) / 1e6
            print(f"  {adata.n_obs:,} cells x {adata.n_vars} genes → "
                  f"{file_mb:.0f} MB ({elapsed:.1f}s)")

        except Exception as e:
            import traceback
            print(f"  ERROR: {e}")
            traceback.print_exc()


# ──────────────────────────────────────────────────────────────────────
# Step 2: Basic QC
# ──────────────────────────────────────────────────────────────────────

def step2_qc(samples):
    """Apply basic QC filtering (n_genes, total_counts thresholds)."""
    print("\n" + "=" * 70)
    print("STEP 2: QC filtering")
    print(f"  min_genes={MIN_GENES}, min_counts={MIN_COUNTS}, "
          f"max_counts_pctl={MAX_COUNTS_PCTL}")
    print("=" * 70)

    for info in samples:
        sid = info["sample_id"]
        h5ad_path = os.path.join(MERSCOPE_H5AD_DIR, f"{sid}_annotated.h5ad")
        if not os.path.exists(h5ad_path):
            print(f"  {sid}: h5ad not found, skipping")
            continue

        adata = ad.read_h5ad(h5ad_path)
        n_total = adata.n_obs

        X = adata.X
        if sparse.issparse(X):
            total_counts = np.array(X.sum(axis=1)).flatten()
            n_genes_detected = np.array((X > 0).sum(axis=1)).flatten()
        else:
            total_counts = X.sum(axis=1)
            n_genes_detected = (X > 0).sum(axis=1)

        adata.obs["total_counts"] = total_counts
        adata.obs["n_genes"] = n_genes_detected

        max_counts_thresh = np.percentile(total_counts, MAX_COUNTS_PCTL)
        qc_pass = (
            (n_genes_detected >= MIN_GENES) &
            (total_counts >= MIN_COUNTS) &
            (total_counts <= max_counts_thresh)
        )
        adata.obs["qc_pass"] = qc_pass

        n_pass = int(qc_pass.sum())
        n_fail = n_total - n_pass
        print(f"  {sid}: {n_total:,} → {n_pass:,} pass, "
              f"{n_fail:,} fail ({100*n_fail/n_total:.1f}%)")

        adata.write_h5ad(h5ad_path)


# ──────────────────────────────────────────────────────────────────────
# Step 3: MapMyCells HANN mapping
# ──────────────────────────────────────────────────────────────────────

def step3_mapmycells(samples):
    """Run MapMyCells HANN label transfer on each MERSCOPE sample."""
    print("\n" + "=" * 70)
    print("STEP 3: MapMyCells HANN mapping")
    print("=" * 70)

    if not os.path.exists(PRECOMPUTED_STATS_PATH):
        print(f"ERROR: SEA-AD MTG precomputed stats not found at:")
        print(f"  {PRECOMPUTED_STATS_PATH}")
        return

    gene_mapping = load_gene_mapping()
    print(f"Comprehensive gene mapping: {len(gene_mapping)} entries")

    # Import MapMyCells runner from existing pipeline step
    pipeline_dir = os.path.dirname(os.path.abspath(__file__))
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "mapmycells_step",
        os.path.join(pipeline_dir, "02_run_mapmycells.py")
    )
    mmc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mmc)

    for i, info in enumerate(samples):
        sid = info["sample_id"]
        panel = info["panel_size"]
        h5ad_path = os.path.join(MERSCOPE_H5AD_DIR, f"{sid}_annotated.h5ad")
        if not os.path.exists(h5ad_path):
            print(f"  {sid}: h5ad not found, skipping")
            continue

        t0 = time.time()
        print(f"\n[{i+1}/{len(samples)}] {sid} ({panel}-gene panel)")

        adata = ad.read_h5ad(h5ad_path)
        n_total = adata.n_obs

        if "qc_pass" in adata.obs.columns:
            qc_mask = adata.obs["qc_pass"].values.astype(bool)
        else:
            qc_mask = np.ones(n_total, dtype=bool)

        n_pass = int(qc_mask.sum())
        print(f"  {n_total:,} total, {n_pass:,} QC pass")

        adata_pass = adata[qc_mask].copy()

        with tempfile.TemporaryDirectory() as tmpdir:
            adata_ensembl = convert_genes_to_ensembl(adata_pass, gene_mapping)
            n_mapped = adata_ensembl.shape[1]
            print(f"  Mapped {n_mapped}/{adata_pass.shape[1]} genes to Ensembl IDs "
                  f"({100*n_mapped/adata_pass.shape[1]:.1f}%)")

            tmp_h5ad = os.path.join(tmpdir, f"{sid}_qcpass.h5ad")
            adata_ensembl.write_h5ad(tmp_h5ad)

            print(f"  Running MapMyCells HANN mapping...")
            t_map = time.time()
            hdf5_path = mmc.run_mapmycells_on_sample(
                tmp_h5ad, tmpdir, n_processors=1
            )
            print(f"  MapMyCells done in {time.time()-t_map:.0f}s")

            labels_df = mmc.parse_mapmycells_output(hdf5_path)

            for level in TAXONOMY_LEVELS:
                col = f"{level}_label"
                conf_col = f"{level}_label_confidence"
                if col in labels_df.columns:
                    n_unique = labels_df[col].nunique()
                    mean_conf = labels_df[conf_col].mean()
                    print(f"    {level}: {n_unique} types, "
                          f"mean confidence={mean_conf:.3f}")

        # Write labels into h5ad
        for level in TAXONOMY_LEVELS:
            adata.obs[f"{level}_label"] = "Unassigned"
            adata.obs[f"{level}_label_confidence"] = np.float32(0.0)

        pass_indices = np.where(qc_mask)[0]
        for col in labels_df.columns:
            if col in adata.obs.columns:
                adata.obs.iloc[pass_indices,
                               adata.obs.columns.get_loc(col)] = \
                    labels_df[col].values

        adata.write_h5ad(h5ad_path)
        print(f"  Done in {time.time()-t0:.0f}s")


# ──────────────────────────────────────────────────────────────────────
# Step 4: Correlation classifier (grouped by panel size)
# ──────────────────────────────────────────────────────────────────────

def step4_correlation_classifier(samples):
    """Run two-stage correlation classifier, grouped by panel size.

    Samples with different panel sizes have different gene sets, so centroids
    must be built independently for each panel size group.
    """
    print("\n" + "=" * 70)
    print("STEP 4: Two-stage correlation classifier")
    print("=" * 70)

    top_n = 100
    qc_pctl = 1.0
    os.makedirs(MERSCOPE_OUTPUT_DIR, exist_ok=True)

    # Group samples by panel size
    panel_groups = defaultdict(list)
    for info in samples:
        panel_groups[info["panel_size"]].append(info)

    for panel_size, group_samples in sorted(panel_groups.items()):
        print(f"\n{'─'*50}")
        print(f"Panel size: {panel_size} genes ({len(group_samples)} samples)")
        print(f"{'─'*50}")

        # Load all samples in this panel group
        adatas = []
        for info in group_samples:
            sid = info["sample_id"]
            fpath = os.path.join(MERSCOPE_H5AD_DIR, f"{sid}_annotated.h5ad")
            if not os.path.exists(fpath):
                print(f"  {sid}: not found, skipping")
                continue

            adata = ad.read_h5ad(fpath)

            if "subclass_label" not in adata.obs.columns:
                print(f"  {sid}: no MapMyCells labels, skipping")
                continue

            if "qc_pass" in adata.obs.columns:
                adata = adata[adata.obs["qc_pass"].values.astype(bool)].copy()

            print(f"  {sid}: {adata.n_obs:,} cells x {adata.n_vars} genes",
                  flush=True)
            adatas.append(adata)

        if not adatas:
            print("  No samples with MapMyCells labels. Run step 3 first.")
            continue

        print("  Concatenating...", flush=True)
        combined = ad.concat(adatas, join="inner")
        print(f"  Combined: {combined.n_obs:,} cells x {combined.n_vars} genes")

        # Build centroids
        print(f"\n  Building subclass centroids (top-{top_n})...")
        sub_centroids, sub_counts, gene_names = build_subclass_centroids(
            combined, top_n=top_n)

        print(f"\n  Building supertype centroids (top-{top_n})...")
        sup_centroids, sup_counts = build_supertype_centroids(
            combined, top_n=top_n)

        # Save centroids per panel size
        centroid_path = os.path.join(
            MERSCOPE_OUTPUT_DIR, f"merscope_{panel_size}_centroids.pkl")
        with open(centroid_path, 'wb') as f:
            pickle.dump({
                'sub_centroids': sub_centroids,
                'sup_centroids': sup_centroids,
                'gene_names': gene_names,
                'sub_counts': sub_counts,
                'sup_counts': sup_counts,
            }, f)
        print(f"  Saved: {centroid_path}")

        # Run classifier
        print("\n  Running two-stage hierarchical classifier...")
        results = run_two_stage_classifier(
            combined, sub_centroids, sup_centroids, gene_names)

        results['corr_class'] = results['corr_subclass'].map(
            lambda x: SUBCLASS_TO_CLASS.get(x, 'Unknown'))

        # QC flagging
        sample_ids_arr = combined.obs['sample_id'].astype(str).values
        margins = results['corr_subclass_margin'].values
        corr_qc_pass, thresholds = flag_low_margin_cells(
            margins, sample_ids_arr, percentile=qc_pctl)

        # Doublet detection
        class_labels = results['corr_class'].values
        doublet_suspect, doublet_type, doublet_stats = flag_doublet_cells(
            combined, class_labels, SUBCLASS_TO_CLASS, panel='xenium')

        corr_qc_pass = corr_qc_pass & (~doublet_suspect)
        results['corr_qc_pass'] = corr_qc_pass
        results['doublet_suspect'] = doublet_suspect
        results['doublet_type'] = doublet_type

        # Summary
        hann_sub = combined.obs['subclass_label'].astype(str).values
        corr_sub = results['corr_subclass'].values
        agree_sub = (hann_sub == corr_sub).mean()
        print(f"\n  Subclass agreement (HANN vs Corr): {100*agree_sub:.1f}%")

        n_hann_supertypes = combined.obs['supertype_label'].astype(str).nunique()
        n_corr_supertypes = results['corr_supertype'].nunique()
        print(f"  HANN supertypes: {n_hann_supertypes}, "
              f"Corr supertypes: {n_corr_supertypes}")

        # Write results back to h5ad files
        cell_sample_ids = combined.obs['sample_id'].astype(str).values

        for sid in sorted(set(cell_sample_ids)):
            fpath = os.path.join(MERSCOPE_H5AD_DIR, f"{sid}_annotated.h5ad")
            adata = ad.read_h5ad(fpath)

            sid_mask = cell_sample_ids == sid

            if "qc_pass" in adata.obs.columns:
                qc_mask = adata.obs["qc_pass"].values.astype(bool)
            else:
                qc_mask = np.ones(adata.n_obs, dtype=bool)

            for col in ['corr_subclass', 'corr_supertype', 'corr_class',
                         'doublet_type']:
                adata.obs[col] = 'Unassigned'
            for col in ['corr_subclass_corr', 'corr_subclass_margin',
                         'corr_supertype_corr']:
                adata.obs[col] = np.float32(0.0)
            adata.obs['corr_qc_pass'] = False
            adata.obs['doublet_suspect'] = False

            sid_results = results.loc[sid_mask]
            qc_indices = np.where(qc_mask)[0]

            for col in sid_results.columns:
                if col not in adata.obs.columns:
                    continue
                vals = sid_results[col].values
                if vals.dtype == object or col in ['corr_subclass', 'corr_supertype',
                                                     'corr_class', 'doublet_type']:
                    adata.obs[col] = adata.obs[col].astype(str)
                    adata.obs.iloc[qc_indices,
                                   adata.obs.columns.get_loc(col)] = vals.astype(str)
                elif vals.dtype == bool or col in ['corr_qc_pass', 'doublet_suspect']:
                    adata.obs.iloc[qc_indices,
                                   adata.obs.columns.get_loc(col)] = vals.astype(bool)
                else:
                    adata.obs[col] = adata.obs[col].astype(np.float32)
                    adata.obs.iloc[qc_indices,
                                   adata.obs.columns.get_loc(col)] = vals.astype(np.float32)

            adata.write_h5ad(fpath)
            print(f"  Written: {sid}")


# ──────────────────────────────────────────────────────────────────────
# Step 5: Depth prediction
# ──────────────────────────────────────────────────────────────────────

def step5_depth_prediction(samples):
    """Predict normalized cortical depth from pia using trained GBR model."""
    print("\n" + "=" * 70)
    print("STEP 5: Depth-from-pia prediction")
    print("=" * 70)

    if not os.path.exists(DEPTH_MODEL_PATH):
        print(f"ERROR: Depth model not found at {DEPTH_MODEL_PATH}")
        print("Train it first via code/pipeline/04_run_depth_prediction.py")
        return

    model_bundle = load_depth_model(DEPTH_MODEL_PATH)
    print(f"  Loaded model: K={model_bundle['K']}, "
          f"R²={model_bundle['test_r2']:.3f}")

    for i, info in enumerate(samples):
        sid = info["sample_id"]
        fpath = os.path.join(MERSCOPE_H5AD_DIR, f"{sid}_annotated.h5ad")
        if not os.path.exists(fpath):
            print(f"  {sid}: not found, skipping")
            continue

        t0 = time.time()
        adata = ad.read_h5ad(fpath)

        if "qc_pass" in adata.obs.columns:
            qc_mask = adata.obs["qc_pass"].values.astype(bool)
        else:
            qc_mask = np.ones(adata.n_obs, dtype=bool)

        n_pass = int(qc_mask.sum())
        adata_pass = adata[qc_mask].copy()

        # Use corr_subclass if available (more accurate), else HANN subclass
        if "corr_subclass" in adata_pass.obs.columns:
            subclass_col = "corr_subclass"
        else:
            subclass_col = "subclass_label"

        pred_depth = predict_depth(adata_pass, model_bundle,
                                   subclass_col=subclass_col)

        # Write back
        adata.obs["predicted_norm_depth"] = np.float32(np.nan)
        qc_indices = np.where(qc_mask)[0]
        adata.obs.iloc[qc_indices,
                       adata.obs.columns.get_loc("predicted_norm_depth")] = \
            pred_depth.astype(np.float32)

        adata.write_h5ad(fpath)

        # Summary stats
        median_d = np.nanmedian(pred_depth)
        in_cortex = ((pred_depth >= 0) & (pred_depth <= 1)).sum()
        pct_cortex = 100 * in_cortex / len(pred_depth)

        elapsed = time.time() - t0
        print(f"  [{i+1}/{len(samples)}] {sid}: {n_pass:,} cells, "
              f"median_depth={median_d:.3f}, cortical={pct_cortex:.0f}% "
              f"({elapsed:.0f}s)")


# ──────────────────────────────────────────────────────────────────────
# Step 6: BANKSY spatial domains + layer assignment
# ──────────────────────────────────────────────────────────────────────

def step6_spatial_domains(samples):
    """Run BANKSY clustering, domain classification, and layer assignment."""
    print("\n" + "=" * 70)
    print("STEP 6: BANKSY spatial domains + layer assignment")
    print("=" * 70)

    for i, info in enumerate(samples):
        sid = info["sample_id"]
        fpath = os.path.join(MERSCOPE_H5AD_DIR, f"{sid}_annotated.h5ad")
        if not os.path.exists(fpath):
            print(f"  {sid}: not found, skipping")
            continue

        t0 = time.time()
        adata = ad.read_h5ad(fpath)
        n_total = adata.n_obs

        if "predicted_norm_depth" not in adata.obs.columns:
            print(f"  {sid}: no depth predictions, skipping (run step 4 first)")
            continue

        if "qc_pass" in adata.obs.columns:
            qc_mask = adata.obs["qc_pass"].values.astype(bool)
        else:
            qc_mask = np.ones(n_total, dtype=bool)

        n_pass = int(qc_mask.sum())
        adata_pass = adata[qc_mask].copy()
        print(f"\n  [{i+1}/{len(samples)}] {sid}: {n_pass:,} QC-pass cells")

        # BANKSY clustering
        t_banksy = time.time()
        adata_b = preprocess_for_banksy(adata_pass)
        banksy_labels = run_banksy(adata_b)
        del adata_b
        n_clusters = len(np.unique(banksy_labels))
        print(f"    BANKSY: {n_clusters} clusters ({time.time()-t_banksy:.0f}s)")

        # Domain classification
        domains, is_l1, cluster_info = classify_banksy_domains(
            adata_pass, banksy_labels)
        for d in ["Cortical", "Vascular", "WM"]:
            n_d = int((domains == d).sum())
            print(f"    {d}: {n_d:,} ({100*n_d/n_pass:.1f}%)", end="  ")
        print(f"L1_border: {int(is_l1.sum()):,}")

        # Layer assignment
        pred_depth = adata_pass.obs["predicted_norm_depth"].values
        depth_layers = assign_discrete_layers(pred_depth)

        # Hybrid: depth bins + Vascular override
        combined_layers = depth_layers.copy()
        combined_layers[domains == "Vascular"] = "Vascular"

        # Spatial smoothing
        t_smooth = time.time()
        smoothed_layers = smooth_layers_spatial(
            coords=adata_pass.obsm["spatial"],
            layers=combined_layers,
            domains=domains,
            is_l1_banksy=is_l1,
            depths=pred_depth,
            verbose=False,
        )
        print(f"    Spatial smoothing: {time.time()-t_smooth:.0f}s")

        # Write back to full adata
        qc_idx = np.where(qc_mask)[0]

        full_banksy = np.full(n_total, -1, dtype=int)
        full_banksy[qc_idx] = banksy_labels
        adata.obs["banksy_cluster"] = full_banksy

        full_domain = np.full(n_total, "", dtype=object)
        full_domain[qc_idx] = domains
        adata.obs["banksy_domain"] = full_domain

        full_is_l1 = np.zeros(n_total, dtype=bool)
        full_is_l1[qc_idx] = is_l1
        adata.obs["banksy_is_l1"] = full_is_l1

        full_layer = np.full(n_total, "Unassigned", dtype=object)
        full_layer[qc_idx] = smoothed_layers
        adata.obs["layer"] = full_layer

        adata.write_h5ad(fpath)

        # Layer distribution
        for lname in ["L1", "L2/3", "L4", "L5", "L6", "WM", "Vascular"]:
            n_l = int((smoothed_layers == lname).sum())
            print(f"    {lname}: {n_l:,} ({100*n_l/n_pass:.1f}%)")

        print(f"    Done in {time.time()-t0:.0f}s")


# ──────────────────────────────────────────────────────────────────────
# Step 7: Validation + panel size comparison
# ──────────────────────────────────────────────────────────────────────

def step7_validate(samples):
    """Compare annotations against MERSCOPE labels and across panel sizes."""
    print("\n" + "=" * 70)
    print("STEP 7: Validation against MERSCOPE labels + panel comparison")
    print("=" * 70)

    os.makedirs(MERSCOPE_OUTPUT_DIR, exist_ok=True)
    all_records = []

    for info in samples:
        sid = info["sample_id"]
        panel = info["panel_size"]
        fpath = os.path.join(MERSCOPE_H5AD_DIR, f"{sid}_annotated.h5ad")
        if not os.path.exists(fpath):
            continue

        adata = ad.read_h5ad(fpath)

        has_hann = "subclass_label" in adata.obs.columns
        has_corr = "corr_subclass" in adata.obs.columns
        has_merscope = "cluster_L2" in adata.obs.columns

        if not has_merscope:
            print(f"  {sid}: no MERSCOPE cluster labels, skipping")
            continue

        if "qc_pass" in adata.obs.columns:
            mask = adata.obs["qc_pass"].values.astype(bool)
        else:
            mask = np.ones(adata.n_obs, dtype=bool)

        n_pass = int(mask.sum())
        print(f"\n  {sid} ({panel}-gene): {n_pass:,} QC-pass cells")

        # Map MERSCOPE L2 to SEA-AD subclass
        merscope_l2 = adata.obs.loc[mask, "cluster_L2"].astype(str).values
        merscope_subclass = np.array([L2_TO_SUBCLASS.get(l, l) for l in merscope_l2])
        merscope_class = np.array([
            SUBCLASS_TO_CLASS.get(L2_TO_SUBCLASS.get(l, ''), 'Unknown')
            for l in merscope_l2
        ])

        # Annotation quality metrics
        if has_hann:
            hann_class_raw = adata.obs.loc[mask, "class_label"].astype(str).values
            # Map HANN class names to standard names
            hann_class = np.array([
                HANN_CLASS_TO_CLASS.get(c, c) for c in hann_class_raw])
            hann_sub = adata.obs.loc[mask, "subclass_label"].astype(str).values
            hann_sup = adata.obs.loc[mask, "supertype_label"].astype(str).values
            hann_conf_sub = adata.obs.loc[mask, "subclass_label_confidence"].astype(float).values
            hann_conf_sup = adata.obs.loc[mask, "supertype_label_confidence"].astype(float).values

            class_agree = (hann_class == merscope_class).mean()
            sub_agree = (hann_sub == merscope_subclass).mean()
            n_supertypes = len(set(hann_sup) - {"Unassigned"})
            mean_conf_sub = hann_conf_sub.mean()
            mean_conf_sup = hann_conf_sup.mean()

            print(f"    HANN class agreement:    {100*class_agree:.1f}%")
            print(f"    HANN subclass agreement: {100*sub_agree:.1f}%")
            print(f"    HANN supertypes found:   {n_supertypes}")
            print(f"    HANN confidence (sub/sup): {mean_conf_sub:.3f} / {mean_conf_sup:.3f}")

            all_records.append({
                "sample_id": sid, "panel_size": panel,
                "donor": info["donor"], "region": info["region"],
                "n_cells": n_pass, "n_genes": adata.n_vars,
                "method": "HANN",
                "class_agreement": class_agree,
                "subclass_agreement": sub_agree,
                "n_supertypes": n_supertypes,
                "mean_subclass_confidence": mean_conf_sub,
                "mean_supertype_confidence": mean_conf_sup,
            })

        if has_corr:
            corr_class = adata.obs.loc[mask, "corr_class"].astype(str).values
            corr_sub = adata.obs.loc[mask, "corr_subclass"].astype(str).values
            corr_sup = adata.obs.loc[mask, "corr_supertype"].astype(str).values
            corr_corr_sub = adata.obs.loc[mask, "corr_subclass_corr"].astype(float).values

            class_agree_c = (corr_class == merscope_class).mean()
            sub_agree_c = (corr_sub == merscope_subclass).mean()
            n_supertypes_c = len(set(corr_sup) - {"Unassigned"})
            mean_corr = corr_corr_sub.mean()

            print(f"    Corr class agreement:    {100*class_agree_c:.1f}%")
            print(f"    Corr subclass agreement: {100*sub_agree_c:.1f}%")
            print(f"    Corr supertypes found:   {n_supertypes_c}")
            print(f"    Corr mean correlation:   {mean_corr:.3f}")

            all_records.append({
                "sample_id": sid, "panel_size": panel,
                "donor": info["donor"], "region": info["region"],
                "n_cells": n_pass, "n_genes": adata.n_vars,
                "method": "Correlation",
                "class_agreement": class_agree_c,
                "subclass_agreement": sub_agree_c,
                "n_supertypes": n_supertypes_c,
                "mean_subclass_confidence": mean_corr,
                "mean_supertype_confidence": np.nan,
            })

        # Per-subclass breakdown (HANN)
        if has_hann:
            for ms_sub in sorted(set(merscope_subclass)):
                ms_mask = merscope_subclass == ms_sub
                n_ms = ms_mask.sum()
                if n_ms < 5:
                    continue
                hann_sub_ms = hann_sub[ms_mask]
                agree = (hann_sub_ms == ms_sub).mean()
                top = pd.Series(hann_sub_ms).value_counts().head(3)
                top_str = ", ".join([f"{v}={c}" for v, c in top.items()])
                all_records.append({
                    "sample_id": sid, "panel_size": panel,
                    "donor": info["donor"], "region": info["region"],
                    "n_cells": int(n_ms), "n_genes": adata.n_vars,
                    "method": "HANN_per_subclass",
                    "merscope_subclass": ms_sub,
                    "class_agreement": np.nan,
                    "subclass_agreement": agree,
                    "n_supertypes": np.nan,
                    "mean_subclass_confidence": np.nan,
                    "mean_supertype_confidence": np.nan,
                    "top_matches": top_str,
                })

    # Save and summarize
    if all_records:
        val_df = pd.DataFrame(all_records)
        val_path = os.path.join(MERSCOPE_OUTPUT_DIR,
                                "merscope_annotation_validation.csv")
        val_df.to_csv(val_path, index=False)
        print(f"\n  Saved: {val_path}")

        # Panel comparison summary
        summary_df = val_df[val_df["method"].isin(["HANN", "Correlation"])]
        if len(summary_df) > 0:
            print(f"\n  {'─'*60}")
            print(f"  PANEL SIZE COMPARISON")
            print(f"  {'─'*60}")
            for method in ["HANN", "Correlation"]:
                mdf = summary_df[summary_df["method"] == method]
                if len(mdf) == 0:
                    continue
                print(f"\n  {method}:")
                for _, row in mdf.iterrows():
                    print(f"    {row['sample_id']:40s} | "
                          f"class={100*row['class_agreement']:.1f}% | "
                          f"subclass={100*row['subclass_agreement']:.1f}% | "
                          f"supertypes={row['n_supertypes']:.0f}")


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="MERSCOPE annotation pipeline (SEA-AD MTG taxonomy)"
    )
    parser.add_argument("--panel-size", type=str, default="all",
                        help="Panel size: 250, 4000, or 'all' (default: all)")
    parser.add_argument("--sample", type=str, default=None,
                        help="Process only this sample ID (repeatable)",
                        action="append")
    parser.add_argument("--skip-load", action="store_true")
    parser.add_argument("--skip-qc", action="store_true")
    parser.add_argument("--skip-mapmycells", action="store_true")
    parser.add_argument("--skip-corr", action="store_true")
    parser.add_argument("--skip-depth", action="store_true")
    parser.add_argument("--skip-domains", action="store_true")
    parser.add_argument("--skip-validate", action="store_true")
    args = parser.parse_args()

    t_start = time.time()

    # Parse panel size
    if args.panel_size.lower() == "all":
        panel_filter = None
    else:
        panel_filter = int(args.panel_size)

    print("=" * 70)
    print("MERSCOPE Annotation Pipeline")
    print(f"  Panel size: {args.panel_size}")
    print(f"  Data dir: {MERSCOPE_DIR}")
    print(f"  Output: {MERSCOPE_H5AD_DIR}")
    print("=" * 70)

    # Discover samples
    samples = discover_merscope_samples(MERSCOPE_DIR, panel_size=panel_filter)
    if not samples:
        print(f"ERROR: No MERSCOPE samples found")
        sys.exit(1)

    print(f"\nFound {len(samples)} samples:")
    for s in samples:
        print(f"  {s['sample_id']:40s}  {s['panel_size']:>4} genes  "
              f"{s['donor']} {s['region']}")

    if args.sample:
        samples = [s for s in samples
                    if s["sample_id"] in args.sample]
        if not samples:
            print(f"ERROR: No matching samples found for {args.sample}")
            sys.exit(1)
        print(f"\nFiltered to {len(samples)} sample(s)")

    # Run pipeline
    if not args.skip_load:
        step1_create_h5ad(samples)

    if not args.skip_qc:
        step2_qc(samples)

    if not args.skip_mapmycells:
        step3_mapmycells(samples)

    if not args.skip_corr:
        step4_correlation_classifier(samples)

    if not args.skip_depth:
        step5_depth_prediction(samples)

    if not args.skip_domains:
        step6_spatial_domains(samples)

    if not args.skip_validate:
        step7_validate(samples)

    total_time = time.time() - t_start
    print(f"\n{'=' * 70}")
    print(f"Pipeline complete in {total_time:.0f}s ({total_time/60:.1f} min)")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
