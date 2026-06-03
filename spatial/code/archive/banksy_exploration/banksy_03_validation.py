#!/usr/bin/env python3
"""
BANKSY Phase 4: Validation against laminar markers and MERFISH reference.

Three validation approaches:
  1. Laminar marker enrichment: compare BANKSY clusters vs depth-bin layers
     for known layer marker genes in the Xenium 300-gene panel
  2. MERFISH benchmark: run BANKSY on selected MERFISH sections with ground
     truth manual layer annotations, compute ARI/NMI/F1
  3. Domain concordance: compare BANKSY domain classification to current pipeline

Usage:
    python3 -u code/analysis/banksy_03_validation.py
    python3 -u code/analysis/banksy_03_validation.py --skip-merfish
"""

import os
import sys
import time
import argparse
import warnings
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from scipy.stats import kendalltau, f_oneway
from sklearn.metrics import (
    adjusted_rand_score, normalized_mutual_info_score,
    confusion_matrix, classification_report,
)

# BANKSY imports
from banksy.initialize_banksy import initialize_banksy
from banksy.embed_banksy import generate_banksy_matrix
from banksy.cluster_methods import run_Leiden_partition
from banksy_utils.umap_pca import pca_umap

# Project imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    H5AD_DIR, MERFISH_PATH, LAYER_COLORS, LAYER_ORDER, CORTICAL_LAYERS,
    SUBCLASS_TO_CLASS, CLASS_COLORS, SAMPLE_TO_DX,
)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "modules"))
from depth_model import LAYER_BINS
from spatial_domains import VASCULAR_TYPES, NON_NEURONAL_TYPES

# ── Constants ─────────────────────────────────────────────────────────

OUT_DIR = os.path.expanduser("~/Github/SCZ_Xenium/output/banksy")
os.makedirs(OUT_DIR, exist_ok=True)

# Best BANKSY parameters from Phase 1
BANKSY_LAMBDA = 0.8
BANKSY_RESOLUTION = 0.3
K_GEOM = 15
PCA_DIMS = [20]

# Xenium laminar marker genes (present in 300-gene panel)
LAYER_MARKERS = {
    "L2/3": ["CUX2", "LAMP5", "RASGRF2"],
    "L4": ["RORB"],
    "L5": ["ADCYAP1", "PCP4", "HTR2C"],
    "L6": ["NR4A2"],
    "WM": ["MBP", "PLP1", "MOG"],
}

# Flatten for easy access
ALL_MARKERS = []
MARKER_TO_LAYER = {}
for layer, genes in LAYER_MARKERS.items():
    for g in genes:
        ALL_MARKERS.append(g)
        MARKER_TO_LAYER[g] = layer

# Domain classification thresholds (matching banksy_02)
VASCULAR_THRESH = 0.50
MENINGEAL_NN_THRESH = 0.50
MENINGEAL_DEPTH_THRESH = 0.20
WM_OLIGO_THRESH = 0.40
WM_DEPTH_THRESH = 0.80


# ── BANKSY helpers ───────────────────────────────────────────────────

def run_banksy(adata_b, lam=BANKSY_LAMBDA, res=BANKSY_RESOLUTION):
    """Run BANKSY, return cluster labels."""
    coord_keys = ("x", "y", "spatial")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        banksy_dict = initialize_banksy(
            adata_b, coord_keys, K_GEOM,
            nbr_weight_decay="scaled_gaussian", max_m=1,
            plt_edge_hist=False, plt_nbr_weights=False, plt_agf_angles=False,
        )
        banksy_dict, _ = generate_banksy_matrix(
            adata_b, banksy_dict, [lam], max_m=1
        )
        pca_umap(banksy_dict, pca_dims=PCA_DIMS, add_umap=False)
        results_df, _ = run_Leiden_partition(
            banksy_dict, [res],
            num_nn=50, num_iterations=-1, partition_seed=42,
            match_labels=False,
        )
    labels = results_df.iloc[0]["labels"]
    if hasattr(labels, "dense"):
        labels = labels.dense
    return np.asarray(labels).astype(int)


def preprocess_for_banksy(adata):
    """Normalize → log1p → z-score for BANKSY."""
    adata_b = adata.copy()
    adata_b.layers["counts"] = adata_b.X.copy()
    sc.pp.normalize_total(adata_b, target_sum=1e4)
    sc.pp.log1p(adata_b)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.scale(adata_b)
    return adata_b


# ── Validation 1: Laminar Marker Enrichment ──────────────────────────

def validate_laminar_markers(sample_id="Br6437"):
    """Compare layer-marker enrichment in BANKSY clusters vs depth-bin layers.

    For each marker gene, computes mean expression per BANKSY cluster and
    per depth-bin layer. Then computes ANOVA F-statistics for both groupings
    to assess which produces better laminar separation.
    """
    print(f"\n{'='*60}")
    print(f"Validation 1: Laminar marker enrichment ({sample_id})")
    print(f"{'='*60}")

    # Load and preprocess
    path = os.path.join(H5AD_DIR, f"{sample_id}_annotated.h5ad")
    adata = ad.read_h5ad(path)
    qc_mask = adata.obs["hybrid_qc_pass"].values.astype(bool)
    adata = adata[qc_mask].copy()
    print(f"  {adata.n_obs:,} QC-pass cells")

    # Get depth-bin layers
    pred_depth = pd.to_numeric(adata.obs["predicted_norm_depth"],
                                errors="coerce").values
    depth_layers = np.full(len(pred_depth), "Unknown", dtype=object)
    for layer_name, (lo, hi) in LAYER_BINS.items():
        mask = (pred_depth >= lo) & (pred_depth < hi)
        depth_layers[mask] = layer_name

    # Run BANKSY
    print("  Running BANKSY...")
    adata_b = preprocess_for_banksy(adata)
    banksy_labels = run_banksy(adata_b)
    n_clusters = len(np.unique(banksy_labels))
    print(f"  BANKSY: {n_clusters} clusters")

    # Order clusters by median depth
    unique_cl = np.unique(banksy_labels)
    cl_depth = {cl: np.nanmedian(pred_depth[banksy_labels == cl]) for cl in unique_cl}
    cl_order = sorted(unique_cl, key=lambda c: cl_depth[c])
    cl_rank = {cl: i for i, cl in enumerate(cl_order)}

    # Get raw counts (not z-scored) for marker expression
    # Use log-normalized expression for interpretability
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)

    # Check which markers are in the panel
    available_markers = [g for g in ALL_MARKERS if g in adata_norm.var_names]
    print(f"  Available markers: {len(available_markers)}/{len(ALL_MARKERS)}")

    # Get expression matrix for markers
    marker_idx = [list(adata_norm.var_names).index(g) for g in available_markers]
    if hasattr(adata_norm.X, "toarray"):
        X_markers = adata_norm.X[:, marker_idx].toarray()
    else:
        X_markers = adata_norm.X[:, marker_idx]

    # Compute mean expression per BANKSY cluster (depth-ordered)
    banksy_means = np.zeros((n_clusters, len(available_markers)))
    for i, cl in enumerate(cl_order):
        mask = banksy_labels == cl
        banksy_means[i] = X_markers[mask].mean(axis=0)

    # Compute mean expression per depth-bin layer
    layer_names = ["L1", "L2/3", "L4", "L5", "L6", "WM"]
    layer_means = np.zeros((len(layer_names), len(available_markers)))
    for i, ln in enumerate(layer_names):
        mask = depth_layers == ln
        if mask.sum() > 0:
            layer_means[i] = X_markers[mask].mean(axis=0)

    # ANOVA F-statistics: marker separation by BANKSY vs depth-bin
    print(f"\n  {'Marker':<12} {'Layer':<6} | {'F(BANKSY)':>10} {'F(DepthBin)':>12} | Better")
    print(f"  {'-'*60}")
    f_results = []
    for j, gene in enumerate(available_markers):
        expr = X_markers[:, j]
        expected_layer = MARKER_TO_LAYER[gene]

        # BANKSY ANOVA
        groups_b = [expr[banksy_labels == cl] for cl in cl_order]
        groups_b = [g for g in groups_b if len(g) > 10]
        if len(groups_b) >= 2:
            f_banksy, p_banksy = f_oneway(*groups_b)
        else:
            f_banksy, p_banksy = 0, 1

        # Depth-bin ANOVA
        groups_d = [expr[depth_layers == ln] for ln in layer_names]
        groups_d = [g for g in groups_d if len(g) > 10]
        if len(groups_d) >= 2:
            f_depth, p_depth = f_oneway(*groups_d)
        else:
            f_depth, p_depth = 0, 1

        better = "BANKSY" if f_banksy > f_depth else "Depth-bin"
        print(f"  {gene:<12} {expected_layer:<6} | {f_banksy:>10.0f} {f_depth:>12.0f} | {better}")
        f_results.append({
            "gene": gene, "expected_layer": expected_layer,
            "f_banksy": f_banksy, "f_depth": f_depth,
            "better": better,
        })

    # Summary
    n_banksy_better = sum(1 for r in f_results if r["better"] == "BANKSY")
    print(f"\n  Summary: BANKSY better for {n_banksy_better}/{len(f_results)} markers")

    # ── Figure: Marker expression heatmaps ──
    fig, axes = plt.subplots(1, 3, figsize=(24, 8), gridspec_kw={"width_ratios": [1, 1, 0.6]})

    # Panel 1: BANKSY cluster × marker heatmap
    ax1 = axes[0]
    im1 = ax1.imshow(banksy_means, aspect="auto", cmap="YlOrRd", interpolation="nearest")
    ax1.set_xticks(range(len(available_markers)))
    ax1.set_xticklabels(available_markers, rotation=45, ha="right", fontsize=12)
    ax1.set_yticks(range(n_clusters))
    ax1.set_yticklabels([f"Cl{cl} (d={cl_depth[cl]:.2f})" for cl in cl_order],
                        fontsize=11)
    ax1.set_title("BANKSY Clusters (depth-ordered)", fontsize=16, fontweight="bold")
    ax1.set_ylabel("Cluster", fontsize=14)
    plt.colorbar(im1, ax=ax1, shrink=0.6, label="Mean log1p expression")

    # Highlight expected peak positions
    for j, gene in enumerate(available_markers):
        expected = MARKER_TO_LAYER[gene]
        # Find which BANKSY cluster has most overlap with expected layer
        best_cl_idx = np.argmax(banksy_means[:, j])
        ax1.plot(j, best_cl_idx, "k*", markersize=8)

    # Panel 2: Depth-bin layer × marker heatmap
    ax2 = axes[1]
    im2 = ax2.imshow(layer_means, aspect="auto", cmap="YlOrRd", interpolation="nearest")
    ax2.set_xticks(range(len(available_markers)))
    ax2.set_xticklabels(available_markers, rotation=45, ha="right", fontsize=12)
    ax2.set_yticks(range(len(layer_names)))
    ax2.set_yticklabels(layer_names, fontsize=13)
    ax2.set_title("Depth-Bin Layers", fontsize=16, fontweight="bold")
    ax2.set_ylabel("Layer", fontsize=14)
    plt.colorbar(im2, ax=ax2, shrink=0.6, label="Mean log1p expression")

    # Panel 3: F-statistic comparison
    ax3 = axes[2]
    x = range(len(f_results))
    f_b = [r["f_banksy"] for r in f_results]
    f_d = [r["f_depth"] for r in f_results]
    ax3.barh([i - 0.2 for i in x], f_b, 0.35, label="BANKSY", color="#4CAF50", alpha=0.8)
    ax3.barh([i + 0.2 for i in x], f_d, 0.35, label="Depth-bin", color="#2196F3", alpha=0.8)
    ax3.set_yticks(x)
    ax3.set_yticklabels([r["gene"] for r in f_results], fontsize=11)
    ax3.set_xlabel("ANOVA F-statistic", fontsize=13)
    ax3.set_title("Marker Separation", fontsize=16, fontweight="bold")
    ax3.legend(fontsize=12)
    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(False)
    ax3.invert_yaxis()

    fig.suptitle(f"{sample_id} — Laminar Marker Validation\n"
                 f"(BANKSY better for {n_banksy_better}/{len(f_results)} markers)",
                 fontsize=20, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])

    fig_path = os.path.join(OUT_DIR, f"validation_markers_{sample_id}.png")
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Saved: {fig_path}")

    # Save results
    df = pd.DataFrame(f_results)
    csv_path = os.path.join(OUT_DIR, f"validation_markers_{sample_id}.csv")
    df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")

    return f_results, banksy_labels, cl_order, cl_depth


# ── Validation 2: MERFISH Benchmark ──────────────────────────────────

def validate_merfish(n_sections=3, max_cells_per_section=30000):
    """Run BANKSY on selected MERFISH sections with ground-truth layer labels.

    Selects sections with the best layer annotation coverage, runs BANKSY,
    compares clusters to manual layers.
    """
    print(f"\n{'='*60}")
    print(f"Validation 2: MERFISH benchmark ({n_sections} sections)")
    print(f"{'='*60}")

    # Load MERFISH
    print("  Loading MERFISH h5ad...")
    t0 = time.time()
    merfish = ad.read_h5ad(MERFISH_PATH)
    print(f"  Loaded: {merfish.n_obs:,} cells × {merfish.n_vars} genes ({time.time()-t0:.0f}s)")

    # Find sections with best layer annotation coverage
    layer_col = "Layer annotation"
    if layer_col not in merfish.obs.columns:
        print(f"  ERROR: '{layer_col}' column not found in MERFISH obs")
        print(f"  Available columns: {list(merfish.obs.columns)[:20]}...")
        return None

    has_layer = ~merfish.obs[layer_col].isna() & (merfish.obs[layer_col] != "")
    print(f"  Cells with layer annotation: {has_layer.sum():,} / {merfish.n_obs:,}")

    # Group by section (or donor) — find sections with most annotated cells
    if "section_id" in merfish.obs.columns:
        section_col = "section_id"
    elif "Donor ID" in merfish.obs.columns:
        section_col = "Donor ID"
    else:
        section_col = merfish.obs.columns[
            merfish.obs.columns.str.contains("section|donor|sample", case=False)
        ][0]

    section_counts = merfish.obs[has_layer].groupby(section_col).size()
    best_sections = section_counts.nlargest(n_sections * 2)
    print(f"\n  Top sections by annotation coverage:")
    for sec, count in best_sections.items():
        total = (merfish.obs[section_col] == sec).sum()
        print(f"    {sec}: {count:,} annotated / {total:,} total ({count/total*100:.0f}%)")

    # Select top n_sections
    selected = list(best_sections.index[:n_sections])
    print(f"\n  Selected: {selected}")

    all_results = []
    for sec_id in selected:
        print(f"\n  Processing section {sec_id}...")

        # Subset to this section, cells with layer annotation
        sec_mask = (merfish.obs[section_col] == sec_id).values
        layer_mask = has_layer.values
        mask = sec_mask & layer_mask
        adata_sec = merfish[mask].copy()

        # Subsample if too large
        if adata_sec.n_obs > max_cells_per_section:
            np.random.seed(42)
            idx = np.random.choice(adata_sec.n_obs, max_cells_per_section, replace=False)
            adata_sec = adata_sec[idx].copy()

        print(f"    {adata_sec.n_obs:,} cells")

        # Get ground truth layers
        gt_layers = adata_sec.obs[layer_col].values.astype(str)
        unique_gt = sorted(set(gt_layers))
        print(f"    Ground truth layers: {unique_gt}")

        # Set up spatial coordinates
        if "X_spatial_raw" in adata_sec.obsm:
            adata_sec.obsm["spatial"] = adata_sec.obsm["X_spatial_raw"]
        elif "spatial" not in adata_sec.obsm:
            print(f"    WARNING: No spatial coordinates found, skipping")
            continue

        # Preprocess for BANKSY
        print(f"    Running BANKSY (λ={BANKSY_LAMBDA}, res={BANKSY_RESOLUTION})...")
        t1 = time.time()
        adata_b = preprocess_for_banksy(adata_sec)
        banksy_labels = run_banksy(adata_b)
        n_cl = len(np.unique(banksy_labels))
        print(f"    BANKSY: {n_cl} clusters ({time.time()-t1:.0f}s)")

        # Map BANKSY clusters to layers by majority vote
        cluster_to_layer = {}
        for cl in np.unique(banksy_labels):
            cl_gt = gt_layers[banksy_labels == cl]
            from collections import Counter
            counts = Counter(cl_gt)
            cluster_to_layer[cl] = counts.most_common(1)[0][0]

        mapped_labels = np.array([cluster_to_layer[cl] for cl in banksy_labels])

        # Compute metrics
        # Filter to standard cortical layers for fair comparison
        cortical = {"L1", "L2/3", "L4", "L5", "L6"}
        cortical_mask = np.array([l in cortical for l in gt_layers])

        if cortical_mask.sum() > 100:
            ari = adjusted_rand_score(gt_layers[cortical_mask],
                                      banksy_labels[cortical_mask])
            nmi = normalized_mutual_info_score(gt_layers[cortical_mask],
                                               banksy_labels[cortical_mask])

            # Per-layer accuracy (majority-vote mapping)
            gt_cortical = gt_layers[cortical_mask]
            mapped_cortical = mapped_labels[cortical_mask]
            accuracy = (gt_cortical == mapped_cortical).mean()

            print(f"    ARI (cortical): {ari:.3f}")
            print(f"    NMI (cortical): {nmi:.3f}")
            print(f"    Accuracy (majority-vote mapping): {accuracy:.3f}")

            # Confusion matrix
            cm_labels = sorted(cortical)
            cm = confusion_matrix(gt_cortical, mapped_cortical, labels=cm_labels)
            print(f"\n    Confusion matrix (rows=truth, cols=predicted):")
            print(f"    {'':>8}", end="")
            for l in cm_labels:
                print(f" {l:>6}", end="")
            print()
            for i, l in enumerate(cm_labels):
                print(f"    {l:>8}", end="")
                for j in range(len(cm_labels)):
                    print(f" {cm[i,j]:>6}", end="")
                print()

            all_results.append({
                "section": sec_id,
                "n_cells": adata_sec.n_obs,
                "n_clusters": n_cl,
                "ari": ari,
                "nmi": nmi,
                "accuracy": accuracy,
            })
        else:
            print(f"    Not enough cortical cells for metrics")

    # Summary
    if all_results:
        df = pd.DataFrame(all_results)
        print(f"\n  MERFISH Benchmark Summary:")
        print(f"  {df.to_string(index=False)}")
        csv_path = os.path.join(OUT_DIR, "validation_merfish_summary.csv")
        df.to_csv(csv_path, index=False)
        print(f"  Saved: {csv_path}")

    return all_results


# ── Validation 3: Domain Classification Concordance ──────────────────

def validate_domain_concordance(sample_ids=None):
    """Compare BANKSY domain classification to current pipeline across samples.

    For cells classified as 'Cortical' by BOTH methods, compare layer
    distributions and cell type proportions.
    """
    if sample_ids is None:
        sample_ids = ["Br2039", "Br6389", "Br8433", "Br6437"]

    print(f"\n{'='*60}")
    print(f"Validation 3: Domain concordance ({len(sample_ids)} samples)")
    print(f"{'='*60}")

    all_concordance = []
    for sample_id in sample_ids:
        path = os.path.join(H5AD_DIR, f"{sample_id}_annotated.h5ad")
        adata = ad.read_h5ad(path)
        qc_mask = adata.obs["hybrid_qc_pass"].values.astype(bool)
        adata = adata[qc_mask].copy()

        # Get current pipeline domains
        current_domains = adata.obs["spatial_domain"].values.astype(str)

        # Run BANKSY
        adata_b = preprocess_for_banksy(adata)
        banksy_labels = run_banksy(adata_b)

        # Classify BANKSY domains
        subclass = adata.obs.get("corr_subclass",
                                  adata.obs.get("subclass_label")).values.astype(str)
        pred_depth = pd.to_numeric(adata.obs["predicted_norm_depth"],
                                    errors="coerce").values

        banksy_domains = np.full(len(banksy_labels), "Cortical", dtype=object)
        for cl in np.unique(banksy_labels):
            mask = banksy_labels == cl
            n_cl = mask.sum()
            cl_sub = subclass[mask]
            cl_depth = pred_depth[mask]

            sub_counts = {}
            for s in cl_sub:
                sub_counts[s] = sub_counts.get(s, 0) + 1

            vasc_frac = sum(sub_counts.get(v, 0) for v in VASCULAR_TYPES) / n_cl
            nn_frac = sum(sub_counts.get(v, 0) for v in NON_NEURONAL_TYPES) / n_cl
            oligo_frac = sub_counts.get("Oligodendrocyte", 0) / n_cl
            mean_d = float(np.nanmean(cl_depth[~np.isnan(cl_depth)])) if (~np.isnan(cl_depth)).sum() > 0 else np.nan

            if vasc_frac > VASCULAR_THRESH:
                domain = "Vascular"
            elif nn_frac > MENINGEAL_NN_THRESH and (np.isnan(mean_d) or mean_d < MENINGEAL_DEPTH_THRESH):
                domain = "Meningeal"
            elif oligo_frac > WM_OLIGO_THRESH and mean_d > WM_DEPTH_THRESH:
                domain = "WM"
            elif (1.0 - nn_frac) > 0.20 and not np.isnan(mean_d) and 0.0 <= mean_d <= 0.90:
                domain = "Cortical"
            elif mean_d > WM_DEPTH_THRESH:
                domain = "WM"
            else:
                domain = "Cortical"

            banksy_domains[mask] = domain

        # Concordance matrix
        # Map current pipeline: Extra-cortical → Meningeal for comparison
        current_mapped = current_domains.copy()
        current_mapped[current_mapped == "Extra-cortical"] = "Meningeal"

        domain_cats = ["Cortical", "Vascular", "Meningeal", "WM"]

        # Compute concordance
        both_cortical = (banksy_domains == "Cortical") & (current_mapped == "Cortical")
        agreement = (banksy_domains == current_mapped).sum() / len(banksy_domains)

        # Cross-tabulation
        ct = pd.crosstab(
            pd.Series(current_mapped, name="Current"),
            pd.Series(banksy_domains, name="BANKSY"),
        )

        print(f"\n  {sample_id} ({SAMPLE_TO_DX.get(sample_id, '?')}):")
        print(f"    Overall agreement: {agreement:.1%}")
        print(f"    Both cortical: {both_cortical.sum():,} / {len(banksy_domains):,} ({both_cortical.mean():.1%})")
        print(f"    Cross-tabulation:")
        print(f"    {ct.to_string()}")

        all_concordance.append({
            "sample_id": sample_id,
            "diagnosis": SAMPLE_TO_DX.get(sample_id, "?"),
            "n_cells": len(banksy_domains),
            "agreement": agreement,
            "both_cortical_n": both_cortical.sum(),
            "both_cortical_pct": both_cortical.mean() * 100,
            "banksy_cortical_pct": (banksy_domains == "Cortical").mean() * 100,
            "current_cortical_pct": (current_mapped == "Cortical").mean() * 100,
            "banksy_vascular_pct": (banksy_domains == "Vascular").mean() * 100,
            "current_vascular_pct": (current_mapped == "Vascular").mean() * 100,
        })

    # Summary
    df = pd.DataFrame(all_concordance)
    print(f"\n  Domain Concordance Summary:")
    print(f"  {df[['sample_id', 'agreement', 'banksy_cortical_pct', 'current_cortical_pct', 'banksy_vascular_pct', 'current_vascular_pct']].to_string(index=False)}")
    csv_path = os.path.join(OUT_DIR, "validation_domain_concordance.csv")
    df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")

    return all_concordance


# ── Main ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-merfish", action="store_true",
                        help="Skip MERFISH benchmark (slow)")
    parser.add_argument("--sample", default="Br6437",
                        help="Sample for marker validation")
    args = parser.parse_args()

    t_start = time.time()

    # Validation 1: Laminar markers
    f_results, banksy_labels, cl_order, cl_depth = validate_laminar_markers(args.sample)

    # Validation 2: MERFISH benchmark
    if not args.skip_merfish:
        merfish_results = validate_merfish(n_sections=3, max_cells_per_section=25000)
    else:
        print("\n  Skipping MERFISH benchmark (--skip-merfish)")

    # Validation 3: Domain concordance
    concordance = validate_domain_concordance()

    total = time.time() - t_start
    print(f"\n{'='*60}")
    print(f"Total validation time: {total:.0f}s ({total/60:.1f} min)")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
