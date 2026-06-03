#!/usr/bin/env python3
"""
BANKSY Phase 1: Single-sample pilot with parameter sweep.

Runs BANKSY on Br2039 (smallest sample) across a grid of lambda and resolution
values. Generates diagnostic figures comparing BANKSY clusters to current
depth-bin layer assignments and cell type composition.

Output: output/banksy/pilot_*.png and output/banksy/pilot_parameter_sweep.csv

Usage:
    python3 -u code/analysis/banksy_01_pilot.py [--sample SAMPLE_ID]
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
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

# BANKSY imports
from banksy.initialize_banksy import initialize_banksy
from banksy.embed_banksy import generate_banksy_matrix
from banksy.cluster_methods import run_Leiden_partition
from banksy_utils.umap_pca import pca_umap

# Project imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    H5AD_DIR, LAYER_COLORS, LAYER_ORDER, CORTICAL_LAYERS,
    SUBCLASS_TO_CLASS, CLASS_COLORS, EXCLUDE_SAMPLES,
)

# ── Constants ─────────────────────────────────────────────────────────

OUT_DIR = os.path.expanduser("~/Github/SCZ_Xenium/output/banksy")
os.makedirs(OUT_DIR, exist_ok=True)

DEFAULT_SAMPLE = "Br2039"  # smallest sample, ~44K cells

# BANKSY parameter grid
LAMBDA_LIST = [0.0, 0.2, 0.5, 0.8]
RESOLUTIONS = [0.3, 0.5, 0.7, 1.0]
K_GEOM = 15
PCA_DIMS = [20]

# Layer depth boundaries for reference lines on figures
LAYER_BOUNDARIES = [0.10, 0.30, 0.45, 0.65, 0.85]
LAYER_TICKS = [0.05, 0.20, 0.375, 0.55, 0.75, 0.925]
LAYER_NAMES = ["L1", "L2/3", "L4", "L5", "L6", "WM"]


# ── Data Loading ──────────────────────────────────────────────────────

def load_sample(sample_id):
    """Load a single sample, subset to QC-pass cells."""
    path = os.path.join(H5AD_DIR, f"{sample_id}_annotated.h5ad")
    print(f"Loading {path}...")
    adata = ad.read_h5ad(path)

    # QC filter
    if "hybrid_qc_pass" in adata.obs.columns:
        qc_mask = adata.obs["hybrid_qc_pass"].values.astype(bool)
    elif "qc_pass" in adata.obs.columns:
        qc_mask = adata.obs["qc_pass"].values.astype(bool)
    else:
        qc_mask = np.ones(adata.n_obs, dtype=bool)

    adata = adata[qc_mask].copy()
    print(f"  {adata.n_obs:,} QC-pass cells, {adata.n_vars} genes")
    return adata


def preprocess_for_banksy(adata):
    """Normalize, log-transform, and z-score for BANKSY.

    Returns a copy with processed .X (original adata unchanged).
    """
    print("Preprocessing for BANKSY...")
    adata_b = adata.copy()

    # Store raw counts
    adata_b.layers["counts"] = adata_b.X.copy()

    # Standard scanpy preprocessing
    sc.pp.normalize_total(adata_b, target_sum=1e4)
    sc.pp.log1p(adata_b)
    sc.pp.scale(adata_b)  # z-score per gene — required for BANKSY

    # Ensure spatial coords are present
    assert "spatial" in adata_b.obsm, "Missing obsm['spatial']"
    print(f"  Spatial range: X=[{adata_b.obsm['spatial'][:,0].min():.0f}, "
          f"{adata_b.obsm['spatial'][:,0].max():.0f}] "
          f"Y=[{adata_b.obsm['spatial'][:,1].min():.0f}, "
          f"{adata_b.obsm['spatial'][:,1].max():.0f}] μm")
    return adata_b


# ── BANKSY Pipeline ───────────────────────────────────────────────────

def run_banksy_sweep(adata_b, lambda_list, resolutions, k_geom=15, pca_dims=[20]):
    """Run BANKSY initialization + parameter sweep. Returns results_df."""

    coord_keys = ("x", "y", "spatial")

    # Step 1: Initialize spatial graph
    print(f"\n{'='*60}")
    print(f"Initializing BANKSY (k_geom={k_geom})...")
    t0 = time.time()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        banksy_dict = initialize_banksy(
            adata_b,
            coord_keys,
            k_geom,
            nbr_weight_decay="scaled_gaussian",
            max_m=1,
            plt_edge_hist=False,
            plt_nbr_weights=False,
            plt_agf_angles=False,
        )
    print(f"  Initialization: {time.time()-t0:.1f}s")

    # Step 2: Generate BANKSY augmented matrices for all lambda values
    print(f"Generating BANKSY matrices for lambdas={lambda_list}...")
    t1 = time.time()
    banksy_dict, banksy_matrix = generate_banksy_matrix(
        adata_b, banksy_dict, lambda_list, max_m=1
    )
    print(f"  Matrix generation: {time.time()-t1:.1f}s")

    # Step 3: PCA + UMAP for each lambda
    print(f"Running PCA (dims={pca_dims})...")
    t2 = time.time()
    pca_umap(banksy_dict, pca_dims=pca_dims, add_umap=True)
    print(f"  PCA+UMAP: {time.time()-t2:.1f}s")

    # Step 4: Leiden clustering across resolutions
    print(f"Running Leiden clustering (resolutions={resolutions})...")
    t3 = time.time()
    results_df, max_num_labels = run_Leiden_partition(
        banksy_dict,
        resolutions,
        num_nn=50,
        num_iterations=-1,
        partition_seed=42,
        match_labels=False,
    )
    print(f"  Leiden clustering: {time.time()-t3:.1f}s")
    print(f"  Total BANKSY time: {time.time()-t0:.1f}s")
    print(f"  Results: {len(results_df)} parameter combinations")

    return banksy_dict, results_df


# ── Evaluation Metrics ────────────────────────────────────────────────

def evaluate_results(results_df, adata_orig):
    """Compute metrics for each parameter combination."""
    # Get reference labels
    depth_layers = adata_orig.obs["layer"].values.astype(str)
    pred_depth = pd.to_numeric(adata_orig.obs["predicted_norm_depth"],
                                errors="coerce").values
    subclass = adata_orig.obs.get("corr_subclass",
                                   adata_orig.obs.get("subclass_label"))

    rows = []
    for idx, row in results_df.iterrows():
        lam = row["lambda_param"]
        res = row["resolution"]
        labels = row["labels"]
        # Handle different label formats
        if hasattr(labels, "dense"):
            labels = labels.dense
        labels = np.asarray(labels).astype(str)
        n_clusters = len(np.unique(labels))

        # ARI and NMI against depth-bin layers
        ari = adjusted_rand_score(depth_layers, labels)
        nmi = normalized_mutual_info_score(depth_layers, labels)

        # Laminar ordering quality: order clusters by median depth,
        # check monotonicity
        cluster_ids = np.unique(labels)
        cluster_median_depth = {}
        for cl in cluster_ids:
            mask = labels == cl
            depths = pred_depth[mask]
            valid = depths[~np.isnan(depths)]
            cluster_median_depth[cl] = np.median(valid) if len(valid) > 0 else 0.5

        ordered_clusters = sorted(cluster_ids, key=lambda c: cluster_median_depth[c])
        ordered_depths = [cluster_median_depth[c] for c in ordered_clusters]
        # Monotonicity: fraction of consecutive pairs in correct order
        n_pairs = len(ordered_depths) - 1
        if n_pairs > 0:
            monotonic = sum(1 for i in range(n_pairs)
                           if ordered_depths[i+1] >= ordered_depths[i]) / n_pairs
        else:
            monotonic = 1.0

        # Vascular detection: find clusters enriched for Endo+VLMC
        vasc_types = {"Endothelial", "VLMC"}
        n_vasc_total = 0
        for cl in cluster_ids:
            mask = labels == cl
            cl_types = subclass[mask].values.astype(str)
            vasc_frac = sum(1 for t in cl_types if t in vasc_types) / len(cl_types)
            if vasc_frac > 0.5:
                n_vasc_total += mask.sum()

        rows.append({
            "lambda": lam,
            "resolution": res,
            "n_clusters": n_clusters,
            "ari_vs_depth_layers": ari,
            "nmi_vs_depth_layers": nmi,
            "laminar_monotonicity": monotonic,
            "n_vascular_cells": n_vasc_total,
            "pct_vascular": 100 * n_vasc_total / len(labels),
        })

    metrics_df = pd.DataFrame(rows)
    return metrics_df


# ── Plotting ──────────────────────────────────────────────────────────

def plot_parameter_heatmap(metrics_df):
    """Heatmap of ARI and n_clusters across lambda x resolution."""
    fig, axes = plt.subplots(1, 3, figsize=(20, 5.5))

    for ax, metric, title, cmap in [
        (axes[0], "ari_vs_depth_layers", "ARI vs Depth-Bin Layers", "YlOrRd"),
        (axes[1], "n_clusters", "Number of Clusters", "viridis"),
        (axes[2], "laminar_monotonicity", "Laminar Ordering Quality", "YlGn"),
    ]:
        pivot = metrics_df.pivot(index="lambda", columns="resolution",
                                  values=metric)
        im = ax.imshow(pivot.values, cmap=cmap, aspect="auto",
                        vmin=pivot.values.min() * 0.9,
                        vmax=pivot.values.max() * 1.05)
        ax.set_xticks(range(len(pivot.columns)))
        ax.set_xticklabels([f"{v:.1f}" for v in pivot.columns], fontsize=14)
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels([f"{v:.1f}" for v in pivot.index], fontsize=14)
        ax.set_xlabel("Resolution", fontsize=16)
        ax.set_ylabel("Lambda", fontsize=16)
        ax.set_title(title, fontsize=18, fontweight="bold")
        # Annotate cells
        for i in range(len(pivot.index)):
            for j in range(len(pivot.columns)):
                val = pivot.values[i, j]
                fmt = f"{val:.2f}" if val < 10 else f"{val:.0f}"
                ax.text(j, i, fmt, ha="center", va="center",
                        fontsize=13, fontweight="bold",
                        color="white" if val < pivot.values.mean() else "black")
        plt.colorbar(im, ax=ax, shrink=0.8)

    fig.suptitle("BANKSY Parameter Sweep — Br2039", fontsize=22, fontweight="bold")
    fig.tight_layout()
    return fig


def plot_spatial_clusters(adata_orig, results_df, metrics_df, n_best=4):
    """Spatial plots for the top-N parameter combos by ARI."""
    top = metrics_df.nlargest(n_best, "ari_vs_depth_layers")
    coords = adata_orig.obsm["spatial"]
    x, y = coords[:, 0], coords[:, 1]

    fig, axes = plt.subplots(1, n_best + 1, figsize=(6 * (n_best + 1), 6))

    # Reference: current depth-bin layers
    ax = axes[0]
    layers = adata_orig.obs["layer"].values.astype(str)
    for layer_name, color in LAYER_COLORS.items():
        mask = layers == layer_name
        if mask.sum() > 0:
            ax.scatter(x[mask], y[mask], c=[color], s=0.3, alpha=0.5,
                       rasterized=True)
    ax.set_title("Current Depth-Bin\nLayers", fontsize=16, fontweight="bold")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_facecolor("#0a0a0a")
    ax.set_xticks([])
    ax.set_yticks([])

    # BANKSY results
    for plot_idx, (_, row_metrics) in enumerate(top.iterrows()):
        ax = axes[plot_idx + 1]
        lam = row_metrics["lambda"]
        res = row_metrics["resolution"]
        ari = row_metrics["ari_vs_depth_layers"]
        n_cl = int(row_metrics["n_clusters"])

        # Find matching result
        match = results_df[
            (results_df["lambda_param"] == lam) &
            (results_df["resolution"] == res)
        ]
        if len(match) == 0:
            continue
        labels = match.iloc[0]["labels"]
        if hasattr(labels, "dense"):
            labels = labels.dense
        labels = np.asarray(labels).astype(int)

        # Color by cluster, ordered by median depth
        unique_cl = np.unique(labels)
        pred_depth = pd.to_numeric(adata_orig.obs["predicted_norm_depth"],
                                    errors="coerce").values
        cl_depth = {cl: np.nanmedian(pred_depth[labels == cl]) for cl in unique_cl}
        ordered = sorted(unique_cl, key=lambda c: cl_depth[c])

        cmap = matplotlib.colormaps.get_cmap("Spectral_r").resampled(len(ordered))
        cl_colors = {cl: cmap(i / max(1, len(ordered) - 1))
                     for i, cl in enumerate(ordered)}

        for cl in ordered:
            mask = labels == cl
            ax.scatter(x[mask], y[mask], c=[cl_colors[cl]], s=0.3, alpha=0.5,
                       rasterized=True)

        ax.set_title(f"λ={lam}, res={res}\n{n_cl} clusters, ARI={ari:.3f}",
                     fontsize=14, fontweight="bold")
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_facecolor("#0a0a0a")
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle("BANKSY Spatial Clusters — Top 4 by ARI", fontsize=20,
                 fontweight="bold")
    fig.tight_layout()
    return fig


def plot_depth_violins(adata_orig, results_df, metrics_df):
    """Violin plot of predicted depth per BANKSY cluster for best combo."""
    best = metrics_df.loc[metrics_df["ari_vs_depth_layers"].idxmax()]
    lam, res = best["lambda"], best["resolution"]

    match = results_df[
        (results_df["lambda_param"] == lam) &
        (results_df["resolution"] == res)
    ]
    labels = match.iloc[0]["labels"]
    if hasattr(labels, "dense"):
        labels = labels.dense
    labels = np.asarray(labels).astype(int)

    pred_depth = pd.to_numeric(adata_orig.obs["predicted_norm_depth"],
                                errors="coerce").values

    # Order clusters by median depth
    unique_cl = np.unique(labels)
    cl_depth = {cl: np.nanmedian(pred_depth[labels == cl]) for cl in unique_cl}
    ordered = sorted(unique_cl, key=lambda c: cl_depth[c])

    fig, ax = plt.subplots(figsize=(max(14, len(ordered) * 0.9), 6))

    data_list = []
    for cl in ordered:
        depths = pred_depth[labels == cl]
        data_list.append(depths[~np.isnan(depths)])

    vp = ax.violinplot(data_list, positions=range(len(ordered)),
                        showmedians=True, showextrema=False)
    cmap = matplotlib.colormaps.get_cmap("Spectral_r").resampled(len(ordered))
    for i, body in enumerate(vp["bodies"]):
        body.set_facecolor(cmap(i / max(1, len(ordered) - 1)))
        body.set_alpha(0.7)
    vp["cmedians"].set_color("white")
    vp["cmedians"].set_linewidth(2)

    # Layer boundary lines
    for boundary in LAYER_BOUNDARIES:
        ax.axhline(boundary, color="#999999", linewidth=0.7, alpha=0.5,
                    linestyle="--")

    ax.set_xticks(range(len(ordered)))
    ax.set_xticklabels([f"C{cl}\n(n={len(data_list[i]):,})"
                         for i, cl in enumerate(ordered)],
                        fontsize=11)
    ax.set_ylim(1.05, -0.05)  # inverted: pia at top
    ax.set_yticks(LAYER_TICKS)
    ax.set_yticklabels(LAYER_NAMES, fontsize=14)
    ax.set_ylabel("Cortical Depth", fontsize=16, fontweight="bold")
    ax.set_xlabel("BANKSY Cluster (ordered by median depth)", fontsize=14)
    ax.set_title(f"BANKSY Clusters vs Cortical Depth — λ={lam}, res={res}",
                 fontsize=20, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return fig


def plot_composition_bars(adata_orig, results_df, metrics_df):
    """Stacked bar of subclass composition per BANKSY cluster for best combo."""
    best = metrics_df.loc[metrics_df["ari_vs_depth_layers"].idxmax()]
    lam, res = best["lambda"], best["resolution"]

    match = results_df[
        (results_df["lambda_param"] == lam) &
        (results_df["resolution"] == res)
    ]
    labels = match.iloc[0]["labels"]
    if hasattr(labels, "dense"):
        labels = labels.dense
    labels = np.asarray(labels).astype(int)

    subclass = adata_orig.obs.get("corr_subclass",
                                   adata_orig.obs.get("subclass_label"))
    subclass = subclass.values.astype(str)
    pred_depth = pd.to_numeric(adata_orig.obs["predicted_norm_depth"],
                                errors="coerce").values

    # Order clusters by median depth
    unique_cl = np.unique(labels)
    cl_depth = {cl: np.nanmedian(pred_depth[labels == cl]) for cl in unique_cl}
    ordered = sorted(unique_cl, key=lambda c: cl_depth[c])

    # Get all subclasses present
    all_subs = sorted(set(subclass))

    # Build composition matrix
    comp = np.zeros((len(ordered), len(all_subs)))
    for i, cl in enumerate(ordered):
        mask = labels == cl
        cl_types = subclass[mask]
        for j, sub in enumerate(all_subs):
            comp[i, j] = (cl_types == sub).sum()

    # Normalize to fractions
    row_sums = comp.sum(axis=1, keepdims=True)
    comp_frac = np.where(row_sums > 0, comp / row_sums, 0)

    # Color by class
    sub_colors = []
    for sub in all_subs:
        cls = SUBCLASS_TO_CLASS.get(sub, "Unknown")
        sub_colors.append(CLASS_COLORS.get(cls, "#888888"))

    fig, ax = plt.subplots(figsize=(max(14, len(ordered) * 0.9), 7))
    bottom = np.zeros(len(ordered))
    bar_width = 0.8

    for j, sub in enumerate(all_subs):
        ax.bar(range(len(ordered)), comp_frac[:, j], bottom=bottom,
               width=bar_width, color=sub_colors[j], label=sub,
               edgecolor="white", linewidth=0.3)
        bottom += comp_frac[:, j]

    ax.set_xticks(range(len(ordered)))
    ax.set_xticklabels([f"C{cl}\n(n={int(row_sums[i,0]):,})"
                         for i, cl in enumerate(ordered)],
                        fontsize=11)
    ax.set_ylabel("Fraction", fontsize=16)
    ax.set_xlabel("BANKSY Cluster (ordered by median depth)", fontsize=14)
    ax.set_title(f"Cell Type Composition per BANKSY Cluster — λ={lam}, res={res}",
                 fontsize=20, fontweight="bold")

    # Legend: show only top-12 most common subclasses
    handles, lbls = ax.get_legend_handles_labels()
    total_per_sub = comp.sum(axis=0)
    top_idx = np.argsort(total_per_sub)[-12:][::-1]
    ax.legend([handles[i] for i in top_idx],
              [lbls[i] for i in top_idx],
              loc="upper right", fontsize=10, ncol=2,
              framealpha=0.9)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return fig


# ── Main ──────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", default=DEFAULT_SAMPLE,
                        help=f"Sample ID (default: {DEFAULT_SAMPLE})")
    args = parser.parse_args()

    sample_id = args.sample
    t_start = time.time()

    # 1. Load and preprocess
    print(f"\n{'='*60}")
    print(f"BANKSY Phase 1: Pilot on {sample_id}")
    print(f"{'='*60}")

    adata_orig = load_sample(sample_id)
    adata_b = preprocess_for_banksy(adata_orig)

    # 2. Run BANKSY sweep
    banksy_dict, results_df = run_banksy_sweep(
        adata_b, LAMBDA_LIST, RESOLUTIONS,
        k_geom=K_GEOM, pca_dims=PCA_DIMS
    )

    # 3. Evaluate
    print(f"\n{'='*60}")
    print("Evaluating parameter combinations...")
    metrics_df = evaluate_results(results_df, adata_orig)
    print(metrics_df.to_string(index=False))

    # Save metrics
    csv_path = os.path.join(OUT_DIR, f"pilot_parameter_sweep_{sample_id}.csv")
    metrics_df.to_csv(csv_path, index=False)
    print(f"\nSaved: {csv_path}")

    # 4. Generate diagnostic figures
    print(f"\n{'='*60}")
    print("Generating diagnostic figures...")

    # 4A. Parameter heatmap
    fig = plot_parameter_heatmap(metrics_df)
    out = os.path.join(OUT_DIR, f"pilot_parameter_heatmap_{sample_id}.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")

    # 4B. Spatial cluster maps
    fig = plot_spatial_clusters(adata_orig, results_df, metrics_df)
    out = os.path.join(OUT_DIR, f"pilot_spatial_clusters_{sample_id}.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")

    # 4C. Depth violins
    fig = plot_depth_violins(adata_orig, results_df, metrics_df)
    out = os.path.join(OUT_DIR, f"pilot_depth_violins_{sample_id}.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")

    # 4D. Composition bars
    fig = plot_composition_bars(adata_orig, results_df, metrics_df)
    out = os.path.join(OUT_DIR, f"pilot_composition_bars_{sample_id}.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")

    # 5. Summary
    total_time = time.time() - t_start
    best = metrics_df.loc[metrics_df["ari_vs_depth_layers"].idxmax()]
    print(f"\n{'='*60}")
    print(f"SUMMARY")
    print(f"{'='*60}")
    print(f"  Sample: {sample_id} ({adata_orig.n_obs:,} cells)")
    print(f"  Total time: {total_time:.0f}s")
    print(f"  Best combo: λ={best['lambda']}, res={best['resolution']}")
    print(f"    ARI vs depth-bin: {best['ari_vs_depth_layers']:.3f}")
    print(f"    NMI vs depth-bin: {best['nmi_vs_depth_layers']:.3f}")
    print(f"    N clusters: {int(best['n_clusters'])}")
    print(f"    Laminar monotonicity: {best['laminar_monotonicity']:.3f}")
    print(f"    Vascular cells: {int(best['n_vascular_cells'])} "
          f"({best['pct_vascular']:.1f}%)")
    print(f"\n  Outputs in: {OUT_DIR}")


if __name__ == "__main__":
    main()
