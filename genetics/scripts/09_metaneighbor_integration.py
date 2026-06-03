#!/usr/bin/env python
"""
09_metaneighbor_integration.py — MetaNeighbor cross-dataset cell-type matching.

Matches SEA-AD supertypes (137) to Siletti clusters (461) using:
  1. Combined mean expression matrix (gene symbols)
  2. 5000 HVGs selected by coefficient of variation
  3. Spearman correlation → neighbor-voting AUROC
  4. Reciprocal best hit identification

Inputs:
  - results/intermediates/seaad_supertype_mean_expression.csv
  - results/intermediates/siletti_cluster_level_stats.npz
  - data/adult_human_20221007.loom (gene ID mapping)

Outputs:
  - results/intermediates/full_seaad_siletti_461_mean_expression.csv
  - results/intermediates/metaneighbor_full461_auroc_matrix.csv
  - results/tables/metaneighbor_full461_reciprocal_best_hits.csv
  - results/tables/metaneighbor_full461_all_matches.csv
  - results/intermediates/rbh_seaad_to_siletti_metaneighbor.csv
  - results/intermediates/rbh_siletti_to_seaad_metaneighbor.csv
  - results/figures/metaneighbor_heatmap.png
  - results/figures/metaneighbor_reciprocal_hits_summary.png
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram

from scz_celltype_enrichment.config import (
    TABLES_DIR, FIGURES_DIR, INTERMEDIATES_DIR,
    FIGURE_DPI, TITLE_FONTSIZE, LABEL_FONTSIZE, TICK_FONTSIZE, LEGEND_FONTSIZE,
    PROJECT_ROOT,
)
from scz_celltype_enrichment.siletti.specificity import build_combined_mean_expression
from scz_celltype_enrichment.siletti.metaneighbor import (
    select_hvgs, compute_auroc_matrix, find_reciprocal_best_hits,
)
from scz_celltype_enrichment.utils import ensure_dirs, Timer


LOOM_PATH = PROJECT_ROOT / "data" / "adult_human_20221007.loom"
CLUSTER_STATS = INTERMEDIATES_DIR / "siletti_cluster_level_stats.npz"
SEAAD_MEAN_EXPR = INTERMEDIATES_DIR / "seaad_supertype_mean_expression.csv"


def build_match_lookups(auroc_matrix, output_dir):
    """Save top-N match lookups for the webapp."""
    # SEA-AD → Siletti (top 3)
    rows = []
    for seaad_type in auroc_matrix.index:
        top3 = auroc_matrix.loc[seaad_type].nlargest(3)
        rows.append({
            "seaad_supertype": seaad_type,
            "best_siletti_match": top3.index[0],
            "best_siletti_auroc": f"{top3.iloc[0]:.4f}",
            "second_siletti_match": top3.index[1],
            "second_siletti_auroc": f"{top3.iloc[1]:.4f}",
            "third_siletti_match": top3.index[2],
            "third_siletti_auroc": f"{top3.iloc[2]:.4f}",
        })
    pd.DataFrame(rows).to_csv(str(output_dir / "rbh_seaad_to_siletti_metaneighbor.csv"), index=False)

    # Siletti → SEA-AD (top 2)
    rows = []
    for sil_cluster in auroc_matrix.columns:
        top2 = auroc_matrix[sil_cluster].nlargest(2)
        rows.append({
            "siletti_cluster": sil_cluster,
            "best_seaad_match": top2.index[0],
            "best_seaad_auroc": f"{top2.iloc[0]:.4f}",
            "second_seaad_match": top2.index[1],
            "second_seaad_auroc": f"{top2.iloc[1]:.4f}",
        })
    pd.DataFrame(rows).to_csv(str(output_dir / "rbh_siletti_to_seaad_metaneighbor.csv"), index=False)
    print("  Saved match lookups for webapp")


def plot_heatmap(auroc_matrix, rbh, output_path):
    """Clustered heatmap of AUROC for reciprocal best hit Siletti clusters."""
    if len(rbh) == 0:
        return

    siletti_show = sorted(rbh["siletti_cluster"].unique())
    sub = auroc_matrix[siletti_show]

    if sub.shape[0] > 2 and sub.shape[1] > 2:
        row_link = linkage(sub.values, method="ward")
        col_link = linkage(sub.values.T, method="ward")
        row_ord = dendrogram(row_link, no_plot=True)["leaves"]
        col_ord = dendrogram(col_link, no_plot=True)["leaves"]
        sub = sub.iloc[row_ord, col_ord]

    fig, ax = plt.subplots(figsize=(max(20, len(siletti_show) * 0.15),
                                     max(12, len(sub) * 0.12)))
    im = ax.imshow(sub.values, aspect="auto", cmap="RdYlBu_r", vmin=0.3, vmax=1.0)
    ax.set_xticks(range(len(sub.columns)))
    ax.set_xticklabels(sub.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(len(sub.index)))
    ax.set_yticklabels(sub.index, fontsize=7)
    ax.set_title("MetaNeighbor AUROC: SEA-AD \u2192 Siletti (5000 HVGs)",
                 fontsize=TITLE_FONTSIZE)
    plt.colorbar(im, ax=ax, label="AUROC", shrink=0.6)
    plt.tight_layout()
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def plot_summary(rbh, all_matches, output_path):
    """Summary plot of reciprocal best hit quality."""
    fig, axes = plt.subplots(1, 2, figsize=(18, 8))

    ax = axes[0]
    ax.hist(rbh["mean_auroc"], bins=30, color="tab:blue", alpha=0.7, edgecolor="black")
    ax.axvline(0.9, color="tab:red", linestyle="--", linewidth=2, label="0.9 threshold")
    ax.set_xlabel("Mean AUROC", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("Count", fontsize=LABEL_FONTSIZE)
    n_high = (rbh["mean_auroc"] > 0.9).sum()
    ax.set_title(f"Reciprocal best hits (n={len(rbh)}, {n_high} with AUROC > 0.9)",
                 fontsize=TITLE_FONTSIZE)
    ax.tick_params(labelsize=TICK_FONTSIZE)
    ax.legend(fontsize=LEGEND_FONTSIZE)

    ax = axes[1]
    recip = all_matches[all_matches["is_reciprocal"]]["best_auroc"]
    non_recip = all_matches[~all_matches["is_reciprocal"]]["best_auroc"]
    ax.hist(recip, bins=30, alpha=0.6, color="tab:blue",
            label=f"Reciprocal (n={len(recip)})", edgecolor="black")
    ax.hist(non_recip, bins=30, alpha=0.6, color="tab:orange",
            label=f"Non-reciprocal (n={len(non_recip)})", edgecolor="black")
    ax.set_xlabel("Best match AUROC", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("Count", fontsize=LABEL_FONTSIZE)
    ax.set_title("Reciprocal vs non-reciprocal", fontsize=TITLE_FONTSIZE)
    ax.tick_params(labelsize=TICK_FONTSIZE)
    ax.legend(fontsize=LEGEND_FONTSIZE)

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def main():
    ensure_dirs(TABLES_DIR, FIGURES_DIR, INTERMEDIATES_DIR)

    print("=" * 60)
    print("Step 09: MetaNeighbor Integration (Full 461 Siletti)")
    print("=" * 60)
    t0 = time.time()

    # --- Build or load combined mean expression ---
    me_path = INTERMEDIATES_DIR / "full_seaad_siletti_461_mean_expression.csv"
    if me_path.exists():
        print("\n1. Loading cached combined mean expression...")
        with Timer("Loading"):
            mean_expr = pd.read_csv(str(me_path), index_col=0)
    else:
        print("\n1. Building combined mean expression...")
        mean_expr = build_combined_mean_expression(SEAAD_MEAN_EXPR, CLUSTER_STATS, LOOM_PATH)
        mean_expr.to_csv(str(me_path))
        print(f"  Saved: {me_path}")

    seaad_cols = [c for c in mean_expr.columns if not c.startswith("Siletti_")]
    siletti_cols = [c for c in mean_expr.columns if c.startswith("Siletti_")]
    print(f"  {len(seaad_cols)} SEA-AD + {len(siletti_cols)} Siletti, {len(mean_expr)} genes")

    mean_seaad = mean_expr[seaad_cols]
    mean_siletti = mean_expr[siletti_cols]
    mean_siletti.columns = [c.replace("Siletti_", "") for c in mean_siletti.columns]

    # --- Select HVGs ---
    print("\n2. Selecting HVGs...")
    hvgs = select_hvgs(mean_seaad, mean_siletti, n_hvgs=5000)

    # --- Compute or load AUROC matrix ---
    auroc_path = INTERMEDIATES_DIR / "metaneighbor_full461_auroc_matrix.csv"
    if auroc_path.exists():
        print(f"\n3. Loading cached AUROC matrix...")
        auroc_matrix = pd.read_csv(str(auroc_path), index_col=0)
    else:
        print("\n3. Computing AUROC matrix...")
        auroc_matrix, _ = compute_auroc_matrix(mean_seaad, mean_siletti, hvgs)
        auroc_matrix.to_csv(str(auroc_path))
        print(f"  Saved: {auroc_path}")

    # --- Reciprocal best hits ---
    print("\n4. Finding reciprocal best hits...")
    rbh, all_matches = find_reciprocal_best_hits(auroc_matrix)

    rbh.to_csv(str(TABLES_DIR / "metaneighbor_full461_reciprocal_best_hits.csv"), index=False)
    all_matches.to_csv(str(TABLES_DIR / "metaneighbor_full461_all_matches.csv"), index=False)

    # --- Match lookups for webapp ---
    print("\n5. Building match lookups...")
    build_match_lookups(auroc_matrix, INTERMEDIATES_DIR)

    # --- Figures ---
    print("\n6. Generating figures...")
    plot_heatmap(auroc_matrix, rbh, FIGURES_DIR / "metaneighbor_heatmap.png")
    plot_summary(rbh, all_matches, FIGURES_DIR / "metaneighbor_reciprocal_hits_summary.png")

    print(f"\nStep 09 complete. Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
