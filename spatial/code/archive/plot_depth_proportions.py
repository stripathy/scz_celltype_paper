#!/usr/bin/env python3
"""
Visualize depth-stratified cell type proportion and density analyses.

Generates three sets of figures:

1. Continuous depth profiles: Smooth curves for each subclass showing
   proportion/density as a function of cortical depth, split by diagnosis.

2. Layer-level heatmap: log2FC (SCZ/Control) for each subclass × layer,
   with significance stars.

3. Summary volcano: Effect size vs significance across all tests.

Requires:
  output/depth_proportions/gam_results.csv
  output/depth_proportions/gam_smooth_predictions.csv
  output/depth_proportions/layer_model_results.csv
  output/depth_proportions/cell_level_data.csv
  output/depth_proportions/layer_counts.csv

Usage:
    python3 -u plot_depth_proportions.py
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.ticker as mticker

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import BASE_DIR

OUTPUT_DIR = os.path.join(BASE_DIR, "output", "depth_proportions")
CORTICAL_LAYERS = ["L1", "L2/3", "L4", "L5", "L6"]

# Consistent colors
CTRL_COLOR = "#4477AA"
SCZ_COLOR = "#CC3311"

# Layer boundaries for shading
LAYER_BOUNDS = {
    "L1": (0, 0.1225), "L2/3": (0.1225, 0.4696), "L4": (0.4696, 0.5443),
    "L5": (0.5443, 0.7079), "L6": (0.7079, 0.9275),
}


def plot_depth_profiles(cells_df, gam_preds, gam_results, outcome="proportion",
                        n_bins=20):
    """Plot continuous depth profiles for each subclass.

    Shows raw binned data (dots + CI) overlaid with GAM-predicted smooth curves.
    One panel per subclass, arranged in a grid.
    """
    # Get subclasses that have predictions
    subclasses = sorted(gam_preds[gam_preds["outcome"] == outcome]["subclass"].unique())
    if len(subclasses) == 0:
        print(f"  No predictions for {outcome}, skipping depth profiles")
        return

    n_sc = len(subclasses)
    ncols = 6
    nrows = int(np.ceil(n_sc / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.5 * nrows),
                              squeeze=False)

    # Compute raw binned means for dots
    cells = cells_df[cells_df["predicted_norm_depth"].between(0, 1)].copy()
    bin_edges = np.linspace(0, 1, n_bins + 1)
    bin_mids = (bin_edges[:-1] + bin_edges[1:]) / 2
    cells["depth_bin"] = pd.cut(cells["predicted_norm_depth"], bins=bin_edges,
                                 labels=bin_mids, include_lowest=True).astype(float)

    for idx, sc in enumerate(subclasses):
        row, col = divmod(idx, ncols)
        ax = axes[row, col]

        # Add layer shading
        for layer, (lo, hi) in LAYER_BOUNDS.items():
            ax.axvspan(lo, hi, alpha=0.06, color="gray")
            ax.text((lo + hi) / 2, 0.98, layer, transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontsize=7, color="gray", alpha=0.7)

        # Raw binned data per diagnosis
        for dx, color in [("Control", CTRL_COLOR), ("SCZ", SCZ_COLOR)]:
            dx_cells = cells[cells["diagnosis"] == dx]

            if outcome == "proportion":
                # Proportion of this subclass per bin per sample
                is_sc = (dx_cells["subclass_label"] == sc).astype(int)
                dx_cells = dx_cells.copy()
                dx_cells["is_type"] = is_sc.values

                bin_props = (dx_cells.groupby(["sample_id", "depth_bin"])
                            .agg(n_type=("is_type", "sum"), n_total=("is_type", "count"))
                            .reset_index())
                bin_props["prop"] = bin_props["n_type"] / bin_props["n_total"]

                # Mean ± SEM across samples
                summary = (bin_props.groupby("depth_bin")["prop"]
                          .agg(["mean", "sem", "count"]).reset_index())
                summary = summary[summary["count"] >= 5]

                ax.errorbar(summary["depth_bin"], summary["mean"],
                           yerr=summary["sem"], fmt="o", color=color,
                           markersize=3, alpha=0.6, capsize=2, linewidth=0.8)
            else:
                # Density: need area estimates (use GAM binned data instead)
                pass  # Smooth predictions carry the story for density

        # GAM smooth curves
        preds = gam_preds[(gam_preds["subclass"] == sc) &
                          (gam_preds["outcome"] == outcome)]
        for dx, color in [("Control", CTRL_COLOR), ("SCZ", SCZ_COLOR)]:
            dx_pred = preds[preds["diagnosis"] == dx].sort_values("depth")
            if len(dx_pred) > 0:
                ax.plot(dx_pred["depth"], dx_pred["predicted"],
                       color=color, linewidth=2, label=dx)

        # Significance annotation
        res = gam_results[(gam_results["subclass"] == sc) &
                          (gam_results["outcome"] == outcome)]
        if len(res) > 0:
            fdr = res.iloc[0].get("interaction_fdr", 1.0)
            pval = res.iloc[0].get("interaction_pval", 1.0)
            if pd.notna(fdr) and fdr < 0.05:
                ax.set_title(f"{sc}\np={pval:.3g} (FDR={fdr:.3g}) ***",
                            fontsize=10, fontweight="bold")
            elif pd.notna(pval) and pval < 0.05:
                ax.set_title(f"{sc}\np={pval:.3g} (FDR={fdr:.3g}) *",
                            fontsize=10)
            else:
                ax.set_title(sc, fontsize=10)
        else:
            ax.set_title(sc, fontsize=10)

        ax.set_xlim(0, 1)
        ax.set_xlabel("Cortical Depth", fontsize=8)
        y_label = "Proportion" if outcome == "proportion" else "Density (cells/mm²)"
        ax.set_ylabel(y_label, fontsize=8)
        ax.tick_params(labelsize=7)

    # Remove empty axes
    for idx in range(n_sc, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row, col].set_visible(False)

    # Legend
    legend_elements = [Line2D([0], [0], color=CTRL_COLOR, lw=2, label="Control"),
                       Line2D([0], [0], color=SCZ_COLOR, lw=2, label="SCZ")]
    fig.legend(handles=legend_elements, loc="upper right", fontsize=12,
              framealpha=0.9)

    fig.suptitle(f"Depth Profiles: {outcome.replace('_', ' ').title()}\n"
                 f"SCZ (n=12) vs Control (n=12)",
                 fontsize=16, fontweight="bold", y=1.02)
    fig.tight_layout()

    outpath = os.path.join(OUTPUT_DIR, f"depth_profiles_{outcome}.png")
    fig.savefig(outpath, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"  Saved: {outpath}")
    plt.close(fig)


def plot_layer_heatmap(layer_results, outcome="proportion"):
    """Heatmap of log2FC per (subclass × layer) with significance markers."""
    res = layer_results[layer_results["outcome"] == outcome].copy()
    if len(res) == 0:
        print(f"  No layer results for {outcome}, skipping heatmap")
        return

    # Pivot to matrix
    pivot = res.pivot_table(index="subclass", columns="layer",
                            values="log2fc", aggfunc="first")
    pivot_p = res.pivot_table(index="subclass", columns="layer",
                              values="fdr", aggfunc="first")

    # Order layers correctly
    layer_order = [l for l in CORTICAL_LAYERS if l in pivot.columns]
    pivot = pivot[layer_order]
    pivot_p = pivot_p.reindex(columns=layer_order)

    # Sort subclasses by mean absolute effect
    sort_order = pivot.abs().mean(axis=1).sort_values(ascending=False).index
    pivot = pivot.loc[sort_order]
    pivot_p = pivot_p.loc[sort_order]

    fig, ax = plt.subplots(figsize=(8, max(6, len(pivot) * 0.35)))

    # Color limits symmetric around 0
    vmax = min(pivot.abs().max().max(), 3.0)
    im = ax.imshow(pivot.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                   aspect="auto")

    # Significance stars
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            val = pivot.iloc[i, j]
            pval = pivot_p.iloc[i, j]
            if pd.notna(pval) and pval < 0.05:
                ax.text(j, i, "***", ha="center", va="center",
                       fontsize=12, fontweight="bold", color="black")
            elif pd.notna(pval) and pval < 0.1:
                ax.text(j, i, "*", ha="center", va="center",
                       fontsize=10, color="black")

    ax.set_xticks(range(len(layer_order)))
    ax.set_xticklabels(layer_order, fontsize=14)
    ax.set_yticks(range(len(pivot)))
    ax.set_yticklabels(pivot.index, fontsize=11)
    ax.set_xlabel("Cortical Layer", fontsize=14)

    cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("log₂FC (SCZ / Control)", fontsize=12)

    outcome_label = outcome.replace("_", " ").title()
    ax.set_title(f"Layer-Specific {outcome_label} Changes\n"
                 f"SCZ vs Control (*** FDR < 0.05, * FDR < 0.10)",
                 fontsize=14, fontweight="bold")

    fig.tight_layout()
    outpath = os.path.join(OUTPUT_DIR, f"layer_heatmap_{outcome}.png")
    fig.savefig(outpath, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"  Saved: {outpath}")
    plt.close(fig)


def plot_summary_volcano(layer_results):
    """Volcano plot: effect size vs significance for all layer × subclass tests."""
    res = layer_results.dropna(subset=["pval", "log2fc"]).copy()
    if len(res) == 0:
        print("  No valid results for volcano plot")
        return

    res["neg_log10_p"] = -np.log10(res["pval"].clip(1e-20))

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    for idx, outcome in enumerate(["proportion", "density_per_mm2"]):
        ax = axes[idx]
        sub = res[res["outcome"] == outcome]

        # Color by layer
        layer_colors = {"L1": "#E69F00", "L2/3": "#56B4E9", "L4": "#009E73",
                        "L5": "#F0E442", "L6": "#CC79A7"}

        for layer in CORTICAL_LAYERS:
            layer_data = sub[sub["layer"] == layer]
            ax.scatter(layer_data["log2fc"], layer_data["neg_log10_p"],
                      c=layer_colors.get(layer, "gray"), s=40, alpha=0.7,
                      label=layer, edgecolors="white", linewidth=0.5)

            # Label significant hits
            sig = layer_data[layer_data["fdr"] < 0.05] if "fdr" in layer_data.columns else pd.DataFrame()
            for _, row in sig.iterrows():
                ax.annotate(row["subclass"], (row["log2fc"], row["neg_log10_p"]),
                           fontsize=7, ha="center", va="bottom",
                           textcoords="offset points", xytext=(0, 5))

        # Significance thresholds
        if len(sub) > 0:
            ax.axhline(-np.log10(0.05), color="gray", linestyle="--",
                      linewidth=0.8, alpha=0.5, label="p=0.05")
        ax.axvline(0, color="gray", linestyle="-", linewidth=0.5, alpha=0.5)

        ax.set_xlabel("log₂FC (SCZ / Control)", fontsize=14)
        ax.set_ylabel("-log₁₀(p-value)", fontsize=14)
        ax.set_title(outcome.replace("_", " ").title(), fontsize=14,
                    fontweight="bold")
        ax.legend(fontsize=10, loc="upper right")
        ax.tick_params(labelsize=11)

    fig.suptitle("Depth-Stratified Differential Analysis: SCZ vs Control\n"
                 "Per-layer models (each dot = one subclass × layer test)",
                 fontsize=16, fontweight="bold")
    fig.tight_layout()

    outpath = os.path.join(OUTPUT_DIR, "layer_volcano.png")
    fig.savefig(outpath, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"  Saved: {outpath}")
    plt.close(fig)


def plot_top_hits_detail(cells_df, layer_df, layer_results, n_top=8):
    """Detailed panels for the top significant hits from layer models."""
    res = layer_results.dropna(subset=["pval"]).sort_values("pval").head(n_top)
    if len(res) == 0:
        return

    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.ravel()

    for idx, (_, hit) in enumerate(res.iterrows()):
        if idx >= len(axes):
            break
        ax = axes[idx]

        sc = hit["subclass"]
        layer = hit["layer"]
        outcome = hit["outcome"]

        # Get data for this subclass × layer
        sub = layer_df[(layer_df["subclass_label"] == sc) &
                       (layer_df["layer"] == layer)].copy()

        y_col = outcome
        for dx, color, offset in [("Control", CTRL_COLOR, -0.15),
                                   ("SCZ", SCZ_COLOR, 0.15)]:
            dx_data = sub[sub["diagnosis"] == dx]
            x_pos = np.ones(len(dx_data)) * (0.5 + offset)
            ax.scatter(x_pos, dx_data[y_col], c=color, s=50, alpha=0.7,
                      edgecolors="white", zorder=3)
            ax.bar(0.5 + offset, dx_data[y_col].mean(), width=0.25,
                  color=color, alpha=0.3, edgecolor=color, linewidth=1.5)

        fdr = hit.get("fdr", np.nan)
        stars = "***" if pd.notna(fdr) and fdr < 0.05 else ("*" if pd.notna(fdr) and fdr < 0.1 else "")
        ax.set_title(f"{sc} in {layer}\n{outcome}: p={hit['pval']:.3g} "
                    f"(FDR={fdr:.3g}) {stars}",
                    fontsize=10, fontweight="bold" if stars else "normal")
        ax.set_xticks([0.35, 0.65])
        ax.set_xticklabels(["Control", "SCZ"], fontsize=11)
        y_label = "Proportion" if outcome == "proportion" else "Density (cells/mm²)"
        ax.set_ylabel(y_label, fontsize=10)

    for idx in range(len(res), len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle("Top Hits: Layer-Specific Proportion/Density Differences",
                 fontsize=16, fontweight="bold")
    fig.tight_layout()

    outpath = os.path.join(OUTPUT_DIR, "top_hits_detail.png")
    fig.savefig(outpath, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"  Saved: {outpath}")
    plt.close(fig)


def main():
    print("Loading results...")
    gam_results = pd.read_csv(os.path.join(OUTPUT_DIR, "gam_results.csv"))
    gam_preds = pd.read_csv(os.path.join(OUTPUT_DIR, "gam_smooth_predictions.csv"))
    layer_results = pd.read_csv(os.path.join(OUTPUT_DIR, "layer_model_results.csv"))
    cells_df = pd.read_csv(os.path.join(OUTPUT_DIR, "cell_level_data.csv"))
    layer_df = pd.read_csv(os.path.join(OUTPUT_DIR, "layer_counts.csv"))

    print(f"  GAM: {len(gam_results)} tests, {len(gam_preds)} predictions")
    print(f"  Layer: {len(layer_results)} tests")
    print(f"  Cells: {len(cells_df):,}")

    # ── Depth profiles ──
    print("\nPlotting depth profiles...")
    for outcome in ["proportion", "density_per_mm2"]:
        plot_depth_profiles(cells_df, gam_preds, gam_results, outcome=outcome)

    # ── Layer heatmaps ──
    print("\nPlotting layer heatmaps...")
    for outcome in ["proportion", "density_per_mm2"]:
        plot_layer_heatmap(layer_results, outcome=outcome)

    # ── Volcano ──
    print("\nPlotting volcano...")
    plot_summary_volcano(layer_results)

    # ── Top hits detail ──
    print("\nPlotting top hits...")
    plot_top_hits_detail(cells_df, layer_df, layer_results)

    print("\nDone!")


if __name__ == "__main__":
    main()
