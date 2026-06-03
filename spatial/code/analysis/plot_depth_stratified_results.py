#!/usr/bin/env python3
"""
Visualize depth-stratified cell type proportion and density analyses.

Generates figures for:
1. Depth profiles (quantile-binned): proportion of each subclass as a function
   of cortical depth, split by SCZ vs Control, with layer boundaries marked.
2. Layer heatmaps: log2FC per (subclass x layer) with significance stars.
3. CLR depth x diagnosis interaction bar chart + significant hits detail.
4. Volcano plot across all per-layer tests.

All depth plots use y-axis = depth with pia at TOP and WM at BOTTOM.

Requires:
  output/crumblr/crumblr_depth_input_subclass.csv
  output/crumblr/crumblr_depth_results_subclass.csv
  output/crumblr/crumblr_depth_interaction_subclass.csv
  output/depth_proportions/layer_model_results.csv
  output/depth_proportions/layer_counts.csv

Usage:
    python3 -u plot_depth_stratified_results.py
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import BASE_DIR, infer_class

CRUMBLR_DIR = os.path.join(BASE_DIR, "output", "crumblr")
DEPTH_DIR = os.path.join(BASE_DIR, "output", "depth_proportions")

# Colors
CTRL_COLOR = "#4477AA"
SCZ_COLOR = "#CC3311"

# Layer boundaries (depth from pia: 0=pia, 1=WM)
LAYER_BOUNDS = {
    "L1": (0.00, 0.1225),
    "L2/3": (0.1225, 0.4696),
    "L4": (0.4696, 0.5443),
    "L5": (0.5443, 0.7079),
    "L6": (0.7079, 0.9275),
}
LAYER_COLORS_LIGHT = {
    "L1": "#FFF3E0",
    "L2/3": "#E3F2FD",
    "L4": "#E8F5E9",
    "L5": "#FFF9C4",
    "L6": "#FCE4EC",
}
CORTICAL_LAYERS = ["L1", "L2/3", "L4", "L5", "L6"]


def add_layer_shading(ax, orientation="horizontal"):
    """Add layer boundary shading to a depth axis.

    orientation: 'horizontal' means depth is on x-axis,
                 'vertical' means depth is on y-axis (pia at top).
    """
    for layer, (lo, hi) in LAYER_BOUNDS.items():
        color = LAYER_COLORS_LIGHT.get(layer, "#f0f0f0")
        if orientation == "vertical":
            ax.axhspan(lo, hi, alpha=0.15, color=color, zorder=0)
            ax.text(0.98, (lo + hi) / 2, layer, transform=ax.get_yaxis_transform(),
                    ha="right", va="center", fontsize=9, color="#666666",
                    fontweight="bold", alpha=0.8)
        else:
            ax.axvspan(lo, hi, alpha=0.15, color=color, zorder=0)
            ax.text((lo + hi) / 2, 0.98, layer, transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontsize=8, color="#666666",
                    fontweight="bold", alpha=0.8)


def add_layer_lines(ax, orientation="vertical"):
    """Add dashed lines at layer boundaries."""
    for layer, (lo, hi) in LAYER_BOUNDS.items():
        if orientation == "vertical":
            ax.axhline(lo, color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
            ax.axhline(hi, color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
        else:
            ax.axvline(lo, color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)
            ax.axvline(hi, color="#999999", linestyle=":", linewidth=0.5, alpha=0.5)


def order_subclasses(subclasses):
    """Order subclasses: Glutamatergic, GABAergic, Non-neuronal."""
    classes = {sc: infer_class(sc) for sc in subclasses}
    order = []
    for cls in ["Glut", "GABA", "NN", "Other"]:
        order.extend(sorted([sc for sc, c in classes.items() if c == cls]))
    return order


def plot_depth_profiles_vertical(depth_df, interaction_results):
    """Plot depth profiles with depth on y-axis (pia at top).

    One panel per subclass. Lines show mean proportion at each depth bin
    for SCZ and Control, with SEM shading.
    """
    subclasses = order_subclasses(depth_df["celltype"].unique())
    n_sc = len(subclasses)
    ncols = 6
    nrows = int(np.ceil(n_sc / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(3.5 * ncols, 4.5 * nrows))
    if nrows == 1:
        axes = axes.reshape(1, -1)

    # Compute proportion = count / total for each (donor, depth_bin, celltype)
    depth_df = depth_df.copy()
    depth_df["proportion"] = depth_df["count"] / depth_df["total"]

    for idx, sc in enumerate(subclasses):
        row, col = divmod(idx, ncols)
        ax = axes[row, col]

        sc_data = depth_df[depth_df["celltype"] == sc]

        add_layer_shading(ax, orientation="vertical")

        for dx, color in [("Control", CTRL_COLOR), ("SCZ", SCZ_COLOR)]:
            dx_data = sc_data[sc_data["diagnosis"] == dx]

            # Mean + SEM per depth_bin across donors
            summary = (dx_data.groupby("depth_midpoint")["proportion"]
                       .agg(["mean", "sem", "count"]).reset_index())
            summary = summary.sort_values("depth_midpoint")

            ax.plot(summary["mean"], summary["depth_midpoint"],
                    color=color, linewidth=2, label=dx, zorder=5)
            ax.fill_betweenx(summary["depth_midpoint"],
                             summary["mean"] - summary["sem"],
                             summary["mean"] + summary["sem"],
                             color=color, alpha=0.2, zorder=4)

        # Significance annotation from interaction results
        if interaction_results is not None:
            hit = interaction_results[interaction_results["celltype"] == sc]
            if len(hit) > 0:
                fdr = hit.iloc[0]["FDR"]
                pval = hit.iloc[0]["P.Value"]
                if fdr < 0.05:
                    ax.set_title(f"{sc}\nFDR={fdr:.3g} ***",
                                fontsize=11, fontweight="bold", color="#CC0000")
                elif pval < 0.05:
                    ax.set_title(f"{sc}\np={pval:.3g} *",
                                fontsize=11, fontweight="bold")
                else:
                    ax.set_title(sc, fontsize=11, fontweight="bold")
            else:
                ax.set_title(sc, fontsize=11, fontweight="bold")
        else:
            ax.set_title(sc, fontsize=11, fontweight="bold")

        ax.set_ylim(1.0, 0.0)  # Pia at top, WM at bottom
        ax.set_ylabel("Depth (pia → WM)", fontsize=9)
        ax.set_xlabel("Proportion", fontsize=9)
        ax.tick_params(labelsize=8)

        # Color class
        cls = infer_class(sc)
        cls_colors = {"Glut": "#4477AA", "GABA": "#EE6677", "NN": "#228833"}
        border_color = cls_colors.get(cls, "#888888")
        for spine in ax.spines.values():
            spine.set_color(border_color)
            spine.set_linewidth(1.5)

    # Remove empty axes
    for idx in range(n_sc, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row, col].set_visible(False)

    # Legend
    legend_elements = [
        Line2D([0], [0], color=CTRL_COLOR, lw=2, label="Control (n=12)"),
        Line2D([0], [0], color=SCZ_COLOR, lw=2, label="SCZ (n=12)"),
    ]
    fig.legend(handles=legend_elements, loc="upper right", fontsize=13,
              framealpha=0.9, bbox_to_anchor=(0.98, 0.99))

    fig.suptitle("Depth Profiles: Cell Type Proportions by Diagnosis\n"
                 "(density-adaptive quantile bins, pia at top)",
                 fontsize=18, fontweight="bold", y=1.01)
    fig.tight_layout()

    outpath = os.path.join(DEPTH_DIR, "depth_profiles_vertical.png")
    fig.savefig(outpath, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"  Saved: {outpath}")
    plt.close(fig)


def plot_layer_heatmap_vertical(layer_results, outcome="proportion"):
    """Heatmap with layers on y-axis (pia at top), subclasses on x-axis."""
    res = layer_results[layer_results["outcome"] == outcome].copy()
    if len(res) == 0:
        print(f"  No layer results for {outcome}")
        return

    # Pivot
    pivot = res.pivot_table(index="layer", columns="subclass",
                            values="log2fc", aggfunc="first")
    pivot_fdr = res.pivot_table(index="layer", columns="subclass",
                                values="fdr", aggfunc="first")

    # Order layers pia to WM
    layer_order = [l for l in CORTICAL_LAYERS if l in pivot.index]
    pivot = pivot.loc[layer_order]
    pivot_fdr = pivot_fdr.loc[layer_order]

    # Order subclasses by class
    sc_order = order_subclasses(pivot.columns)
    sc_order = [s for s in sc_order if s in pivot.columns]
    pivot = pivot[sc_order]
    pivot_fdr = pivot_fdr[sc_order]

    fig, ax = plt.subplots(figsize=(max(10, len(sc_order) * 0.6), 5))

    vmax = min(pivot.abs().max().max(), 2.0)
    im = ax.imshow(pivot.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                   aspect="auto")

    # Significance stars
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            fdr = pivot_fdr.iloc[i, j]
            if pd.notna(fdr) and fdr < 0.05:
                ax.text(j, i, "***", ha="center", va="center",
                       fontsize=14, fontweight="bold", color="black")
            elif pd.notna(fdr) and fdr < 0.1:
                ax.text(j, i, "**", ha="center", va="center",
                       fontsize=11, color="black")

    ax.set_xticks(range(len(sc_order)))
    ax.set_xticklabels(sc_order, fontsize=10, rotation=45, ha="right")
    ax.set_yticks(range(len(layer_order)))
    ax.set_yticklabels(layer_order, fontsize=13)
    ax.set_ylabel("Layer (pia → WM)", fontsize=14)

    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label("log₂FC (SCZ / Control)", fontsize=12)

    outcome_label = outcome.replace("_", " ").replace("per mm2", "/mm²").title()
    ax.set_title(f"Layer-Specific {outcome_label} Changes\n"
                 f"SCZ vs Control (*** FDR<0.05, ** FDR<0.10)",
                 fontsize=15, fontweight="bold")

    fig.tight_layout()
    outpath = os.path.join(DEPTH_DIR, f"layer_heatmap_{outcome}_vertical.png")
    fig.savefig(outpath, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"  Saved: {outpath}")
    plt.close(fig)


def plot_interaction_bar(all_results, interaction_results):
    """Bar chart of depth × diagnosis interaction effect sizes."""
    if interaction_results is None or len(interaction_results) == 0:
        return

    df = interaction_results.copy()
    df["class"] = df["celltype"].apply(infer_class)
    df = df.sort_values("logFC", ascending=True)

    fig, ax = plt.subplots(figsize=(10, 8))

    class_colors = {"Glut": "#4477AA", "GABA": "#EE6677", "NN": "#228833", "Other": "#888888"}
    colors = [class_colors.get(c, "#888888") for c in df["class"]]

    bars = ax.barh(range(len(df)), df["logFC"], color=colors,
                   edgecolor="white", linewidth=0.5)

    # Add significance stars
    for i, (_, row) in enumerate(df.iterrows()):
        x = row["logFC"]
        star = ""
        if row["FDR"] < 0.05:
            star = " ***"
        elif row["FDR"] < 0.1:
            star = " **"
        elif row["P.Value"] < 0.05:
            star = " *"

        offset = 0.02 if x >= 0 else -0.02
        ha = "left" if x >= 0 else "right"
        ax.text(x + offset, i, f"{row['logFC']:+.3f}{star}",
               va="center", ha=ha, fontsize=9,
               fontweight="bold" if row["FDR"] < 0.05 else "normal")

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["celltype"], fontsize=11)
    ax.set_xlabel("Depth × Diagnosis Interaction (logFC)", fontsize=14)
    ax.axvline(0, color="black", linewidth=0.8)

    # Legend for class colors
    legend_elements = [
        mpatches.Patch(facecolor=class_colors["Glut"], label="Glutamatergic"),
        mpatches.Patch(facecolor=class_colors["GABA"], label="GABAergic"),
        mpatches.Patch(facecolor=class_colors["NN"], label="Non-neuronal"),
    ]
    ax.legend(handles=legend_elements, fontsize=11, loc="lower right")

    ax.set_title("Depth × Diagnosis Interaction (CLR mixed models)\n"
                 "Positive = SCZ depth profile flatter than Control\n"
                 "*** FDR<0.05  ** FDR<0.10  * p<0.05",
                 fontsize=14, fontweight="bold")
    fig.tight_layout()

    outpath = os.path.join(DEPTH_DIR, "depth_interaction_bar.png")
    fig.savefig(outpath, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"  Saved: {outpath}")
    plt.close(fig)


def plot_significant_hits_detail(depth_df, interaction_results, n_top=4):
    """Detailed depth profiles for top interaction hits with individual donor traces."""
    if interaction_results is None or len(interaction_results) == 0:
        return

    top = interaction_results.sort_values("P.Value").head(n_top)

    fig, axes = plt.subplots(1, n_top, figsize=(5 * n_top, 7))
    if n_top == 1:
        axes = [axes]

    depth_df = depth_df.copy()
    depth_df["proportion"] = depth_df["count"] / depth_df["total"]

    for idx, (_, hit) in enumerate(top.iterrows()):
        ax = axes[idx]
        sc = hit["celltype"]

        sc_data = depth_df[depth_df["celltype"] == sc]

        add_layer_shading(ax, orientation="vertical")

        # Individual donor traces (thin lines)
        for dx, color in [("Control", CTRL_COLOR), ("SCZ", SCZ_COLOR)]:
            dx_data = sc_data[sc_data["diagnosis"] == dx]
            for donor in dx_data["donor"].unique():
                d_data = dx_data[dx_data["donor"] == donor].sort_values("depth_midpoint")
                ax.plot(d_data["proportion"], d_data["depth_midpoint"],
                       color=color, linewidth=0.5, alpha=0.3, zorder=3)

        # Group means (thick lines)
        for dx, color in [("Control", CTRL_COLOR), ("SCZ", SCZ_COLOR)]:
            dx_data = sc_data[sc_data["diagnosis"] == dx]
            summary = (dx_data.groupby("depth_midpoint")["proportion"]
                       .agg(["mean", "sem"]).reset_index().sort_values("depth_midpoint"))
            ax.plot(summary["mean"], summary["depth_midpoint"],
                   color=color, linewidth=3, label=dx, zorder=5)
            ax.fill_betweenx(summary["depth_midpoint"],
                             summary["mean"] - summary["sem"],
                             summary["mean"] + summary["sem"],
                             color=color, alpha=0.25, zorder=4)

        ax.set_ylim(1.0, 0.0)  # Pia at top
        ax.set_ylabel("Depth (pia → WM)", fontsize=12)
        ax.set_xlabel("Proportion", fontsize=12)

        fdr = hit["FDR"]
        pval = hit["P.Value"]
        star = "***" if fdr < 0.05 else ("**" if fdr < 0.1 else ("*" if pval < 0.05 else ""))
        ax.set_title(f"{sc}\nInteraction logFC={hit['logFC']:+.3f}\n"
                     f"p={pval:.2g}, FDR={fdr:.2g} {star}",
                     fontsize=13, fontweight="bold")
        ax.tick_params(labelsize=10)

        if idx == 0:
            ax.legend(fontsize=11, loc="lower right")

    fig.suptitle("Top Depth × Diagnosis Interactions\n"
                 "(thin lines = individual donors, thick = group mean ± SEM)",
                 fontsize=16, fontweight="bold", y=1.02)
    fig.tight_layout()

    outpath = os.path.join(DEPTH_DIR, "depth_interaction_significant_hits.png")
    fig.savefig(outpath, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"  Saved: {outpath}")
    plt.close(fig)


def plot_layer_volcano(layer_results):
    """Volcano plot across all per-layer tests."""
    res = layer_results.dropna(subset=["pval", "log2fc"]).copy()
    if len(res) == 0:
        return

    res["neg_log10_p"] = -np.log10(res["pval"].clip(1e-20))

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    layer_colors = {"L1": "#E69F00", "L2/3": "#56B4E9", "L4": "#009E73",
                    "L5": "#F0E442", "L6": "#CC79A7"}

    for idx, outcome in enumerate(["proportion", "density_per_mm2"]):
        ax = axes[idx]
        sub = res[res["outcome"] == outcome]

        for layer in CORTICAL_LAYERS:
            layer_data = sub[sub["layer"] == layer]
            ax.scatter(layer_data["log2fc"], layer_data["neg_log10_p"],
                      c=layer_colors.get(layer, "gray"), s=50, alpha=0.7,
                      label=layer, edgecolors="white", linewidth=0.5, zorder=3)

            # Label FDR-significant hits
            sig = layer_data[layer_data.get("fdr", pd.Series(dtype=float)) < 0.05]
            for _, row in sig.iterrows():
                ax.annotate(row["subclass"],
                           (row["log2fc"], row["neg_log10_p"]),
                           fontsize=9, ha="center", va="bottom",
                           textcoords="offset points", xytext=(0, 6),
                           fontweight="bold")

        ax.axhline(-np.log10(0.05), color="gray", linestyle="--",
                  linewidth=0.8, alpha=0.5)
        ax.axvline(0, color="gray", linestyle="-", linewidth=0.5, alpha=0.5)

        ax.set_xlabel("log₂FC (SCZ / Control)", fontsize=14)
        ax.set_ylabel("-log₁₀(p-value)", fontsize=14)
        outcome_label = outcome.replace("_", " ").replace("per mm2", "/mm²").title()
        ax.set_title(outcome_label, fontsize=15, fontweight="bold")
        ax.legend(fontsize=11, loc="upper right")
        ax.tick_params(labelsize=11)

    fig.suptitle("Per-Layer Differential Analysis: SCZ vs Control\n"
                 "(each dot = one subclass × layer test)",
                 fontsize=16, fontweight="bold")
    fig.tight_layout()

    outpath = os.path.join(DEPTH_DIR, "layer_volcano.png")
    fig.savefig(outpath, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"  Saved: {outpath}")
    plt.close(fig)


def main():
    print("Loading data...")

    # Crumblr depth-binned data
    depth_df = pd.read_csv(os.path.join(CRUMBLR_DIR,
                                         "crumblr_depth_input_subclass.csv"))
    print(f"  Depth data: {len(depth_df):,} rows, "
          f"{depth_df['donor'].nunique()} donors, "
          f"{depth_df['celltype'].nunique()} types, "
          f"{depth_df['depth_bin'].nunique()} bins")

    # CLR interaction results
    try:
        all_results = pd.read_csv(os.path.join(CRUMBLR_DIR,
                                                "crumblr_depth_results_subclass.csv"))
        interaction_results = pd.read_csv(os.path.join(CRUMBLR_DIR,
                                                        "crumblr_depth_interaction_subclass.csv"))
        print(f"  CLR results: {len(all_results)} rows, "
              f"{len(interaction_results)} interaction tests")
    except FileNotFoundError:
        print("  WARNING: CLR results not found, skipping interaction plots")
        all_results = None
        interaction_results = None

    # Layer model results
    try:
        layer_results = pd.read_csv(os.path.join(DEPTH_DIR,
                                                   "layer_model_results.csv"))
        print(f"  Layer results: {len(layer_results)} tests")
    except FileNotFoundError:
        print("  WARNING: Layer results not found, skipping layer plots")
        layer_results = None

    # ── 1. Depth profiles (vertical, pia at top) ──
    print("\nPlotting depth profiles...")
    plot_depth_profiles_vertical(depth_df, interaction_results)

    # ── 2. Layer heatmaps ──
    if layer_results is not None:
        print("\nPlotting layer heatmaps...")
        for outcome in ["proportion", "density_per_mm2"]:
            plot_layer_heatmap_vertical(layer_results, outcome=outcome)

    # ── 3. CLR interaction bar chart ──
    if interaction_results is not None:
        print("\nPlotting interaction bar chart...")
        plot_interaction_bar(all_results, interaction_results)

    # ── 4. Significant hits detail ──
    if interaction_results is not None:
        print("\nPlotting significant hits detail...")
        plot_significant_hits_detail(depth_df, interaction_results, n_top=4)

    # ── 5. Layer volcano ──
    if layer_results is not None:
        print("\nPlotting layer volcano...")
        plot_layer_volcano(layer_results)

    print("\nDone!")


if __name__ == "__main__":
    main()
