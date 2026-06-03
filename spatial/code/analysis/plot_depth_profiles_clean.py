#!/usr/bin/env python3
"""
Publication-quality depth profiles: cell type proportions and density by diagnosis.

Two figure pairs (proportion + density):
  1. Neuronal subclasses (Glutamatergic + GABAergic)
  2. Non-neuronal subclasses

Cowplot-style aesthetics: shared axis labels, no redundant labeling,
large readable text, clean grid layout.

Proportions are shown as percentages within each class (neuronal or
non-neuronal). Density is shown as cells/mm² (convex hull area per bin).

Per-bin significance: at each depth bin, tests SCZ vs Control via
Wilcoxon rank-sum, marks nominally significant bins (p<0.05) with
a star marker on the right margin.

Output:
  output/depth_proportions/depth_profiles_neuronal.png
  output/depth_proportions/depth_profiles_nonneuronal.png
  output/depth_proportions/depth_density_neuronal.png
  output/depth_proportions/depth_density_nonneuronal.png
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import BASE_DIR, infer_class

CRUMBLR_DIR = os.path.join(BASE_DIR, "output", "crumblr")
DEPTH_DIR = os.path.join(BASE_DIR, "output", "depth_proportions")

CTRL_COLOR = "#4477AA"
SCZ_COLOR = "#CC3311"
SIG_COLOR = "#333333"

LAYER_BOUNDS = {
    "L1": (0.00, 0.1225),
    "L2/3": (0.1225, 0.4696),
    "L4": (0.4696, 0.5443),
    "L5": (0.5443, 0.7079),
    "L6": (0.7079, 0.9275),
}


def compute_per_bin_significance(depth_df, celltype, metric="pct",
                                 class_filter=None):
    """Test SCZ vs Control metric at each depth bin via Wilcoxon.

    metric: 'pct' for within-class percentage, 'density' for cells/mm².
    Returns dict: depth_midpoint -> p-value
    """
    df = depth_df.copy()

    if metric == "pct":
        if class_filter is not None:
            class_totals = (df.groupby(["donor", "depth_bin"])["count"]
                            .sum().reset_index(name="class_total"))
            df = df.merge(class_totals, on=["donor", "depth_bin"])
            df["value"] = 100.0 * df["count"] / df["class_total"]
        else:
            df["value"] = 100.0 * df["count"] / df["total"]
    else:  # density
        df["value"] = df["density_per_mm2"]

    ct_data = df[df["celltype"] == celltype]

    results = {}
    for mid in sorted(ct_data["depth_midpoint"].unique()):
        bin_data = ct_data[ct_data["depth_midpoint"] == mid]
        ctrl = bin_data[bin_data["diagnosis"] == "Control"]["value"].dropna().values
        scz = bin_data[bin_data["diagnosis"] == "SCZ"]["value"].dropna().values
        if len(ctrl) >= 3 and len(scz) >= 3:
            try:
                _, pval = mannwhitneyu(ctrl, scz, alternative="two-sided")
                results[mid] = pval
            except ValueError:
                results[mid] = 1.0
        else:
            results[mid] = 1.0
    return results


def plot_depth_panel(ax, depth_df, celltype, bin_pvals,
                     show_ylabel=False, show_xlabel=False,
                     show_layer_labels=False, class_filter=None,
                     metric="pct"):
    """Plot a single depth profile panel."""
    df = depth_df.copy()

    if metric == "pct":
        if class_filter is not None:
            class_totals = (df.groupby(["donor", "depth_bin"])["count"]
                            .sum().reset_index(name="class_total"))
            df = df.merge(class_totals, on=["donor", "depth_bin"])
            df["value"] = 100.0 * df["count"] / df["class_total"]
        else:
            df["value"] = 100.0 * df["count"] / df["total"]
    else:
        df["value"] = df["density_per_mm2"]

    ct_data = df[df["celltype"] == celltype]

    # Alternating layer shading
    shaded = False
    for layer, (lo, hi) in LAYER_BOUNDS.items():
        if shaded:
            ax.axhspan(lo, hi, alpha=0.045, color="#000000", zorder=0)
        shaded = not shaded

    # Layer boundary lines
    for _, (lo, hi) in LAYER_BOUNDS.items():
        ax.axhline(lo, color="#d0d0d0", linewidth=0.5, zorder=1)
    ax.axhline(0.90, color="#d0d0d0", linewidth=0.5, zorder=1)

    # Layer labels on leftmost column — placed far left to avoid y-tick overlap
    if show_layer_labels:
        for layer, (lo, hi) in LAYER_BOUNDS.items():
            ax.text(-0.22, (lo + hi) / 2, layer,
                    transform=ax.get_yaxis_transform(),
                    ha="center", va="center", fontsize=11, color="#777777",
                    fontweight="bold", clip_on=False)

    # Plot group mean + SEM
    for dx, color in [("Control", CTRL_COLOR), ("SCZ", SCZ_COLOR)]:
        dx_data = ct_data[ct_data["diagnosis"] == dx]
        summary = (dx_data.groupby("depth_midpoint")["value"]
                   .agg(["mean", "sem"]).reset_index().sort_values("depth_midpoint"))

        ax.plot(summary["mean"], summary["depth_midpoint"],
                color=color, linewidth=2.0, zorder=5)
        ax.fill_betweenx(summary["depth_midpoint"],
                         summary["mean"] - summary["sem"],
                         summary["mean"] + summary["sem"],
                         color=color, alpha=0.18, zorder=4)

    # Auto-scale x, then add right padding for stars
    ax.set_xlim(left=0)
    ax.autoscale(axis="x")
    xlim = ax.get_xlim()
    x_range = xlim[1] - xlim[0]
    ax.set_xlim(xlim[0], xlim[1] + x_range * 0.18)

    # Mark significant bins with stars at right margin
    # * p<0.05, ** p<0.01, *** p<0.005
    star_x = xlim[1] + x_range * 0.06
    for mid, pval in bin_pvals.items():
        if pval < 0.005:
            label = "***"
        elif pval < 0.01:
            label = "**"
        elif pval < 0.05:
            label = "*"
        else:
            continue
        ax.text(star_x, mid, label, fontsize=10, fontweight="bold",
               color=SIG_COLOR, ha="left", va="center",
               clip_on=False, zorder=10)

    ax.set_ylim(1.0, 0.0)  # Pia at top
    ax.tick_params(labelsize=12, length=3, width=0.5)

    if show_ylabel:
        ax.set_ylabel("Cortical Depth", fontsize=16, labelpad=45)
    else:
        ax.tick_params(axis="y", labelleft=False)

    # Always show x-tick labels (scales differ across panels)
    xlabel_text = "% of class" if metric == "pct" else "cells / mm²"
    if show_xlabel:
        ax.set_xlabel(xlabel_text, fontsize=14)
    # x-tick labels always visible since scales vary per panel

    # Clean spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.spines["left"].set_color("#888888")
    ax.spines["bottom"].set_color("#888888")

    ax.set_title(celltype, fontsize=15, fontweight="bold", pad=6)


def make_figure(depth_df, celltypes, class_label, class_filter, filename,
                ncols=None, metric="pct"):
    """Make a cowplot-style multi-panel figure for one class of cell types."""
    n_ct = len(celltypes)
    if n_ct == 0:
        return

    if ncols is None:
        ncols = min(n_ct, 6)
    nrows = int(np.ceil(n_ct / ncols))

    panel_w = 3.0
    panel_h = 4.2 if nrows == 1 else 3.8
    fig_w = panel_w * ncols + 2.2
    fig_h = panel_h * nrows + 2.5

    fig = plt.figure(figsize=(fig_w, fig_h))

    gs = gridspec.GridSpec(nrows, ncols, figure=fig,
                           wspace=0.06, hspace=0.28,
                           left=0.13, right=0.97,
                           top=0.91, bottom=0.08)

    # Filter to class celltypes
    class_celltypes_all = [ct for ct in depth_df["celltype"].unique()
                           if _matches_class(ct, class_filter)]
    df_class = depth_df[depth_df["celltype"].isin(class_celltypes_all)].copy()

    metric_label = "proportion" if metric == "pct" else "density"
    print(f"  Computing per-bin significance ({metric_label}) for {class_label}...")
    all_pvals = {}
    n_sig_total = 0
    for ct in celltypes:
        pvals = compute_per_bin_significance(
            df_class, ct, metric=metric, class_filter=class_filter)
        all_pvals[ct] = pvals
        n_sig = sum(1 for p in pvals.values() if p < 0.05)
        n_sig_total += n_sig
        if n_sig > 0:
            print(f"    {ct}: {n_sig}/{len(pvals)} bins nom. sig.")
    print(f"    Total: {n_sig_total} significant bins across {n_ct} cell types")

    for idx, ct in enumerate(celltypes):
        row, col = divmod(idx, ncols)
        ax = fig.add_subplot(gs[row, col])

        is_left = (col == 0)
        is_bottom = (row == nrows - 1) or (idx + ncols >= n_ct)

        plot_depth_panel(ax, df_class, ct, all_pvals[ct],
                        show_ylabel=is_left,
                        show_xlabel=is_bottom,
                        show_layer_labels=is_left,
                        class_filter=class_filter,
                        metric=metric)

    # Legend + significance key combined
    legend_elements = [
        Line2D([0], [0], color=CTRL_COLOR, lw=2.5, label="Control (n=12)"),
        Line2D([0], [0], color=SCZ_COLOR, lw=2.5, label="SCZ (n=12)"),
        Line2D([0], [0], color="none", label=""),  # spacer
        Line2D([0], [0], color="none",
               label="Wilcoxon per bin:"),
        Line2D([0], [0], color="none",
               label="* p<.05  ** p<.01  *** p<.005"),
    ]
    fig.legend(handles=legend_elements, loc="upper right",
              fontsize=12, framealpha=0.95, edgecolor="#cccccc",
              bbox_to_anchor=(0.97, 0.99))

    title_suffix = "Proportion" if metric == "pct" else "Density"
    fig.suptitle(f"{class_label} — {title_suffix}",
                 fontsize=20, fontweight="bold", x=0.13, ha="left", y=1.0)

    outpath = os.path.join(DEPTH_DIR, filename)
    fig.savefig(outpath, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"  Saved: {outpath}")
    plt.close(fig)


def _matches_class(celltype, class_filter):
    """Check if celltype belongs to the class filter."""
    cls = infer_class(celltype)
    if class_filter == "neuronal":
        return cls in ("Glut", "GABA")
    elif class_filter == "nonneuronal":
        return cls in ("NN", "Other")
    return True


def main():
    print("Loading data...")
    depth_df = pd.read_csv(os.path.join(CRUMBLR_DIR,
                                         "crumblr_depth_input_subclass.csv"))
    print(f"  {len(depth_df):,} rows, {depth_df['donor'].nunique()} donors, "
          f"{depth_df['celltype'].nunique()} types")

    has_density = "density_per_mm2" in depth_df.columns
    if has_density:
        print(f"  Density data available ({depth_df['density_per_mm2'].notna().sum():,} non-null)")
    else:
        print("  WARNING: no density data — run build_crumblr_depth_input.py first")

    all_cts = depth_df["celltype"].unique()

    glut = sorted([ct for ct in all_cts if infer_class(ct) == "Glut"])
    gaba = sorted([ct for ct in all_cts if infer_class(ct) == "GABA"])
    nn = sorted([ct for ct in all_cts if infer_class(ct) in ("NN", "Other")])

    neuronal = glut + gaba
    nonneuronal = nn

    print(f"  Neuronal: {len(neuronal)} ({len(glut)} Glut + {len(gaba)} GABA)")
    print(f"  Non-neuronal: {len(nonneuronal)}")

    # ── Proportion figures ──
    print("\n── Proportion figures ──")
    make_figure(depth_df, neuronal, "Neuronal Subclasses",
                class_filter="neuronal", ncols=6, metric="pct",
                filename="depth_profiles_neuronal.png")

    make_figure(depth_df, nonneuronal, "Non-Neuronal Subclasses",
                class_filter="nonneuronal", ncols=6, metric="pct",
                filename="depth_profiles_nonneuronal.png")

    # ── Density figures ──
    if has_density:
        print("\n── Density figures ──")
        make_figure(depth_df, neuronal, "Neuronal Subclasses",
                    class_filter=None, ncols=6, metric="density",
                    filename="depth_density_neuronal.png")

        make_figure(depth_df, nonneuronal, "Non-Neuronal Subclasses",
                    class_filter=None, ncols=6, metric="density",
                    filename="depth_density_nonneuronal.png")

    print("\nDone!")


if __name__ == "__main__":
    main()
