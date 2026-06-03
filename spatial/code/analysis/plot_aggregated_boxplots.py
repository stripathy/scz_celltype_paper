#!/usr/bin/env python3
"""
Aggregated boxplots: summed proportions and densities for vulnerable
Sst subtypes and all L6b cells (subclass), SCZ vs Control.

Layout (2x2):
  Row 1: Vulnerable Sst subtypes (Sst_2 + Sst_22 + Sst_25 + Sst_20 + Sst_3)
  Row 2: All L6b cells (subclass level)
  Col 1: Proportion (% of cortical cells)
  Col 2: Density (cells / mm²)

T-test p-values reported on each panel.

Input:
  output/presentation/xenium_composition_by_sample.csv  (Sst supertypes)
  output/density_analysis/density_per_sample_supertype.csv  (all L6b supertypes)

Output: output/presentation/slide_aggregated_boxplots.png
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    BASE_DIR, PRESENTATION_DIR, DX_COLORS, BG_COLOR, SST_TYPES,
    EXCLUDE_SAMPLES, format_pval,
)

OUT_DIR = PRESENTATION_DIR
DATA_CSV = os.path.join(OUT_DIR, "xenium_composition_by_sample.csv")
DENSITY_RAW = os.path.join(BASE_DIR, "output", "density_analysis",
                            "density_per_sample_supertype.csv")

BG = BG_COLOR


def aggregate_group(data, cell_types, metric):
    """Sum metric across cell types for each sample, return per-sample Series."""
    sub = data[data["celltype"].isin(cell_types)]
    agg = sub.groupby(["sample_id", "diagnosis"])[metric].sum().reset_index()
    return agg


def aggregate_by_prefix(data, prefix, metric):
    """Sum metric across all supertypes matching prefix for each sample."""
    sub = data[data["celltype"].str.startswith(prefix + "_")]
    agg = sub.groupby(["sample_id", "diagnosis"])[metric].sum().reset_index()
    return agg


def make_panel(ax, agg, metric_col, ylabel, title, subtitle):
    """Draw boxplot + strip for aggregated data."""
    ctrl = agg[agg["diagnosis"] == "Control"][metric_col].values
    scz = agg[agg["diagnosis"] == "SCZ"][metric_col].values

    positions = [0, 1]
    bp = ax.boxplot([ctrl, scz], positions=positions,
                    widths=0.5, patch_artist=True, showfliers=False,
                    medianprops=dict(color="white", linewidth=2),
                    whiskerprops=dict(color="#888888", linewidth=1.2),
                    capprops=dict(color="#888888", linewidth=1.2))

    for patch, color in zip(bp['boxes'], [DX_COLORS["Control"], DX_COLORS["SCZ"]]):
        patch.set_facecolor(color)
        patch.set_alpha(0.4)
        patch.set_edgecolor("white")
        patch.set_linewidth(1.0)

    # Individual points
    for i, (vals, dx) in enumerate([(ctrl, "Control"), (scz, "SCZ")]):
        jitter = np.random.default_rng(42).uniform(-0.15, 0.15, len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals,
                   c=DX_COLORS[dx], s=55, alpha=0.85, edgecolors="white",
                   linewidths=0.7, zorder=5)

    n_ctrl = len(ctrl)
    n_scz = len(scz)
    ax.set_xticks(positions)
    ax.set_xticklabels([f"Control\n(n={n_ctrl})", f"SCZ\n(n={n_scz})"],
                       fontsize=14, color="white")

    # T-test
    _, pval = ttest_ind(ctrl, scz, equal_var=False)
    pval_str = format_pval(pval)
    pval_color = "white" if pval < 0.05 else ("#cccccc" if pval < 0.1 else "#888888")
    pval_weight = "bold" if pval < 0.1 else "normal"

    # Report means too
    ctrl_mean = np.mean(ctrl)
    scz_mean = np.mean(scz)
    fc_str = f"Ctrl: {ctrl_mean:.2f}, SCZ: {scz_mean:.2f}"

    ax.text(0.5, 0.96, f"Welch's t-test {pval_str}\n{fc_str}",
            transform=ax.transAxes, ha="center", va="top",
            fontsize=12, color=pval_color, fontweight=pval_weight,
            fontstyle="italic",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#1a1a1a",
                      edgecolor="#444444", alpha=0.85))

    # Subtitle (cell type list) above the panel
    ax.text(0.5, 1.01, subtitle, transform=ax.transAxes, ha="center", va="bottom",
            fontsize=11, color="#aaaaaa", fontstyle="italic")

    ax.set_ylabel(ylabel, fontsize=15, color="white")
    ax.set_facecolor(BG)
    ax.tick_params(colors="white", labelsize=13)
    for spine in ax.spines.values():
        spine.set_color("#555555")
    ax.grid(axis="y", alpha=0.15, color="#555555")

    return pval


def main():
    # ── Sst data from composition CSV ──
    data = pd.read_csv(DATA_CSV)
    data = data[~data["sample_id"].isin(EXCLUDE_SAMPLES)]
    print(f"Loaded {len(data)} rows from composition CSV")

    sst_agg_types = SST_TYPES  # ["Sst_2", "Sst_22", "Sst_25", "Sst_20", "Sst_3"]
    sst_prop = aggregate_group(data, sst_agg_types, "proportion_pct")
    sst_dens = aggregate_group(data, sst_agg_types, "density_per_mm2")
    sst_label = " + ".join(sst_agg_types)

    # ── L6b data from density CSV (all L6b supertypes at subclass level) ──
    # Exclude samples with very low L6 representation (<3% L6 cells),
    # which lack sufficient deep cortex tissue coverage for meaningful
    # L6b comparisons: Br2039 (2.0%), Br5973 (2.4%), Br2719 (3.4%), Br5314 (2.3%)
    L6B_EXCLUDE = {"Br2039", "Br5973", "Br2719", "Br5314"}

    raw = pd.read_csv(DENSITY_RAW)
    raw = raw[raw["region"] == "cortical"]
    raw = raw[~raw["sample_id"].isin(EXCLUDE_SAMPLES | L6B_EXCLUDE)]

    # Rename to match expected columns
    raw = raw.rename(columns={"supertype": "celltype"})

    # Compute proportion_pct from count / total
    raw["proportion_pct"] = (raw["count"] / raw["total_cells"]) * 100

    l6b_types_found = sorted(raw[raw["celltype"].str.startswith("L6b_")]["celltype"].unique())
    print(f"L6b supertypes found: {l6b_types_found}")
    print(f"L6b samples excluded (low L6 coverage): {sorted(L6B_EXCLUDE)}")

    l6b_prop = aggregate_by_prefix(raw, "L6b", "proportion_pct")
    l6b_dens = aggregate_by_prefix(raw, "L6b", "density_per_mm2")
    l6b_label = "All L6b cells (subclass, excl. 4 low-L6 samples)"

    # --- Figure ---
    fig, axes = plt.subplots(2, 2, figsize=(12, 11), facecolor=BG)

    # Row 1: Sst
    p1 = make_panel(axes[0, 0], sst_prop, "proportion_pct",
                     "% of cortical cells", "", sst_label)
    p2 = make_panel(axes[0, 1], sst_dens, "density_per_mm2",
                     "cells / mm²", "", sst_label)

    # Row 2: L6b (all subtypes)
    p3 = make_panel(axes[1, 0], l6b_prop, "proportion_pct",
                     "% of cortical cells", "", l6b_label)
    p4 = make_panel(axes[1, 1], l6b_dens, "density_per_mm2",
                     "cells / mm²", "", l6b_label)

    print(f"\nSst proportion: {format_pval(p1)}")
    print(f"Sst density:    {format_pval(p2)}")
    print(f"L6b proportion: {format_pval(p3)}")
    print(f"L6b density:    {format_pval(p4)}")

    plt.tight_layout(pad=2.0, rect=[0.08, 0.0, 1.0, 0.94])

    # Row labels (rotated 90°, left side)
    fig.text(0.02, 0.73, "Vulnerable\nSst subtypes",
             ha="center", va="center", fontsize=18, fontweight="bold",
             color="white", rotation=90, transform=fig.transFigure)
    fig.text(0.02, 0.28, "L6b\n(subclass)",
             ha="center", va="center", fontsize=18, fontweight="bold",
             color="white", rotation=90, transform=fig.transFigure)

    # Column headers
    fig.text(0.32, 0.96, "Proportion", ha="center", va="bottom",
             fontsize=18, fontweight="bold", color="#dddddd",
             transform=fig.transFigure)
    fig.text(0.74, 0.96, "Density", ha="center", va="bottom",
             fontsize=18, fontweight="bold", color="#dddddd",
             transform=fig.transFigure)

    outpath = os.path.join(OUT_DIR, "slide_aggregated_boxplots.png")
    plt.savefig(outpath, dpi=200, facecolor=BG, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {outpath}")


if __name__ == "__main__":
    main()
