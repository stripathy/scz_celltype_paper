#!/usr/bin/env python3
"""
Boxplots of cell type proportions (%) and densities for selected cell types
in Xenium data, SCZ vs Control.

Organized in 2 rows:
  Row 1: Sst subtypes + L2/3 IT_10
  Row 2: Layer 6 subtypes + L5/6 NP_4

Cell types are ordered left-to-right by increasing median depth from pia
(most superficial first).

Proportions and densities are raw (unadjusted).
P-values are from the crumblr compositional model
(which covaries for age + sex + diagnosis on CLR-transformed proportions).

Output: output/presentation/slide_xenium_proportion_boxplots.png
        output/presentation/slide_xenium_density_boxplots.png
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    SAMPLE_TO_DX, DX_COLORS, BG_COLOR, H5AD_DIR, PRESENTATION_DIR,
    CRUMBLR_DIR, SST_TYPES, CORTICAL_LAYERS, EXCLUDE_SAMPLES,
    format_pval, load_cells,
)

OUT_DIR = PRESENTATION_DIR
CRUMBLR_PATH = os.path.join(CRUMBLR_DIR, "crumblr_results_supertype.csv")
DEPTH_PATH = os.path.join(PRESENTATION_DIR, "median_depth_supertype.csv")

BG = BG_COLOR

# Row 1: Sst subtypes + L2/3 IT_10
ROW1_TYPES = SST_TYPES + ["L2/3 IT_10"]
ROW1_LABEL = "Sst subtypes + L2/3 IT"

# Row 2: Layer 6 subtypes + L5/6 NP_4
L6_TYPES = ["L6b_1", "L6b_2", "L6b_4", "L6 CT_1", "L6 CT_3"]
ROW2_TYPES = L6_TYPES + ["L5/6 NP_4"]
ROW2_LABEL = "Layer 6 + L5/6 NP subtypes"

ALL_TYPES = ROW1_TYPES + ROW2_TYPES


def sort_types_by_depth(type_list):
    """Sort cell types by increasing median depth from pia (most superficial first).

    Uses Xenium median depth from the cross-platform comparison output.
    """
    depth_df = pd.read_csv(DEPTH_PATH)
    depth_map = dict(zip(depth_df["supertype"], depth_df["median_depth_xenium"]))
    # Sort by increasing depth (low depth = superficial = left side of plot)
    return sorted(type_list, key=lambda ct: depth_map.get(ct, 0.5))


def compute_cortical_area_mm2(coords):
    """Compute convex hull area of cortical cells. Coords in um, returns mm2."""
    if len(coords) < 3:
        return np.nan
    try:
        hull = ConvexHull(coords)
        return hull.volume / 1e6  # um2 -> mm2
    except Exception:
        return np.nan


def load_all_samples():
    """Load per-sample data: cortical proportions and densities for target cell types."""
    records = []

    h5ad_files = sorted(glob.glob(os.path.join(H5AD_DIR, "*_annotated.h5ad")))

    for fpath in h5ad_files:
        sample_id = os.path.basename(fpath).replace("_annotated.h5ad", "")
        if sample_id in EXCLUDE_SAMPLES:
            print(f"  {sample_id}... SKIPPED (excluded)")
            continue
        print(f"  {sample_id}...", end="", flush=True)

        obs = load_cells(sample_id, cortical_only=True)
        n_cortical = len(obs)

        # Cortical area from convex hull
        cortical_coords = obs[["x", "y"]].values
        area_mm2 = compute_cortical_area_mm2(cortical_coords)

        print(f" {n_cortical:,} cortical, {area_mm2:.2f} mm2")

        dx = SAMPLE_TO_DX.get(sample_id, "Unknown")

        for ct in ALL_TYPES:
            ct_count = (obs["supertype_label"] == ct).sum()
            prop = ct_count / n_cortical if n_cortical > 0 else np.nan
            density = ct_count / area_mm2 if area_mm2 and area_mm2 > 0 else np.nan

            records.append({
                "sample_id": sample_id,
                "diagnosis": dx,
                "celltype": ct,
                "count": ct_count,
                "n_cortical": n_cortical,
                "cortical_area_mm2": area_mm2,
                "proportion_pct": prop * 100,  # as percent
                "density_per_mm2": density,
            })

    return pd.DataFrame(records)


def make_boxplot(ax, data, ct, crumblr_pval):
    """Draw a single boxplot panel for one cell type, SCZ vs Control."""
    ctrl = data[(data["celltype"] == ct) & (data["diagnosis"] == "Control")]["proportion_pct"].dropna()
    scz = data[(data["celltype"] == ct) & (data["diagnosis"] == "SCZ")]["proportion_pct"].dropna()

    positions = [0, 1]
    bp = ax.boxplot([ctrl.values, scz.values], positions=positions,
                    widths=0.5, patch_artist=True, showfliers=False,
                    medianprops=dict(color="white", linewidth=1.5),
                    whiskerprops=dict(color="#888888"),
                    capprops=dict(color="#888888"))

    for patch, color in zip(bp['boxes'], [DX_COLORS["Control"], DX_COLORS["SCZ"]]):
        patch.set_facecolor(color)
        patch.set_alpha(0.4)
        patch.set_edgecolor("white")
        patch.set_linewidth(0.8)

    # Overlay individual points with jitter
    n_ctrl = len(ctrl)
    n_scz = len(scz)
    for i, (vals, dx) in enumerate([(ctrl, "Control"), (scz, "SCZ")]):
        jitter = np.random.default_rng(42).uniform(-0.12, 0.12, len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals,
                   c=DX_COLORS[dx], s=40, alpha=0.85, edgecolors="white",
                   linewidths=0.5, zorder=5)

    ax.set_xticks(positions)
    ax.set_xticklabels([f"Ctrl\n(n={n_ctrl})", f"SCZ\n(n={n_scz})"],
                       fontsize=12, color="white")

    # Title = cell type name
    ax.set_title(ct, fontsize=15, fontweight="bold", color="white", pad=6)

    # Crumblr p-value annotation
    pval_str = format_pval(crumblr_pval)
    pval_color = "white" if crumblr_pval < 0.05 else ("#cccccc" if crumblr_pval < 0.1 else "#888888")
    pval_weight = "bold" if crumblr_pval < 0.1 else "normal"
    ax.text(0.5, 0.02, f"crumblr {pval_str}",
            transform=ax.transAxes, ha="center", va="bottom",
            fontsize=10, color=pval_color, fontweight=pval_weight,
            fontstyle="italic")

    ax.set_facecolor(BG)
    ax.tick_params(colors="white", labelsize=11)
    for spine in ax.spines.values():
        spine.set_color("#555555")
    ax.grid(axis="y", alpha=0.15, color="#555555")


def make_density_boxplot(ax, data, ct, crumblr_pval):
    """Draw a density boxplot panel with crumblr p-value."""
    ctrl = data[(data["celltype"] == ct) & (data["diagnosis"] == "Control")]["density_per_mm2"].dropna()
    scz = data[(data["celltype"] == ct) & (data["diagnosis"] == "SCZ")]["density_per_mm2"].dropna()

    positions = [0, 1]
    bp = ax.boxplot([ctrl.values, scz.values], positions=positions,
                    widths=0.5, patch_artist=True, showfliers=False,
                    medianprops=dict(color="white", linewidth=1.5),
                    whiskerprops=dict(color="#888888"),
                    capprops=dict(color="#888888"))

    for patch, color in zip(bp['boxes'], [DX_COLORS["Control"], DX_COLORS["SCZ"]]):
        patch.set_facecolor(color)
        patch.set_alpha(0.4)
        patch.set_edgecolor("white")
        patch.set_linewidth(0.8)

    # Overlay individual points with jitter
    n_ctrl = len(ctrl)
    n_scz = len(scz)
    for i, (vals, dx) in enumerate([(ctrl, "Control"), (scz, "SCZ")]):
        jitter = np.random.default_rng(42).uniform(-0.12, 0.12, len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals,
                   c=DX_COLORS[dx], s=40, alpha=0.85, edgecolors="white",
                   linewidths=0.5, zorder=5)

    ax.set_xticks(positions)
    ax.set_xticklabels([f"Ctrl\n(n={n_ctrl})", f"SCZ\n(n={n_scz})"],
                       fontsize=12, color="white")

    # Title = cell type name
    ax.set_title(ct, fontsize=15, fontweight="bold", color="white", pad=6)

    # Crumblr p-value
    pval_str = format_pval(crumblr_pval)
    pval_color = "white" if crumblr_pval < 0.05 else ("#cccccc" if crumblr_pval < 0.1 else "#888888")
    pval_weight = "bold" if crumblr_pval < 0.1 else "normal"
    ax.text(0.5, 0.02, f"crumblr {pval_str}",
            transform=ax.transAxes, ha="center", va="bottom",
            fontsize=10, color=pval_color, fontweight=pval_weight,
            fontstyle="italic")

    ax.set_facecolor(BG)
    ax.tick_params(colors="white", labelsize=11)
    for spine in ax.spines.values():
        spine.set_color("#555555")
    ax.grid(axis="y", alpha=0.15, color="#555555")


def main():
    print("Loading Xenium samples...")
    data = load_all_samples()

    # Save raw data
    csv_path = os.path.join(OUT_DIR, "xenium_composition_by_sample.csv")
    data.to_csv(csv_path, index=False)
    print(f"\nSaved data: {csv_path}")

    # Load crumblr results for p-values
    crumblr = pd.read_csv(CRUMBLR_PATH)
    crumblr_pvals = dict(zip(crumblr["celltype"], crumblr["P.Value"]))
    print("\nCrumblr p-values:")
    for ct in ALL_TYPES:
        p = crumblr_pvals.get(ct, np.nan)
        print(f"  {ct:25s} p = {p:.4f}")

    # Sort cell types by increasing median depth from pia (most superficial first)
    row1_sorted = sort_types_by_depth(ROW1_TYPES)
    row2_sorted = sort_types_by_depth(ROW2_TYPES)
    print(f"\nRow 1 order (by depth): {row1_sorted}")
    print(f"Row 2 order (by depth): {row2_sorted}")

    row_groups = [
        (ROW1_LABEL, row1_sorted),
        (ROW2_LABEL, row2_sorted),
    ]

    n_cols = max(len(ROW1_TYPES), len(ROW2_TYPES))
    n_rows = 2

    # ===== Proportion figure (raw) =====
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.5 * n_cols, 4.5 * n_rows),
                              facecolor=BG)

    for row_idx, (row_label, types) in enumerate(row_groups):
        for j, ct in enumerate(types):
            pval = crumblr_pvals.get(ct, np.nan)
            make_boxplot(axes[row_idx, j], data, ct, pval)

        # Hide unused columns
        for j in range(len(types), n_cols):
            axes[row_idx, j].set_visible(False)

        # Row label
        axes[row_idx, 0].set_ylabel("% of cortical cells",
                                     fontsize=12, color="white")

    # Row group labels on the left side
    row_positions = [0.75, 0.25]
    for (label, _), ypos in zip(row_groups, row_positions):
        fig.text(0.01, ypos, label, ha="left", va="center",
                 fontsize=16, fontweight="bold", color="#dddddd",
                 rotation=90, transform=fig.transFigure)

    fig.suptitle("Xenium SCZ vs Control: cell type proportions\n"
                 "(crumblr p-values, adj. for age + sex)",
                 fontsize=22, fontweight="bold", color="white", y=1.01)

    plt.tight_layout(pad=1.5, rect=[0.04, 0.0, 1.0, 0.97])

    outpath = os.path.join(OUT_DIR, "slide_xenium_proportion_boxplots.png")
    plt.savefig(outpath, dpi=200, facecolor=BG, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {outpath}")

    # ===== Density figure (raw) =====
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.5 * n_cols, 4.5 * n_rows),
                              facecolor=BG)

    print("\nCrumblr p-values for density:")
    for row_idx, (row_label, types) in enumerate(row_groups):
        for j, ct in enumerate(types):
            pval = crumblr_pvals.get(ct, np.nan)
            make_density_boxplot(axes[row_idx, j], data, ct, pval)
            print(f"  {ct:25s} crumblr p = {pval:.4f}")

        # Hide unused columns
        for j in range(len(types), n_cols):
            axes[row_idx, j].set_visible(False)

        # Row label
        axes[row_idx, 0].set_ylabel("cells / mm²",
                                     fontsize=12, color="white")

    # Row group labels
    for (label, _), ypos in zip(row_groups, row_positions):
        fig.text(0.01, ypos, label, ha="left", va="center",
                 fontsize=16, fontweight="bold", color="#dddddd",
                 rotation=90, transform=fig.transFigure)

    fig.suptitle("Xenium SCZ vs Control: cell type density\n"
                 "(crumblr p-values, adj. for age + sex)",
                 fontsize=22, fontweight="bold", color="white", y=1.01)

    plt.tight_layout(pad=1.5, rect=[0.04, 0.0, 1.0, 0.97])

    outpath = os.path.join(OUT_DIR, "slide_xenium_density_boxplots.png")
    plt.savefig(outpath, dpi=200, facecolor=BG, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {outpath}")


if __name__ == "__main__":
    main()
