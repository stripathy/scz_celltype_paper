#!/usr/bin/env python3
"""
Combined figure:
  Top: Stacked bar plot of layer cell counts across Xenium samples
       (sorted by ascending WM proportion, x-axis colored by diagnosis)
  Bottom: Spatial plot of the WM outlier Br2039 showing layer assignments,
          plus a 'typical' sample for comparison.

Output: output/presentation/slide_layer_bar_with_outlier.png
"""

import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    SAMPLE_TO_DX, DX_COLORS, BG_COLOR, H5AD_DIR, PRESENTATION_DIR, LAYER_COLORS,
)

OUT_DIR = PRESENTATION_DIR

BG = BG_COLOR

# Bottom-to-top stacking order: WM at bottom, L1 at top
LAYER_STACK_ORDER = ["Vascular", "WM", "L6", "L5", "L4", "L2/3", "L1"]
# Legend order: L1 first
LAYER_LEGEND_ORDER = ["L1", "L2/3", "L4", "L5", "L6", "WM", "Vascular"]


def load_sample(sample_id):
    """Load obs + spatial coords for a sample."""
    fpath = os.path.join(H5AD_DIR, f"{sample_id}_annotated.h5ad")
    adata = ad.read_h5ad(fpath)
    obs = adata.obs[["layer", "qc_pass"]].copy()
    coords = adata.obsm["spatial"]
    obs["x"] = coords[:, 0]
    obs["y"] = coords[:, 1]
    obs = obs[obs["qc_pass"] == True]
    obs["layer"] = obs["layer"].astype(str)
    return obs


def plot_spatial(ax, obs, title, s=0.6):
    """Plot cells colored by layer on a spatial axis."""
    x = obs["x"].values
    y = obs["y"].values
    layers = obs["layer"].values

    # Plot each layer
    for lname in LAYER_LEGEND_ORDER:
        mask = layers == lname
        if mask.sum() > 0:
            ax.scatter(x[mask], y[mask], s=s, color=[LAYER_COLORS[lname]],
                       alpha=0.7, rasterized=True, linewidths=0, zorder=3,
                       label=lname)

    # WM proportion annotation
    n_wm = (layers == "WM").sum()
    n_total = len(layers)
    wm_pct = n_wm / n_total * 100

    ax.text(0.03, 0.97,
            f"{n_total:,} cells\nWM: {n_wm:,} ({wm_pct:.0f}%)",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=14, color="#dddddd",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#333333",
                      edgecolor="#555555", alpha=0.85))

    ax.set_title(title, fontsize=20, fontweight="bold", color="white", pad=8)
    ax.set_facecolor(BG)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def main():
    # --- Collect per-sample layer counts ---
    print("Loading all samples for bar plot...")
    records = []
    h5ad_files = sorted([f for f in os.listdir(H5AD_DIR) if f.endswith("_annotated.h5ad")])

    for fname in h5ad_files:
        sample_id = fname.replace("_annotated.h5ad", "")
        fpath = os.path.join(H5AD_DIR, fname)
        adata = ad.read_h5ad(fpath, backed="r")
        obs = adata.obs[["layer", "qc_pass"]].copy()
        obs = obs[obs["qc_pass"] == True]
        obs["layer"] = obs["layer"].astype(str)
        total = len(obs)
        layer_counts = obs["layer"].value_counts()
        row = {"sample_id": sample_id, "total": total}
        for layer in LAYER_LEGEND_ORDER:
            row[layer] = layer_counts.get(layer, 0)
        records.append(row)

    df = pd.DataFrame(records)
    df["diagnosis"] = df["sample_id"].map(SAMPLE_TO_DX)
    df["wm_prop"] = df["WM"] / df["total"]
    df = df.sort_values("wm_prop").reset_index(drop=True)

    # --- Pick a typical sample (median WM) for comparison ---
    median_idx = len(df) // 2
    typical_sample = df.iloc[median_idx]["sample_id"]
    typical_wm = df.iloc[median_idx]["wm_prop"]
    print(f"Typical sample: {typical_sample} (WM = {typical_wm:.1%})")
    print(f"Outlier sample: Br2039 (WM = {df[df['sample_id']=='Br2039']['wm_prop'].values[0]:.1%})")

    # --- Create figure ---
    fig = plt.figure(figsize=(20, 16), facecolor=BG)
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1], hspace=0.25, wspace=0.15,
                           left=0.06, right=0.96, top=0.95, bottom=0.06)

    # Top row: stacked bar (spans both columns)
    ax_bar = fig.add_subplot(gs[0, :])
    ax_bar.set_facecolor(BG)

    x_pos = np.arange(len(df))
    bar_width = 0.75
    bottoms = np.zeros(len(df))

    # Normalize each sample to proportions (sum to 1.0)
    for layer in LAYER_STACK_ORDER:
        props = df[layer].values.astype(float) / df["total"].values.astype(float)
        ax_bar.bar(x_pos, props, bar_width, bottom=bottoms,
                   color=LAYER_COLORS[layer], label=layer,
                   edgecolor="#222222", linewidth=0.5)
        bottoms += props

    # Highlight Br2039 bar
    outlier_idx = df[df["sample_id"] == "Br2039"].index[0]
    ax_bar.annotate("", xy=(outlier_idx, 1.02),
                    xytext=(outlier_idx, 1.10),
                    arrowprops=dict(arrowstyle="->", color="#ff6666", lw=2.5),
                    annotation_clip=False)
    ax_bar.text(outlier_idx, 1.12, "54% WM",
                ha="center", va="bottom", fontsize=13, color="#ff6666",
                fontweight="bold", clip_on=False)

    ax_bar.set_xticks(x_pos)

    ax_bar.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, p: f"{v*100:.0f}%"))
    ax_bar.set_ylim(0, 1.0)
    ax_bar.set_ylabel("Proportion of cells", fontsize=18, color="white")
    ax_bar.set_title("Layer composition across Xenium samples",
                     fontsize=22, fontweight="bold", color="white", pad=12)
    ax_bar.tick_params(axis="y", colors="white", labelsize=13)
    ax_bar.tick_params(axis="x", labelsize=14)
    for spine in ax_bar.spines.values():
        spine.set_color("#555555")
    ax_bar.grid(axis="y", alpha=0.2, color="#555555")

    # Color x-tick labels by diagnosis using manual text objects
    # (matplotlib's set_xticklabels colors get reset by bbox_inches="tight")
    ax_bar.set_xticklabels([])  # hide default labels
    for i, (_, row) in enumerate(df.iterrows()):
        ax_bar.text(i, -0.02, row["sample_id"],
                    transform=ax_bar.get_xaxis_transform(),
                    rotation=45, ha="right", va="top",
                    fontsize=14, fontweight="bold",
                    color=DX_COLORS[row["diagnosis"]])

    # Layer legend inside bar plot (L1 first in legend)
    layer_handles = [
        Patch(facecolor=LAYER_COLORS[l], edgecolor="#222222", label=l)
        for l in LAYER_LEGEND_ORDER
    ]
    layer_legend = ax_bar.legend(handles=layer_handles, loc="upper left",
                                 fontsize=12, frameon=True,
                                 fancybox=True, framealpha=0.85,
                                 edgecolor="#555555", title="Layer", title_fontsize=13,
                                 ncol=1)
    layer_legend.get_frame().set_facecolor("#222222")
    for text in layer_legend.get_texts():
        text.set_color("white")
    layer_legend.get_title().set_color("white")

    # Diagnosis legend
    dx_patches = [
        Patch(facecolor=DX_COLORS["Control"], label="Control"),
        Patch(facecolor=DX_COLORS["SCZ"], label="SCZ"),
    ]
    dx_legend = ax_bar.legend(handles=dx_patches, loc="upper right",
                              fontsize=13, frameon=True, fancybox=True,
                              framealpha=0.85, edgecolor="#555555",
                              title="Diagnosis (x-axis)", title_fontsize=13)
    dx_legend.get_frame().set_facecolor("#222222")
    for text in dx_legend.get_texts():
        text.set_color("white")
    dx_legend.get_title().set_color("white")
    ax_bar.add_artist(layer_legend)  # re-add layer legend (overwritten by second legend call)

    # --- Bottom row: spatial plots ---
    print(f"\nLoading {typical_sample} for spatial plot...")
    obs_typical = load_sample(typical_sample)

    print(f"Loading Br2039 for spatial plot...")
    obs_outlier = load_sample("Br2039")

    ax_typical = fig.add_subplot(gs[1, 0])
    dx_typical = SAMPLE_TO_DX[typical_sample]
    plot_spatial(ax_typical, obs_typical,
                 f"{typical_sample} ({dx_typical}) -- typical section", s=1.5)

    ax_outlier = fig.add_subplot(gs[1, 1])
    plot_spatial(ax_outlier, obs_outlier,
                 "Br2039 (SCZ) -- WM outlier", s=1.5)

    # Save
    outpath = os.path.join(OUT_DIR, "slide_layer_bar_with_outlier.png")
    plt.savefig(outpath, dpi=200, facecolor=BG, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {outpath}")


if __name__ == "__main__":
    main()
