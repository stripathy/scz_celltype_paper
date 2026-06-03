#!/usr/bin/env python3
"""
Unified layer bar plot script.

Replaces:
  - plot_layer_stacked_bar.py      (stacked bar of raw counts)
  - plot_layer_bar_with_outlier.py (proportion bar + outlier spatial panels)

Usage:
  # Replicate plot_layer_stacked_bar.py output (raw counts, no outlier):
  python plot_layer_bars.py --mode counts

  # Replicate plot_layer_bar_with_outlier.py output (proportions + outlier panels):
  python plot_layer_bars.py --mode proportions --outlier

  # Customize output path or samples:
  python plot_layer_bars.py --mode counts --output my_plot.png
  python plot_layer_bars.py --mode proportions --outlier --samples Br6389 Br2039 Br8433

CLI arguments:
  --mode          {counts, proportions}  Y-axis: raw cell counts or proportions (default: counts)
  --outlier       If set, add bottom-row spatial panels for WM outlier vs typical sample
  --outlier-id    Sample ID to annotate as outlier (default: Br2039)
  --samples       Restrict to these sample IDs (default: all in H5AD_DIR)
  --output        Output PNG path (default: auto-named under PRESENTATION_DIR)
  --save-csv      Also save layer counts CSV alongside the plot
  --dpi           DPI for saved figure (default: 200)
  --color-scheme  {config, hex}  Layer color palette (default: config)
                  'config' uses LAYER_COLORS from config.py (RGB tuples),
                  'hex' uses the standalone hex palette from plot_layer_stacked_bar.py
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    SAMPLE_TO_DX, DX_COLORS, BG_COLOR, H5AD_DIR, PRESENTATION_DIR,
    LAYER_COLORS as CONFIG_LAYER_COLORS,
)

# Standalone hex palette (from the former plot_layer_stacked_bar.py)
HEX_LAYER_COLORS = {
    "L1":       "#e6194b",
    "L2/3":     "#f58231",
    "L4":       "#ffe119",
    "L5":       "#3cb44b",
    "L6":       "#4363d8",
    "WM":       "#911eb4",
    "Vascular": "#a9a9a9",
}

# Bottom-to-top stacking order: WM at bottom, L1 at top
LAYER_STACK_ORDER = ["Vascular", "WM", "L6", "L5", "L4", "L2/3", "L1"]
# Legend / display order: L1 first
LAYER_LEGEND_ORDER = ["L1", "L2/3", "L4", "L5", "L6", "WM", "Vascular"]

BG = BG_COLOR


# ──────────────────────────────────────────────────────────────────────
# Data loading
# ──────────────────────────────────────────────────────────────────────

def collect_layer_counts(sample_ids=None):
    """Load per-sample layer counts from annotated h5ad files.

    Parameters
    ----------
    sample_ids : list of str, optional
        Restrict to these samples. If None, use all *_annotated.h5ad in H5AD_DIR.

    Returns
    -------
    DataFrame with columns: sample_id, total, L1, L2/3, ..., WM, Vascular,
    diagnosis, wm_prop.  Sorted by ascending WM proportion.
    """
    if sample_ids is None:
        h5ad_files = sorted(
            f for f in os.listdir(H5AD_DIR) if f.endswith("_annotated.h5ad")
        )
        sample_ids = [f.replace("_annotated.h5ad", "") for f in h5ad_files]

    records = []
    for sid in sample_ids:
        fpath = os.path.join(H5AD_DIR, f"{sid}_annotated.h5ad")
        if not os.path.exists(fpath):
            print(f"  WARNING: {fpath} not found, skipping")
            continue
        print(f"  Loading {sid}...")
        adata = ad.read_h5ad(fpath, backed="r")
        obs = adata.obs[["layer", "qc_pass"]].copy()
        obs = obs[obs["qc_pass"] == True]
        obs["layer"] = obs["layer"].astype(str)
        total = len(obs)
        layer_counts = obs["layer"].value_counts()
        row = {"sample_id": sid, "total": total}
        for layer in LAYER_LEGEND_ORDER:
            row[layer] = layer_counts.get(layer, 0)
        records.append(row)

    df = pd.DataFrame(records)
    df["diagnosis"] = df["sample_id"].map(SAMPLE_TO_DX)
    df["wm_prop"] = df["WM"] / df["total"]
    df = df.sort_values("wm_prop").reset_index(drop=True)
    return df


def load_sample_obs(sample_id):
    """Load obs + spatial coords for one sample (for spatial panels)."""
    fpath = os.path.join(H5AD_DIR, f"{sample_id}_annotated.h5ad")
    adata = ad.read_h5ad(fpath)
    obs = adata.obs[["layer", "qc_pass"]].copy()
    coords = adata.obsm["spatial"]
    obs["x"] = coords[:, 0]
    obs["y"] = coords[:, 1]
    obs = obs[obs["qc_pass"] == True]
    obs["layer"] = obs["layer"].astype(str)
    return obs


# ──────────────────────────────────────────────────────────────────────
# Plotting helpers
# ──────────────────────────────────────────────────────────────────────

def _style_bar_axis(ax, layer_colors):
    """Common dark-theme styling for the bar axis."""
    ax.set_facecolor(BG)
    ax.tick_params(axis="y", colors="white", labelsize=13)
    ax.tick_params(axis="x", labelsize=14)
    for spine in ax.spines.values():
        spine.set_color("#555555")
    ax.grid(axis="y", alpha=0.2, color="#555555")


def _add_layer_legend(ax, layer_colors):
    """Add layer legend (upper-left)."""
    handles = [
        Patch(facecolor=layer_colors[l], edgecolor="#222222", label=l)
        for l in LAYER_LEGEND_ORDER
    ]
    leg = ax.legend(
        handles=handles, loc="upper left", fontsize=12, frameon=True,
        fancybox=True, framealpha=0.85, edgecolor="#555555",
        title="Layer", title_fontsize=13, ncol=1,
    )
    leg.get_frame().set_facecolor("#222222")
    for text in leg.get_texts():
        text.set_color("white")
    leg.get_title().set_color("white")
    return leg


def _add_dx_legend(ax_or_fig, layer_legend, use_fig=False):
    """Add diagnosis legend (upper-right).

    If use_fig is True, place legend on the figure rather than the axes
    (needed for the raw-counts layout to avoid clipping).
    """
    dx_patches = [
        Patch(facecolor=DX_COLORS["Control"], label="Control"),
        Patch(facecolor=DX_COLORS["SCZ"], label="SCZ"),
    ]
    if use_fig:
        dx_leg = ax_or_fig.legend(
            handles=dx_patches, loc="upper right", fontsize=14,
            frameon=True, fancybox=True, framealpha=0.85,
            edgecolor="#555555", title="Diagnosis", title_fontsize=15,
            bbox_to_anchor=(0.98, 0.95),
        )
    else:
        dx_leg = ax_or_fig.legend(
            handles=dx_patches, loc="upper right", fontsize=13,
            frameon=True, fancybox=True, framealpha=0.85,
            edgecolor="#555555", title="Diagnosis (x-axis)",
            title_fontsize=13,
        )
    dx_leg.get_frame().set_facecolor("#222222")
    for text in dx_leg.get_texts():
        text.set_color("white")
    dx_leg.get_title().set_color("white")

    # Re-add layer legend if it was placed on the same axes (second legend call
    # replaces the first)
    if not use_fig and layer_legend is not None:
        ax_or_fig.add_artist(layer_legend)

    return dx_leg


def plot_spatial(ax, obs, title, layer_colors, s=1.5):
    """Plot cells colored by layer on a spatial axis."""
    x = obs["x"].values
    y = obs["y"].values
    layers = obs["layer"].values

    for lname in LAYER_LEGEND_ORDER:
        mask = layers == lname
        if mask.sum() > 0:
            ax.scatter(
                x[mask], y[mask], s=s, color=[layer_colors[lname]],
                alpha=0.7, rasterized=True, linewidths=0, zorder=3,
                label=lname,
            )

    n_wm = (layers == "WM").sum()
    n_total = len(layers)
    wm_pct = n_wm / n_total * 100
    ax.text(
        0.03, 0.97,
        f"{n_total:,} cells\nWM: {n_wm:,} ({wm_pct:.0f}%)",
        transform=ax.transAxes, ha="left", va="top", fontsize=14,
        color="#dddddd",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#333333",
                  edgecolor="#555555", alpha=0.85),
    )

    ax.set_title(title, fontsize=20, fontweight="bold", color="white", pad=8)
    ax.set_facecolor(BG)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


# ──────────────────────────────────────────────────────────────────────
# Main plot functions
# ──────────────────────────────────────────────────────────────────────

def plot_counts(df, layer_colors, outpath, dpi, save_csv):
    """Raw-count stacked bar (replicates plot_layer_stacked_bar.py)."""
    fig, ax = plt.subplots(figsize=(18, 8), facecolor=BG)
    ax.set_facecolor(BG)

    x_pos = np.arange(len(df))
    bar_width = 0.75
    bottoms = np.zeros(len(df))

    # Stack in L1..Vascular order (same as LAYER_LEGEND_ORDER used in original)
    for layer in LAYER_LEGEND_ORDER:
        counts = df[layer].values.astype(float)
        ax.bar(
            x_pos, counts, bar_width, bottom=bottoms,
            color=layer_colors[layer], label=layer,
            edgecolor="#222222", linewidth=0.5,
        )
        bottoms += counts

    # X-axis: sample IDs colored by diagnosis
    ax.set_xticks(x_pos)
    ax.set_xticklabels(df["sample_id"], rotation=45, ha="right", fontsize=14)
    for i, (_, row) in enumerate(df.iterrows()):
        color = DX_COLORS[row["diagnosis"]]
        ax.get_xticklabels()[i].set_color(color)
        ax.get_xticklabels()[i].set_fontweight("bold")

    ax.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda v, p: f"{int(v):,}")
    )
    ax.set_ylabel("Number of cells (QC-pass)", fontsize=18, color="white")
    ax.set_xlabel(
        "Sample (sorted by ascending WM proportion)", fontsize=16, color="white"
    )
    ax.set_title(
        "Layer composition across Xenium samples",
        fontsize=22, fontweight="bold", color="white", pad=15,
    )
    _style_bar_axis(ax, layer_colors)

    # Legends
    layer_leg = _add_layer_legend(ax, layer_colors)
    _add_dx_legend(fig, layer_leg, use_fig=True)

    plt.tight_layout(pad=2.0)
    plt.savefig(outpath, dpi=dpi, facecolor=BG, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {outpath}")

    if save_csv:
        csv_path = outpath.replace(".png", ".csv")
        df.to_csv(csv_path, index=False)
        print(f"Saved: {csv_path}")


def plot_proportions(df, layer_colors, outpath, dpi, save_csv,
                     outlier=False, outlier_id="Br2039"):
    """Proportion stacked bar, optionally with outlier spatial panels.

    Replicates plot_layer_bar_with_outlier.py when outlier=True.
    """
    if outlier:
        fig = plt.figure(figsize=(20, 16), facecolor=BG)
        gs = gridspec.GridSpec(
            2, 2, height_ratios=[1, 1], hspace=0.25, wspace=0.15,
            left=0.06, right=0.96, top=0.95, bottom=0.06,
        )
        ax_bar = fig.add_subplot(gs[0, :])
    else:
        fig, ax_bar = plt.subplots(figsize=(18, 8), facecolor=BG)

    ax_bar.set_facecolor(BG)

    x_pos = np.arange(len(df))
    bar_width = 0.75
    bottoms = np.zeros(len(df))

    for layer in LAYER_STACK_ORDER:
        props = df[layer].values.astype(float) / df["total"].values.astype(float)
        ax_bar.bar(
            x_pos, props, bar_width, bottom=bottoms,
            color=layer_colors[layer], label=layer,
            edgecolor="#222222", linewidth=0.5,
        )
        bottoms += props

    # Outlier annotation arrow
    if outlier and outlier_id in df["sample_id"].values:
        outlier_idx = df[df["sample_id"] == outlier_id].index[0]
        outlier_wm_pct = df.loc[outlier_idx, "wm_prop"] * 100
        ax_bar.annotate(
            "", xy=(outlier_idx, 1.02), xytext=(outlier_idx, 1.10),
            arrowprops=dict(arrowstyle="->", color="#ff6666", lw=2.5),
            annotation_clip=False,
        )
        ax_bar.text(
            outlier_idx, 1.12, f"{outlier_wm_pct:.0f}% WM",
            ha="center", va="bottom", fontsize=13, color="#ff6666",
            fontweight="bold", clip_on=False,
        )

    ax_bar.set_xticks(x_pos)
    ax_bar.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda v, p: f"{v*100:.0f}%")
    )
    ax_bar.set_ylim(0, 1.0)
    ax_bar.set_ylabel("Proportion of cells", fontsize=18, color="white")
    ax_bar.set_title(
        "Layer composition across Xenium samples",
        fontsize=22, fontweight="bold", color="white", pad=12,
    )
    _style_bar_axis(ax_bar, layer_colors)

    # X-tick labels colored by diagnosis (manual text objects to survive
    # bbox_inches="tight")
    ax_bar.set_xticklabels([])
    for i, (_, row) in enumerate(df.iterrows()):
        ax_bar.text(
            i, -0.02, row["sample_id"],
            transform=ax_bar.get_xaxis_transform(),
            rotation=45, ha="right", va="top", fontsize=14,
            fontweight="bold", color=DX_COLORS[row["diagnosis"]],
        )

    # Legends
    layer_leg = _add_layer_legend(ax_bar, layer_colors)
    _add_dx_legend(ax_bar, layer_leg, use_fig=False)

    # --- Optional spatial panels ---
    if outlier:
        # Pick typical sample (median WM proportion)
        median_idx = len(df) // 2
        typical_sample = df.iloc[median_idx]["sample_id"]
        typical_wm = df.iloc[median_idx]["wm_prop"]
        print(f"Typical sample: {typical_sample} (WM = {typical_wm:.1%})")
        if outlier_id in df["sample_id"].values:
            print(
                f"Outlier sample: {outlier_id} "
                f"(WM = {df[df['sample_id']==outlier_id]['wm_prop'].values[0]:.1%})"
            )

        print(f"\nLoading {typical_sample} for spatial plot...")
        obs_typical = load_sample_obs(typical_sample)
        dx_typical = SAMPLE_TO_DX.get(typical_sample, "?")

        ax_typical = fig.add_subplot(gs[1, 0])
        plot_spatial(
            ax_typical, obs_typical,
            f"{typical_sample} ({dx_typical}) -- typical section",
            layer_colors,
        )

        if outlier_id in df["sample_id"].values:
            print(f"Loading {outlier_id} for spatial plot...")
            obs_outlier = load_sample_obs(outlier_id)
            ax_outlier = fig.add_subplot(gs[1, 1])
            dx_outlier = SAMPLE_TO_DX.get(outlier_id, "?")
            plot_spatial(
                ax_outlier, obs_outlier,
                f"{outlier_id} ({dx_outlier}) -- WM outlier",
                layer_colors,
            )
    else:
        plt.tight_layout(pad=2.0)

    plt.savefig(outpath, dpi=dpi, facecolor=BG, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {outpath}")

    if save_csv:
        csv_path = outpath.replace(".png", ".csv")
        df.to_csv(csv_path, index=False)
        print(f"Saved: {csv_path}")


# ──────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="Unified layer bar plot for Xenium samples.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--mode", choices=["counts", "proportions"], default="counts",
        help="Y-axis: raw cell counts or proportions (default: counts)",
    )
    parser.add_argument(
        "--outlier", action="store_true",
        help="Add bottom-row spatial panels for WM outlier vs typical sample "
             "(only applies to proportions mode)",
    )
    parser.add_argument(
        "--outlier-id", default="Br2039",
        help="Sample ID to annotate as WM outlier (default: Br2039)",
    )
    parser.add_argument(
        "--samples", nargs="+", default=None,
        help="Restrict to these sample IDs (default: all in H5AD_DIR)",
    )
    parser.add_argument(
        "--output", default=None,
        help="Output PNG path (default: auto-named under PRESENTATION_DIR)",
    )
    parser.add_argument(
        "--save-csv", action="store_true",
        help="Also save layer counts CSV alongside the plot",
    )
    parser.add_argument(
        "--dpi", type=int, default=200,
        help="DPI for saved figure (default: 200)",
    )
    parser.add_argument(
        "--color-scheme", choices=["config", "hex"], default="config",
        help="Layer color palette: 'config' (RGB tuples from config.py) or "
             "'hex' (standalone hex palette). Default: config",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Select color palette
    if args.color_scheme == "hex":
        layer_colors = HEX_LAYER_COLORS
    else:
        layer_colors = CONFIG_LAYER_COLORS

    # Default output paths
    if args.output is None:
        if args.mode == "counts":
            args.output = os.path.join(
                PRESENTATION_DIR, "slide_layer_stacked_bar.png"
            )
        elif args.outlier:
            args.output = os.path.join(
                PRESENTATION_DIR, "slide_layer_bar_with_outlier.png"
            )
        else:
            args.output = os.path.join(
                PRESENTATION_DIR, "slide_layer_proportions.png"
            )

    # Always save CSV in counts mode (matches original behavior)
    if args.mode == "counts":
        args.save_csv = True

    print("Loading layer counts...")
    df = collect_layer_counts(sample_ids=args.samples)

    print(f"\n{'Sample':<10} {'Dx':<10} {'Total':>8} {'WM%':>8}")
    print("-" * 40)
    for _, row in df.iterrows():
        print(
            f"{row['sample_id']:<10} {row['diagnosis']:<10} "
            f"{row['total']:>8,} {row['wm_prop']:>7.1%}"
        )

    if args.mode == "counts":
        plot_counts(df, layer_colors, args.output, args.dpi, args.save_csv)
    else:
        plot_proportions(
            df, layer_colors, args.output, args.dpi, args.save_csv,
            outlier=args.outlier, outlier_id=args.outlier_id,
        )

    print("Done.")


if __name__ == "__main__":
    main()
