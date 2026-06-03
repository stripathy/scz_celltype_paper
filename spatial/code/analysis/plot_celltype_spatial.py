#!/usr/bin/env python3
"""
Spatial visualization of vulnerable cell type subtypes in SCZ vs Control.

Unified script replacing 6 near-identical scripts:
  plot_sst_spatial_representative.py  →  --group sst --background plain
  plot_sst_spatial_inset.py           →  --group sst --background inset
  plot_sst_spatial_layershading.py    →  --group sst --background layers
  plot_l6b_spatial_representative.py  →  --group l6b --background plain
  plot_l6b_spatial_inset.py           →  --group l6b --background inset
  plot_l6b_spatial_layershading.py    →  --group l6b --background layers

Shows median-representative Xenium tissue sections with highlighted cell subtypes
against configurable backgrounds (plain cortical, full-tissue inset, layer shading).

Examples:
  python plot_celltype_spatial.py --group sst --background plain
  python plot_celltype_spatial.py --group l6b --background layers
  python plot_celltype_spatial.py --group sst --background inset
  python plot_celltype_spatial.py --group sst --background inset --types Sst_25 Sst_22 Sst_2
  python plot_celltype_spatial.py --group sst --background inset --samples Br6389:Control Br8433:Control Br6437:SCZ Br6032:SCZ

Output: output/presentation/slide_{group}_spatial_{background}.png
"""

import os
import sys
import argparse
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    PRESENTATION_DIR, BG_COLOR, DX_COLORS, ALL_CELL_COLOR, MARKER_SIZE_BG,
    SST_TYPES, SST_COLORS, L6B_TYPES, L6B_COLORS,
    REPRESENTATIVE_SAMPLES, CORTICAL_LAYERS,
    LAYER_COLORS, LAYER_ORDER,
    load_cortical, load_all_cells, draw_inset, draw_layer_shading,
    style_dark_axis,
)

# ──────────────────────────────────────────────────────────────────────
# Group-specific defaults
# ──────────────────────────────────────────────────────────────────────

GROUP_DEFAULTS = {
    "sst": {
        "types": SST_TYPES,
        "colors": SST_COLORS,
        "marker_size": 18,
        "alpha": 0.9,
        "edge_width": 0.5,
        "label": "Vulnerable Sst",
        "title_prefix": "Vulnerable Sst subtypes",
    },
    "l6b": {
        "types": L6B_TYPES,
        "colors": L6B_COLORS,
        "marker_size": 14,
        "alpha": 0.85,
        "edge_width": 0.4,
        "label": "L6b subtypes",
        "title_prefix": "L6b subtypes",
    },
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Spatial plot of vulnerable cell types in SCZ vs Control")
    parser.add_argument("--group", choices=list(GROUP_DEFAULTS.keys()),
                        required=True,
                        help="Cell type group to highlight (sst or l6b)")
    parser.add_argument("--background", choices=["plain", "inset", "layers"],
                        default="plain",
                        help="Background mode: plain (cortical scatter), "
                             "inset (full-tissue inset), layers (layer shading)")
    parser.add_argument("--types", nargs="+", default=None,
                        help="Override which subtypes to plot (e.g., Sst_25 Sst_22 Sst_2)")
    parser.add_argument("--samples", nargs="+", default=None,
                        help="Override representative samples as ID:Dx pairs "
                             "(e.g., Br6389:Control Br8433:Control Br6437:SCZ Br6032:SCZ)")
    parser.add_argument("--output", default=None,
                        help="Override output filename (default: auto-generated)")
    return parser.parse_args()


def parse_sample_list(sample_args):
    """Parse 'Br6389:Control Br8433:Control' into list of (id, dx) tuples."""
    samples = []
    for s in sample_args:
        parts = s.split(":")
        if len(parts) != 2:
            raise ValueError(f"Invalid sample format '{s}', expected 'SampleID:Diagnosis'")
        samples.append((parts[0], parts[1]))
    return samples


def main():
    args = parse_args()

    # Resolve group settings
    grp = GROUP_DEFAULTS[args.group]
    cell_types = args.types if args.types else grp["types"]
    colors = grp["colors"]
    marker_size = grp["marker_size"]
    alpha = grp["alpha"]
    edge_width = grp["edge_width"]
    group_label = grp["label"]
    title_prefix = grp["title_prefix"]

    # Validate that all requested types have colors defined
    for t in cell_types:
        if t not in colors:
            print(f"Warning: no color defined for '{t}', using grey")
            colors[t] = "#999999"

    # Resolve sample list
    samples = (parse_sample_list(args.samples) if args.samples
               else REPRESENTATIVE_SAMPLES)

    # Determine if we need full tissue data (inset/layers) or cortical only
    needs_all_cells = args.background in ("inset", "layers")

    # Output filename
    bg_suffix = {"plain": "representative", "inset": "inset",
                 "layers": "layershading"}[args.background]
    outname = args.output or f"slide_{args.group}_spatial_{bg_suffix}.png"
    outpath = os.path.join(PRESENTATION_DIR, outname)

    # ── Build figure ─────────────────────────────────────────────────
    bottom_margin = 0.08 if args.background == "layers" else 0.06
    fig = plt.figure(figsize=(16, 14), facecolor=BG_COLOR)
    gs = gridspec.GridSpec(2, 2, figure=fig,
                           wspace=0.05, hspace=0.22,
                           left=0.02, right=0.98, top=0.89,
                           bottom=bottom_margin)

    for idx, (sample_id, diagnosis) in enumerate(samples):
        row, col = divmod(idx, 2)
        ax = fig.add_subplot(gs[row, col])
        ax.set_facecolor(BG_COLOR)

        print(f"Loading {sample_id} ({diagnosis})...")

        if needs_all_cells:
            df = load_all_cells(sample_id)
            cortical = df[df["layer"].isin(CORTICAL_LAYERS)]
        else:
            df = load_cortical(sample_id)
            cortical = df

        total_cortical = len(cortical)

        # ── Draw background ──────────────────────────────────────────
        if args.background == "plain":
            ax.scatter(cortical["x"], cortical["y"], s=MARKER_SIZE_BG,
                       c=ALL_CELL_COLOR, alpha=0.25, rasterized=True,
                       linewidths=0)

        elif args.background == "inset":
            ax.scatter(cortical["x"], cortical["y"], s=MARKER_SIZE_BG,
                       c=ALL_CELL_COLOR, alpha=0.25, rasterized=True,
                       linewidths=0)

        elif args.background == "layers":
            shading_layers = {"L1", "L2/3", "L4", "L5", "L6", "WM"}
            shading_cells = df[df["layer"].isin(shading_layers)]
            draw_layer_shading(ax, shading_cells, alpha=0.12)

        # ── Overlay highlighted cell types ────────────────────────────
        type_counts = {}
        total_highlighted = 0
        for ct in cell_types:
            mask = cortical["supertype_label"] == ct
            sub = cortical[mask]
            n = len(sub)
            type_counts[ct] = n
            total_highlighted += n
            if n > 0:
                ax.scatter(sub["x"], sub["y"], s=marker_size,
                           c=colors[ct], alpha=alpha,
                           linewidths=edge_width, edgecolors="white",
                           zorder=5, rasterized=True)

        pct = total_highlighted / total_cortical * 100 if total_cortical > 0 else 0

        # Title
        title_color = DX_COLORS[diagnosis]
        ax.set_title(f"{sample_id} ({diagnosis})", fontsize=16,
                     fontweight="bold", color=title_color, pad=4)

        # Counts annotation
        count_lines = [f"{group_label}: {total_highlighted} ({pct:.2f}%)"]
        for t in cell_types:
            count_lines.append(f"  {t}: {type_counts[t]}")

        ax.text(0.02, 0.98, "\n".join(count_lines), transform=ax.transAxes,
                ha="left", va="top", fontsize=11, color="#dddddd",
                fontfamily="monospace",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="#1a1a1a",
                          edgecolor="#444444", alpha=0.85))

        style_dark_axis(ax)

        # Draw inset if requested
        if args.background == "inset":
            draw_inset(ax, df, loc="lower right")

        print(f"  {group_label}: {total_highlighted} ({pct:.2f}%)")

    # ── Row labels ───────────────────────────────────────────────────
    fig.text(0.50, 0.905, "Control", ha="center", fontsize=20,
             fontweight="bold", color=DX_COLORS["Control"])
    fig.text(0.50, 0.455, "SCZ", ha="center", fontsize=20,
             fontweight="bold", color=DX_COLORS["SCZ"])

    # ── Legends ──────────────────────────────────────────────────────
    if args.background == "layers":
        # Dual legend: cell types bottom-left, layers bottom-right
        ct_legend = [
            Line2D([0], [0], marker='o', color=BG_COLOR,
                   markerfacecolor=colors[t],
                   markeredgecolor="white", markeredgewidth=0.5,
                   markersize=12, label=t, linewidth=0)
            for t in cell_types
        ]
        fig.legend(handles=ct_legend, loc="lower left",
                   ncol=max(len(cell_types), 3),
                   fontsize=13, frameon=False, labelcolor="white",
                   bbox_to_anchor=(0.02, 0.005))

        layer_legend = [
            Line2D([0], [0], marker='s', color=BG_COLOR,
                   markerfacecolor=LAYER_COLORS[layer],
                   markersize=10, label=layer, linewidth=0, alpha=0.6)
            for layer in LAYER_ORDER
        ]
        fig.legend(handles=layer_legend, loc="lower right",
                   ncol=len(LAYER_ORDER),
                   fontsize=12, frameon=False, labelcolor="white",
                   bbox_to_anchor=(0.98, 0.005))
    else:
        # Single centered legend: cell types + background marker
        legend_elements = [
            Line2D([0], [0], marker='o', color=BG_COLOR,
                   markerfacecolor=colors[t],
                   markeredgecolor="white", markeredgewidth=0.5,
                   markersize=12, label=t, linewidth=0)
            for t in cell_types
        ]
        legend_elements.append(
            Line2D([0], [0], marker='o', color=BG_COLOR,
                   markerfacecolor=ALL_CELL_COLOR,
                   markersize=6, label="All cortical cells",
                   linewidth=0, alpha=0.5)
        )
        fig.legend(handles=legend_elements, loc="lower center",
                   ncol=len(cell_types) + 1,
                   fontsize=13, frameon=False, labelcolor="white",
                   bbox_to_anchor=(0.5, 0.005))

    fig.suptitle(
        f"{title_prefix} in SCZ vs Control (Xenium)",
        fontsize=22, fontweight="bold", color="white", y=0.96
    )

    plt.savefig(outpath, dpi=150, facecolor=BG_COLOR, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {outpath}")


if __name__ == "__main__":
    main()
