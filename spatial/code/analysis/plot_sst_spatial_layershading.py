#!/usr/bin/env python3
"""
Vulnerable Sst spatial plot with very transparent layer shading for cortical context.

Each cortical layer is shown as a faint colored background, giving a sense of
tissue shape and laminar organization without obscuring the highlighted Sst cells.

Output: output/presentation/slide_sst_spatial_layershading.png
"""

import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    PRESENTATION_DIR, BG_COLOR, DX_COLORS,
    SST_TYPES, SST_COLORS, REPRESENTATIVE_SAMPLES, CORTICAL_LAYERS,
    LAYER_COLORS, LAYER_ORDER,
    load_all_cells, draw_layer_shading, style_dark_axis,
)

MARKER_SIZE_SST = 18
LAYER_ALPHA = 0.12


def main():
    fig = plt.figure(figsize=(16, 14), facecolor=BG_COLOR)
    gs = gridspec.GridSpec(2, 2, figure=fig,
                           wspace=0.05, hspace=0.22,
                           left=0.02, right=0.98, top=0.89, bottom=0.08)

    for idx, (sample_id, diagnosis) in enumerate(REPRESENTATIVE_SAMPLES):
        row, col = divmod(idx, 2)
        ax = fig.add_subplot(gs[row, col])
        ax.set_facecolor(BG_COLOR)

        print(f"Loading {sample_id} ({diagnosis})...")
        df = load_all_cells(sample_id)

        # All cells with layers for shading (including WM for context)
        shading_layers = {"L1", "L2/3", "L4", "L5", "L6", "WM"}
        shading_cells = df[df["layer"].isin(shading_layers)]

        # Cortical cells only (L1-L6) for counting
        cortical = df[df["layer"].isin(CORTICAL_LAYERS)]
        total_cortical = len(cortical)

        # Draw layer shading first (background)
        draw_layer_shading(ax, shading_cells, alpha=LAYER_ALPHA)

        # Overlay vulnerable Sst cells (from cortical L1-L6 only)
        sst_counts = {}
        total_sst_vuln = 0
        for sst_type in SST_TYPES:
            mask = cortical["supertype_label"] == sst_type
            sub = cortical[mask]
            n = len(sub)
            sst_counts[sst_type] = n
            total_sst_vuln += n
            if n > 0:
                ax.scatter(sub["x"], sub["y"], s=MARKER_SIZE_SST,
                           c=SST_COLORS[sst_type], alpha=0.9,
                           linewidths=0.5, edgecolors="white",
                           zorder=5, rasterized=True)

        vuln_pct = total_sst_vuln / total_cortical * 100 if total_cortical > 0 else 0

        # Title
        title_color = DX_COLORS[diagnosis]
        ax.set_title(f"{sample_id} ({diagnosis})", fontsize=16,
                     fontweight="bold", color=title_color, pad=4)

        # Counts annotation
        count_lines = [f"Vulnerable Sst: {total_sst_vuln} ({vuln_pct:.2f}%)"]
        for t in SST_TYPES:
            count_lines.append(f"  {t}: {sst_counts[t]}")
        count_text = "\n".join(count_lines)

        ax.text(0.02, 0.98, count_text, transform=ax.transAxes,
                ha="left", va="top", fontsize=11, color="#dddddd",
                fontfamily="monospace",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="#1a1a1a",
                          edgecolor="#444444", alpha=0.85))

        style_dark_axis(ax)

        print(f"  Vulnerable Sst: {total_sst_vuln} ({vuln_pct:.2f}%)")

    # Row labels
    fig.text(0.50, 0.905, "Control", ha="center", fontsize=20,
             fontweight="bold", color=DX_COLORS["Control"])
    fig.text(0.50, 0.455, "SCZ", ha="center", fontsize=20,
             fontweight="bold", color=DX_COLORS["SCZ"])

    # Cell type legend (bottom left)
    sst_legend = [
        Line2D([0], [0], marker='o', color=BG_COLOR,
               markerfacecolor=SST_COLORS[t],
               markeredgecolor="white", markeredgewidth=0.5,
               markersize=12, label=t, linewidth=0)
        for t in SST_TYPES
    ]
    fig.legend(handles=sst_legend, loc="lower left", ncol=4,
               fontsize=13, frameon=False, labelcolor="white",
               bbox_to_anchor=(0.02, 0.005))

    # Layer legend (bottom right)
    layer_legend = [
        Line2D([0], [0], marker='s', color=BG_COLOR,
               markerfacecolor=LAYER_COLORS[layer],
               markersize=10, label=layer, linewidth=0, alpha=0.6)
        for layer in LAYER_ORDER
    ]
    fig.legend(handles=layer_legend, loc="lower right", ncol=len(LAYER_ORDER),
               fontsize=12, frameon=False, labelcolor="white",
               bbox_to_anchor=(0.98, 0.005))

    fig.suptitle(
        "Vulnerable Sst subtypes in SCZ vs Control (Xenium)",
        fontsize=22, fontweight="bold", color="white", y=0.96
    )

    outpath = os.path.join(PRESENTATION_DIR, "slide_sst_spatial_layershading.png")
    plt.savefig(outpath, dpi=150, facecolor=BG_COLOR, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {outpath}")


if __name__ == "__main__":
    main()
