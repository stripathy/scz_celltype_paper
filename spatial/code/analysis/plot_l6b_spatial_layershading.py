#!/usr/bin/env python3
"""
L6b spatial plot with very transparent layer shading to provide cortical context.

Each cortical layer is shown as a faint colored background using a 2D histogram
(hexbin-style density), giving a sense of tissue shape and laminar organization
without obscuring the highlighted L6b cells.

Output: output/presentation/slide_l6b_spatial_layershading.png
"""

import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    PRESENTATION_DIR, BG_COLOR, DX_COLORS,
    L6B_TYPES, L6B_COLORS, REPRESENTATIVE_SAMPLES, CORTICAL_LAYERS,
    LAYER_COLORS, LAYER_ORDER,
    load_all_cells, draw_layer_shading, style_dark_axis,
)

MARKER_SIZE_L6B = 14
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

        # Cortical cells only (L1-L6) for counting and L6b overlay
        cortical = df[df["layer"].isin(CORTICAL_LAYERS)]
        total_cortical = len(cortical)

        # Draw layer shading first (background)
        draw_layer_shading(ax, shading_cells, alpha=LAYER_ALPHA)

        # Overlay L6b cells (from cortical L1-L6 only)
        l6b_counts = {}
        total_l6b = 0
        for l6b_type in L6B_TYPES:
            mask = cortical["supertype_label"] == l6b_type
            sub = cortical[mask]
            n = len(sub)
            l6b_counts[l6b_type] = n
            total_l6b += n
            if n > 0:
                ax.scatter(sub["x"], sub["y"], s=MARKER_SIZE_L6B,
                           c=L6B_COLORS[l6b_type], alpha=0.85,
                           linewidths=0.4, edgecolors="white",
                           zorder=5, rasterized=True)

        l6b_pct = total_l6b / total_cortical * 100 if total_cortical > 0 else 0

        # Title
        title_color = DX_COLORS[diagnosis]
        ax.set_title(f"{sample_id} ({diagnosis})", fontsize=16,
                     fontweight="bold", color=title_color, pad=4)

        # Counts annotation
        count_lines = [f"L6b subtypes: {total_l6b} ({l6b_pct:.2f}%)"]
        for t in L6B_TYPES:
            count_lines.append(f"  {t}: {l6b_counts[t]}")
        count_text = "\n".join(count_lines)

        ax.text(0.02, 0.98, count_text, transform=ax.transAxes,
                ha="left", va="top", fontsize=11, color="#dddddd",
                fontfamily="monospace",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="#1a1a1a",
                          edgecolor="#444444", alpha=0.85))

        style_dark_axis(ax)

        print(f"  L6b: {total_l6b} ({l6b_pct:.2f}%)")

    # Row labels
    fig.text(0.50, 0.905, "Control", ha="center", fontsize=20,
             fontweight="bold", color=DX_COLORS["Control"])
    fig.text(0.50, 0.455, "SCZ", ha="center", fontsize=20,
             fontweight="bold", color=DX_COLORS["SCZ"])

    # Cell type legend (bottom left area)
    l6b_legend = [
        Line2D([0], [0], marker='o', color=BG_COLOR,
               markerfacecolor=L6B_COLORS[t],
               markeredgecolor="white", markeredgewidth=0.5,
               markersize=12, label=t, linewidth=0)
        for t in L6B_TYPES
    ]
    fig.legend(handles=l6b_legend, loc="lower left", ncol=3,
               fontsize=13, frameon=False, labelcolor="white",
               bbox_to_anchor=(0.02, 0.005))

    # Layer legend (bottom right area)
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
        "L6b subtypes in SCZ vs Control (Xenium)",
        fontsize=22, fontweight="bold", color="white", y=0.96
    )

    outpath = os.path.join(PRESENTATION_DIR, "slide_l6b_spatial_layershading.png")
    plt.savefig(outpath, dpi=150, facecolor=BG_COLOR, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {outpath}")


if __name__ == "__main__":
    main()
