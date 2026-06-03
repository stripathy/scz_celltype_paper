#!/usr/bin/env python3
"""
L6b spatial plot with full-tissue inset showing where the cortical region sits.

Each main panel shows L6b cells highlighted on cortical background.
A small inset in the corner shows ALL cells in the full tissue, with cortical
cells highlighted and extra-cortical/WM cells dimmed, so you can see the
section geometry and where cortex is.

Output: output/presentation/slide_l6b_spatial_inset.png
"""

import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    PRESENTATION_DIR, BG_COLOR, DX_COLORS, ALL_CELL_COLOR, MARKER_SIZE_BG,
    L6B_TYPES, L6B_COLORS, REPRESENTATIVE_SAMPLES, CORTICAL_LAYERS,
    load_all_cells, draw_inset, style_dark_axis,
)

MARKER_SIZE_L6B = 14


def main():
    fig = plt.figure(figsize=(16, 14), facecolor=BG_COLOR)
    gs = gridspec.GridSpec(2, 2, figure=fig,
                           wspace=0.05, hspace=0.22,
                           left=0.02, right=0.98, top=0.89, bottom=0.06)

    for idx, (sample_id, diagnosis) in enumerate(REPRESENTATIVE_SAMPLES):
        row, col = divmod(idx, 2)
        ax = fig.add_subplot(gs[row, col])
        ax.set_facecolor(BG_COLOR)

        print(f"Loading {sample_id} ({diagnosis})...")
        df = load_all_cells(sample_id)

        cortical = df[df["layer"].isin(CORTICAL_LAYERS)]
        total_cortical = len(cortical)

        # Background: all cortical cells
        ax.scatter(cortical["x"], cortical["y"], s=MARKER_SIZE_BG,
                   c=ALL_CELL_COLOR, alpha=0.25, rasterized=True, linewidths=0)

        # Overlay L6b cells
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

        # Counts annotation (top left)
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
        draw_inset(ax, df, loc="lower right")

        print(f"  L6b: {total_l6b} ({l6b_pct:.2f}%)")

    # Row labels
    fig.text(0.50, 0.905, "Control", ha="center", fontsize=20,
             fontweight="bold", color=DX_COLORS["Control"])
    fig.text(0.50, 0.455, "SCZ", ha="center", fontsize=20,
             fontweight="bold", color=DX_COLORS["SCZ"])

    # Legend at bottom
    legend_elements = [
        Line2D([0], [0], marker='o', color=BG_COLOR,
               markerfacecolor=L6B_COLORS[t],
               markeredgecolor="white", markeredgewidth=0.5,
               markersize=12, label=t, linewidth=0)
        for t in L6B_TYPES
    ]
    legend_elements.append(
        Line2D([0], [0], marker='o', color=BG_COLOR,
               markerfacecolor=ALL_CELL_COLOR,
               markersize=6, label="All cortical cells", linewidth=0, alpha=0.5)
    )
    fig.legend(handles=legend_elements, loc="lower center", ncol=4,
               fontsize=13, frameon=False, labelcolor="white",
               bbox_to_anchor=(0.5, 0.005))

    fig.suptitle(
        "L6b subtypes in SCZ vs Control (Xenium)",
        fontsize=22, fontweight="bold", color="white", y=0.96
    )

    outpath = os.path.join(PRESENTATION_DIR, "slide_l6b_spatial_inset.png")
    plt.savefig(outpath, dpi=150, facecolor=BG_COLOR, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {outpath}")


if __name__ == "__main__":
    main()
