#!/usr/bin/env python3
"""
Sst_25 spatial plot with layer shading for cortical context.

Shows only Sst_25 cells in median-representative samples to illustrate
the depletion of this specific supertype in SCZ.

Representative samples are chosen as closest to group median Sst_25
proportion:
  Control median = 0.199%: Br5400 (0.203%, d=0.005), Br5314 (0.194%, d=0.005)
  SCZ median     = 0.083%: Br6496 (0.058%, d=0.025), Br6437 (0.107%, d=0.025)

Output: output/presentation/slide_sst25_spatial_layershading.png
"""

import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    PRESENTATION_DIR, BG_COLOR, DX_COLORS,
    SST_COLORS, CORTICAL_LAYERS,
    LAYER_COLORS, LAYER_ORDER,
    load_all_cells, draw_layer_shading, style_dark_axis,
)

# Sst_25-specific representative samples (closest to group medians)
SST25_REPRESENTATIVE = [
    ("Br5400", "Control"),
    ("Br5314", "Control"),
    ("Br6496", "SCZ"),
    ("Br6437", "SCZ"),
]

SST25_COLOR = SST_COLORS["Sst_25"]
MARKER_SIZE = 22
LAYER_ALPHA = 0.12


def main():
    fig = plt.figure(figsize=(16, 14), facecolor=BG_COLOR)
    gs = gridspec.GridSpec(2, 2, figure=fig,
                           wspace=0.05, hspace=0.22,
                           left=0.02, right=0.98, top=0.89, bottom=0.08)

    for idx, (sample_id, diagnosis) in enumerate(SST25_REPRESENTATIVE):
        row, col = divmod(idx, 2)
        ax = fig.add_subplot(gs[row, col])
        ax.set_facecolor(BG_COLOR)

        print(f"Loading {sample_id} ({diagnosis})...")
        df = load_all_cells(sample_id)

        # All cells with layers for shading
        shading_layers = {"L1", "L2/3", "L4", "L5", "L6", "WM"}
        shading_cells = df[df["layer"].isin(shading_layers)]

        # Cortical cells only (L1-L6)
        cortical = df[df["layer"].isin(CORTICAL_LAYERS)]
        total_cortical = len(cortical)

        # Draw layer shading (background)
        draw_layer_shading(ax, shading_cells, alpha=LAYER_ALPHA)

        # Overlay Sst_25 cells
        mask = cortical["supertype_label"] == "Sst_25"
        sst25 = cortical[mask]
        n_sst25 = len(sst25)

        if n_sst25 > 0:
            ax.scatter(sst25["x"], sst25["y"], s=MARKER_SIZE,
                       c=SST25_COLOR, alpha=0.9,
                       linewidths=0.5, edgecolors="white",
                       zorder=5, rasterized=True)

        pct = n_sst25 / total_cortical * 100 if total_cortical > 0 else 0

        # Title
        title_color = DX_COLORS[diagnosis]
        ax.set_title(f"{sample_id} ({diagnosis})", fontsize=16,
                     fontweight="bold", color=title_color, pad=4)

        # Count annotation
        count_text = f"Sst_25: {n_sst25} ({pct:.3f}%)\nCortical: {total_cortical:,}"
        ax.text(0.02, 0.98, count_text, transform=ax.transAxes,
                ha="left", va="top", fontsize=12, color="#dddddd",
                fontfamily="monospace",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="#1a1a1a",
                          edgecolor="#444444", alpha=0.85))

        style_dark_axis(ax)
        print(f"  Sst_25: {n_sst25} ({pct:.3f}%)")

    # Row labels
    fig.text(0.50, 0.905, "Control", ha="center", fontsize=20,
             fontweight="bold", color=DX_COLORS["Control"])
    fig.text(0.50, 0.455, "SCZ", ha="center", fontsize=20,
             fontweight="bold", color=DX_COLORS["SCZ"])

    # Cell type legend (bottom left)
    sst25_legend = [
        Line2D([0], [0], marker='o', color=BG_COLOR,
               markerfacecolor=SST25_COLOR,
               markeredgecolor="white", markeredgewidth=0.5,
               markersize=12, label="Sst_25", linewidth=0)
    ]
    fig.legend(handles=sst25_legend, loc="lower left", ncol=1,
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
        "Sst_25 in SCZ vs Control (Xenium)",
        fontsize=22, fontweight="bold", color="white", y=0.96
    )

    outpath = os.path.join(PRESENTATION_DIR, "slide_sst25_spatial_layershading.png")
    plt.savefig(outpath, dpi=150, facecolor=BG_COLOR, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {outpath}")


if __name__ == "__main__":
    main()
