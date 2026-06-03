#!/usr/bin/env python3
"""
Vulnerable Sst spatial plot with full-tissue inset showing section context.

Each main panel shows Sst cells highlighted on cortical background.
A small inset in the corner shows ALL cells in the full tissue, with cortical
layers colored so you can see the section geometry.

Output: output/presentation/slide_sst_spatial_inset.png
"""

import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    PRESENTATION_DIR, BG_COLOR, DX_COLORS, ALL_CELL_COLOR, MARKER_SIZE_BG,
    SST_COLORS, REPRESENTATIVE_SAMPLES, CORTICAL_LAYERS,
    load_all_cells, draw_inset, style_dark_axis,
)

# Subset of SST types for this figure
SST_PLOT_TYPES = ["Sst_25", "Sst_22", "Sst_2"]

# Override representative samples for SST figure:
# Br6032 replaces Br2421 — it is the exact median for summed Sst_25+22+2
# proportion among SCZ samples (0.53%) and has all cortical layers represented.
SST_REPRESENTATIVE = [
    ("Br6389", "Control"),
    ("Br8433", "Control"),
    ("Br6437", "SCZ"),
    ("Br6032", "SCZ"),
]

MARKER_SIZE_SST = 18


def main():
    fig = plt.figure(figsize=(16, 14), facecolor=BG_COLOR)
    gs = gridspec.GridSpec(2, 2, figure=fig,
                           wspace=0.05, hspace=0.22,
                           left=0.02, right=0.98, top=0.89, bottom=0.06)

    for idx, (sample_id, diagnosis) in enumerate(SST_REPRESENTATIVE):
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

        # Overlay vulnerable Sst cells
        sst_counts = {}
        total_sst_vuln = 0
        for sst_type in SST_PLOT_TYPES:
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

        title_color = DX_COLORS[diagnosis]
        ax.set_title(f"{sample_id} ({diagnosis})", fontsize=16,
                     fontweight="bold", color=title_color, pad=4)

        count_lines = [f"Vulnerable Sst: {total_sst_vuln} ({vuln_pct:.2f}%)"]
        for t in SST_PLOT_TYPES:
            count_lines.append(f"  {t}: {sst_counts[t]}")

        ax.text(0.02, 0.98, "\n".join(count_lines), transform=ax.transAxes,
                ha="left", va="top", fontsize=11, color="#dddddd",
                fontfamily="monospace",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="#1a1a1a",
                          edgecolor="#444444", alpha=0.85))

        style_dark_axis(ax)
        draw_inset(ax, df, loc="lower right")

        print(f"  Vulnerable Sst: {total_sst_vuln} ({vuln_pct:.2f}%)")

    # Row labels
    fig.text(0.50, 0.905, "Control", ha="center", fontsize=20,
             fontweight="bold", color=DX_COLORS["Control"])
    fig.text(0.50, 0.455, "SCZ", ha="center", fontsize=20,
             fontweight="bold", color=DX_COLORS["SCZ"])

    # Legend
    legend_elements = [
        Line2D([0], [0], marker='o', color=BG_COLOR,
               markerfacecolor=SST_COLORS[t],
               markeredgecolor="white", markeredgewidth=0.5,
               markersize=12, label=t, linewidth=0)
        for t in SST_PLOT_TYPES
    ]
    legend_elements.append(
        Line2D([0], [0], marker='o', color=BG_COLOR,
               markerfacecolor=ALL_CELL_COLOR,
               markersize=6, label="All cortical cells", linewidth=0, alpha=0.5)
    )
    fig.legend(handles=legend_elements, loc="lower center", ncol=4,
               fontsize=13, frameon=False, labelcolor="white",
               bbox_to_anchor=(0.5, 0.005))

    fig.suptitle("Vulnerable Sst subtypes in SCZ vs Control (Xenium)",
                 fontsize=22, fontweight="bold", color="white", y=0.96)

    outpath = os.path.join(PRESENTATION_DIR, "slide_sst_spatial_inset.png")
    plt.savefig(outpath, dpi=150, facecolor=BG_COLOR, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {outpath}")


if __name__ == "__main__":
    main()
