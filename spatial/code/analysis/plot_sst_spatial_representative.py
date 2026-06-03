#!/usr/bin/env python3
"""
Spatial visualization of vulnerable Sst subtypes in SCZ vs Control.

Shows median-representative Xenium tissue sections with Sst_2, Sst_25, Sst_20,
Sst_3 cells highlighted against a background of all cortical cells.

Samples chosen to be near the group median for BOTH SST proportion and L6b
proportion, with good cortical layer distributions.

  Control: Br6389, Br8433
  SCZ:     Br6437, Br2421

Output: output/presentation/slide_sst_spatial_representative.png
"""

import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    PRESENTATION_DIR, BG_COLOR, DX_COLORS, ALL_CELL_COLOR, MARKER_SIZE_BG,
    SST_TYPES, SST_COLORS, REPRESENTATIVE_SAMPLES,
    load_cortical, style_dark_axis,
)

MARKER_SIZE_SST = 18


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
        df = load_cortical(sample_id)
        total_cells = len(df)
        print(f"  {total_cells:,} cortical cells")

        # Background: all cortical cells
        ax.scatter(df["x"], df["y"], s=MARKER_SIZE_BG, c=ALL_CELL_COLOR,
                   alpha=0.25, rasterized=True, linewidths=0)

        # Overlay each Sst subtype
        sst_counts = {}
        total_sst_vuln = 0
        for sst_type in SST_TYPES:
            sst_mask = df["supertype_label"] == sst_type
            sst_df = df[sst_mask]
            n = len(sst_df)
            sst_counts[sst_type] = n
            total_sst_vuln += n
            if n > 0:
                ax.scatter(sst_df["x"], sst_df["y"], s=MARKER_SIZE_SST,
                           c=SST_COLORS[sst_type], alpha=0.9,
                           linewidths=0.5, edgecolors="white",
                           zorder=5, rasterized=True)

        # Proportion
        vuln_pct = total_sst_vuln / total_cells * 100 if total_cells > 0 else 0

        # Title: sample ID and diagnosis
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

    # Shared legend at bottom
    legend_elements = [
        Line2D([0], [0], marker='o', color=BG_COLOR,
               markerfacecolor=SST_COLORS[t],
               markeredgecolor="white", markeredgewidth=0.5,
               markersize=12, label=t, linewidth=0)
        for t in SST_TYPES
    ]
    legend_elements.append(
        Line2D([0], [0], marker='o', color=BG_COLOR,
               markerfacecolor=ALL_CELL_COLOR,
               markersize=6, label="All cortical cells", linewidth=0, alpha=0.5)
    )
    fig.legend(handles=legend_elements, loc="lower center", ncol=5,
               fontsize=13, frameon=False, labelcolor="white",
               bbox_to_anchor=(0.5, 0.005))

    fig.suptitle(
        "Vulnerable Sst subtypes in SCZ vs Control (Xenium)",
        fontsize=22, fontweight="bold", color="white", y=0.96
    )

    outpath = os.path.join(PRESENTATION_DIR, "slide_sst_spatial_representative.png")
    plt.savefig(outpath, dpi=150, facecolor=BG_COLOR, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {outpath}")


if __name__ == "__main__":
    main()
