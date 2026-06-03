#!/usr/bin/env python3
"""
Geometry gallery: show 4 Xenium sections illustrating tissue shape variability.

Cells colored by class (Glutamatergic, GABAergic, Non-neuronal) on dark
background. Shows that sections differ dramatically in shape, size, and
proportion of cortical vs non-cortical tissue.

Output: output/presentation/slide4_geometry_gallery.png
"""

import os
import sys
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import H5AD_DIR, PRESENTATION_DIR, BG_COLOR, SAMPLE_TO_DX, DX_COLORS

OUT_DIR = PRESENTATION_DIR
BG = BG_COLOR

# Pick 4 samples that show range of tissue geometries
# Mix of Control/SCZ, varied sizes and non-cortical %
SAMPLES = [
    "Br5622",   # Control - largest (69k cells), 7% non-cortical
    "Br6432",   # Control - smallest (40k cells), 14% non-cortical
    "Br8772",   # SCZ - 55k cells, 16% non-cortical (most)
    "Br5746",   # SCZ - 44k cells, 4% non-cortical (least)
]

# Class colors matching Allen Institute / SEA-AD reference
CLASS_COLORS = {
    "Neuronal: Glutamatergic": "#00ADF8",
    "Neuronal: GABAergic": "#F05A28",
    "Non-neuronal and Non-neural": "#808080",
    "Glutamatergic": "#00ADF8",
    "GABAergic": "#F05A28",
    "Non-neuronal": "#808080",
}

# Short labels for legend
CLASS_SHORT = {
    "Neuronal: Glutamatergic": "Glutamatergic",
    "Neuronal: GABAergic": "GABAergic",
    "Non-neuronal and Non-neural": "Non-neuronal",
    "Glutamatergic": "Glutamatergic",
    "GABAergic": "GABAergic",
    "Non-neuronal": "Non-neuronal",
}


def load_sample(sample_id):
    """Load all QC-pass cells with spatial coords and class labels."""
    fpath = os.path.join(H5AD_DIR, f"{sample_id}_annotated.h5ad")
    adata = ad.read_h5ad(fpath, backed="r")
    obs = adata.obs[["sample_id", "class_label", "spatial_domain", "qc_pass"]].copy()
    coords = adata.obsm["spatial"]
    obs["x"] = coords[:, 0]
    obs["y"] = coords[:, 1]
    # QC-pass only
    obs = obs[obs["qc_pass"] == True].copy()
    return obs


DIAGNOSIS = SAMPLE_TO_DX


def main():
    fig = plt.figure(figsize=(26, 8), facecolor=BG)
    gs = gridspec.GridSpec(1, 4, figure=fig,
                           wspace=0.08,
                           left=0.02, right=0.98, top=0.82, bottom=0.14)

    for idx, sample_id in enumerate(SAMPLES):
        ax = fig.add_subplot(gs[0, idx])
        ax.set_facecolor(BG)

        print(f"Loading {sample_id}...")
        df = load_sample(sample_id)
        total = len(df)
        diagnosis = DIAGNOSIS.get(sample_id, "Unknown")

        # Count cortical vs non-cortical
        n_cortical = (df["spatial_domain"] == "Cortical").sum()
        n_noncortical = total - n_cortical
        pct_noncortical = n_noncortical / total * 100

        # Plot all cells colored by class
        for cls_name, color in CLASS_COLORS.items():
            mask = df["class_label"] == cls_name
            sub = df[mask]
            if len(sub) > 0:
                ax.scatter(sub["x"], sub["y"], s=1.2, c=color,
                           alpha=0.6, rasterized=True, linewidths=0)

        # Panel label
        diag_color = DX_COLORS[diagnosis]
        ax.set_title(f"{sample_id} ({diagnosis})", fontsize=18,
                     fontweight="bold", color=diag_color, pad=8)

        # Stats annotation
        stats_text = f"{total:,} cells  |  {pct_noncortical:.0f}% non-cortical"
        ax.text(0.5, -0.04, stats_text, transform=ax.transAxes,
                ha="center", va="top", fontsize=14, color="#cccccc",
                fontfamily="monospace")

        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

    # Shared legend at bottom
    legend_labels = ["Glutamatergic", "GABAergic", "Non-neuronal"]
    legend_colors = ["#3399dd", "#ee4433", "#44bb44"]
    legend_elements = [
        Line2D([0], [0], marker='o', color=BG,
               markerfacecolor=c, markersize=12,
               label=l, linewidth=0)
        for l, c in zip(legend_labels, legend_colors)
    ]
    fig.legend(handles=legend_elements, loc="lower center", ncol=3,
               fontsize=14, frameon=False, labelcolor="white",
               bbox_to_anchor=(0.5, 0.0))

    outpath = os.path.join(OUT_DIR, "slide4_geometry_gallery.png")
    plt.savefig(outpath, dpi=200, facecolor=BG)
    plt.close()
    print(f"\nSaved: {outpath}")


if __name__ == "__main__":
    main()
