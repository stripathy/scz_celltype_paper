#!/usr/bin/env python3
"""
Stacked bar plot: number of cells per layer across Xenium samples.

Samples sorted by ascending WM proportion (lowest WM on left).
X-axis: sample IDs, colored by disease status (SCZ vs Control).
Bars stacked by layer (L1, L2/3, L4, L5, L6, WM, Vascular).

Output: output/presentation/slide_layer_stacked_bar.png
"""

import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    SAMPLE_TO_DX, DX_COLORS, BG_COLOR, H5AD_DIR, PRESENTATION_DIR,
)

OUT_DIR = PRESENTATION_DIR

BG = BG_COLOR

LAYER_ORDER = ["L1", "L2/3", "L4", "L5", "L6", "WM", "Vascular"]

# Script-specific hex-based layer colors (different from config's RGB tuple palette)
LAYER_COLORS = {
    "L1":        "#e6194b",
    "L2/3":      "#f58231",
    "L4":        "#ffe119",
    "L5":        "#3cb44b",
    "L6":        "#4363d8",
    "WM":        "#911eb4",
    "Vascular":  "#a9a9a9",
}


def main():
    # --- Collect per-sample layer counts ---
    records = []
    h5ad_files = sorted([f for f in os.listdir(H5AD_DIR) if f.endswith("_annotated.h5ad")])

    for fname in h5ad_files:
        sample_id = fname.replace("_annotated.h5ad", "")
        fpath = os.path.join(H5AD_DIR, fname)
        print(f"  Loading {sample_id}...")
        adata = ad.read_h5ad(fpath, backed="r")
        obs = adata.obs[["layer", "qc_pass"]].copy()
        obs = obs[obs["qc_pass"] == True]
        obs["layer"] = obs["layer"].astype(str)

        total = len(obs)
        layer_counts = obs["layer"].value_counts()

        row = {"sample_id": sample_id, "total": total}
        for layer in LAYER_ORDER:
            row[layer] = layer_counts.get(layer, 0)
        records.append(row)

    df = pd.DataFrame(records)

    # Add diagnosis
    df["diagnosis"] = df["sample_id"].map(SAMPLE_TO_DX)

    # Compute WM proportion for sorting
    df["wm_prop"] = df["WM"] / df["total"]

    # Sort by ascending WM proportion
    df = df.sort_values("wm_prop").reset_index(drop=True)

    print(f"\n{'Sample':<10} {'Dx':<10} {'Total':>8} {'WM%':>8}")
    print("-" * 40)
    for _, row in df.iterrows():
        print(f"{row['sample_id']:<10} {row['diagnosis']:<10} {row['total']:>8,} {row['wm_prop']:>7.1%}")

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(18, 8), facecolor=BG)
    ax.set_facecolor(BG)

    x = np.arange(len(df))
    bar_width = 0.75
    bottoms = np.zeros(len(df))

    for layer in LAYER_ORDER:
        counts = df[layer].values.astype(float)
        ax.bar(x, counts, bar_width, bottom=bottoms,
               color=LAYER_COLORS[layer], label=layer,
               edgecolor="#222222", linewidth=0.5)
        bottoms += counts

    # X-axis: sample IDs colored by diagnosis
    ax.set_xticks(x)
    ax.set_xticklabels(df["sample_id"], rotation=45, ha="right", fontsize=14)

    # Color x-tick labels by diagnosis
    for i, (_, row) in enumerate(df.iterrows()):
        color = DX_COLORS[row["diagnosis"]]
        ax.get_xticklabels()[i].set_color(color)
        ax.get_xticklabels()[i].set_fontweight("bold")

    # Y-axis formatting
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, p: f"{int(v):,}"))
    ax.set_ylabel("Number of cells (QC-pass)", fontsize=18, color="white")
    ax.set_xlabel("Sample (sorted by ascending WM proportion)", fontsize=16, color="white")
    ax.set_title("Layer composition across Xenium samples",
                 fontsize=22, fontweight="bold", color="white", pad=15)

    ax.tick_params(colors="white", labelsize=13)
    for spine in ax.spines.values():
        spine.set_color("#555555")
    ax.grid(axis="y", alpha=0.2, color="#555555")

    # Layer legend (top-right)
    layer_legend = ax.legend(loc="upper left", fontsize=13, frameon=True,
                             fancybox=True, framealpha=0.85,
                             edgecolor="#555555", title="Layer", title_fontsize=14)
    layer_legend.get_frame().set_facecolor("#222222")
    for text in layer_legend.get_texts():
        text.set_color("white")
    layer_legend.get_title().set_color("white")

    # Diagnosis legend (manually placed)
    from matplotlib.patches import Patch
    dx_patches = [
        Patch(facecolor=DX_COLORS["Control"], label="Control"),
        Patch(facecolor=DX_COLORS["SCZ"], label="SCZ"),
    ]
    dx_legend = fig.legend(handles=dx_patches, loc="upper right",
                           fontsize=14, frameon=True, fancybox=True,
                           framealpha=0.85, edgecolor="#555555",
                           title="Diagnosis", title_fontsize=15,
                           bbox_to_anchor=(0.98, 0.95))
    dx_legend.get_frame().set_facecolor("#222222")
    for text in dx_legend.get_texts():
        text.set_color("white")
    dx_legend.get_title().set_color("white")

    plt.tight_layout(pad=2.0)

    outpath = os.path.join(OUT_DIR, "slide_layer_stacked_bar.png")
    plt.savefig(outpath, dpi=200, facecolor=BG, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {outpath}")

    # Save data
    csv_path = os.path.join(OUT_DIR, "layer_counts_by_sample.csv")
    df.to_csv(csv_path, index=False)
    print(f"Saved: {csv_path}")


if __name__ == "__main__":
    main()
