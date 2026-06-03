#!/usr/bin/env python3
"""
Plot median depth-from-pia for each cell type at subclass and supertype levels,
comparing MERFISH vs Xenium (cortical-cropped datasets).

For MERFISH: uses manually annotated "Normalized depth from pia" for cells
that have it, restricted to cortical layers (L1-L6).

For Xenium: uses predicted_norm_depth for cortical cells (layer in L1-L6),
aggregated across all 24 samples.

Produces two scatter plots:
  1) Subclass-level median depth: MERFISH vs Xenium
  2) Supertype-level median depth: MERFISH vs Xenium

Output: output/presentation/slide_median_depth_by_celltype.png
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
from adjustText import adjust_text

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    BG_COLOR, CLASS_COLORS, SUBCLASS_TO_CLASS, H5AD_DIR,
    PRESENTATION_DIR, EXCLUDE_SAMPLES,
    load_cells, load_merfish_cortical as _load_merfish_cortical,
)

OUT_DIR = PRESENTATION_DIR

# Minimum number of cells to include a cell type
MIN_CELLS = 20


def load_merfish():
    """Load MERFISH cortical cells with manual depth annotation."""
    print("Loading MERFISH reference...")
    obs = _load_merfish_cortical()
    print(f"  MERFISH cortical with manual depth: {len(obs):,} cells")
    return obs


def load_xenium_cortical():
    """Load all Xenium cortical cells (L1-L6) across all samples."""
    print("Loading Xenium samples...")
    h5ad_files = sorted(glob.glob(os.path.join(H5AD_DIR, "*_annotated.h5ad")))

    dfs = []
    for fpath in h5ad_files:
        sample_id = os.path.basename(fpath).replace("_annotated.h5ad", "")
        if sample_id in EXCLUDE_SAMPLES:
            print(f"  {sample_id}: SKIPPED (excluded)")
            continue

        obs = load_cells(sample_id, cortical_only=True,
                         extra_obs_columns=["predicted_norm_depth"])
        obs = obs.rename(columns={
            "subclass_label": "subclass",
            "supertype_label": "supertype",
            "predicted_norm_depth": "depth",
        })
        dfs.append(obs[["subclass", "supertype", "depth"]])
        print(f"  {sample_id}: {len(obs):,} cortical cells")

    df = pd.concat(dfs, ignore_index=True)
    print(f"  Xenium total cortical: {len(df):,} cells")
    return df


def compute_median_depth(df, level="subclass"):
    """Compute median depth per cell type at given level."""
    grouped = df.groupby(level)["depth"].agg(["median", "count"]).reset_index()
    grouped.columns = [level, "median_depth", "n_cells"]
    grouped = grouped[grouped["n_cells"] >= MIN_CELLS]
    return grouped


def get_class_for_supertype(supertype_name):
    """Infer class from supertype prefix."""
    for subclass, cls in SUBCLASS_TO_CLASS.items():
        if supertype_name.startswith(subclass.replace(" ", "_")) or \
           supertype_name.startswith(subclass.split()[0]):
            return cls
    # Try matching by common prefixes
    prefix = supertype_name.split("_")[0]
    for subclass, cls in SUBCLASS_TO_CLASS.items():
        if prefix == subclass.split()[0].split("/")[0]:
            return cls
    return "Non-neuronal"  # fallback


def make_scatter(ax, merfish_df, xenium_df, level, label_all=True,
                 label_threshold=0.1):
    """Make a depth scatter plot for one level (subclass or supertype)."""
    merged = merfish_df.merge(xenium_df, on=level, suffixes=("_merfish", "_xenium"))

    if len(merged) == 0:
        ax.text(0.5, 0.5, "No shared cell types", transform=ax.transAxes,
                ha="center", color="white", fontsize=14)
        return

    x = merged["median_depth_merfish"].values
    y = merged["median_depth_xenium"].values
    names = merged[level].values

    # Assign colors by class
    if level == "subclass":
        colors = [CLASS_COLORS.get(SUBCLASS_TO_CLASS.get(n, "Non-neuronal"), "#888888")
                  for n in names]
    else:
        colors = [CLASS_COLORS.get(get_class_for_supertype(n), "#888888")
                  for n in names]

    # Size by min cell count
    min_n = np.minimum(merged["n_cells_merfish"].values,
                       merged["n_cells_xenium"].values)
    sizes = np.clip(np.log10(min_n) * 25, 15, 120)

    ax.scatter(x, y, c=colors, s=sizes, alpha=0.8, edgecolors="white",
               linewidths=0.5, zorder=3)

    # Identity line
    lims = [min(x.min(), y.min()) - 0.02, max(x.max(), y.max()) + 0.02]
    ax.plot(lims, lims, "--", color="#666666", linewidth=1, zorder=1)

    # Correlation
    r, p = pearsonr(x, y)
    rho, p_s = spearmanr(x, y)
    ax.text(0.03, 0.97, f"r = {r:.3f}\nrho = {rho:.3f}\nn = {len(merged)}",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=13, color="white",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#1a1a1a",
                      edgecolor="#444444", alpha=0.85))

    # Labels
    texts = []
    if level == "subclass" or label_all:
        for i, name in enumerate(names):
            if level == "supertype" and not label_all:
                # Only label if far from identity line
                if abs(x[i] - y[i]) < label_threshold:
                    continue
            texts.append(ax.text(x[i], y[i], name, fontsize=8 if level == "supertype" else 11,
                                  color="#dddddd", ha="center", va="bottom"))
    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="#555555",
                                                    lw=0.5),
                    force_text=(0.3, 0.3), force_points=(0.3, 0.3),
                    expand_text=(1.1, 1.2), expand_points=(1.1, 1.2))

    ax.set_xlabel("MERFISH median depth from pia", fontsize=14, color="white")
    ax.set_ylabel("Xenium median depth from pia", fontsize=14, color="white")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.tick_params(colors="white", labelsize=12)
    ax.set_facecolor(BG_COLOR)
    for spine in ax.spines.values():
        spine.set_color("#555555")


def main():
    # Load data
    merfish_df = load_merfish()
    xenium_df = load_xenium_cortical()

    # Compute median depths at both levels
    merfish_subclass = compute_median_depth(merfish_df, "subclass")
    xenium_subclass = compute_median_depth(xenium_df, "subclass")

    merfish_supertype = compute_median_depth(merfish_df, "supertype")
    xenium_supertype = compute_median_depth(xenium_df, "supertype")

    print(f"\nMERFISH subclasses: {len(merfish_subclass)}")
    print(f"Xenium subclasses:  {len(xenium_subclass)}")
    print(f"MERFISH supertypes: {len(merfish_supertype)}")
    print(f"Xenium supertypes:  {len(xenium_supertype)}")

    # Make figure
    fig, axes = plt.subplots(1, 2, figsize=(18, 8), facecolor=BG_COLOR)

    # Left: subclass
    make_scatter(axes[0], merfish_subclass, xenium_subclass, "subclass")
    axes[0].set_title("Subclass level", fontsize=18, fontweight="bold",
                      color="white", pad=10)

    # Right: supertype — only label outliers (far from identity)
    make_scatter(axes[1], merfish_supertype, xenium_supertype, "supertype",
                 label_all=False, label_threshold=0.08)
    axes[1].set_title("Supertype level", fontsize=18, fontweight="bold",
                      color="white", pad=10)

    # Class legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color=BG_COLOR,
               markerfacecolor=CLASS_COLORS[cls],
               markersize=10, label=cls, linewidth=0)
        for cls in ["Glutamatergic", "GABAergic", "Non-neuronal"]
    ]
    fig.legend(handles=legend_elements, loc="lower center", ncol=3,
               fontsize=13, frameon=False, labelcolor="white",
               bbox_to_anchor=(0.5, 0.01))

    fig.suptitle("Median depth from pia: MERFISH vs Xenium (cortical cells)",
                 fontsize=20, fontweight="bold", color="white", y=0.98)

    plt.tight_layout(rect=[0, 0.06, 1, 0.94])

    outpath = os.path.join(OUT_DIR, "slide_median_depth_by_celltype.png")
    plt.savefig(outpath, dpi=150, facecolor=BG_COLOR, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {outpath}")

    # Save data
    # Merge and save for reference
    sub_merged = merfish_subclass.merge(xenium_subclass, on="subclass",
                                         suffixes=("_merfish", "_xenium"))
    sub_merged.to_csv(os.path.join(OUT_DIR, "median_depth_subclass.csv"), index=False)

    sup_merged = merfish_supertype.merge(xenium_supertype, on="supertype",
                                          suffixes=("_merfish", "_xenium"))
    sup_merged.to_csv(os.path.join(OUT_DIR, "median_depth_supertype.csv"), index=False)
    print(f"Saved CSVs to {OUT_DIR}")


if __name__ == "__main__":
    main()
