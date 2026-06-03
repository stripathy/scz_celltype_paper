#!/usr/bin/env python3
"""
Side-by-side comparison: manually annotated MERFISH vs Xenium uncropped/cropped.

Left panel: MERFISH (manually annotated cortical cells) vs Xenium (uncropped)
Right panel: MERFISH (manually annotated cortical cells) vs Xenium (cropped L1-L6)

The MERFISH x-axis is fixed: only cells with manual depth annotations (curated
cortical regions, ~369k cells). The Xenium y-axis changes from all QC-pass cells
(left) to L1-L6 cropped (right), showing how cropping improves agreement.

Output: output/presentation/slide_cropped_proportion_scatter.png
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    BG_COLOR, PRESENTATION_DIR, CONTROL_SAMPLES, CORTICAL_LAYERS,
    classify_celltype, load_cells, load_merfish_cortical,
)

OUT_DIR = PRESENTATION_DIR

BG = BG_COLOR

# Xenium control samples
CONTROLS = CONTROL_SAMPLES


def compute_merfish_manual_proportions(level="Supertype"):
    """Compute mean proportions across MERFISH donors using only manually
    annotated cells (those with 'Normalized depth from pia')."""
    print(f"  MERFISH manually annotated ({level})...")
    obs = load_merfish_cortical()

    # Map level to column name
    col_map = {"Supertype": "supertype", "Subclass": "subclass"}
    celltype_col = col_map.get(level, level.lower())
    obs = obs.rename(columns={celltype_col: "celltype"})

    n_cells = len(obs)
    n_donors = obs["donor"].nunique()
    print(f"    {n_cells:,} manually annotated cells from {n_donors} donors")

    # Per-donor proportions
    records = []
    for donor in obs["donor"].unique():
        donor_df = obs[obs["donor"] == donor]
        total = len(donor_df)
        counts = donor_df["celltype"].value_counts()
        for ct, n in counts.items():
            records.append({"donor": donor, "celltype": ct,
                            "proportion": n / total})

    df = pd.DataFrame(records)
    mean_props = df.groupby("celltype")["proportion"].mean().reset_index()
    mean_props.columns = ["celltype", "merfish_prop"]
    return mean_props


def compute_xenium_proportions(level="supertype_label", crop_layers=None):
    """Compute mean proportions across Xenium controls.

    If crop_layers is provided, filter to those layer values first.
    """
    tag = "L1-L6 cropped" if crop_layers else "uncropped"
    print(f"  Xenium controls ({level}, {tag})...")

    # Map level to standard column name
    col_map = {"subclass_label": "subclass_label",
               "supertype_label": "supertype_label"}
    use_col = col_map.get(level, level)

    records = []
    total_cells = 0
    for sample_id in CONTROLS:
        cortical_only = crop_layers is not None
        obs = load_cells(sample_id, cortical_only=cortical_only)

        if crop_layers:
            obs = obs[obs["layer"].isin(crop_layers)]

        total = len(obs)
        total_cells += total
        counts = obs[use_col].value_counts()
        for ct, n in counts.items():
            records.append({"donor": sample_id, "celltype": str(ct),
                            "proportion": n / total})

    print(f"    {total_cells:,} cells from {len(CONTROLS)} samples")
    df = pd.DataFrame(records)
    mean_props = df.groupby("celltype")["proportion"].mean().reset_index()
    mean_props.columns = ["celltype", "xenium_prop"]
    return mean_props


def plot_scatter(ax, merged, title, xlabel, ylabel, max_labels=15):
    """Plot MERFISH vs Xenium proportion scatter on semilog axes."""
    x = merged["merfish_prop"].values
    y = merged["xenium_prop"].values

    colors = [classify_celltype(ct)[0] for ct in merged["celltype"]]

    ax.scatter(x, y, c=colors, s=70, alpha=0.8, edgecolors="white",
               linewidths=0.5, zorder=5)

    # Label most deviant points
    merged = merged.copy()
    merged["log_ratio"] = np.log2(
        (merged["xenium_prop"] + 1e-7) / (merged["merfish_prop"] + 1e-7)
    )
    merged["abs_log_ratio"] = merged["log_ratio"].abs()
    top = merged.nlargest(max_labels, "abs_log_ratio")

    for _, row in top.iterrows():
        ax.annotate(row["celltype"],
                    (row["merfish_prop"], row["xenium_prop"]),
                    fontsize=9, color="#dddddd", alpha=0.9,
                    xytext=(5, 5), textcoords="offset points")

    # Diagonal line
    lo = min(x[x > 0].min(), y[y > 0].min()) * 0.3
    hi = max(x.max(), y.max()) * 3
    ax.plot([lo, hi], [lo, hi], "--", color="#888888", linewidth=1.5,
            alpha=0.6, zorder=1)

    # Correlation (log-transformed)
    valid = (x > 0) & (y > 0)
    if valid.sum() > 3:
        r, p = pearsonr(np.log10(x[valid]), np.log10(y[valid]))
        rho, _ = spearmanr(x[valid], y[valid])
        ax.text(0.04, 0.96,
                f"r = {r:.2f} (log-scale)\n\u03c1 = {rho:.2f}\nn = {valid.sum()} types",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=14, color="#dddddd",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="#333333",
                          edgecolor="#555555", alpha=0.85))

    ax.set_xscale("log")
    ax.set_yscale("log")

    # Human-readable percentage ticks
    def pct_formatter(val, pos):
        pct = val * 100
        if pct >= 1:
            return f"{pct:.0f}%"
        elif pct >= 0.1:
            return f"{pct:.1f}%"
        elif pct >= 0.01:
            return f"{pct:.2f}%"
        else:
            return f"{pct:.3f}%"

    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_major_formatter(mticker.FuncFormatter(pct_formatter))
        axis.set_minor_formatter(mticker.NullFormatter())

    ax.set_xlabel(xlabel, fontsize=16, color="white")
    ax.set_ylabel(ylabel, fontsize=16, color="white")
    ax.set_title(title, fontsize=20, fontweight="bold", color="white", pad=10)
    ax.set_facecolor(BG)
    ax.tick_params(colors="white", labelsize=13)
    for spine in ax.spines.values():
        spine.set_color("#555555")
    ax.grid(True, alpha=0.2, color="#555555")


def main():
    # --- MERFISH: manually annotated cells (fixed x-axis for both panels) ---
    print("Computing MERFISH manually annotated proportions:")
    merfish_manual = compute_merfish_manual_proportions("Supertype")

    # --- Left panel: Xenium uncropped ---
    print("\nComputing Xenium UNCROPPED proportions:")
    xenium_uncrop = compute_xenium_proportions("supertype_label", crop_layers=None)
    merged_uncrop = merfish_manual.merge(xenium_uncrop, on="celltype", how="inner")
    print(f"  Matched: {len(merged_uncrop)} supertypes")

    # --- Right panel: Xenium cropped to L1-L6 ---
    print("\nComputing Xenium CROPPED (L1-L6) proportions:")
    xenium_crop = compute_xenium_proportions("supertype_label",
                                              crop_layers=CORTICAL_LAYERS)
    merged_crop = merfish_manual.merge(xenium_crop, on="celltype", how="inner")
    print(f"  Matched: {len(merged_crop)} supertypes")

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(20, 9), facecolor=BG)

    # Shared x-axis range for visual comparison
    plot_scatter(axes[0], merged_uncrop,
                 "Xenium uncropped",
                 "MERFISH (manually annotated, cortical)",
                 "Xenium (all QC-pass cells)")

    plot_scatter(axes[1], merged_crop,
                 "Xenium cropped to L1-L6",
                 "MERFISH (manually annotated, cortical)",
                 "Xenium (L1-L6 only)")

    # Shared legend
    legend_elements = [
        Line2D([0], [0], marker="o", color=BG, markerfacecolor="#3399dd",
               markersize=10, label="Glutamatergic", linewidth=0),
        Line2D([0], [0], marker="o", color=BG, markerfacecolor="#ee4433",
               markersize=10, label="GABAergic", linewidth=0),
        Line2D([0], [0], marker="o", color=BG, markerfacecolor="#44bb44",
               markersize=10, label="Non-neuronal", linewidth=0),
    ]
    fig.legend(handles=legend_elements, loc="lower center", ncol=3,
               fontsize=13, frameon=False, labelcolor="white",
               bbox_to_anchor=(0.5, -0.01))

    plt.tight_layout(pad=2.0)

    outpath = os.path.join(OUT_DIR, "slide_cropped_proportion_scatter.png")
    plt.savefig(outpath, dpi=200, facecolor=BG, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {outpath}")

    # Save data
    merged_uncrop.to_csv(
        os.path.join(OUT_DIR, "manual_merfish_vs_xenium_uncropped.csv"),
        index=False)
    merged_crop.to_csv(
        os.path.join(OUT_DIR, "manual_merfish_vs_xenium_cropped.csv"),
        index=False)


if __name__ == "__main__":
    main()
