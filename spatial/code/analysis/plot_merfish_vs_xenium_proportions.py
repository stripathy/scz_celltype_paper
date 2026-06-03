#!/usr/bin/env python3
"""
Compare naive (uncropped) Xenium cell type proportions vs MERFISH reference.

The MERFISH data is from SEA-AD (cropped to a specific cortical region with
manual depth annotations). The Xenium data uses entire tissue sections
(uncropped). This mismatch in tissue sampling is a key confound.

Produces 2-panel figure: subclass-level and supertype/cluster-level scatter
plots on semilog axes.

Output: output/presentation/slide5_proportion_scatter.png
"""

import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    BG_COLOR, H5AD_DIR, MERFISH_PATH, PRESENTATION_DIR,
    CONTROL_SAMPLES, SUBCLASS_CONF_THRESH, classify_celltype,
)

OUT_DIR = PRESENTATION_DIR

BG = BG_COLOR

# Xenium control samples only (for fair comparison to MERFISH controls)
CONTROLS = CONTROL_SAMPLES


def compute_merfish_proportions(level="Subclass"):
    """Compute mean proportions across MERFISH donors (all depth-annotated cells)."""
    print(f"Loading MERFISH ({level})...")
    merfish = ad.read_h5ad(MERFISH_PATH, backed="r")
    obs = merfish.obs[["Donor ID", level]].copy()
    obs.columns = ["donor", "celltype"]

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
    # Mean across donors
    mean_props = df.groupby("celltype")["proportion"].mean().reset_index()
    mean_props.columns = ["celltype", "merfish_prop"]
    return mean_props


def compute_xenium_proportions(level="subclass_label"):
    """Compute mean proportions across Xenium control samples (all QC-pass cells, uncropped)."""
    print(f"Loading Xenium controls ({level})...")
    records = []
    for sample_id in CONTROLS:
        fpath = os.path.join(H5AD_DIR, f"{sample_id}_annotated.h5ad")
        adata = ad.read_h5ad(fpath, backed="r")

        # Use correlation classifier labels if available
        has_corr = "corr_subclass" in adata.obs.columns
        corr_level_map = {"subclass_label": "corr_subclass",
                          "supertype_label": "corr_supertype"}
        use_col = corr_level_map.get(level, level) if has_corr else level

        cols = [use_col, "qc_pass"]
        if has_corr:
            cols.append("corr_qc_pass")
        else:
            cols.append("subclass_label_confidence")
        obs = adata.obs[cols].copy()
        obs = obs[obs["qc_pass"] == True]

        # Correlation QC or legacy confidence filter
        if has_corr:
            obs = obs[obs["corr_qc_pass"] == True]
        else:
            obs = obs[obs["subclass_label_confidence"].astype(float) >= SUBCLASS_CONF_THRESH]

        total = len(obs)
        counts = obs[use_col].value_counts()
        for ct, n in counts.items():
            records.append({"donor": sample_id, "celltype": str(ct),
                            "proportion": n / total})

    df = pd.DataFrame(records)
    mean_props = df.groupby("celltype")["proportion"].mean().reset_index()
    mean_props.columns = ["celltype", "xenium_prop"]
    return mean_props


def build_supertype_mapping():
    """Build mapping from Xenium supertype_label to MERFISH Supertype.

    The Xenium labels were transferred from the SEA-AD taxonomy via
    MapMyCells, so they match the MERFISH taxonomy directly.
    """
    # Just return the label columns -- we'll merge on matching names
    pass


def plot_scatter(ax, merged, title, max_labels=15):
    """Plot MERFISH vs Xenium proportion scatter on semilog axes."""
    x = merged["merfish_prop"].values
    y = merged["xenium_prop"].values

    colors = [classify_celltype(ct)[0] for ct in merged["celltype"]]

    ax.scatter(x, y, c=colors, s=70, alpha=0.8, edgecolors="white",
               linewidths=0.5, zorder=5)

    # Label only the most deviant points (by log-ratio distance from diagonal)
    merged = merged.copy()
    merged["log_ratio"] = np.log2(
        (merged["xenium_prop"] + 1e-7) / (merged["merfish_prop"] + 1e-7)
    )
    merged["abs_log_ratio"] = merged["log_ratio"].abs()
    # Also label abundant types
    merged["score"] = merged["abs_log_ratio"] + np.log10(
        merged[["merfish_prop", "xenium_prop"]].max(axis=1) + 1e-7
    ) * 0.5
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

    ax.set_xlabel("MERFISH (SEA-AD) proportion", fontsize=16, color="white")
    ax.set_ylabel("Xenium (uncropped) proportion", fontsize=16, color="white")
    ax.set_title(title, fontsize=20, fontweight="bold", color="white", pad=10)
    ax.set_facecolor(BG)
    ax.tick_params(colors="white", labelsize=13)
    for spine in ax.spines.values():
        spine.set_color("#555555")
    ax.grid(True, alpha=0.2, color="#555555")


def main():
    # --- Subclass level ---
    merfish_sub = compute_merfish_proportions("Subclass")
    xenium_sub = compute_xenium_proportions("subclass_label")
    merged_sub = merfish_sub.merge(xenium_sub, on="celltype", how="inner")
    print(f"Subclass: {len(merged_sub)} cell types matched")

    # --- Supertype / cluster level ---
    merfish_sup = compute_merfish_proportions("Supertype")
    xenium_clust = compute_xenium_proportions("supertype_label")
    merged_sup = merfish_sup.merge(xenium_clust, on="celltype", how="inner")
    print(f"Supertype: {len(merged_sup)} cell types matched")

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(20, 9), facecolor=BG)

    plot_scatter(axes[0], merged_sub, "Subclass level")
    plot_scatter(axes[1], merged_sup, "Supertype level")

    plt.tight_layout(pad=2.0)

    outpath = os.path.join(OUT_DIR, "slide5_proportion_scatter.png")
    plt.savefig(outpath, dpi=200, facecolor=BG)
    plt.close()
    print(f"\nSaved: {outpath}")

    # Also save the data
    merged_sub.to_csv(os.path.join(OUT_DIR, "slide5_subclass_proportions.csv"), index=False)
    merged_sup.to_csv(os.path.join(OUT_DIR, "slide5_supertype_proportions.csv"), index=False)


if __name__ == "__main__":
    main()
