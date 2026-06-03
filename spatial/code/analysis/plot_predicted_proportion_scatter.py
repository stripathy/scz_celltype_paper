#!/usr/bin/env python3
"""
Neurons-only proportion scatter: MERFISH (all cells, L1-L6) vs Xenium.

Computes per-donor/sample proportions from scratch, then plots:
  Left panel: MERFISH vs Xenium uncropped
  Right panel: MERFISH vs Xenium cropped (L1-L6)

Features:
  - Neurons only (Glutamatergic + GABAergic)
  - +/-SD error bars across donors/samples
  - Frequency-weighted Pearson r (log-scale)
  - Dot size scaled by geometric mean frequency
  - Mean CV annotation

Output: output/presentation/slide_neurons_only_proportion_scatter.png
        output/presentation/proportion_stats_predicted_merfish.csv
"""

import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    BG_COLOR, H5AD_DIR, MERFISH_PATH, PRESENTATION_DIR,
    CONTROL_SAMPLES, CORTICAL_LAYERS, EXCLUDE_SAMPLES,
    classify_celltype, load_cells,
)

OUT_DIR = PRESENTATION_DIR
STATS_CSV = os.path.join(OUT_DIR, "proportion_stats_predicted_merfish.csv")

BG = BG_COLOR


def is_neuron(ct):
    _, cls = classify_celltype(ct)
    return cls in ("Glutamatergic", "GABAergic")


def weighted_corr(x, y, w):
    """Weighted Pearson correlation."""
    w = w / w.sum()
    mx = np.sum(w * x)
    my = np.sum(w * y)
    dx = x - mx
    dy = y - my
    cov = np.sum(w * dx * dy)
    sx = np.sqrt(np.sum(w * dx**2))
    sy = np.sqrt(np.sum(w * dy**2))
    if sx == 0 or sy == 0:
        return np.nan
    return cov / (sx * sy)


def compute_stats_csv():
    """Compute per-donor proportion stats for MERFISH and Xenium, save to CSV.

    MERFISH: all cells (cortical layers L1-L6, using Layer annotation),
             per-donor proportions at supertype level.
    Xenium: control samples only, using correlation classifier labels if available.
             Both uncropped (all QC-pass) and cropped (L1-L6) versions.
    """
    print("=== Computing proportion stats from scratch ===\n")

    # --- MERFISH: all cells with layer annotation in L1-L6 ---
    print("Loading MERFISH reference (all cortical cells)...")
    merfish = ad.read_h5ad(MERFISH_PATH, backed="r")
    obs_m = merfish.obs[["Donor ID", "Supertype", "Layer annotation"]].copy()
    obs_m["Layer annotation"] = obs_m["Layer annotation"].astype(str)
    obs_m = obs_m[obs_m["Layer annotation"].isin(CORTICAL_LAYERS)]
    obs_m = obs_m.rename(columns={"Donor ID": "donor", "Supertype": "celltype"})
    obs_m["celltype"] = obs_m["celltype"].astype(str)

    n_merfish = len(obs_m)
    n_donors = obs_m["donor"].nunique()
    print(f"  {n_merfish:,} cortical cells from {n_donors} donors")

    # Per-donor proportions
    merfish_records = []
    for donor in obs_m["donor"].unique():
        donor_df = obs_m[obs_m["donor"] == donor]
        total = len(donor_df)
        counts = donor_df["celltype"].value_counts()
        for ct, n in counts.items():
            merfish_records.append({"donor": donor, "celltype": ct,
                                    "proportion": n / total})
    merfish_props = pd.DataFrame(merfish_records)
    merfish_stats = merfish_props.groupby("celltype")["proportion"].agg(
        ["mean", "std"]).reset_index()
    merfish_stats.columns = ["celltype", "merfish_mean", "merfish_std"]
    print(f"  {len(merfish_stats)} MERFISH cell types")

    # --- Xenium: control samples ---
    print("\nLoading Xenium control samples...")
    controls = [s for s in CONTROL_SAMPLES if s not in EXCLUDE_SAMPLES]

    xenium_uncrop_records = []
    xenium_crop_records = []

    for sample_id in controls:
        # Uncropped: all QC-pass cells
        obs_all = load_cells(sample_id, cortical_only=False)
        obs_all = obs_all.rename(columns={"supertype_label": "celltype"})

        total_uncrop = len(obs_all)
        counts_uncrop = obs_all["celltype"].value_counts()
        for ct, n in counts_uncrop.items():
            xenium_uncrop_records.append({"donor": sample_id, "celltype": ct,
                                          "proportion": n / total_uncrop})

        # Cropped: L1-L6 only
        obs_crop = obs_all[obs_all["layer"].isin(CORTICAL_LAYERS)]
        total_crop = len(obs_crop)
        counts_crop = obs_crop["celltype"].value_counts()
        for ct, n in counts_crop.items():
            xenium_crop_records.append({"donor": sample_id, "celltype": ct,
                                        "proportion": n / total_crop})

        print(f"  {sample_id}: {total_uncrop:,} uncrop, {total_crop:,} crop")

    # Aggregate Xenium stats
    xenium_uncrop = pd.DataFrame(xenium_uncrop_records)
    xenium_uncrop_stats = xenium_uncrop.groupby("celltype")["proportion"].agg(
        ["mean", "std"]).reset_index()
    xenium_uncrop_stats.columns = ["celltype", "xenium_mean_uncrop", "xenium_std_uncrop"]

    xenium_crop = pd.DataFrame(xenium_crop_records)
    xenium_crop_stats = xenium_crop.groupby("celltype")["proportion"].agg(
        ["mean", "std"]).reset_index()
    xenium_crop_stats.columns = ["celltype", "xenium_mean_crop", "xenium_std_crop"]

    # Merge all
    stats = merfish_stats.merge(xenium_uncrop_stats, on="celltype", how="outer")
    stats = stats.merge(xenium_crop_stats, on="celltype", how="outer")
    stats = stats.sort_values("celltype").reset_index(drop=True)

    stats.to_csv(STATS_CSV, index=False)
    print(f"\nSaved: {STATS_CSV} ({len(stats)} cell types)")

    return stats


def plot_panel(ax, df, x_col, y_col, x_sd_col, y_sd_col, title, xlabel, ylabel,
               max_labels=12):
    """Plot neurons-only scatter with error bars, weighted corr, scaled dots."""
    # Filter to rows with valid data in both columns
    valid_mask = df[x_col].notna() & df[y_col].notna()
    valid_mask &= (df[x_col] > 0) & (df[y_col] > 0)
    d = df[valid_mask].copy()

    x = d[x_col].values
    y = d[y_col].values
    x_sd = d[x_sd_col].fillna(0).values
    y_sd = d[y_sd_col].fillna(0).values

    colors = [classify_celltype(ct)[0] for ct in d["celltype"]]

    # Size scaled by geometric mean frequency
    geo_mean = np.sqrt(x * y)
    sizes = 30 + 800 * (geo_mean / geo_mean.max())

    # Error bars
    for i in range(len(x)):
        ax.errorbar(x[i], y[i], xerr=x_sd[i], yerr=y_sd[i],
                    fmt="none", ecolor=colors[i], alpha=0.3, linewidth=1.0,
                    capsize=0, zorder=2)

    ax.scatter(x, y, c=colors, s=sizes, alpha=0.8, edgecolors="white",
               linewidths=0.5, zorder=5)

    # Label most deviant points
    d["log_ratio"] = np.log2((y + 1e-7) / (x + 1e-7))
    d["abs_log_ratio"] = d["log_ratio"].abs()
    top = d.nlargest(max_labels, "abs_log_ratio")

    for _, row in top.iterrows():
        ax.annotate(row["celltype"],
                    (row[x_col], row[y_col]),
                    fontsize=9, color="#dddddd", alpha=0.9,
                    xytext=(5, 5), textcoords="offset points")

    # Diagonal line
    lo = min(x.min(), y.min()) * 0.3
    hi = max(x.max(), y.max()) * 3
    ax.plot([lo, hi], [lo, hi], "--", color="#888888", linewidth=1.5,
            alpha=0.6, zorder=1)

    # Correlations
    log_x = np.log10(x)
    log_y = np.log10(y)

    r_log, _ = pearsonr(log_x, log_y)
    rho, _ = spearmanr(x, y)

    # Frequency-weighted correlation (on log scale)
    weights = geo_mean
    r_weighted = weighted_corr(log_x, log_y, weights)

    # Mean CV
    cv_x = x_sd / (x + 1e-10)
    cv_y = y_sd / (y + 1e-10)
    mean_cv_x = np.mean(cv_x[x_sd > 0])
    mean_cv_y = np.mean(cv_y[y_sd > 0])

    ax.text(0.04, 0.96,
            f"r = {r_log:.2f} (log)\n"
            f"r_w = {r_weighted:.2f} (freq-weighted)\n"
            f"\u03c1 = {rho:.2f}\n"
            f"n = {len(x)} types\n"
            f"CV: MERFISH={mean_cv_x:.2f}, Xenium={mean_cv_y:.2f}",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=13, color="#dddddd",
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
    # Regenerate stats CSV from scratch
    stats = compute_stats_csv()

    # Filter to neurons only
    neuron_mask = stats["celltype"].apply(is_neuron)
    neurons = stats[neuron_mask].copy()
    print(f"\nNeurons only: {neuron_mask.sum()} types")

    n_glut = sum(1 for ct in neurons["celltype"] if classify_celltype(ct)[1] == "Glutamatergic")
    n_gaba = sum(1 for ct in neurons["celltype"] if classify_celltype(ct)[1] == "GABAergic")
    print(f"  Glutamatergic: {n_glut}, GABAergic: {n_gaba}")

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(22, 10), facecolor=BG)

    # Left: MERFISH vs Xenium uncropped
    plot_panel(axes[0], neurons,
               "merfish_mean", "xenium_mean_uncrop",
               "merfish_std", "xenium_std_uncrop",
               "Xenium uncropped (neurons only)",
               "MERFISH (all cells, L1-L6)",
               "Xenium (all QC-pass cells)")

    # Right: MERFISH vs Xenium cropped
    plot_panel(axes[1], neurons,
               "merfish_mean", "xenium_mean_crop",
               "merfish_std", "xenium_std_crop",
               "Xenium cropped to L1-L6 (neurons only)",
               "MERFISH (all cells, L1-L6)",
               "Xenium (L1-L6 only)")

    # Shared legend
    from config import CLASS_COLORS
    legend_elements = [
        Line2D([0], [0], marker="o", color=BG, markerfacecolor=CLASS_COLORS["Glutamatergic"],
               markersize=12, label="Glutamatergic", linewidth=0),
        Line2D([0], [0], marker="o", color=BG, markerfacecolor=CLASS_COLORS["GABAergic"],
               markersize=12, label="GABAergic", linewidth=0),
    ]
    fig.legend(handles=legend_elements, loc="lower center", ncol=2,
               fontsize=14, frameon=False, labelcolor="white",
               bbox_to_anchor=(0.5, -0.01))

    plt.tight_layout(pad=2.0)

    outpath = os.path.join(OUT_DIR, "slide_neurons_only_proportion_scatter.png")
    plt.savefig(outpath, dpi=200, facecolor=BG, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {outpath}")


if __name__ == "__main__":
    main()
