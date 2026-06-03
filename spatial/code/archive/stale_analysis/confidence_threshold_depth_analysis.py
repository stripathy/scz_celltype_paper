#!/usr/bin/env python3
"""
Assess how confidence filtering affects Xenium-MERFISH depth consistency.

1. Cell loss per subclass and supertype at various confidence thresholds
2. Median depth comparison (MERFISH vs Xenium) at default (0.28) vs stringent (0.5)
3. Correlation improvement at both subclass and supertype levels

Output: output/plots/confidence_threshold_depth_analysis.png
        output/plots/confidence_threshold_cell_loss.csv
"""

import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
from adjustText import adjust_text

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    BG_COLOR, CLASS_COLORS, SUBCLASS_TO_CLASS, H5AD_DIR, MERFISH_PATH,
    CORTICAL_LAYERS, SUBCLASS_CONF_THRESH,
)

BASE_DIR = os.path.expanduser("~/Github/SCZ_Xenium")
OUT_DIR = os.path.join(BASE_DIR, "output", "plots")
os.makedirs(OUT_DIR, exist_ok=True)

# Confidence thresholds to evaluate
THRESHOLDS = [0.28, 0.3, 0.4, 0.5, 0.6, 0.7]
# Highlighted thresholds for the depth scatter comparison
COMPARE_THRESHOLDS = [0.28, 0.5]

MIN_CELLS = 20


# ── Data loading ──────────────────────────────────────────────────────

def load_merfish_cortical():
    """Load MERFISH cortical cells with manual depth annotation."""
    print("Loading MERFISH reference...")
    adata = ad.read_h5ad(MERFISH_PATH, backed="r")
    obs = adata.obs[["Subclass", "Supertype", "Normalized depth from pia",
                      "Layer annotation"]].copy()
    obs = obs.dropna(subset=["Normalized depth from pia"])
    obs = obs[obs["Layer annotation"].astype(str).isin(CORTICAL_LAYERS)]
    obs = obs.rename(columns={
        "Subclass": "subclass",
        "Supertype": "supertype",
        "Normalized depth from pia": "depth",
    })
    print(f"  MERFISH cortical with manual depth: {len(obs):,} cells")
    return obs


def load_all_xenium_cortical():
    """Load all Xenium cortical cells with confidence scores (no filtering yet)."""
    print("Loading all Xenium cortical cells...")
    h5ad_files = sorted([f for f in os.listdir(H5AD_DIR)
                         if f.endswith("_annotated.h5ad")])
    dfs = []
    for fname in h5ad_files:
        sample_id = fname.replace("_annotated.h5ad", "")
        fpath = os.path.join(H5AD_DIR, fname)
        adata = ad.read_h5ad(fpath, backed="r")
        cols = ["subclass_label", "supertype_label", "predicted_norm_depth",
                "layer", "qc_pass", "subclass_label_confidence",
                "supertype_label_confidence"]
        obs = adata.obs[cols].copy()
        # QC pass + cortical only
        mask = (obs["qc_pass"] == True) & (obs["layer"].astype(str).isin(CORTICAL_LAYERS))
        obs = obs[mask].copy()
        obs["sample_id"] = sample_id
        obs = obs.rename(columns={
            "subclass_label": "subclass",
            "supertype_label": "supertype",
            "predicted_norm_depth": "depth",
            "subclass_label_confidence": "subclass_conf",
            "supertype_label_confidence": "supertype_conf",
        })
        dfs.append(obs[["sample_id", "subclass", "supertype", "depth",
                        "subclass_conf", "supertype_conf"]])
    df = pd.concat(dfs, ignore_index=True)
    df["subclass_conf"] = df["subclass_conf"].astype(float)
    df["supertype_conf"] = df["supertype_conf"].astype(float)
    print(f"  Total cortical QC-pass: {len(df):,} cells")
    return df


# ── Analysis functions ────────────────────────────────────────────────

def compute_cell_loss(xenium_df, thresholds):
    """Compute number and fraction of cells lost per subclass at each threshold."""
    records = []
    for thresh in thresholds:
        for subclass in sorted(xenium_df["subclass"].unique()):
            sub_df = xenium_df[xenium_df["subclass"] == subclass]
            n_total = len(sub_df)
            n_pass = (sub_df["subclass_conf"] >= thresh).sum()
            n_lost = n_total - n_pass
            frac_lost = n_lost / n_total if n_total > 0 else 0
            records.append({
                "subclass": subclass,
                "threshold": thresh,
                "n_total": n_total,
                "n_pass": n_pass,
                "n_lost": n_lost,
                "frac_lost": frac_lost,
            })
    return pd.DataFrame(records)


def compute_median_depth(df, level="subclass"):
    """Compute median depth per cell type at given level."""
    grouped = df.groupby(level)["depth"].agg(["median", "count"]).reset_index()
    grouped.columns = [level, "median_depth", "n_cells"]
    grouped = grouped[grouped["n_cells"] >= MIN_CELLS]
    return grouped


def get_class(name, level="subclass"):
    """Get class for a cell type name."""
    if level == "subclass":
        return SUBCLASS_TO_CLASS.get(name, "Non-neuronal")
    # Supertype: infer from prefix
    for subclass, cls in SUBCLASS_TO_CLASS.items():
        if name.startswith(subclass.replace(" ", "_")) or \
           name.startswith(subclass.split()[0]):
            return cls
    prefix = name.split("_")[0]
    for subclass, cls in SUBCLASS_TO_CLASS.items():
        if prefix == subclass.split()[0].split("/")[0]:
            return cls
    return "Non-neuronal"


# ── Plotting ──────────────────────────────────────────────────────────

def plot_cell_loss_heatmap(ax, loss_df):
    """Heatmap of fraction lost by subclass x threshold."""
    pivot = loss_df.pivot(index="subclass", columns="threshold", values="frac_lost")
    # Sort by total cells (largest first)
    totals = loss_df[loss_df["threshold"] == loss_df["threshold"].min()].set_index("subclass")["n_total"]
    pivot = pivot.loc[totals.sort_values(ascending=False).index]

    im = ax.imshow(pivot.values, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels([f"{t:.2f}" for t in pivot.columns], fontsize=11, color="white")
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=10, color="white")
    ax.set_xlabel("Subclass confidence threshold", fontsize=14, color="white")
    ax.set_title("Fraction of cells lost per subclass", fontsize=16,
                 fontweight="bold", color="white", pad=10)

    # Annotate cells
    for i in range(len(pivot.index)):
        for j in range(len(pivot.columns)):
            val = pivot.values[i, j]
            n_pass = loss_df[(loss_df["subclass"] == pivot.index[i]) &
                             (loss_df["threshold"] == pivot.columns[j])]["n_pass"].values[0]
            text_color = "white" if val > 0.5 else "black"
            ax.text(j, i, f"{val:.0%}\n({n_pass:,})",
                    ha="center", va="center", fontsize=7, color=text_color)

    ax.set_facecolor(BG_COLOR)
    ax.tick_params(colors="white")
    return im


def plot_depth_scatter(ax, merfish_med, xenium_med, level, title,
                       label_all_subclass=True):
    """Scatter of median depth: MERFISH vs Xenium."""
    merged = merfish_med.merge(xenium_med, on=level, suffixes=("_merfish", "_xenium"))
    if len(merged) == 0:
        ax.text(0.5, 0.5, "No shared types", transform=ax.transAxes,
                ha="center", color="white", fontsize=14)
        return None, None

    x = merged["median_depth_merfish"].values
    y = merged["median_depth_xenium"].values
    names = merged[level].values

    colors = [CLASS_COLORS.get(get_class(n, level), "#888888") for n in names]
    min_n = np.minimum(merged["n_cells_merfish"].values, merged["n_cells_xenium"].values)
    sizes = np.clip(np.log10(min_n + 1) * 25, 15, 120)

    ax.scatter(x, y, c=colors, s=sizes, alpha=0.8, edgecolors="white",
               linewidths=0.5, zorder=3)

    lims = [min(x.min(), y.min()) - 0.02, max(x.max(), y.max()) + 0.02]
    ax.plot(lims, lims, "--", color="#666666", linewidth=1, zorder=1)

    r, p_r = pearsonr(x, y)
    rho, p_s = spearmanr(x, y)
    # Mean absolute deviation from identity
    mad = np.mean(np.abs(x - y))

    ax.text(0.03, 0.97,
            f"r = {r:.3f}\n\u03c1 = {rho:.3f}\nMAD = {mad:.3f}\nn = {len(merged)}",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=12, color="white",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#1a1a1a",
                      edgecolor="#444444", alpha=0.85))

    # Labels
    texts = []
    if level == "subclass" and label_all_subclass:
        for i, name in enumerate(names):
            texts.append(ax.text(x[i], y[i], name, fontsize=10, color="#dddddd"))
    elif level == "supertype":
        # Only label outliers
        deviations = np.abs(x - y)
        threshold = np.percentile(deviations, 90)
        for i, name in enumerate(names):
            if deviations[i] >= threshold:
                texts.append(ax.text(x[i], y[i], name, fontsize=7, color="#dddddd"))

    if texts:
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="#555555", lw=0.5),
                    force_text=(0.3, 0.3), force_points=(0.3, 0.3))

    ax.set_xlabel("MERFISH median depth", fontsize=13, color="white")
    ax.set_ylabel("Xenium median depth", fontsize=13, color="white")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_title(title, fontsize=14, fontweight="bold", color="white", pad=8)
    ax.tick_params(colors="white", labelsize=11)
    ax.set_facecolor(BG_COLOR)
    for spine in ax.spines.values():
        spine.set_color("#555555")

    return r, rho


def plot_correlation_vs_threshold(ax, xenium_df, merfish_df, level):
    """Line plot of correlation vs confidence threshold."""
    rs, rhos, mads, ns = [], [], [], []
    thresholds_fine = np.arange(0.2, 0.81, 0.05)

    merfish_med = compute_median_depth(merfish_df, level)

    for thresh in thresholds_fine:
        filtered = xenium_df[xenium_df["subclass_conf"] >= thresh]
        if len(filtered) < 100:
            rs.append(np.nan)
            rhos.append(np.nan)
            mads.append(np.nan)
            ns.append(0)
            continue

        xen_med = compute_median_depth(filtered, level)
        merged = merfish_med.merge(xen_med, on=level, suffixes=("_m", "_x"))
        if len(merged) < 5:
            rs.append(np.nan)
            rhos.append(np.nan)
            mads.append(np.nan)
            ns.append(len(merged))
            continue

        x = merged["median_depth_m"].values
        y = merged["median_depth_x"].values
        r, _ = pearsonr(x, y)
        rho, _ = spearmanr(x, y)
        mad = np.mean(np.abs(x - y))
        rs.append(r)
        rhos.append(rho)
        mads.append(mad)
        ns.append(len(merged))

    ax.plot(thresholds_fine, rs, "o-", color="#4fc3f7", label="Pearson r",
            markersize=5, linewidth=2)
    ax.plot(thresholds_fine, rhos, "s-", color="#ef5350", label="Spearman \u03c1",
            markersize=5, linewidth=2)

    # MAD on secondary axis
    ax2 = ax.twinx()
    ax2.plot(thresholds_fine, mads, "^-", color="#66bb6a", label="MAD",
             markersize=5, linewidth=2, alpha=0.8)
    ax2.set_ylabel("Mean abs deviation from identity", fontsize=12, color="#66bb6a")
    ax2.tick_params(axis="y", colors="#66bb6a", labelsize=11)
    ax2.spines["right"].set_color("#66bb6a")

    # Mark key thresholds
    for t in [0.28, 0.5]:
        ax.axvline(t, color="#888888", linestyle=":", alpha=0.6)
        ax.text(t, ax.get_ylim()[0] + 0.001, f"{t}", color="#aaaaaa",
                fontsize=9, ha="center", va="bottom")

    ax.set_xlabel("Subclass confidence threshold", fontsize=13, color="white")
    ax.set_ylabel("Correlation with MERFISH", fontsize=13, color="white")
    ax.set_title(f"Depth correlation vs threshold ({level})",
                 fontsize=14, fontweight="bold", color="white", pad=8)
    ax.legend(loc="lower left", fontsize=11, framealpha=0.7)
    ax2.legend(loc="lower right", fontsize=11, framealpha=0.7)
    ax.set_facecolor(BG_COLOR)
    ax.tick_params(colors="white", labelsize=11)
    for spine in ax.spines.values():
        spine.set_color("#555555")
    ax2.spines["right"].set_color("#66bb6a")
    ax2.spines["left"].set_color("#555555")
    ax2.spines["top"].set_color("#555555")
    ax2.spines["bottom"].set_color("#555555")


# ── Main ──────────────────────────────────────────────────────────────

def main():
    merfish_df = load_merfish_cortical()
    xenium_df = load_all_xenium_cortical()

    # ── Panel 1: Cell loss heatmap ──
    print("\nComputing cell loss by threshold...")
    loss_df = compute_cell_loss(xenium_df, THRESHOLDS)
    loss_df.to_csv(os.path.join(OUT_DIR, "confidence_threshold_cell_loss.csv"),
                   index=False)
    print("  Saved cell loss CSV")

    # Print summary table
    print("\n── Cell loss summary (total cortical cells) ──")
    for thresh in THRESHOLDS:
        sub = loss_df[loss_df["threshold"] == thresh]
        total_pass = sub["n_pass"].sum()
        total_all = sub["n_total"].sum()
        frac = 1 - total_pass / total_all
        print(f"  threshold={thresh:.2f}: {total_pass:>9,} / {total_all:,} "
              f"kept ({frac:.1%} lost)")

    # ── MERFISH reference medians ──
    merfish_sub = compute_median_depth(merfish_df, "subclass")
    merfish_sup = compute_median_depth(merfish_df, "supertype")

    # ── Create figure ──
    # Layout: 3 rows x 2 cols
    #   Row 0: cell loss heatmap (spans both cols)
    #   Row 1: subclass depth scatter at 0.28 vs 0.5
    #   Row 2: supertype depth scatter at 0.28 vs 0.5
    #   + correlation vs threshold plots

    fig = plt.figure(figsize=(24, 28), facecolor=BG_COLOR)

    # Row 0: heatmap
    ax_heatmap = fig.add_axes([0.05, 0.78, 0.88, 0.18])
    im = plot_cell_loss_heatmap(ax_heatmap, loss_df)
    cbar_ax = fig.add_axes([0.94, 0.78, 0.015, 0.18])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.ax.tick_params(colors="white")
    cbar.set_label("Fraction lost", color="white", fontsize=12)

    # Row 1: subclass depth scatter (0.28 vs 0.5) + correlation curve
    print("\nComputing depth comparisons...")
    axes_row1 = [
        fig.add_axes([0.05, 0.52, 0.28, 0.22]),   # subclass @ 0.28
        fig.add_axes([0.37, 0.52, 0.28, 0.22]),   # subclass @ 0.50
        fig.add_axes([0.69, 0.52, 0.28, 0.22]),   # corr vs thresh (subclass)
    ]

    for i, thresh in enumerate(COMPARE_THRESHOLDS):
        filtered = xenium_df[xenium_df["subclass_conf"] >= thresh]
        n_cells = len(filtered)
        xen_med = compute_median_depth(filtered, "subclass")
        r, rho = plot_depth_scatter(
            axes_row1[i], merfish_sub, xen_med, "subclass",
            f"Subclass (thresh={thresh:.2f}, n={n_cells:,} cells)")
        if r is not None:
            print(f"  Subclass @ {thresh}: r={r:.3f}, rho={rho:.3f}")

    plot_correlation_vs_threshold(axes_row1[2], xenium_df, merfish_df, "subclass")

    # Row 2: supertype depth scatter (0.28 vs 0.5) + correlation curve
    axes_row2 = [
        fig.add_axes([0.05, 0.26, 0.28, 0.22]),
        fig.add_axes([0.37, 0.26, 0.28, 0.22]),
        fig.add_axes([0.69, 0.26, 0.28, 0.22]),
    ]

    for i, thresh in enumerate(COMPARE_THRESHOLDS):
        filtered = xenium_df[xenium_df["subclass_conf"] >= thresh]
        n_cells = len(filtered)
        xen_med = compute_median_depth(filtered, "supertype")
        r, rho = plot_depth_scatter(
            axes_row2[i], merfish_sup, xen_med, "supertype",
            f"Supertype (thresh={thresh:.2f}, n={n_cells:,} cells)")
        if r is not None:
            print(f"  Supertype @ {thresh}: r={r:.3f}, rho={rho:.3f}")

    plot_correlation_vs_threshold(axes_row2[2], xenium_df, merfish_df, "supertype")

    # Row 3: per-subclass depth distribution comparison (selected types)
    # Show depth distributions for types that improve most with filtering
    ax_bottom = fig.add_axes([0.05, 0.03, 0.9, 0.19])
    plot_depth_shift_by_type(ax_bottom, xenium_df, merfish_df)

    fig.suptitle("Effect of confidence filtering on Xenium–MERFISH depth consistency",
                 fontsize=22, fontweight="bold", color="white", y=0.99)

    outpath = os.path.join(OUT_DIR, "confidence_threshold_depth_analysis.png")
    fig.savefig(outpath, dpi=150, facecolor=BG_COLOR, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {outpath}")


def plot_depth_shift_by_type(ax, xenium_df, merfish_df):
    """Show how median depth changes for each subclass with filtering."""
    merfish_med = compute_median_depth(merfish_df, "subclass").set_index("subclass")

    records = []
    for thresh in [0.28, 0.5]:
        filtered = xenium_df[xenium_df["subclass_conf"] >= thresh]
        xen_med = compute_median_depth(filtered, "subclass").set_index("subclass")

        for sub in xen_med.index:
            if sub not in merfish_med.index:
                continue
            records.append({
                "subclass": sub,
                "threshold": thresh,
                "xenium_depth": xen_med.loc[sub, "median_depth"],
                "merfish_depth": merfish_med.loc[sub, "median_depth"],
                "deviation": abs(xen_med.loc[sub, "median_depth"] -
                                 merfish_med.loc[sub, "median_depth"]),
                "n_cells": xen_med.loc[sub, "n_cells"],
            })

    rdf = pd.DataFrame(records)
    if len(rdf) == 0:
        return

    # Compare deviation at 0.28 vs 0.5
    dev_028 = rdf[rdf["threshold"] == 0.28].set_index("subclass")["deviation"]
    dev_050 = rdf[rdf["threshold"] == 0.50].set_index("subclass")["deviation"]
    shared = dev_028.index.intersection(dev_050.index)
    improvement = (dev_028[shared] - dev_050[shared]).sort_values(ascending=False)

    # Plot: grouped bar chart of deviation from MERFISH
    subclasses = improvement.index.tolist()
    x = np.arange(len(subclasses))
    width = 0.35

    colors_028 = [CLASS_COLORS.get(get_class(s, "subclass"), "#888") for s in subclasses]
    bars1 = ax.bar(x - width / 2, dev_028[subclasses].values, width,
                   color=colors_028, alpha=0.5, edgecolor="white", linewidth=0.5,
                   label="thresh = 0.28")
    bars2 = ax.bar(x + width / 2, dev_050[subclasses].values, width,
                   color=colors_028, alpha=1.0, edgecolor="white", linewidth=0.5,
                   label="thresh = 0.50")

    ax.set_xticks(x)
    ax.set_xticklabels(subclasses, rotation=45, ha="right", fontsize=10, color="white")
    ax.set_ylabel("|Xenium − MERFISH| median depth", fontsize=13, color="white")
    ax.set_title("Per-subclass depth deviation from MERFISH (sorted by improvement with filtering)",
                 fontsize=14, fontweight="bold", color="white", pad=8)
    ax.legend(fontsize=12, loc="upper right", framealpha=0.7)
    ax.set_facecolor(BG_COLOR)
    ax.tick_params(colors="white", labelsize=11)
    for spine in ax.spines.values():
        spine.set_color("#555555")
    ax.axhline(0, color="#555555", linewidth=0.5)


if __name__ == "__main__":
    main()
