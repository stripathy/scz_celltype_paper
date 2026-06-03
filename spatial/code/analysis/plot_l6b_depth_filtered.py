#!/usr/bin/env python3
"""
Regenerate L6b depth Xenium vs MERFISH comparison at two confidence thresholds.

Row 1: current threshold (0.28) — same as original figure
Row 2: stringent threshold (0.50)
Row 3: summary of improvement

Output: output/plots/l6b_depth_xenium_vs_merfish.png (overwrites original)
"""

import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import H5AD_DIR, MERFISH_PATH, CORTICAL_LAYERS, EXCLUDE_SAMPLES

BASE_DIR = os.path.expanduser("~/Github/SCZ_Xenium")
OUT_DIR = os.path.join(BASE_DIR, "output", "plots")

# Layer boundaries for vertical lines on depth plots
LAYER_BOUNDS = [0.1225, 0.4696, 0.5443, 0.7079, 0.9275]
LAYER_NAMES = ["L1", "L2/3", "L4", "L5", "L6", "WM"]


def load_merfish_l6b():
    """Load MERFISH L6b cells with manual depth annotations."""
    print("Loading MERFISH L6b...")
    adata = ad.read_h5ad(MERFISH_PATH, backed="r")
    obs = adata.obs[["Subclass", "Normalized depth from pia",
                      "Layer annotation"]].copy()
    obs = obs[obs["Subclass"].astype(str) == "L6b"]
    obs = obs.dropna(subset=["Normalized depth from pia"])
    obs["depth"] = obs["Normalized depth from pia"].astype(float)
    obs["layer"] = obs["Layer annotation"].astype(str)
    print(f"  MERFISH L6b with depth: {len(obs):,}")
    return obs


def load_xenium_l6b():
    """Load all Xenium L6b cells (QC-pass, using correlation labels if available)."""
    print("Loading Xenium L6b...")
    h5ad_files = sorted([f for f in os.listdir(H5AD_DIR)
                         if f.endswith("_annotated.h5ad")])
    dfs = []
    for fname in h5ad_files:
        sample_id = fname.replace("_annotated.h5ad", "")
        if sample_id in EXCLUDE_SAMPLES:
            print(f"  {sample_id}: SKIPPED (excluded)")
            continue
        fpath = os.path.join(H5AD_DIR, fname)
        adata = ad.read_h5ad(fpath, backed="r")

        # Use correlation classifier labels if available
        has_corr = "corr_subclass" in adata.obs.columns
        subclass_col = "corr_subclass" if has_corr else "subclass_label"

        cols = [subclass_col, "predicted_norm_depth", "layer", "qc_pass"]
        if has_corr:
            cols.append("corr_qc_pass")
            cols.append("corr_subclass_margin")
        else:
            cols.append("subclass_label_confidence")
        obs = adata.obs[cols].copy()

        # Filter to QC-pass L6b cells
        obs = obs[(obs["qc_pass"] == True) &
                  (obs[subclass_col].astype(str) == "L6b")]

        # Apply correlation QC filter
        if has_corr:
            obs = obs[obs["corr_qc_pass"] == True]

        obs["depth"] = obs["predicted_norm_depth"].astype(float)
        obs["layer"] = obs["layer"].astype(str)
        # Use margin as confidence proxy for filtering in plot_row
        if has_corr:
            obs["conf"] = obs["corr_subclass_margin"].astype(float)
        else:
            obs["conf"] = obs["subclass_label_confidence"].astype(float)
        dfs.append(obs[["depth", "layer", "conf"]])
    df = pd.concat(dfs, ignore_index=True)
    print(f"  Xenium L6b (QC-pass): {len(df):,}")
    return df


def assign_layer_from_depth(depth):
    """Convert continuous depth to discrete layer."""
    if depth < 0.1225:
        return "L1"
    elif depth < 0.4696:
        return "L2/3"
    elif depth < 0.5443:
        return "L4"
    elif depth < 0.7079:
        return "L5"
    elif depth < 0.9275:
        return "L6"
    else:
        return "WM"


def plot_row(axes, merfish, xenium, thresh, row_label):
    """Plot one row: depth histogram (left) + layer bar chart (right)."""
    ax_hist, ax_bar = axes

    # Filter Xenium by threshold
    xen = xenium[xenium["conf"] >= thresh].copy()

    if len(xen) == 0:
        for ax in axes:
            ax.text(0.5, 0.5, f"No cells pass threshold {thresh}",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=14, color="#999")
            ax.set_title(f"L6b — {row_label}", fontsize=15, fontweight="bold")
        return {
            "thresh": thresh, "n_merfish": len(merfish), "n_xenium": 0,
            "merfish_median": merfish["depth"].median(), "xenium_median": np.nan,
            "merfish_std": merfish["depth"].std(), "xenium_std": np.nan,
            "upper_merfish_pct": 0, "upper_xenium_pct": 0,
        }

    # --- Left: depth histogram ---
    bins = np.linspace(0, 1, 60)
    ax_hist.hist(merfish["depth"], bins=bins, density=True, alpha=0.55,
                 color="#4fc3f7", label=f"MERFISH (n={len(merfish):,})",
                 edgecolor="none")
    ax_hist.hist(xen["depth"], bins=bins, density=True, alpha=0.55,
                 color="#ef5350", label=f"Xenium (n={len(xen):,})",
                 edgecolor="none")

    for b in LAYER_BOUNDS:
        ax_hist.axvline(b, color="#888888", linestyle="--", linewidth=0.8, alpha=0.6)

    ax_hist.set_xlabel("Normalized depth from pia (0=pia, 1=WM)", fontsize=13)
    ax_hist.set_ylabel("Density", fontsize=13)
    ax_hist.set_title(f"L6b Depth Distribution\n{row_label}",
                      fontsize=15, fontweight="bold")
    ax_hist.legend(fontsize=12, loc="upper left")
    ax_hist.set_xlim(0, 1.05)

    # Stats annotation
    merfish_median = merfish["depth"].median()
    xen_median = xen["depth"].median()
    merfish_std = merfish["depth"].std()
    xen_std = xen["depth"].std()
    ax_hist.text(0.97, 0.97,
                 f"MERFISH: med={merfish_median:.3f}, sd={merfish_std:.3f}\n"
                 f"Xenium: med={xen_median:.3f}, sd={xen_std:.3f}",
                 transform=ax_hist.transAxes, ha="right", va="top",
                 fontsize=10, bbox=dict(boxstyle="round,pad=0.3",
                                        facecolor="white", alpha=0.8))

    # --- Right: layer assignment bar chart ---
    # MERFISH layers from manual annotation
    merfish_layers = merfish["layer"].value_counts()
    # Xenium layers from predicted layer column
    xen_layers = xen["layer"].value_counts()

    layers = LAYER_NAMES
    x_pos = np.arange(len(layers))
    width = 0.35

    merfish_pcts = [merfish_layers.get(l, 0) / len(merfish) * 100 for l in layers]
    xen_pcts = [xen_layers.get(l, 0) / len(xen) * 100 for l in layers]

    bars1 = ax_bar.bar(x_pos - width / 2, merfish_pcts, width,
                        color="#4fc3f7", label="MERFISH", edgecolor="white",
                        linewidth=0.5)
    bars2 = ax_bar.bar(x_pos + width / 2, xen_pcts, width,
                        color="#ef5350", label="Xenium", edgecolor="white",
                        linewidth=0.5)

    # Add percentage labels
    for bar, pct in zip(bars1, merfish_pcts):
        if pct > 1:
            ax_bar.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                        f"{pct:.0f}%", ha="center", va="bottom", fontsize=9,
                        color="#4fc3f7", fontweight="bold")
    for bar, pct in zip(bars2, xen_pcts):
        if pct > 1:
            ax_bar.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                        f"{pct:.0f}%", ha="center", va="bottom", fontsize=9,
                        color="#ef5350", fontweight="bold")

    ax_bar.set_xticks(x_pos)
    ax_bar.set_xticklabels(layers, fontsize=12)
    ax_bar.set_ylabel("% of L6b cells", fontsize=13)
    ax_bar.set_title(f"L6b Layer Assignment\n{row_label}",
                     fontsize=15, fontweight="bold")
    ax_bar.legend(fontsize=12)

    # Compute upper-cortex fraction
    upper_merfish = sum(merfish_pcts[:3])  # L1 + L2/3 + L4
    upper_xen = sum(xen_pcts[:3])

    return {
        "thresh": thresh,
        "n_merfish": len(merfish),
        "n_xenium": len(xen),
        "merfish_median": merfish_median,
        "xenium_median": xen_median,
        "merfish_std": merfish_std,
        "xenium_std": xen_std,
        "upper_merfish_pct": upper_merfish,
        "upper_xenium_pct": upper_xen,
    }


def main():
    merfish = load_merfish_l6b()
    xenium = load_xenium_l6b()

    fig, axes = plt.subplots(2, 2, figsize=(18, 14))

    stats_028 = plot_row(axes[0], merfish, xenium, 0.28,
                          "Confidence threshold = 0.28 (current)")
    stats_050 = plot_row(axes[1], merfish, xenium, 0.50,
                          "Confidence threshold = 0.50 (stringent)")

    # Summary suptitle
    improvement = stats_028["upper_xenium_pct"] - stats_050["upper_xenium_pct"]
    cell_loss = stats_028["n_xenium"] - stats_050["n_xenium"]
    cell_loss_pct = (cell_loss / stats_028["n_xenium"] * 100
                     if stats_028["n_xenium"] > 0 else 0)

    if stats_028["n_xenium"] > 0:
        fig.suptitle(
            f"L6b depth: Xenium vs MERFISH — effect of confidence filtering\n"
            f"Stringent filter reduces upper-cortex L6b from "
            f"{stats_028['upper_xenium_pct']:.1f}% → {stats_050['upper_xenium_pct']:.1f}% "
            f"(MERFISH: {stats_028['upper_merfish_pct']:.1f}%) | "
            f"Cost: {cell_loss:,} cells lost ({cell_loss_pct:.0f}%)",
            fontsize=16, fontweight="bold", color="#cc0000", y=1.02)
    else:
        fig.suptitle(
            "L6b depth: Xenium vs MERFISH — confidence filtering\n"
            "(No Xenium cells pass either threshold — likely using correlation margins)",
            fontsize=16, fontweight="bold", color="#cc0000", y=1.02)

    plt.tight_layout()

    outpath = os.path.join(OUT_DIR, "l6b_depth_xenium_vs_merfish.png")
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {outpath}")

    # Print summary
    print("\n── Summary ──")
    for s in [stats_028, stats_050]:
        print(f"\n  Threshold = {s['thresh']:.2f}:")
        print(f"    Xenium L6b cells: {s['n_xenium']:,}")
        if s['n_xenium'] > 0:
            print(f"    Xenium median depth: {s['xenium_median']:.3f} "
                  f"(MERFISH: {s['merfish_median']:.3f})")
            print(f"    Xenium depth std: {s['xenium_std']:.3f} "
                  f"(MERFISH: {s['merfish_std']:.3f})")
            print(f"    Upper cortex (L1-L4): Xenium {s['upper_xenium_pct']:.1f}% "
                  f"vs MERFISH {s['upper_merfish_pct']:.1f}%")
        else:
            print("    (No cells pass this threshold)")


if __name__ == "__main__":
    main()
