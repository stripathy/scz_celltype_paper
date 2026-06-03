#!/usr/bin/env python3
"""
Compare proportions of upper-layer vs deep-layer excitatory cell types
between Xenium and MERFISH. This is a depth-model-independent check:
if Xenium has inflated deep-layer types (L6b, L6 CT, L6 IT, L5 ET) and
deflated upper-layer types (L2/3 IT, L4 IT), it suggests label transfer
bias rather than depth model issues.

Uses all cells (not restricted to cortical layers) to avoid circularity
with the depth model.
"""

import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    H5AD_DIR, MERFISH_PATH, SUBCLASS_TO_CLASS, CLASS_COLORS,
    SUBCLASS_CONF_THRESH,
)

BASE_DIR = os.path.expanduser("~/Github/SCZ_Xenium")
OUT_DIR = os.path.join(BASE_DIR, "output", "plots")

# Group excitatory subclasses by expected laminar position
UPPER_EXCITATORY = ["L2/3 IT", "L4 IT"]
DEEP_EXCITATORY = ["L5 IT", "L5 ET", "L5/6 NP", "L6 IT", "L6 IT Car3", "L6 CT", "L6b"]
ALL_EXCITATORY = UPPER_EXCITATORY + DEEP_EXCITATORY


def load_merfish_proportions():
    """Per-donor subclass proportions from MERFISH."""
    print("Loading MERFISH...")
    adata = ad.read_h5ad(MERFISH_PATH, backed="r")
    obs = adata.obs[["Donor ID", "Subclass"]].copy()
    obs.columns = ["donor", "subclass"]
    obs["subclass"] = obs["subclass"].astype(str)

    records = []
    for donor in obs["donor"].unique():
        d = obs[obs["donor"] == donor]
        total = len(d)
        counts = d["subclass"].value_counts()
        for ct, n in counts.items():
            records.append({"donor": str(donor), "subclass": ct,
                            "proportion": n / total, "n": n, "total": total})
    return pd.DataFrame(records)


def load_xenium_proportions(conf_thresh=SUBCLASS_CONF_THRESH):
    """Per-sample subclass proportions from Xenium (all QC-pass cells)."""
    print(f"Loading Xenium (conf >= {conf_thresh})...")
    h5ad_files = sorted([f for f in os.listdir(H5AD_DIR)
                         if f.endswith("_annotated.h5ad")])
    records = []
    for fname in h5ad_files:
        sample_id = fname.replace("_annotated.h5ad", "")
        fpath = os.path.join(H5AD_DIR, fname)
        adata = ad.read_h5ad(fpath, backed="r")
        obs = adata.obs[["subclass_label", "qc_pass",
                          "subclass_label_confidence"]].copy()
        obs = obs[obs["qc_pass"] == True]
        obs = obs[obs["subclass_label_confidence"].astype(float) >= conf_thresh]
        total = len(obs)
        counts = obs["subclass_label"].value_counts()
        for ct, n in counts.items():
            records.append({"donor": sample_id, "subclass": str(ct),
                            "proportion": n / total, "n": n, "total": total})
    return pd.DataFrame(records)


def summarize(df, label):
    """Compute mean proportion per subclass across donors."""
    mean = df.groupby("subclass")["proportion"].mean()
    print(f"\n── {label} ──")
    print(f"  Upper excitatory (L2/3 IT + L4 IT):")
    for s in UPPER_EXCITATORY:
        print(f"    {s:12s}: {mean.get(s, 0)*100:5.2f}%")
    print(f"    TOTAL:       {sum(mean.get(s, 0) for s in UPPER_EXCITATORY)*100:5.2f}%")
    print(f"  Deep excitatory (L5-L6):")
    for s in DEEP_EXCITATORY:
        print(f"    {s:12s}: {mean.get(s, 0)*100:5.2f}%")
    print(f"    TOTAL:       {sum(mean.get(s, 0) for s in DEEP_EXCITATORY)*100:5.2f}%")
    return mean


def main():
    merfish = load_merfish_proportions()
    xenium_028 = load_xenium_proportions(0.28)
    xenium_050 = load_xenium_proportions(0.50)

    merfish_mean = summarize(merfish, "MERFISH (all donors)")
    xenium_028_mean = summarize(xenium_028, "Xenium (thresh=0.28)")
    xenium_050_mean = summarize(xenium_050, "Xenium (thresh=0.50)")

    # ── Figure ──
    fig, axes = plt.subplots(1, 3, figsize=(22, 8))

    # Panel 1: All excitatory subclass proportions side by side
    ax = axes[0]
    subclasses = ALL_EXCITATORY
    x = np.arange(len(subclasses))
    w = 0.25

    m_vals = [merfish_mean.get(s, 0) * 100 for s in subclasses]
    x028_vals = [xenium_028_mean.get(s, 0) * 100 for s in subclasses]
    x050_vals = [xenium_050_mean.get(s, 0) * 100 for s in subclasses]

    bars_m = ax.bar(x - w, m_vals, w, color="#4fc3f7", label="MERFISH",
                     edgecolor="white", linewidth=0.5)
    bars_028 = ax.bar(x, x028_vals, w, color="#ef5350", alpha=0.6,
                       label="Xenium (0.28)", edgecolor="white", linewidth=0.5)
    bars_050 = ax.bar(x + w, x050_vals, w, color="#ef5350",
                       label="Xenium (0.50)", edgecolor="white", linewidth=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(subclasses, rotation=45, ha="right", fontsize=11)
    ax.set_ylabel("Mean proportion across donors (%)", fontsize=13)
    ax.set_title("Excitatory subclass proportions", fontsize=16, fontweight="bold")
    ax.legend(fontsize=11)

    # Add dividing line between upper and deep
    ax.axvline(len(UPPER_EXCITATORY) - 0.5, color="#888888", linestyle="--",
               linewidth=1, alpha=0.7)
    ax.text(0.5, 0.95, "Upper", transform=ax.transAxes, ha="left", va="top",
            fontsize=10, color="#888888", style="italic")
    ax.text(len(UPPER_EXCITATORY) / len(subclasses) + 0.05, 0.95, "Deep",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=10, color="#888888", style="italic")

    # Panel 2: Upper vs deep aggregate proportions
    ax = axes[1]
    upper_m = sum(merfish_mean.get(s, 0) for s in UPPER_EXCITATORY) * 100
    upper_028 = sum(xenium_028_mean.get(s, 0) for s in UPPER_EXCITATORY) * 100
    upper_050 = sum(xenium_050_mean.get(s, 0) for s in UPPER_EXCITATORY) * 100
    deep_m = sum(merfish_mean.get(s, 0) for s in DEEP_EXCITATORY) * 100
    deep_028 = sum(xenium_028_mean.get(s, 0) for s in DEEP_EXCITATORY) * 100
    deep_050 = sum(xenium_050_mean.get(s, 0) for s in DEEP_EXCITATORY) * 100

    groups = ["Upper excitatory\n(L2/3 IT + L4 IT)", "Deep excitatory\n(L5-L6b)"]
    x2 = np.arange(2)
    w2 = 0.25
    ax.bar(x2 - w2, [upper_m, deep_m], w2, color="#4fc3f7", label="MERFISH",
           edgecolor="white", linewidth=0.5)
    ax.bar(x2, [upper_028, deep_028], w2, color="#ef5350", alpha=0.6,
           label="Xenium (0.28)", edgecolor="white", linewidth=0.5)
    ax.bar(x2 + w2, [upper_050, deep_050], w2, color="#ef5350",
           label="Xenium (0.50)", edgecolor="white", linewidth=0.5)

    # Add value labels
    for i, (vm, v028, v050) in enumerate([(upper_m, upper_028, upper_050),
                                           (deep_m, deep_028, deep_050)]):
        ax.text(i - w2, vm + 0.3, f"{vm:.1f}%", ha="center", va="bottom",
                fontsize=10, color="#4fc3f7", fontweight="bold")
        ax.text(i, v028 + 0.3, f"{v028:.1f}%", ha="center", va="bottom",
                fontsize=10, color="#ef5350", fontweight="bold", alpha=0.6)
        ax.text(i + w2, v050 + 0.3, f"{v050:.1f}%", ha="center", va="bottom",
                fontsize=10, color="#ef5350", fontweight="bold")

    ax.set_xticks(x2)
    ax.set_xticklabels(groups, fontsize=12)
    ax.set_ylabel("Mean proportion across donors (%)", fontsize=13)
    ax.set_title("Upper vs deep excitatory\n(aggregate)", fontsize=16,
                 fontweight="bold")
    ax.legend(fontsize=11)

    # Panel 3: Per-donor distributions (boxplot)
    ax = axes[2]

    # Compute per-donor upper and deep fractions
    def get_donor_layer_fracs(df):
        records = []
        for donor in df["donor"].unique():
            d = df[df["donor"] == donor]
            total_prop = d.set_index("subclass")["proportion"]
            upper = sum(total_prop.get(s, 0) for s in UPPER_EXCITATORY) * 100
            deep = sum(total_prop.get(s, 0) for s in DEEP_EXCITATORY) * 100
            ratio = deep / upper if upper > 0 else np.nan
            records.append({"donor": donor, "upper": upper, "deep": deep,
                            "deep_to_upper_ratio": ratio})
        return pd.DataFrame(records)

    m_fracs = get_donor_layer_fracs(merfish)
    x028_fracs = get_donor_layer_fracs(xenium_028)
    x050_fracs = get_donor_layer_fracs(xenium_050)

    positions = [0, 1, 2, 4, 5, 6]
    data = [m_fracs["deep_to_upper_ratio"], x028_fracs["deep_to_upper_ratio"],
            x050_fracs["deep_to_upper_ratio"],
            m_fracs["upper"], x028_fracs["upper"], x050_fracs["upper"]]
    colors = ["#4fc3f7", "#ef5350", "#ef5350", "#4fc3f7", "#ef5350", "#ef5350"]
    alphas = [1.0, 0.5, 1.0, 1.0, 0.5, 1.0]

    bp = ax.boxplot(data, positions=positions, widths=0.6, patch_artist=True,
                     showfliers=True, flierprops=dict(markersize=4))
    for patch, color, alpha in zip(bp['boxes'], colors, alphas):
        patch.set_facecolor(color)
        patch.set_alpha(alpha)
    for element in ['whiskers', 'caps', 'medians']:
        for line in bp[element]:
            line.set_color('#333333')

    ax.set_xticks([1, 5])
    ax.set_xticklabels(["Deep/Upper ratio", "Upper excitatory %"], fontsize=12)
    ax.set_title("Per-donor distributions", fontsize=16, fontweight="bold")

    # Manual legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#4fc3f7", label="MERFISH"),
        Patch(facecolor="#ef5350", alpha=0.5, label="Xenium (0.28)"),
        Patch(facecolor="#ef5350", label="Xenium (0.50)"),
    ]
    ax.legend(handles=legend_elements, fontsize=11)

    # Print ratio stats
    print(f"\n── Deep/Upper excitatory ratio ──")
    print(f"  MERFISH:       {m_fracs['deep_to_upper_ratio'].median():.2f} "
          f"(median, IQR: {m_fracs['deep_to_upper_ratio'].quantile(0.25):.2f}-"
          f"{m_fracs['deep_to_upper_ratio'].quantile(0.75):.2f})")
    print(f"  Xenium (0.28): {x028_fracs['deep_to_upper_ratio'].median():.2f} "
          f"(median, IQR: {x028_fracs['deep_to_upper_ratio'].quantile(0.25):.2f}-"
          f"{x028_fracs['deep_to_upper_ratio'].quantile(0.75):.2f})")
    print(f"  Xenium (0.50): {x050_fracs['deep_to_upper_ratio'].median():.2f} "
          f"(median, IQR: {x050_fracs['deep_to_upper_ratio'].quantile(0.25):.2f}-"
          f"{x050_fracs['deep_to_upper_ratio'].quantile(0.75):.2f})")

    fig.suptitle("Are deep-layer excitatory types over-represented in Xenium vs MERFISH?",
                 fontsize=18, fontweight="bold", y=1.02)
    plt.tight_layout()

    outpath = os.path.join(OUT_DIR, "upper_vs_deep_excitatory_proportions.png")
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {outpath}")


if __name__ == "__main__":
    main()
