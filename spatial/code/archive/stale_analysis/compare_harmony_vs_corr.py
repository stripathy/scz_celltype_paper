#!/usr/bin/env python3
"""
Compare Harmony label transfer vs Correlation Classifier vs HANN.

Generates:
  1. Subclass agreement heatmap (confusion matrix style)
  2. Per-type proportion comparison (Harmony vs Corr vs HANN)
  3. Spatial coherence check (layer distributions by method)
  4. Confidence distribution comparison
  5. Summary CSV with per-type agreement rates

Usage:
    python3 -u compare_harmony_vs_corr.py
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    H5AD_DIR, BASE_DIR, SAMPLE_TO_DX, EXCLUDE_SAMPLES,
    SUBCLASS_TO_CLASS, PRESENTATION_DIR, load_cells,
)

OUT_DIR = os.path.join(BASE_DIR, "output", "plots", "harmony_validation")
os.makedirs(OUT_DIR, exist_ok=True)


def load_harmony_data():
    """Load all cells with both harmony and corr labels."""
    import anndata as ad

    h5ad_files = sorted(
        f for f in os.listdir(H5AD_DIR)
        if f.endswith("_annotated.h5ad")
    )

    dfs = []
    for fname in h5ad_files:
        sample_id = fname.replace("_annotated.h5ad", "")
        if sample_id in EXCLUDE_SAMPLES:
            continue
        adata = ad.read_h5ad(os.path.join(H5AD_DIR, fname))
        obs = adata.obs.copy()
        obs["sample_id"] = sample_id
        # Only keep QC-pass cells
        if "qc_pass" in obs.columns:
            obs = obs[obs["qc_pass"] == True]
        dfs.append(obs)

    df = pd.concat(dfs, ignore_index=True)
    print(f"Loaded {len(df):,} QC-pass cells from {len(dfs)} samples")

    # Check which columns exist
    for col in ["harmony_subclass", "corr_subclass", "subclass_label",
                "harmony_supertype", "corr_supertype", "supertype_label",
                "harmony_subclass_confidence", "harmony_supertype_confidence",
                "corr_subclass_corr", "subclass_label_confidence",
                "predicted_norm_depth", "layer"]:
        has = col in df.columns and df[col].notna().sum() > 0
        print(f"  {col}: {'present' if has else 'MISSING'}")

    return df


def plot_subclass_confusion(df):
    """Confusion matrix: Harmony subclass vs Corr subclass."""
    print("\n=== Subclass Confusion Matrix ===")
    valid = df["harmony_subclass"].notna() & df["corr_subclass"].notna()
    sub = df.loc[valid, ["harmony_subclass", "corr_subclass"]].copy()

    # Build confusion matrix
    all_types = sorted(set(sub["harmony_subclass"]) | set(sub["corr_subclass"]))
    ct = pd.crosstab(sub["corr_subclass"], sub["harmony_subclass"])
    ct = ct.reindex(index=all_types, columns=all_types, fill_value=0)

    # Normalize by row (corr classifier counts)
    ct_norm = ct.div(ct.sum(axis=1), axis=0)

    fig, ax = plt.subplots(figsize=(16, 14))
    im = ax.imshow(ct_norm.values, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(all_types)))
    ax.set_yticks(range(len(all_types)))
    ax.set_xticklabels(all_types, rotation=90, fontsize=11)
    ax.set_yticklabels(all_types, fontsize=11)
    ax.set_xlabel("Harmony Subclass", fontsize=14)
    ax.set_ylabel("Correlation Classifier Subclass", fontsize=14)
    ax.set_title("Subclass Agreement: Correlation Classifier vs Harmony\n"
                 "(Row-normalized: fraction of Corr cells assigned to each Harmony type)",
                 fontsize=16)

    # Add diagonal agreement values
    for i in range(len(all_types)):
        val = ct_norm.values[i, i]
        color = "white" if val > 0.5 else "black"
        ax.text(i, i, f"{val:.2f}", ha="center", va="center",
                fontsize=9, fontweight="bold", color=color)

    plt.colorbar(im, ax=ax, label="Fraction", shrink=0.8)
    plt.tight_layout()
    path = os.path.join(OUT_DIR, "subclass_confusion_corr_vs_harmony.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {path}")

    # Print per-type agreement
    diag = pd.Series(np.diag(ct_norm.values), index=all_types)
    print("\n  Per-type agreement (fraction of Corr cells matching Harmony):")
    for t in diag.sort_values().index:
        n = ct.loc[t].sum()
        print(f"    {t:20s}: {diag[t]:.3f}  (n={n:,})")

    return ct, ct_norm


def plot_proportion_comparison(df):
    """Compare subclass proportions across methods."""
    print("\n=== Subclass Proportion Comparison ===")

    methods = {}
    for col, name in [("harmony_subclass", "Harmony"),
                       ("corr_subclass", "Corr Classifier"),
                       ("subclass_label", "HANN")]:
        if col in df.columns:
            valid = df[col].notna()
            if valid.sum() > 0:
                methods[name] = df.loc[valid, col].value_counts(normalize=True)

    if len(methods) < 2:
        print("  Not enough methods with data, skipping")
        return

    prop_df = pd.DataFrame(methods).fillna(0)
    prop_df = prop_df.sort_values(list(methods.keys())[0], ascending=True)

    fig, ax = plt.subplots(figsize=(12, 10))
    x = np.arange(len(prop_df))
    width = 0.8 / len(methods)

    colors = {"Harmony": "#2196F3", "Corr Classifier": "#FF9800", "HANN": "#4CAF50"}
    for i, (name, props) in enumerate(methods.items()):
        vals = prop_df[name].values
        ax.barh(x + i * width, vals, width, label=name,
                color=colors.get(name, f"C{i}"), alpha=0.8)

    ax.set_yticks(x + width * (len(methods) - 1) / 2)
    ax.set_yticklabels(prop_df.index, fontsize=11)
    ax.set_xlabel("Proportion of Cells", fontsize=14)
    ax.set_title("Subclass Proportions by Classification Method", fontsize=16)
    ax.legend(fontsize=12)
    plt.tight_layout()

    path = os.path.join(OUT_DIR, "subclass_proportions_by_method.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {path}")

    # Print notable differences
    print("\n  Largest proportion differences (Harmony vs Corr):")
    if "Harmony" in prop_df.columns and "Corr Classifier" in prop_df.columns:
        diff = (prop_df["Harmony"] - prop_df["Corr Classifier"]).sort_values()
        for t in list(diff.index[:5]) + list(diff.index[-5:]):
            print(f"    {t:20s}: Harmony={prop_df.loc[t, 'Harmony']:.3f}  "
                  f"Corr={prop_df.loc[t, 'Corr Classifier']:.3f}  "
                  f"diff={diff[t]:+.3f}")


def plot_spatial_coherence(df):
    """Check whether Harmony labels respect cortical depth."""
    print("\n=== Spatial Coherence Check ===")

    if "predicted_norm_depth" not in df.columns:
        print("  No depth predictions, skipping")
        return

    valid = df["predicted_norm_depth"].notna() & df["harmony_subclass"].notna()
    sub = df.loc[valid].copy()

    # Expected depth ordering for excitatory types
    exc_types = ["L2/3 IT", "L4 IT", "L5 IT", "L5 ET", "L5/6 NP",
                 "L6 CT", "L6 IT", "L6 IT Car3", "L6b"]
    exc_present = [t for t in exc_types if t in sub["harmony_subclass"].values]

    if not exc_present:
        print("  No excitatory types found, skipping")
        return

    fig, axes = plt.subplots(1, 2, figsize=(20, 8))

    for ax, (col, title) in zip(axes, [
        ("harmony_subclass", "Harmony"),
        ("corr_subclass", "Corr Classifier"),
    ]):
        if col not in sub.columns:
            ax.set_visible(False)
            continue

        data = []
        labels = []
        for t in exc_present:
            mask = sub[col] == t
            if mask.sum() > 0:
                depths = sub.loc[mask, "predicted_norm_depth"].values
                data.append(depths)
                labels.append(f"{t}\n(n={mask.sum():,})")

        if data:
            bp = ax.boxplot(data, labels=labels, vert=True, patch_artist=True,
                           showfliers=False)
            for patch in bp["boxes"]:
                patch.set_facecolor("#90CAF9" if "Harmony" in title else "#FFCC80")
            ax.set_ylabel("Predicted Normalized Depth", fontsize=14)
            ax.set_title(f"{title}: Depth Distribution by Excitatory Subclass",
                        fontsize=14)
            ax.tick_params(axis="x", rotation=45, labelsize=10)
            # Depth should increase from L2/3 → L6b
            ax.set_ylim(-0.05, 1.05)

    plt.tight_layout()
    path = os.path.join(OUT_DIR, "spatial_coherence_depth_by_subclass.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {path}")

    # Compute median depth per type for each method
    print("\n  Median depth by excitatory type:")
    print(f"  {'Type':15s} {'Harmony':>10s} {'Corr':>10s} {'Expected':>10s}")
    for i, t in enumerate(exc_present):
        expected = (i + 0.5) / len(exc_present)
        harm_depth = sub.loc[sub["harmony_subclass"] == t, "predicted_norm_depth"].median()
        corr_depth = np.nan
        if "corr_subclass" in sub.columns:
            corr_mask = sub["corr_subclass"] == t
            if corr_mask.sum() > 0:
                corr_depth = sub.loc[corr_mask, "predicted_norm_depth"].median()
        print(f"  {t:15s} {harm_depth:10.3f} {corr_depth:10.3f} {expected:10.3f}")


def plot_confidence_comparison(df):
    """Compare confidence distributions across methods."""
    print("\n=== Confidence Distributions ===")

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Subclass confidence
    ax = axes[0]
    if "harmony_subclass_confidence" in df.columns:
        vals = df["harmony_subclass_confidence"].dropna()
        ax.hist(vals, bins=50, alpha=0.6, label=f"Harmony (med={vals.median():.3f})",
                color="#2196F3")
    if "corr_subclass_corr" in df.columns:
        vals = df["corr_subclass_corr"].dropna()
        ax.hist(vals, bins=50, alpha=0.6, label=f"Corr (med={vals.median():.3f})",
                color="#FF9800")
    ax.set_xlabel("Confidence / Correlation Score", fontsize=14)
    ax.set_ylabel("Count", fontsize=14)
    ax.set_title("Subclass Confidence Distribution", fontsize=16)
    ax.legend(fontsize=11)

    # Supertype confidence
    ax = axes[1]
    if "harmony_supertype_confidence" in df.columns:
        vals = df["harmony_supertype_confidence"].dropna()
        ax.hist(vals, bins=50, alpha=0.6, label=f"Harmony (med={vals.median():.3f})",
                color="#2196F3")
    if "corr_supertype_corr" in df.columns:
        vals = df["corr_supertype_corr"].dropna()
        ax.hist(vals, bins=50, alpha=0.6, label=f"Corr (med={vals.median():.3f})",
                color="#FF9800")
    ax.set_xlabel("Confidence / Correlation Score", fontsize=14)
    ax.set_ylabel("Count", fontsize=14)
    ax.set_title("Supertype Confidence Distribution", fontsize=16)
    ax.legend(fontsize=11)

    plt.tight_layout()
    path = os.path.join(OUT_DIR, "confidence_distributions.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {path}")


def plot_per_sample_agreement(df):
    """Per-sample agreement rates between Harmony and Corr."""
    print("\n=== Per-Sample Agreement ===")

    valid = df["harmony_subclass"].notna() & df["corr_subclass"].notna()
    sub = df.loc[valid].copy()
    sub["agree"] = sub["harmony_subclass"] == sub["corr_subclass"]

    by_sample = sub.groupby("sample_id")["agree"].mean().sort_values()

    fig, ax = plt.subplots(figsize=(14, 6))
    colors = [("#FF9800" if SAMPLE_TO_DX.get(s) == "SCZ" else "#4CAF50")
              for s in by_sample.index]
    ax.bar(range(len(by_sample)), by_sample.values, color=colors, alpha=0.8)
    ax.set_xticks(range(len(by_sample)))
    ax.set_xticklabels(by_sample.index, rotation=45, ha="right", fontsize=11)
    ax.set_ylabel("Fraction Agreement (Subclass)", fontsize=14)
    ax.set_title("Harmony vs Corr Classifier Agreement by Sample\n"
                 "(Orange = SCZ, Green = Control)", fontsize=16)
    ax.axhline(by_sample.mean(), ls="--", color="gray", lw=1)
    ax.text(len(by_sample) - 1, by_sample.mean() + 0.01,
            f"mean={by_sample.mean():.3f}", fontsize=11, color="gray")
    ax.set_ylim(0, 1)
    plt.tight_layout()

    path = os.path.join(OUT_DIR, "per_sample_agreement.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {path}")

    for s, v in by_sample.items():
        dx = SAMPLE_TO_DX.get(s, "?")
        print(f"  {s} ({dx:7s}): {v:.3f}")


def save_summary_csv(df):
    """Save agreement summary to CSV."""
    print("\n=== Saving Summary ===")

    records = []
    valid = df["harmony_subclass"].notna() & df["corr_subclass"].notna()
    sub = df.loc[valid]

    for t in sorted(sub["corr_subclass"].unique()):
        mask = sub["corr_subclass"] == t
        n_corr = mask.sum()
        harm_labels = sub.loc[mask, "harmony_subclass"]
        agree = (harm_labels == t).mean()
        top_harm = harm_labels.value_counts().head(3)
        records.append({
            "subclass": t,
            "n_corr": n_corr,
            "agree_rate": agree,
            "top1_harmony": top_harm.index[0] if len(top_harm) > 0 else "",
            "top1_frac": top_harm.values[0] / n_corr if len(top_harm) > 0 else 0,
            "top2_harmony": top_harm.index[1] if len(top_harm) > 1 else "",
            "top2_frac": top_harm.values[1] / n_corr if len(top_harm) > 1 else 0,
        })

    summary = pd.DataFrame(records)
    path = os.path.join(PRESENTATION_DIR, "harmony_vs_corr_agreement.csv")
    summary.to_csv(path, index=False)
    print(f"  Saved: {path}")
    print(summary.to_string(index=False))


def main():
    print("=" * 70)
    print("Harmony vs Correlation Classifier Comparison")
    print("=" * 70)

    df = load_harmony_data()

    ct, ct_norm = plot_subclass_confusion(df)
    plot_proportion_comparison(df)
    plot_spatial_coherence(df)
    plot_confidence_comparison(df)
    plot_per_sample_agreement(df)
    save_summary_csv(df)

    # Overall summary
    print("\n" + "=" * 70)
    print("OVERALL SUMMARY")
    print("=" * 70)
    valid = df["harmony_subclass"].notna() & df["corr_subclass"].notna()
    agree = (df.loc[valid, "harmony_subclass"] == df.loc[valid, "corr_subclass"]).mean()
    print(f"  Subclass agreement (Harmony vs Corr): {100*agree:.1f}%")

    if "subclass_label" in df.columns:
        valid2 = df["harmony_subclass"].notna() & df["subclass_label"].notna()
        harm_str = df.loc[valid2, "harmony_subclass"].astype(str).values
        hann_str = df.loc[valid2, "subclass_label"].astype(str).values
        agree2 = (harm_str == hann_str).mean()
        print(f"  Subclass agreement (Harmony vs HANN): {100*agree2:.1f}%")

    # Endothelial inflation check
    harm_endo = (df["harmony_subclass"] == "Endothelial").mean()
    corr_endo = (df.loc[df["corr_subclass"].notna(), "corr_subclass"] == "Endothelial").mean()
    print(f"\n  Endothelial proportion: Harmony={100*harm_endo:.1f}%, Corr={100*corr_endo:.1f}%")

    print(f"\nAll figures saved to: {OUT_DIR}")


if __name__ == "__main__":
    main()
