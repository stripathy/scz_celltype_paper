#!/usr/bin/env python3
"""
Compare SCZ vs Control results between corr and hybrid QC modes.

Side-by-side comparison of crumblr compositional analysis results,
cell count changes, and (if available) DE results.

Generates:
  output/presentation/qc_mode_comparison_crumblr.png
  output/presentation/qc_mode_comparison_counts.png
  output/presentation/qc_mode_comparison_de.png (if DE results exist)
  output/presentation/qc_mode_comparison_summary.csv

Usage:
    python3 -u compare_qc_modes.py
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    CRUMBLR_DIR, PRESENTATION_DIR, BG_COLOR, DX_COLORS,
    SUBCLASS_TO_CLASS, CLASS_COLORS,
)

OUT_DIR = PRESENTATION_DIR
DE_DIR = os.path.join(os.path.dirname(PRESENTATION_DIR), "de")
BG = BG_COLOR


# ──────────────────────────────────────────────────────────────────────
# 1. Crumblr comparison at subclass level
# ──────────────────────────────────────────────────────────────────────

def compare_crumblr_subclass():
    """Compare crumblr results between corr and hybrid QC at subclass level."""

    corr_path = os.path.join(CRUMBLR_DIR, "crumblr_results_subclass.csv")
    hybrid_path = os.path.join(CRUMBLR_DIR, "crumblr_results_subclass_hybrid.csv")

    if not os.path.exists(corr_path) or not os.path.exists(hybrid_path):
        print("Missing crumblr subclass results — skipping")
        return None

    corr = pd.read_csv(corr_path)
    hybrid = pd.read_csv(hybrid_path)

    # Merge on celltype
    merged = corr.merge(hybrid, on="celltype", suffixes=("_corr", "_hybrid"))

    # Add class info
    merged["class"] = merged["celltype"].map(
        lambda x: SUBCLASS_TO_CLASS.get(x, "Other"))

    return merged


def compare_crumblr_supertype():
    """Compare crumblr results between corr and hybrid QC at supertype level."""

    corr_path = os.path.join(CRUMBLR_DIR, "crumblr_results_supertype.csv")
    hybrid_path = os.path.join(CRUMBLR_DIR, "crumblr_results_supertype_hybrid.csv")

    if not os.path.exists(corr_path) or not os.path.exists(hybrid_path):
        print("Missing crumblr supertype results — skipping")
        return None

    corr = pd.read_csv(corr_path)
    hybrid = pd.read_csv(hybrid_path)

    merged = corr.merge(hybrid, on="celltype", suffixes=("_corr", "_hybrid"))
    return merged


def plot_crumblr_comparison(sub_merged, sup_merged):
    """Plot side-by-side crumblr logFC comparison."""

    fig = plt.figure(figsize=(24, 14), facecolor=BG)
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.3)

    # ── Panel A: Subclass logFC scatter ──
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_facecolor(BG)

    for cls, color in CLASS_COLORS.items():
        mask = sub_merged["class"] == cls
        if mask.sum() == 0:
            continue
        ax1.scatter(sub_merged.loc[mask, "logFC_corr"],
                    sub_merged.loc[mask, "logFC_hybrid"],
                    c=color, s=100, alpha=0.9, edgecolors="white",
                    linewidths=0.5, label=cls, zorder=5)

    # Label significant points
    for _, row in sub_merged.iterrows():
        if row["FDR_corr"] < 0.1 or row["FDR_hybrid"] < 0.1:
            ax1.annotate(row["celltype"],
                         (row["logFC_corr"], row["logFC_hybrid"]),
                         fontsize=10, color="white", fontweight="bold",
                         xytext=(5, 5), textcoords="offset points")

    lim = max(abs(sub_merged["logFC_corr"].max()),
              abs(sub_merged["logFC_hybrid"].max())) * 1.1
    ax1.plot([-lim, lim], [-lim, lim], "w--", alpha=0.3, linewidth=1)
    ax1.axhline(0, color="#555555", linewidth=0.5)
    ax1.axvline(0, color="#555555", linewidth=0.5)
    ax1.set_xlabel("logFC (corr QC)", fontsize=14, color="white")
    ax1.set_ylabel("logFC (hybrid QC)", fontsize=14, color="white")
    ax1.set_title("Subclass: logFC comparison", fontsize=18,
                  fontweight="bold", color="white")
    ax1.legend(fontsize=11, loc="upper left", facecolor="#222222",
               edgecolor="#555555", labelcolor="white")
    ax1.tick_params(colors="white")
    for spine in ax1.spines.values():
        spine.set_color("#555555")

    # Correlation
    r = np.corrcoef(sub_merged["logFC_corr"], sub_merged["logFC_hybrid"])[0, 1]
    ax1.text(0.95, 0.05, f"r = {r:.4f}", transform=ax1.transAxes,
             fontsize=14, color="white", ha="right", va="bottom",
             bbox=dict(boxstyle="round", facecolor="#333333", alpha=0.8))

    # ── Panel B: Subclass -log10(p) scatter ──
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_facecolor(BG)

    neglog_corr = -np.log10(sub_merged["P.Value_corr"].clip(lower=1e-10))
    neglog_hybrid = -np.log10(sub_merged["P.Value_hybrid"].clip(lower=1e-10))

    for cls, color in CLASS_COLORS.items():
        mask = sub_merged["class"] == cls
        if mask.sum() == 0:
            continue
        ax2.scatter(neglog_corr[mask], neglog_hybrid[mask],
                    c=color, s=100, alpha=0.9, edgecolors="white",
                    linewidths=0.5, label=cls, zorder=5)

    for _, row in sub_merged.iterrows():
        if row["FDR_corr"] < 0.1 or row["FDR_hybrid"] < 0.1:
            ax2.annotate(row["celltype"],
                         (-np.log10(max(row["P.Value_corr"], 1e-10)),
                          -np.log10(max(row["P.Value_hybrid"], 1e-10))),
                         fontsize=10, color="white", fontweight="bold",
                         xytext=(5, 5), textcoords="offset points")

    pmax = max(neglog_corr.max(), neglog_hybrid.max()) * 1.1
    ax2.plot([0, pmax], [0, pmax], "w--", alpha=0.3, linewidth=1)
    ax2.axhline(-np.log10(0.05), color="#888888", linestyle=":", linewidth=0.8)
    ax2.axvline(-np.log10(0.05), color="#888888", linestyle=":", linewidth=0.8)
    ax2.set_xlabel("-log10(p) corr QC", fontsize=14, color="white")
    ax2.set_ylabel("-log10(p) hybrid QC", fontsize=14, color="white")
    ax2.set_title("Subclass: significance comparison", fontsize=18,
                  fontweight="bold", color="white")
    ax2.tick_params(colors="white")
    for spine in ax2.spines.values():
        spine.set_color("#555555")

    # ── Panel C: Supertype logFC comparison (top hits) ──
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.set_facecolor(BG)

    # Show the top 30 most significant cell types (either mode)
    sup_merged_sorted = sup_merged.copy()
    sup_merged_sorted["min_p"] = sup_merged_sorted[
        ["P.Value_corr", "P.Value_hybrid"]].min(axis=1)
    top30 = sup_merged_sorted.nsmallest(30, "min_p")

    # Bar plot: logFC side by side
    y_pos = np.arange(len(top30))
    bar_height = 0.35

    bars_corr = ax3.barh(y_pos - bar_height/2, top30["logFC_corr"],
                          height=bar_height, color="#4fc3f7", alpha=0.7,
                          label="corr QC")
    bars_hybrid = ax3.barh(y_pos + bar_height/2, top30["logFC_hybrid"],
                            height=bar_height, color="#ef5350", alpha=0.7,
                            label="hybrid QC")

    ax3.set_yticks(y_pos)
    ax3.set_yticklabels(top30["celltype"], fontsize=9, color="white")
    ax3.axvline(0, color="white", linewidth=0.5)
    ax3.set_xlabel("logFC (SCZ vs Control)", fontsize=14, color="white")
    ax3.set_title("Supertype: top 30 most significant\n(logFC comparison)",
                  fontsize=16, fontweight="bold", color="white")
    ax3.legend(fontsize=11, loc="lower right", facecolor="#222222",
               edgecolor="#555555", labelcolor="white")
    ax3.tick_params(colors="white")
    for spine in ax3.spines.values():
        spine.set_color("#555555")
    ax3.invert_yaxis()

    # Add FDR stars
    for i, (_, row) in enumerate(top30.iterrows()):
        for mode, x_offset, color in [("_corr", -0.02, "#4fc3f7"),
                                        ("_hybrid", 0.02, "#ef5350")]:
            fdr = row[f"FDR{mode}"]
            if fdr < 0.01:
                star = "***"
            elif fdr < 0.05:
                star = "**"
            elif fdr < 0.1:
                star = "*"
            else:
                star = ""
            if star:
                lfc = row[f"logFC{mode}"]
                x_pos = lfc + (0.05 if lfc > 0 else -0.05)
                y = i - 0.175 if mode == "_corr" else i + 0.175
                ax3.text(x_pos, y, star, fontsize=8, color=color,
                         ha="left" if lfc > 0 else "right", va="center",
                         fontweight="bold")

    # ── Panel D: Summary table ──
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.set_facecolor(BG)
    ax4.axis("off")

    # Build summary text
    lines = []
    lines.append("CRUMBLR COMPOSITIONAL ANALYSIS")
    lines.append("=" * 50)
    lines.append(f"")
    lines.append(f"{'Metric':<35s} {'corr':>10s} {'hybrid':>10s}")
    lines.append(f"{'-'*35} {'-'*10} {'-'*10}")

    # Subclass
    n_fdr05_corr = (sub_merged["FDR_corr"] < 0.05).sum()
    n_fdr05_hybrid = (sub_merged["FDR_hybrid"] < 0.05).sum()
    n_fdr10_corr = (sub_merged["FDR_corr"] < 0.10).sum()
    n_fdr10_hybrid = (sub_merged["FDR_hybrid"] < 0.10).sum()
    n_nom_corr = (sub_merged["P.Value_corr"] < 0.05).sum()
    n_nom_hybrid = (sub_merged["P.Value_hybrid"] < 0.05).sum()

    lines.append(f"{'Subclass FDR < 0.05':<35s} {n_fdr05_corr:>10d} {n_fdr05_hybrid:>10d}")
    lines.append(f"{'Subclass FDR < 0.10':<35s} {n_fdr10_corr:>10d} {n_fdr10_hybrid:>10d}")
    lines.append(f"{'Subclass nom p < 0.05':<35s} {n_nom_corr:>10d} {n_nom_hybrid:>10d}")

    # Supertype
    n_fdr05_corr_s = (sup_merged["FDR_corr"] < 0.05).sum()
    n_fdr05_hybrid_s = (sup_merged["FDR_hybrid"] < 0.05).sum()
    n_fdr10_corr_s = (sup_merged["FDR_corr"] < 0.10).sum()
    n_fdr10_hybrid_s = (sup_merged["FDR_hybrid"] < 0.10).sum()
    n_nom_corr_s = (sup_merged["P.Value_corr"] < 0.05).sum()
    n_nom_hybrid_s = (sup_merged["P.Value_hybrid"] < 0.05).sum()

    lines.append(f"")
    lines.append(f"{'Supertype FDR < 0.05':<35s} {n_fdr05_corr_s:>10d} {n_fdr05_hybrid_s:>10d}")
    lines.append(f"{'Supertype FDR < 0.10':<35s} {n_fdr10_corr_s:>10d} {n_fdr10_hybrid_s:>10d}")
    lines.append(f"{'Supertype nom p < 0.05':<35s} {n_nom_corr_s:>10d} {n_nom_hybrid_s:>10d}")

    # Key hits
    lines.append(f"")
    lines.append(f"KEY SUBCLASS HITS (FDR < 0.10):")
    for _, row in sub_merged.sort_values("P.Value_corr").iterrows():
        if row["FDR_corr"] < 0.1 or row["FDR_hybrid"] < 0.1:
            d_corr = "↑SCZ" if row["logFC_corr"] > 0 else "↓SCZ"
            d_hybrid = "↑SCZ" if row["logFC_hybrid"] > 0 else "↓SCZ"
            lines.append(
                f"  {row['celltype']:20s} "
                f"corr: {d_corr} {row['logFC_corr']:+.3f} FDR={row['FDR_corr']:.4f} | "
                f"hybrid: {d_hybrid} {row['logFC_hybrid']:+.3f} FDR={row['FDR_hybrid']:.4f}"
            )

    text = "\n".join(lines)
    ax4.text(0.02, 0.98, text, transform=ax4.transAxes, fontsize=11,
             color="white", family="monospace", va="top", ha="left")

    fig.suptitle("QC Mode Comparison: corr vs hybrid\n"
                 "(crumblr compositional analysis, SCZ vs Control)",
                 fontsize=22, fontweight="bold", color="white", y=1.02)

    outpath = os.path.join(OUT_DIR, "qc_mode_comparison_crumblr.png")
    plt.savefig(outpath, dpi=200, facecolor=BG, bbox_inches="tight")
    plt.close()
    print(f"Saved: {outpath}")

    return sub_merged, sup_merged


# ──────────────────────────────────────────────────────────────────────
# 2. Cell count comparison
# ──────────────────────────────────────────────────────────────────────

def compare_cell_counts():
    """Compare per-subclass cell counts between corr and hybrid QC."""

    corr_path = os.path.join(CRUMBLR_DIR, "crumblr_input_subclass.csv")
    hybrid_path = os.path.join(CRUMBLR_DIR, "crumblr_input_subclass_hybrid.csv")

    if not os.path.exists(corr_path) or not os.path.exists(hybrid_path):
        print("Missing crumblr input files — skipping count comparison")
        return

    corr = pd.read_csv(corr_path)
    hybrid = pd.read_csv(hybrid_path)

    # Aggregate across all samples
    corr_totals = corr.groupby("celltype")["count"].sum().reset_index()
    corr_totals.columns = ["celltype", "count_corr"]
    hybrid_totals = hybrid.groupby("celltype")["count"].sum().reset_index()
    hybrid_totals.columns = ["celltype", "count_hybrid"]

    counts = corr_totals.merge(hybrid_totals, on="celltype", how="outer").fillna(0)
    counts["delta"] = counts["count_hybrid"] - counts["count_corr"]
    counts["pct_change"] = 100 * counts["delta"] / counts["count_corr"].clip(lower=1)
    counts["class"] = counts["celltype"].map(
        lambda x: SUBCLASS_TO_CLASS.get(x, "Other"))
    counts = counts.sort_values("delta", ascending=True)

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 10), facecolor=BG)

    # Panel A: Absolute change
    ax1.set_facecolor(BG)
    y_pos = np.arange(len(counts))
    colors = [CLASS_COLORS.get(c, "#888888") for c in counts["class"]]
    ax1.barh(y_pos, counts["delta"], color=colors, alpha=0.8,
             edgecolor="white", linewidth=0.3)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(counts["celltype"], fontsize=11, color="white")
    ax1.axvline(0, color="white", linewidth=0.5)
    ax1.set_xlabel("Cell count change (hybrid - corr)", fontsize=14, color="white")
    ax1.set_title("Absolute change in cell counts",
                  fontsize=18, fontweight="bold", color="white")
    ax1.tick_params(colors="white")
    for spine in ax1.spines.values():
        spine.set_color("#555555")

    # Add count labels
    for i, (_, row) in enumerate(counts.iterrows()):
        d = row["delta"]
        ax1.text(d + (50 if d >= 0 else -50), i,
                 f"{d:+,.0f}", fontsize=9, color="white",
                 ha="left" if d >= 0 else "right", va="center")

    # Panel B: Percent change
    ax2.set_facecolor(BG)
    ax2.barh(y_pos, counts["pct_change"], color=colors, alpha=0.8,
             edgecolor="white", linewidth=0.3)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(counts["celltype"], fontsize=11, color="white")
    ax2.axvline(0, color="white", linewidth=0.5)
    ax2.set_xlabel("% change (hybrid vs corr)", fontsize=14, color="white")
    ax2.set_title("Percent change in cell counts",
                  fontsize=18, fontweight="bold", color="white")
    ax2.tick_params(colors="white")
    for spine in ax2.spines.values():
        spine.set_color("#555555")

    # Totals
    total_corr = counts["count_corr"].sum()
    total_hybrid = counts["count_hybrid"].sum()
    total_delta = total_hybrid - total_corr
    fig.suptitle(f"Cell count impact of hybrid QC\n"
                 f"Total: {total_corr:,.0f} (corr) → {total_hybrid:,.0f} (hybrid) = "
                 f"{total_delta:+,.0f} ({100*total_delta/total_corr:.1f}%)",
                 fontsize=20, fontweight="bold", color="white", y=1.02)

    plt.tight_layout()
    outpath = os.path.join(OUT_DIR, "qc_mode_comparison_counts.png")
    plt.savefig(outpath, dpi=200, facecolor=BG, bbox_inches="tight")
    plt.close()
    print(f"Saved: {outpath}")

    return counts


# ──────────────────────────────────────────────────────────────────────
# 3. DE comparison
# ──────────────────────────────────────────────────────────────────────

def compare_de():
    """Compare DE results between corr and hybrid QC (if available)."""

    corr_path = os.path.join(DE_DIR, "de_summary_subclass.csv")
    hybrid_path = os.path.join(DE_DIR, "de_summary_subclass_hybrid.csv")

    if not os.path.exists(corr_path) or not os.path.exists(hybrid_path):
        print("DE results not yet available for both modes — skipping")
        return None

    corr = pd.read_csv(corr_path)
    hybrid = pd.read_csv(hybrid_path)

    merged = corr.merge(hybrid, on="celltype", suffixes=("_corr", "_hybrid"))
    merged["class"] = merged["celltype"].map(
        lambda x: SUBCLASS_TO_CLASS.get(x, "Other"))

    # Plot: number of DE genes comparison
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 10), facecolor=BG)

    # Panel A: Scatter of n_sig
    ax1.set_facecolor(BG)
    for cls, color in CLASS_COLORS.items():
        mask = merged["class"] == cls
        if mask.sum() == 0:
            continue
        ax1.scatter(merged.loc[mask, "n_sig_fdr05_corr"],
                    merged.loc[mask, "n_sig_fdr05_hybrid"],
                    c=color, s=100, alpha=0.9, edgecolors="white",
                    linewidths=0.5, label=cls, zorder=5)

    for _, row in merged.iterrows():
        if row["n_sig_fdr05_corr"] > 0 or row["n_sig_fdr05_hybrid"] > 0:
            ax1.annotate(row["celltype"],
                         (row["n_sig_fdr05_corr"], row["n_sig_fdr05_hybrid"]),
                         fontsize=9, color="white",
                         xytext=(5, 5), textcoords="offset points")

    maxval = max(merged["n_sig_fdr05_corr"].max(),
                 merged["n_sig_fdr05_hybrid"].max()) * 1.1
    ax1.plot([0, maxval], [0, maxval], "w--", alpha=0.3, linewidth=1)
    ax1.set_xlabel("# DE genes (corr QC)", fontsize=14, color="white")
    ax1.set_ylabel("# DE genes (hybrid QC)", fontsize=14, color="white")
    ax1.set_title("DE genes per subclass (FDR < 0.20)",
                  fontsize=16, fontweight="bold", color="white")
    ax1.legend(fontsize=11, loc="upper left", facecolor="#222222",
               edgecolor="#555555", labelcolor="white")
    ax1.tick_params(colors="white")
    for spine in ax1.spines.values():
        spine.set_color("#555555")

    # Panel B: Bar comparison
    ax2.set_facecolor(BG)
    merged_sorted = merged.sort_values("n_sig_fdr05_corr", ascending=True)
    y_pos = np.arange(len(merged_sorted))
    bar_height = 0.35

    ax2.barh(y_pos - bar_height/2, merged_sorted["n_sig_fdr05_corr"],
             height=bar_height, color="#4fc3f7", alpha=0.7, label="corr QC")
    ax2.barh(y_pos + bar_height/2, merged_sorted["n_sig_fdr05_hybrid"],
             height=bar_height, color="#ef5350", alpha=0.7, label="hybrid QC")

    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(merged_sorted["celltype"], fontsize=10, color="white")
    ax2.set_xlabel("# DE genes (FDR < 0.20)", fontsize=14, color="white")
    ax2.set_title("DE genes by subclass",
                  fontsize=16, fontweight="bold", color="white")
    ax2.legend(fontsize=11, loc="lower right", facecolor="#222222",
               edgecolor="#555555", labelcolor="white")
    ax2.tick_params(colors="white")
    for spine in ax2.spines.values():
        spine.set_color("#555555")

    total_corr = merged["n_sig_fdr05_corr"].sum()
    total_hybrid = merged["n_sig_fdr05_hybrid"].sum()
    fig.suptitle(f"Differential Expression: corr vs hybrid QC\n"
                 f"Total DE genes: {total_corr} (corr) vs {total_hybrid} (hybrid)",
                 fontsize=20, fontweight="bold", color="white", y=1.02)

    plt.tight_layout()
    outpath = os.path.join(OUT_DIR, "qc_mode_comparison_de.png")
    plt.savefig(outpath, dpi=200, facecolor=BG, bbox_inches="tight")
    plt.close()
    print(f"Saved: {outpath}")

    return merged


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("QC Mode Comparison: corr vs hybrid")
    print("=" * 60)

    # 1. Crumblr comparison
    print("\n1. Crumblr compositional analysis...")
    sub_merged = compare_crumblr_subclass()
    sup_merged = compare_crumblr_supertype()

    if sub_merged is not None and sup_merged is not None:
        plot_crumblr_comparison(sub_merged, sup_merged)

        # Print key comparison
        print("\n  SUBCLASS-LEVEL COMPARISON:")
        print(f"  {'Cell type':<20s} {'logFC_corr':>10s} {'logFC_hyb':>10s} "
              f"{'FDR_corr':>10s} {'FDR_hyb':>10s} {'Direction':>12s}")
        print(f"  {'-'*20} {'-'*10} {'-'*10} {'-'*10} {'-'*10} {'-'*12}")
        for _, row in sub_merged.sort_values("P.Value_corr").head(10).iterrows():
            direction = "SAME" if (np.sign(row["logFC_corr"]) ==
                                   np.sign(row["logFC_hybrid"])) else "FLIPPED"
            print(f"  {row['celltype']:<20s} {row['logFC_corr']:>+10.4f} "
                  f"{row['logFC_hybrid']:>+10.4f} {row['FDR_corr']:>10.4f} "
                  f"{row['FDR_hybrid']:>10.4f} {direction:>12s}")

    # 2. Cell count comparison
    print("\n2. Cell count changes...")
    counts = compare_cell_counts()
    if counts is not None:
        print(f"\n  Top gainers:")
        top = counts.nlargest(5, "delta")
        for _, row in top.iterrows():
            print(f"    {row['celltype']:<20s}: {row['delta']:>+6,.0f} "
                  f"({row['pct_change']:>+5.1f}%)")

    # 3. DE comparison
    print("\n3. Differential expression...")
    de_merged = compare_de()

    # Save summary CSV
    print("\n4. Saving summary...")
    summary_rows = []
    if sub_merged is not None:
        for _, row in sub_merged.iterrows():
            summary_rows.append({
                "level": "subclass",
                "celltype": row["celltype"],
                "class": row["class"],
                "logFC_corr": row["logFC_corr"],
                "logFC_hybrid": row["logFC_hybrid"],
                "logFC_delta": row["logFC_hybrid"] - row["logFC_corr"],
                "pval_corr": row["P.Value_corr"],
                "pval_hybrid": row["P.Value_hybrid"],
                "FDR_corr": row["FDR_corr"],
                "FDR_hybrid": row["FDR_hybrid"],
            })

    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        csv_path = os.path.join(OUT_DIR, "qc_mode_comparison_summary.csv")
        summary_df.to_csv(csv_path, index=False)
        print(f"  Saved: {csv_path}")

    print("\nDone!")


if __name__ == "__main__":
    main()
