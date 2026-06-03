#!/usr/bin/env python3
"""
Presentation-quality scatter: snRNAseq meta-analysis betas vs Xenium density logFC.

Two-panel layout: Neuronal supertypes (left) | Non-neuronal supertypes (right).

Input:
  output/density_analysis/density_results_supertype_cortical.csv
  data/nicole_scz_snrnaseq_betas/final_results_crumblr_7_cohorts.csv (neuronal)
  data/nicole_scz_snrnaseq_betas/final_results_crumblr_7_nonN_cohorts.csv (non-neuronal)

Output:
  output/density_analysis/snrnaseq_vs_density_supertype.png
  output/density_analysis/snrnaseq_vs_density_comparison.csv
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import pearsonr
from adjustText import adjust_text

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "analysis"))
from config import BASE_DIR, infer_class, BG_COLOR

DENSITY_PATH = os.path.join(BASE_DIR, "output", "density_analysis",
                             "density_results_supertype_cortical.csv")
NICOLE_NEURONAL_PATH = os.path.join(BASE_DIR, "data", "nicole_scz_snrnaseq_betas",
                                     "final_results_crumblr_7_cohorts.csv")
NICOLE_NONNEURONAL_PATH = os.path.join(BASE_DIR, "data", "nicole_scz_snrnaseq_betas",
                                        "final_results_crumblr_7_nonN_cohorts.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "density_analysis")

BG = BG_COLOR
CLASS_COLORS = {'Glut': '#00ADF8', 'GABA': '#F05A28', 'NN': '#808080', 'Other': '#999999'}


def plot_panel(ax, panel_df, snrna_sig, snrna_nom, title, show_ylabel=True, label_all=False):
    """Plot a single density scatter panel."""
    ax.set_facecolor(BG)

    for _, row in panel_df.iterrows():
        c = CLASS_COLORS.get(row["class"], '#999999')
        if pd.notna(row.get("se")):
            ax.plot([row["beta_snrnaseq"] - row["se"], row["beta_snrnaseq"] + row["se"]],
                    [row["logFC_density"], row["logFC_density"]],
                    color=c, alpha=0.15, linewidth=1.2, zorder=2, solid_capstyle='round')
        if pd.notna(row.get("se_density")) and pd.notna(row["logFC_density"]):
            ax.plot([row["beta_snrnaseq"], row["beta_snrnaseq"]],
                    [row["logFC_density"] - row["se_density"],
                     row["logFC_density"] + row["se_density"]],
                    color=c, alpha=0.15, linewidth=1.2, zorder=2, solid_capstyle='round')

    for cls in ['Glut', 'GABA', 'NN', 'Other']:
        mask = panel_df['class'] == cls
        if mask.sum() == 0:
            continue
        sub = panel_df[mask]
        ax.scatter(sub["beta_snrnaseq"], sub["logFC_density"],
                   c=CLASS_COLORS[cls], s=80, alpha=0.8,
                   edgecolors="white", linewidth=0.5, zorder=5,
                   label=f"{cls} (n={mask.sum()})")

    valid = panel_df["beta_snrnaseq"].notna() & panel_df["logFC_density"].notna()
    m = panel_df[valid]
    if len(m) > 2:
        z = np.polyfit(m["beta_snrnaseq"], m["logFC_density"], 1)
        lim_x = max(abs(m["beta_snrnaseq"]).max(), 0.3) * 1.3
        x_line = np.linspace(-lim_x, lim_x, 100)
        ax.plot(x_line, np.polyval(z, x_line), color="#888888", alpha=0.5,
                linewidth=1.5, linestyle="-", zorder=1)

    ax.axhline(0, color="#555555", alpha=0.4, linewidth=0.8, zorder=1)
    ax.axvline(0, color="#555555", alpha=0.4, linewidth=0.8, zorder=1)

    # Labels
    texts = []
    for _, row in panel_df.iterrows():
        ct = row["celltype"]
        if label_all or ct in snrna_sig:
            txt = ax.text(row["beta_snrnaseq"], row["logFC_density"],
                          f"  {ct}", fontsize=14, fontweight="bold",
                          color="white", alpha=0.95, zorder=10)
            texts.append(txt)
        elif ct in snrna_nom:
            txt = ax.text(row["beta_snrnaseq"], row["logFC_density"],
                          f"  {ct}", fontsize=12, fontweight="normal",
                          color="#bbbbbb", alpha=0.85, zorder=10)
            texts.append(txt)

    if texts:
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="#aaaaaa",
                                   alpha=0.5, linewidth=0.7),
                    force_text=(0.9, 0.9), force_points=(0.6, 0.6),
                    expand_text=(1.3, 1.3), expand_points=(1.5, 1.5))

    if len(m) > 2:
        r_val, p_val = pearsonr(m["beta_snrnaseq"], m["logFC_density"])
        ax.text(0.97, 0.04,
                f"r = {r_val:.2f} (p = {p_val:.1e})\nn = {len(m)}",
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=16, color="#dddddd",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="#222222",
                          edgecolor="#555555", alpha=0.9))

    ax.set_xlabel("snRNAseq meta-analysis beta", fontsize=20, color="white")
    if show_ylabel:
        ax.set_ylabel("Xenium density logFC (cells/mm²)", fontsize=20, color="white")
    ax.set_title(title, fontsize=22, fontweight="bold", color="white", pad=10)

    ax.tick_params(colors="white", labelsize=16)
    for spine in ax.spines.values():
        spine.set_color("#555555")
    ax.grid(True, alpha=0.15, color="#555555")

    ax.text(0.03, 0.04,
            "Bold = FDR < 0.1  |  Light = nom. p < 0.05 or |logFC| > 0.8",
            transform=ax.transAxes, ha="left", va="bottom",
            fontsize=13, color="#aaaaaa", fontstyle="italic")


def main():
    # Load density results
    density = pd.read_csv(DENSITY_PATH)
    density = density.rename(columns={
        "supertype": "celltype", "logFC": "logFC_density",
        "pval": "pval_density", "fdr": "fdr_density", "se": "se_density",
    })

    # Load Nicole snRNAseq betas (stratified)
    dfs = []
    for path, stratum in [(NICOLE_NEURONAL_PATH, "neuronal"),
                           (NICOLE_NONNEURONAL_PATH, "non-neuronal")]:
        df = pd.read_csv(path)
        df = df[~df["CellType"].str.contains("SEAAD", na=False)]
        df["analysis_stratum"] = stratum
        dfs.append(df)
    nicole = pd.concat(dfs, ignore_index=True)
    nicole = nicole.rename(columns={
        "CellType": "celltype", "estimate": "beta_snrnaseq",
        "pval": "pval_snrnaseq", "padj": "padj_snrnaseq",
    })
    print(f"Nicole: {len(nicole)} types")

    merged = pd.merge(
        nicole[["celltype", "beta_snrnaseq", "se", "pval_snrnaseq", "padj_snrnaseq"]],
        density[["celltype", "logFC_density", "pval_density", "fdr_density", "se_density"]],
        on="celltype", how="inner"
    )
    merged["class"] = merged["celltype"].apply(infer_class)
    print(f"Merged: {len(merged)} shared supertypes")

    neuronal = merged[merged["class"].isin(["Glut", "GABA"])].copy()
    nonneuronal = merged[merged["class"] == "NN"].copy()
    print(f"Neuronal: {len(neuronal)}, Non-neuronal: {len(nonneuronal)}")

    valid = merged["logFC_density"].notna() & merged["beta_snrnaseq"].notna()
    m = merged[valid]
    r_all, p_all = pearsonr(m["beta_snrnaseq"], m["logFC_density"])
    r_neur, p_neur = pearsonr(neuronal["beta_snrnaseq"], neuronal["logFC_density"])
    r_nn, p_nn = pearsonr(nonneuronal["beta_snrnaseq"], nonneuronal["logFC_density"]) if len(nonneuronal) > 2 else (np.nan, np.nan)
    print(f"All: r={r_all:.3f}, Neuronal: r={r_neur:.3f}, Non-neuronal: r={r_nn:.3f}")

    snrna_sig = set(merged[merged["padj_snrnaseq"] < 0.1]["celltype"].values)
    snrna_nom = set(merged[(merged["pval_snrnaseq"] < 0.05) &
                            (merged["padj_snrnaseq"] >= 0.1)]["celltype"].values)
    large_density = set(merged[abs(merged["logFC_density"]) > 0.8]["celltype"].values)
    snrna_nom = snrna_nom | large_density

    # Two-panel figure
    fig, (ax_n, ax_nn) = plt.subplots(1, 2, figsize=(22, 10), facecolor=BG)

    plot_panel(ax_n, neuronal, snrna_sig, snrna_nom,
              title=f'Neuronal supertypes (n={len(neuronal)})', show_ylabel=True)

    plot_panel(ax_nn, nonneuronal, snrna_sig, snrna_nom,
              title=f'Non-neuronal supertypes (n={len(nonneuronal)})',
              show_ylabel=False, label_all=True)

    legend_elements = [
        Line2D([0], [0], marker="o", color=BG, markerfacecolor=CLASS_COLORS["Glut"],
               markersize=12, label=f"Glutamatergic (n={(neuronal['class']=='Glut').sum()})", linewidth=0),
        Line2D([0], [0], marker="o", color=BG, markerfacecolor=CLASS_COLORS["GABA"],
               markersize=12, label=f"GABAergic (n={(neuronal['class']=='GABA').sum()})", linewidth=0),
    ]
    leg = ax_n.legend(handles=legend_elements, loc="upper left", fontsize=15,
                      frameon=True, fancybox=True, framealpha=0.85,
                      edgecolor="#555555", labelcolor="white")
    leg.get_frame().set_facecolor("#222222")

    fig.suptitle("SCZ density effects: snRNAseq meta-analysis vs Xenium spatial",
                 fontsize=26, fontweight="bold", color="white", y=1.02)

    plt.tight_layout(pad=2.0)
    out_png = os.path.join(OUTPUT_DIR, "snrnaseq_vs_density_supertype.png")
    plt.savefig(out_png, dpi=200, facecolor=BG, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_png}")

    merged.to_csv(os.path.join(OUTPUT_DIR, "snrnaseq_vs_density_comparison.csv"), index=False)

    # Print FDR-significant
    sig_either = merged[(merged["padj_snrnaseq"] < 0.05) | (merged.get("fdr_density", pd.Series(dtype=float)) < 0.05)]
    if len(sig_either) > 0:
        print(f"\nFDR-significant in either platform:")
        for _, row in sig_either.sort_values("celltype").iterrows():
            concordant = "Y" if (row["beta_snrnaseq"] * row["logFC_density"]) > 0 else "N"
            print(f"  {row['celltype']:30s} snRNA={row['beta_snrnaseq']:+.3f} "
                  f"(padj={row['padj_snrnaseq']:.3f})  "
                  f"density={row['logFC_density']:+.3f}  concordant={concordant}")


if __name__ == "__main__":
    main()
