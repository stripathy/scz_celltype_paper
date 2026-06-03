#!/usr/bin/env python
"""
13_gwas_vs_composition.py — Correlate GWAS enrichment with case-control composition.

Compares GWAS cell-type enrichment (from RBH combined taxonomy) with
case-control compositional changes in SCZ (from Endresz et al., in prep;
7-cohort snRNA-seq meta-analysis, crumblr model).

Tests whether cell types with stronger genetic risk also show altered
proportions in SCZ brains. Performs the analysis across all matched
cell types and within SST interneurons specifically.

Inputs:
  - results/tables/rbh_combined_enrichment.csv (GWAS enrichment)
  - External: SCZ_Xenium/data/nicole_scz_snrnaseq_betas/final_results_crumblr_7_cohorts.csv

Outputs:
  - results/tables/gwas_vs_casecontrol_composition.csv
  - results/figures/gwas_vs_casecontrol_composition.png
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats

from scz_celltype_enrichment.config import (
    TABLES_DIR, FIGURES_DIR, PROJECT_ROOT,
    FIGURE_DPI, FDR_THRESHOLD,
)
from scz_celltype_enrichment.utils import ensure_dirs, classify_supertypes_by_class


# Path to Endresz et al. composition betas (7-cohort crumblr meta-analysis)
COMPOSITION_CSV = (   # snRNA-seq composition export, via the shared interface (was SCZ_Xenium/data/...)
    Path(__file__).resolve().parents[2] / "shared" / "snrnaseq_de"
    / "nicole_scz_snrnaseq_betas" / "final_results_crumblr_7_cohorts.csv"
)


def load_and_merge(enrich_path, comp_path):
    """Load GWAS enrichment and composition data, merge on cell type name."""
    enrich = pd.read_csv(str(enrich_path))
    comp = pd.read_csv(str(comp_path))

    # GWAS enrichment uses "supertype" column
    enrich["gwas_logp"] = -np.log10(np.clip(enrich["p_value"], 1e-300, 1))

    # Merge on cell type name (SEA-AD types overlap between datasets)
    merged = enrich.merge(
        comp[["CellType", "estimate", "se", "zval", "pval", "padj"]],
        left_on="supertype", right_on="CellType", how="inner",
    ).drop(columns=["CellType"])

    merged.rename(columns={
        "estimate": "comp_beta",
        "se": "comp_se",
        "zval": "comp_zval",
        "pval": "comp_pval",
        "padj": "comp_padj",
    }, inplace=True)

    merged["abs_comp_beta"] = merged["comp_beta"].abs()
    merged["comp_logp"] = -np.log10(np.clip(merged["comp_pval"], 1e-300, 1))

    # Cell class annotation
    cls_map = classify_supertypes_by_class(merged["supertype"].tolist())
    merged["cell_class"] = [cls_map.get(ct, "Other") for ct in merged["supertype"]]

    # SST flag
    merged["is_sst"] = merged["supertype"].str.startswith("Sst")

    return merged


def compute_correlations(merged):
    """Compute Spearman correlations for all types and SST subset."""
    results = {}

    # All types: beta vs beta
    r, p = stats.spearmanr(merged["beta"], merged["comp_beta"])
    results["all_beta_vs_comp_beta"] = {"r": r, "p": p, "n": len(merged)}

    # All types: |comp_beta| vs GWAS -log10(p)
    r, p = stats.spearmanr(merged["abs_comp_beta"], merged["gwas_logp"])
    results["all_abscomp_vs_gwas_logp"] = {"r": r, "p": p, "n": len(merged)}

    # SST subset
    sst = merged[merged["is_sst"]]
    if len(sst) >= 5:
        r, p = stats.spearmanr(sst["abs_comp_beta"], sst["gwas_logp"])
        results["sst_abscomp_vs_gwas_logp"] = {"r": r, "p": p, "n": len(sst)}

        r, p = stats.spearmanr(sst["comp_logp"], sst["gwas_logp"])
        results["sst_complogp_vs_gwas_logp"] = {"r": r, "p": p, "n": len(sst)}

        r, p = stats.spearmanr(sst["comp_beta"], sst["beta"])
        results["sst_comp_beta_vs_gwas_beta"] = {"r": r, "p": p, "n": len(sst)}

    return results


def plot_gwas_vs_composition(merged, corr_results, output_path):
    """
    4-panel figure comparing GWAS enrichment with case-control composition.

    Cowplot style: white background, no grid, half-open axes, large text.
    """
    mpl.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 1.5,
        "font.size": 14,
    })

    LS = 18  # label size
    TS = 20  # title size
    TKS = 14  # tick size
    ANNOT = 10
    class_colors = {
        "GABAergic": "#D55E00",
        "Glutamatergic": "#0072B2",
        "Non-neuronal": "#009E73",
    }

    sst = merged[merged["is_sst"]]

    fig, axes = plt.subplots(2, 2, figsize=(20, 18))

    # --- Panel A: All types, |composition beta| vs GWAS -log10(p) ---
    ax = axes[0, 0]
    for cls in ["Non-neuronal", "Glutamatergic", "GABAergic"]:
        mask = merged["cell_class"] == cls
        ax.scatter(
            merged.loc[mask, "abs_comp_beta"],
            merged.loc[mask, "gwas_logp"],
            c=class_colors.get(cls, "gray"), s=50, alpha=0.7,
            label=cls, edgecolors="black", linewidth=0.3,
        )

    # Annotate top GWAS enriched
    for _, row in merged.nlargest(8, "gwas_logp").iterrows():
        ax.annotate(row["supertype"], (row["abs_comp_beta"], row["gwas_logp"]),
                    fontsize=ANNOT, xytext=(4, 4), textcoords="offset points")

    info = corr_results.get("all_abscomp_vs_gwas_logp", {})
    ax.text(0.03, 0.97, f"r = {info.get('r', 0):.2f}, p = {info.get('p', 1):.2e}\nn = {info.get('n', 0)}",
            transform=ax.transAxes, fontsize=13, va="top", fontstyle="italic", color="#555")
    ax.set_xlabel("|Composition change| (case − control)", fontsize=LS)
    ax.set_ylabel(r"GWAS enrichment $-\log_{10}(p)$", fontsize=LS)
    ax.set_title("A. All cell types", fontsize=TS, fontweight="bold")
    ax.tick_params(labelsize=TKS)
    ax.legend(fontsize=12, frameon=False, loc="upper right")

    # --- Panel B: All types, composition beta vs GWAS beta ---
    ax = axes[0, 1]
    for cls in ["Non-neuronal", "Glutamatergic", "GABAergic"]:
        mask = merged["cell_class"] == cls
        ax.scatter(
            merged.loc[mask, "comp_beta"],
            merged.loc[mask, "beta"],
            c=class_colors.get(cls, "gray"), s=50, alpha=0.7,
            label=cls, edgecolors="black", linewidth=0.3,
        )

    info = corr_results.get("all_beta_vs_comp_beta", {})
    ax.text(0.03, 0.97, f"r = {info.get('r', 0):.2f}, p = {info.get('p', 1):.2e}\nn = {info.get('n', 0)}",
            transform=ax.transAxes, fontsize=13, va="top", fontstyle="italic", color="#555")
    ax.axhline(0, color="gray", linestyle=":", linewidth=0.8, alpha=0.5)
    ax.axvline(0, color="gray", linestyle=":", linewidth=0.8, alpha=0.5)
    ax.set_xlabel("Composition change beta (case − control)", fontsize=LS)
    ax.set_ylabel("GWAS enrichment beta", fontsize=LS)
    ax.set_title("B. All cell types (signed)", fontsize=TS, fontweight="bold")
    ax.tick_params(labelsize=TKS)

    # --- Panel C: SST only, |composition beta| vs GWAS -log10(p) ---
    ax = axes[1, 0]
    ax.scatter(
        sst["abs_comp_beta"], sst["gwas_logp"],
        c="#D55E00", s=80, alpha=0.8,
        edgecolors="black", linewidth=0.5,
    )
    for _, row in sst.iterrows():
        ax.annotate(row["supertype"], (row["abs_comp_beta"], row["gwas_logp"]),
                    fontsize=ANNOT, xytext=(5, 3), textcoords="offset points",
                    fontweight="bold")

    info = corr_results.get("sst_abscomp_vs_gwas_logp", {})
    if info:
        ax.text(0.03, 0.97, f"r = {info.get('r', 0):.2f}, p = {info.get('p', 1):.2e}\nn = {info.get('n', 0)}",
                transform=ax.transAxes, fontsize=13, va="top", fontstyle="italic", color="#555")
    ax.set_xlabel("|Composition change| (case − control)", fontsize=LS)
    ax.set_ylabel(r"GWAS enrichment $-\log_{10}(p)$", fontsize=LS)
    ax.set_title("C. SST interneurons", fontsize=TS, fontweight="bold")
    ax.tick_params(labelsize=TKS)

    # --- Panel D: SST only, significance vs significance ---
    ax = axes[1, 1]
    ax.scatter(
        sst["comp_logp"], sst["gwas_logp"],
        c="#D55E00", s=80, alpha=0.8,
        edgecolors="black", linewidth=0.5,
    )
    for _, row in sst.iterrows():
        ax.annotate(row["supertype"], (row["comp_logp"], row["gwas_logp"]),
                    fontsize=ANNOT, xytext=(5, 3), textcoords="offset points",
                    fontweight="bold")

    info = corr_results.get("sst_complogp_vs_gwas_logp", {})
    if info:
        ax.text(0.03, 0.97, f"r = {info.get('r', 0):.2f}, p = {info.get('p', 1):.2e}\nn = {info.get('n', 0)}",
                transform=ax.transAxes, fontsize=13, va="top", fontstyle="italic", color="#555")
    ax.set_xlabel(r"Composition change $-\log_{10}(p)$", fontsize=LS)
    ax.set_ylabel(r"GWAS enrichment $-\log_{10}(p)$", fontsize=LS)
    ax.set_title("D. SST: significance vs significance", fontsize=TS, fontweight="bold")
    ax.tick_params(labelsize=TKS)

    plt.tight_layout(h_pad=3, w_pad=3)
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    mpl.rcParams.update(mpl.rcParamsDefault)
    print(f"  Saved: {output_path}")


def main():
    ensure_dirs(TABLES_DIR, FIGURES_DIR)

    print("=" * 60)
    print("Step 13: GWAS vs Case-Control Composition")
    print("=" * 60)
    t0 = time.time()

    # Check inputs
    enrich_path = TABLES_DIR / "rbh_combined_enrichment.csv"
    if not enrich_path.exists():
        print(f"ERROR: {enrich_path} not found. Run script 10 first.")
        return
    if not COMPOSITION_CSV.exists():
        print(f"ERROR: Composition data not found at {COMPOSITION_CSV}")
        return

    # Load and merge
    print("\n1. Loading and merging data...")
    merged = load_and_merge(enrich_path, COMPOSITION_CSV)
    n_sst = merged["is_sst"].sum()
    print(f"  Matched cell types: {len(merged)}")
    print(f"  SST subtypes: {n_sst}")
    print(f"  Cell classes: {merged['cell_class'].value_counts().to_dict()}")

    # Compute correlations
    print("\n2. Computing correlations...")
    corr_results = compute_correlations(merged)
    for name, vals in corr_results.items():
        print(f"  {name}: r={vals['r']:.3f}, p={vals['p']:.2e} (n={vals['n']})")

    # Save merged table
    output_csv = TABLES_DIR / "gwas_vs_casecontrol_composition.csv"
    merged.to_csv(str(output_csv), index=False)
    print(f"\n  Saved: {output_csv}")

    # Generate figure
    print("\n3. Generating figure...")
    plot_gwas_vs_composition(
        merged, corr_results,
        FIGURES_DIR / "gwas_vs_casecontrol_composition.png",
    )

    # Print key SST findings
    sst = merged[merged["is_sst"]].sort_values("gwas_logp", ascending=False)
    print(f"\n  SST subtypes ranked by GWAS enrichment:")
    for _, row in sst.head(10).iterrows():
        sig = "***" if row["p_fdr"] < 0.001 else "**" if row["p_fdr"] < 0.01 else "*" if row["p_fdr"] < 0.05 else ""
        print(f"    {row['supertype']:20s}  GWAS_logp={row['gwas_logp']:5.1f}{sig:4s}  "
              f"comp_beta={row['comp_beta']:+.3f}  comp_padj={row['comp_padj']:.3f}")

    print(f"\nStep 13 complete. Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
