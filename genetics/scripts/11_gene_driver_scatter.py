#!/usr/bin/env python
"""
11_gene_driver_scatter.py — Gene driver scatter plots for cell-type enrichment.

For a given cell type, plots:
  X-axis: cell-type specificity score
  Y-axis: -log10(GWAS p-value)

Genes in the top-right quadrant (high specificity + FDR-significant GWAS) are
the "driver genes" — the genes that make this cell type enriched for SCZ risk.

Based on the approach recommended by Laramie Duncan and colleagues.

Usage:
    python scripts/11_gene_driver_scatter.py                     # All FDR-significant cell types
    python scripts/11_gene_driver_scatter.py --type Sst_2        # Single cell type
    python scripts/11_gene_driver_scatter.py --type Sst_2 Lamp5_5 Pvalb_6  # Multiple types
    python scripts/11_gene_driver_scatter.py --top 10            # Top 10 by enrichment p-value

Outputs (per cell type):
    results/figures/gene_drivers/gene_driver_scatter_{TYPE}.png
    results/tables/gene_drivers/gene_drivers_{TYPE}.csv
"""
import sys
import time
import argparse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests

from scz_celltype_enrichment.config import (
    TABLES_DIR, FIGURES_DIR, INTERMEDIATES_DIR,
    FIGURE_DPI, TITLE_FONTSIZE, LABEL_FONTSIZE, TICK_FONTSIZE, LEGEND_FONTSIZE,
    FDR_THRESHOLD,
)
from scz_celltype_enrichment.utils import ensure_dirs


def load_gene_data():
    """Load specificity matrix and GWAS gene-level results."""
    spec = pd.read_csv(str(INTERMEDIATES_DIR / "seaad_supertype_specificity.csv"), index_col=0)
    gwas = pd.read_csv(str(INTERMEDIATES_DIR / "gwas_matched.csv"))

    # Compute FDR on GWAS gene-level p-values
    _, gwas_fdr, _, _ = multipletests(gwas["P"].values, method="fdr_bh")
    gwas["P_FDR"] = gwas_fdr

    # Align specificity to GWAS gene order
    overlap = gwas["SYMBOL"].values
    spec_matched = spec.loc[overlap]

    return spec_matched, gwas


def identify_driver_genes(specificity_values, gwas_df, spec_percentile=90):
    """
    Identify driver genes for a cell type.

    A driver gene has:
      - Cell-type specificity in the top `spec_percentile`% of nonzero values
      - GWAS gene-level FDR < 0.05

    Parameters
    ----------
    specificity_values : np.ndarray
        Specificity scores for this cell type (aligned with gwas_df).
    gwas_df : pd.DataFrame
        GWAS results with columns P, P_FDR, ZSTAT, SYMBOL.
    spec_percentile : int
        Percentile threshold for specificity (default: 90 = top 10%).

    Returns
    -------
    pd.DataFrame
        Driver genes sorted by specificity × -log10(p) product.
    """
    nonzero = specificity_values > 0
    spec_threshold = np.percentile(specificity_values[nonzero], spec_percentile) if nonzero.sum() > 0 else 0

    logp = -np.log10(np.clip(gwas_df["P"].values, 1e-300, 1))

    is_high_spec = specificity_values > spec_threshold
    is_gwas_fdr = gwas_df["P_FDR"].values < FDR_THRESHOLD
    is_driver = is_high_spec & is_gwas_fdr

    result = pd.DataFrame({
        "gene": gwas_df["SYMBOL"].values,
        "specificity": specificity_values,
        "gwas_z": gwas_df["ZSTAT"].values,
        "gwas_p": gwas_df["P"].values,
        "gwas_fdr": gwas_df["P_FDR"].values,
        "neglog10p": logp,
        "is_high_spec": is_high_spec,
        "is_gwas_fdr": is_gwas_fdr,
        "is_driver": is_driver,
        "driver_score": specificity_values * logp,
    })

    return result.sort_values("driver_score", ascending=False).reset_index(drop=True), spec_threshold


def plot_gene_driver_scatter(gene_df, cell_type, spec_threshold, output_path):
    """
    Generate the gene driver scatter plot.

    X: specificity, Y: -log10(GWAS p), colored by quadrant.
    """
    fig, ax = plt.subplots(figsize=(12, 10))

    spec = gene_df["specificity"].values
    logp = gene_df["neglog10p"].values
    is_driver = gene_df["is_driver"].values
    is_high_spec = gene_df["is_high_spec"].values
    is_gwas_fdr = gene_df["is_gwas_fdr"].values

    # Assign colors
    colors = np.full(len(gene_df), "#dddddd")
    sizes = np.full(len(gene_df), 6.0)
    alphas = np.full(len(gene_df), 0.3)

    # High specificity only (blue)
    mask_hs = is_high_spec & ~is_gwas_fdr
    colors[mask_hs] = "#3498db"
    sizes[mask_hs] = 15
    alphas[mask_hs] = 0.5

    # FDR-significant GWAS only (orange)
    mask_hg = is_gwas_fdr & ~is_high_spec
    colors[mask_hg] = "#e67e22"
    sizes[mask_hg] = 15
    alphas[mask_hg] = 0.5

    # Driver genes: both (red)
    colors[is_driver] = "#e74c3c"
    sizes[is_driver] = 45
    alphas[is_driver] = 0.8

    # Plot in layers (gray first, then colored, then drivers on top)
    for mask, zorder in [
        (~is_high_spec & ~is_gwas_fdr, 1),
        (mask_hs, 2),
        (mask_hg, 2),
        (is_driver, 3),
    ]:
        idx = np.where(mask)[0]
        if len(idx) == 0:
            continue
        ax.scatter(spec[idx], logp[idx], c=colors[idx], s=sizes[idx],
                   alpha=alphas[idx], edgecolors="none", zorder=zorder)

    # Label top driver genes
    drivers = gene_df[gene_df["is_driver"]].nlargest(25, "driver_score")
    for _, row in drivers.iterrows():
        ax.annotate(
            row["gene"], (row["specificity"], row["neglog10p"]),
            fontsize=9, fontweight="bold", color="#c0392b",
            xytext=(5, 5), textcoords="offset points",
            zorder=4,
        )

    # Also label a few high-specificity nominally significant genes
    nominal = gene_df[
        gene_df["is_high_spec"] & ~gene_df["is_gwas_fdr"] & (gene_df["gwas_p"] < 0.05)
    ].nlargest(8, "driver_score")
    for _, row in nominal.iterrows():
        ax.annotate(
            row["gene"], (row["specificity"], row["neglog10p"]),
            fontsize=8, color="#2980b9",
            xytext=(5, 3), textcoords="offset points",
            zorder=4,
        )

    # Threshold lines
    fdr_logp = -np.log10(gene_df[gene_df["is_gwas_fdr"]]["gwas_p"].max()) if is_gwas_fdr.any() else 2
    ax.axhline(fdr_logp, color="#e74c3c", linestyle="--", linewidth=1, alpha=0.5,
               label=f"GWAS FDR < 0.05 (-log10p ≈ {fdr_logp:.1f})")
    ax.axhline(-np.log10(0.05), color="#f39c12", linestyle=":", linewidth=1, alpha=0.4,
               label="Nominal (p < 0.05)")
    ax.axvline(spec_threshold, color="gray", linestyle="--", linewidth=1, alpha=0.3,
               label=f"Top 10% specificity (> {spec_threshold:.4f})")

    # Counts
    n_driver = is_driver.sum()
    n_hs = mask_hs.sum()
    n_hg = mask_hg.sum()
    ax.text(0.98, 0.98,
            f"Driver genes (top-right): {n_driver}\n"
            f"High specificity only: {n_hs}\n"
            f"GWAS FDR-sig only: {n_hg}",
            transform=ax.transAxes, fontsize=11, va="top", ha="right",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

    ax.set_xlabel(f"{cell_type} specificity score", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("-log$_{10}$(GWAS p-value)", fontsize=LABEL_FONTSIZE)
    ax.set_title(f"Gene drivers of SCZ enrichment: {cell_type}\n"
                 f"Red = high specificity + GWAS FDR < 0.05",
                 fontsize=TITLE_FONTSIZE)
    ax.tick_params(labelsize=TICK_FONTSIZE)
    ax.legend(fontsize=LEGEND_FONTSIZE - 1, loc="upper left")

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def process_cell_type(cell_type, spec_matched, gwas_df, fig_dir, table_dir):
    """Generate driver scatter plot and CSV for one cell type."""
    if cell_type not in spec_matched.columns:
        print(f"  WARNING: {cell_type} not found in specificity matrix, skipping.")
        return None

    specificity = spec_matched[cell_type].values
    gene_df, spec_threshold = identify_driver_genes(specificity, gwas_df)

    n_drivers = gene_df["is_driver"].sum()

    # Save figure
    fig_path = fig_dir / f"gene_driver_scatter_{cell_type}.png"
    plot_gene_driver_scatter(gene_df, cell_type, spec_threshold, fig_path)

    # Save driver genes CSV
    drivers_csv = gene_df[gene_df["is_driver"]].copy()
    drivers_csv = drivers_csv[["gene", "specificity", "gwas_z", "gwas_p", "gwas_fdr",
                                "driver_score"]].reset_index(drop=True)
    csv_path = table_dir / f"gene_drivers_{cell_type}.csv"
    drivers_csv.to_csv(str(csv_path), index=False)

    print(f"  {cell_type}: {n_drivers} driver genes → {fig_path.name}")
    return n_drivers


def main():
    parser = argparse.ArgumentParser(description="Gene driver scatter plots")
    parser.add_argument("--type", nargs="+", help="Cell type name(s)")
    parser.add_argument("--top", type=int, help="Top N enriched types")
    parser.add_argument("--all-sig", action="store_true",
                        help="All FDR-significant types from enrichment")
    args = parser.parse_args()

    fig_dir = FIGURES_DIR / "gene_drivers"
    table_dir = TABLES_DIR / "gene_drivers"
    ensure_dirs(fig_dir, table_dir)

    print("=" * 60)
    print("Step 11: Gene Driver Scatter Plots")
    print("=" * 60)
    t0 = time.time()

    print("\n1. Loading data...")
    spec_matched, gwas_df = load_gene_data()
    print(f"  Specificity: {spec_matched.shape[0]:,} genes × {spec_matched.shape[1]} types")
    print(f"  GWAS FDR < 0.05: {(gwas_df['P_FDR'] < FDR_THRESHOLD).sum():,} genes")

    # Determine which cell types to process
    if args.type:
        cell_types = args.type
    elif args.top or args.all_sig:
        enrich = pd.read_csv(str(TABLES_DIR / "seaad_magma_scz_enrichment_all_supertypes.csv"))
        if args.all_sig:
            cell_types = enrich[enrich["p_fdr"] < FDR_THRESHOLD]["supertype"].tolist()
        else:
            cell_types = enrich.nsmallest(args.top, "p_value")["supertype"].tolist()
    else:
        # Default: top 10
        enrich = pd.read_csv(str(TABLES_DIR / "seaad_magma_scz_enrichment_all_supertypes.csv"))
        cell_types = enrich.nsmallest(10, "p_value")["supertype"].tolist()

    print(f"\n2. Processing {len(cell_types)} cell types...")
    summary = []
    for ct in cell_types:
        n = process_cell_type(ct, spec_matched, gwas_df, fig_dir, table_dir)
        if n is not None:
            summary.append({"cell_type": ct, "n_driver_genes": n})

    # Save summary
    if summary:
        summary_df = pd.DataFrame(summary).sort_values("n_driver_genes", ascending=False)
        summary_df.to_csv(str(table_dir / "driver_gene_summary.csv"), index=False)
        print(f"\n  Summary saved: {table_dir / 'driver_gene_summary.csv'}")
        print(f"\n  Driver gene counts:")
        for _, row in summary_df.iterrows():
            print(f"    {row['cell_type']:25s} {row['n_driver_genes']:4d} drivers")

    print(f"\nStep 11 complete. Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
