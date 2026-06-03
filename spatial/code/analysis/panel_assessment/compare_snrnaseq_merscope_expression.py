#!/usr/bin/env python3
"""
Compare snRNAseq and MERSCOPE gene expression for matched cell types.

Computes pseudobulk mean expression per subclass for both platforms,
finds shared genes and subclasses, and computes per-gene, per-cell-type,
and overall correlations. Produces publication-quality figures and CSV outputs.

Usage:
    python compare_snrnaseq_merscope_expression.py
"""

import os
import sys
import time
import warnings
from pathlib import Path
from glob import glob

import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import pearsonr, spearmanr

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ── Shared module imports ─────────────────────────────────────────────────
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "modules"))
from constants import SUBCLASS_TO_CLASS, CLASS_COLORS
from pseudobulk import compute_pseudobulk_mean
from reference_utils import load_and_normalize_reference

# ── Paths ──────────────────────────────────────────────────────────────
PROJECT_DIR = Path("/Users/shreejoy/Github/SCZ_Xenium")
SNRNASEQ_PATH = PROJECT_DIR / "data/reference/nicole_sea_ad_snrnaseq_reference.h5ad"
MERSCOPE_DIR = PROJECT_DIR / "output/merscope_h5ad"
OUTPUT_DIR = PROJECT_DIR / "output/presentation"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ══════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════

def main():
    t0 = time.time()

    # ------------------------------------------------------------------
    # 1. Load & normalize snRNAseq reference
    # ------------------------------------------------------------------
    print("=" * 70)
    print("STEP 1: Loading snRNAseq reference")
    print("=" * 70)

    adata_sn = load_and_normalize_reference(str(SNRNASEQ_PATH),
                                              normalize=True, min_cells=0)
    print(f"  X dtype: {adata_sn.X.dtype}")

    # Map column names: snRNAseq uses 'Subclass' (capital S)
    sn_subclass_col = "Subclass"
    print(f"\n  snRNAseq subclass column: '{sn_subclass_col}'")
    print(f"  Unique subclasses: {adata_sn.obs[sn_subclass_col].nunique()}")

    # ------------------------------------------------------------------
    # 2. Load & normalize MERSCOPE 4K files
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STEP 2: Loading MERSCOPE 4K annotated files")
    print("=" * 70)

    merscope_files = sorted(glob(str(MERSCOPE_DIR / "*_4000_*_annotated.h5ad")))
    print(f"  Found {len(merscope_files)} MERSCOPE 4K files")

    merscope_adatas = []
    for f in merscope_files:
        fname = Path(f).name
        adata_m = ad.read_h5ad(f)

        # Filter to QC-pass cells (use corr_qc_pass)
        if "corr_qc_pass" in adata_m.obs.columns:
            n_before = adata_m.shape[0]
            adata_m = adata_m[adata_m.obs["corr_qc_pass"] == True].copy()
            n_after = adata_m.shape[0]
            print(f"  {fname}: {n_before} -> {n_after} cells (corr_qc_pass filter)")
        else:
            print(f"  {fname}: {adata_m.shape[0]} cells (no corr_qc_pass column)")

        # Filter out Unassigned
        if "corr_subclass" in adata_m.obs.columns:
            adata_m = adata_m[adata_m.obs["corr_subclass"] != "Unassigned"].copy()

        merscope_adatas.append(adata_m)

    # Concatenate all MERSCOPE files
    adata_merscope = ad.concat(merscope_adatas, join="inner")
    print(f"\n  Combined MERSCOPE shape: {adata_merscope.shape}")
    print(f"  MERSCOPE X dtype: {adata_merscope.X.dtype}")

    # Normalize MERSCOPE (auto-detect raw counts)
    X_sample_m = adata_merscope.X[:100, :100]
    if sp.issparse(X_sample_m):
        X_sample_m = X_sample_m.toarray()
    if np.all(X_sample_m == X_sample_m.astype(int)):
        print("  MERSCOPE data appears to be raw counts. Normalizing...")
        sc.pp.normalize_total(adata_merscope, target_sum=1e4)
        sc.pp.log1p(adata_merscope)
    else:
        print("  MERSCOPE data appears already normalized.")

    merscope_subclass_col = "corr_subclass"
    print(f"\n  MERSCOPE subclass column: '{merscope_subclass_col}'")
    print(f"  Unique subclasses: {adata_merscope.obs[merscope_subclass_col].nunique()}")
    print(f"  MERSCOPE subclass distribution:")
    print(adata_merscope.obs[merscope_subclass_col].value_counts().to_string())

    # ------------------------------------------------------------------
    # 3. Find shared genes and subclasses
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STEP 3: Finding shared genes and subclasses")
    print("=" * 70)

    sn_genes = set(adata_sn.var_names)
    merscope_genes = set(adata_merscope.var_names)
    shared_genes = sorted(sn_genes & merscope_genes)
    print(f"  snRNAseq genes: {len(sn_genes)}")
    print(f"  MERSCOPE genes: {len(merscope_genes)}")
    print(f"  Shared genes: {len(shared_genes)}")

    sn_subclasses = set(adata_sn.obs[sn_subclass_col].unique())
    merscope_subclasses = set(adata_merscope.obs[merscope_subclass_col].unique())
    shared_subclasses = sorted(sn_subclasses & merscope_subclasses)
    print(f"\n  snRNAseq subclasses: {sn_subclasses}")
    print(f"  MERSCOPE subclasses: {merscope_subclasses}")
    print(f"  Shared subclasses ({len(shared_subclasses)}): {shared_subclasses}")

    # ------------------------------------------------------------------
    # 4. Compute pseudobulk means
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STEP 4: Computing pseudobulk mean expression")
    print("=" * 70)

    # Subset snRNAseq to shared subclasses
    adata_sn_sub = adata_sn[adata_sn.obs[sn_subclass_col].isin(shared_subclasses)].copy()

    print("\n  snRNAseq pseudobulk (per subclass):")
    pb_sn = compute_pseudobulk_mean(adata_sn_sub, sn_subclass_col, gene_subset=shared_genes)

    # Subset MERSCOPE to shared subclasses
    adata_merscope_sub = adata_merscope[
        adata_merscope.obs[merscope_subclass_col].isin(shared_subclasses)
    ].copy()

    print("\n  MERSCOPE pseudobulk (per subclass):")
    pb_merscope = compute_pseudobulk_mean(
        adata_merscope_sub, merscope_subclass_col, gene_subset=shared_genes
    )

    # Align: ensure same subclasses in both
    common_subclasses = sorted(set(pb_sn.index) & set(pb_merscope.index))
    pb_sn = pb_sn.loc[common_subclasses]
    pb_merscope = pb_merscope.loc[common_subclasses]

    print(f"\n  Final shared subclasses in pseudobulk: {len(common_subclasses)}")
    print(f"  Final shared genes in pseudobulk: {len(shared_genes)}")
    print(f"  snRNAseq pseudobulk shape: {pb_sn.shape}")
    print(f"  MERSCOPE pseudobulk shape: {pb_merscope.shape}")

    # ------------------------------------------------------------------
    # 5. Compute correlations
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STEP 5: Computing correlations")
    print("=" * 70)

    # 5a. Per-cell-type correlation: for each subclass, correlate expression
    #     profile across genes between platforms
    print("\n  5a. Per-cell-type correlations (across genes):")
    celltype_corrs = {}
    for ct in common_subclasses:
        sn_vals = pb_sn.loc[ct].values
        mer_vals = pb_merscope.loc[ct].values
        r, p = pearsonr(sn_vals, mer_vals)
        rho, p_sp = spearmanr(sn_vals, mer_vals)
        celltype_corrs[ct] = {
            "pearson_r": r, "pearson_p": p,
            "spearman_rho": rho, "spearman_p": p_sp,
            "class": SUBCLASS_TO_CLASS.get(ct, "Unknown"),
            "n_cells_snrnaseq": int((adata_sn.obs[sn_subclass_col] == ct).sum()),
            "n_cells_merscope": int((adata_merscope.obs[merscope_subclass_col] == ct).sum()),
        }
        print(f"    {ct:20s}  r={r:.3f}  rho={rho:.3f}  "
              f"(sn={celltype_corrs[ct]['n_cells_snrnaseq']}, "
              f"mer={celltype_corrs[ct]['n_cells_merscope']})")

    df_ct_corr = pd.DataFrame(celltype_corrs).T
    df_ct_corr.index.name = "subclass"

    # 5b. Per-gene correlation: for each gene, correlate mean expression
    #     across subclasses between platforms
    print("\n  5b. Per-gene correlations (across cell types):")
    gene_corrs = {}
    for gene in shared_genes:
        sn_vals = pb_sn[gene].values
        mer_vals = pb_merscope[gene].values
        # Skip genes with zero variance in either platform
        if np.std(sn_vals) < 1e-10 or np.std(mer_vals) < 1e-10:
            continue
        r, p = pearsonr(sn_vals, mer_vals)
        rho, p_sp = spearmanr(sn_vals, mer_vals)
        gene_corrs[gene] = {
            "pearson_r": r, "pearson_p": p,
            "spearman_rho": rho, "spearman_p": p_sp,
            "mean_sn": np.mean(sn_vals),
            "mean_merscope": np.mean(mer_vals),
        }

    df_gene_corr = pd.DataFrame(gene_corrs).T
    df_gene_corr.index.name = "gene"
    print(f"    Genes with valid correlations: {len(df_gene_corr)}")
    print(f"    Median Pearson r: {df_gene_corr['pearson_r'].median():.3f}")
    print(f"    Mean Pearson r: {df_gene_corr['pearson_r'].mean():.3f}")
    print(f"    Genes with r > 0.5: {(df_gene_corr['pearson_r'] > 0.5).sum()}")
    print(f"    Genes with r > 0.8: {(df_gene_corr['pearson_r'] > 0.8).sum()}")

    # 5c. Overall correlation of the full pseudobulk matrix
    print("\n  5c. Overall correlation (full pseudobulk matrix):")
    sn_flat = pb_sn.values.flatten()
    mer_flat = pb_merscope.values.flatten()
    r_overall, p_overall = pearsonr(sn_flat, mer_flat)
    rho_overall, p_sp_overall = spearmanr(sn_flat, mer_flat)
    print(f"    Overall Pearson r: {r_overall:.4f} (p={p_overall:.2e})")
    print(f"    Overall Spearman rho: {rho_overall:.4f} (p={p_sp_overall:.2e})")

    # ------------------------------------------------------------------
    # 6. Create figures
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STEP 6: Creating figures")
    print("=" * 70)

    # ── Figure style setup ─────────────────────────────────────────────
    plt.rcParams.update({
        "font.size": 16,
        "axes.titlesize": 20,
        "axes.labelsize": 18,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "legend.fontsize": 14,
        "figure.dpi": 150,
    })

    # ── 6a. Scatter: overall pseudobulk (colored by cell class) ────────
    print("\n  6a. Overall pseudobulk scatter (snRNAseq vs MERSCOPE)...")
    fig, ax = plt.subplots(figsize=(10, 9))

    for ct in common_subclasses:
        cell_class = SUBCLASS_TO_CLASS.get(ct, "Unknown")
        color = CLASS_COLORS.get(cell_class, "#999999")
        sn_vals = pb_sn.loc[ct].values
        mer_vals = pb_merscope.loc[ct].values
        ax.scatter(sn_vals, mer_vals, c=color, alpha=0.15, s=8, rasterized=True)

    # Overlay per-subclass means as larger labeled points
    for ct in common_subclasses:
        cell_class = SUBCLASS_TO_CLASS.get(ct, "Unknown")
        color = CLASS_COLORS.get(cell_class, "#999999")
        sn_mean = pb_sn.loc[ct].mean()
        mer_mean = pb_merscope.loc[ct].mean()
        ax.scatter(sn_mean, mer_mean, c=color, s=120, edgecolors="black",
                   linewidths=1.0, zorder=10)
        ax.annotate(ct, (sn_mean, mer_mean), fontsize=9, ha="left",
                    xytext=(5, 5), textcoords="offset points")

    # Legend for classes
    for cls_name, cls_color in CLASS_COLORS.items():
        ax.scatter([], [], c=cls_color, s=80, label=cls_name, edgecolors="black")
    ax.legend(loc="upper left", fontsize=14, framealpha=0.9)

    # Diagonal
    lims = [0, max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, "k--", alpha=0.3, linewidth=1)

    ax.set_xlabel("snRNAseq mean expression (log1p)", fontsize=18)
    ax.set_ylabel("MERSCOPE mean expression (log1p)", fontsize=18)
    ax.set_title(
        f"Pseudobulk Expression: snRNAseq vs MERSCOPE\n"
        f"Pearson r = {r_overall:.3f}, Spearman rho = {rho_overall:.3f}\n"
        f"({len(shared_genes)} genes, {len(common_subclasses)} subclasses)",
        fontsize=18,
    )
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "snrnaseq_vs_merscope_pseudobulk_scatter.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("    Saved: snrnaseq_vs_merscope_pseudobulk_scatter.png")

    # ── 6b. Per-cell-type correlation bar chart ────────────────────────
    print("  6b. Per-cell-type correlation bar chart...")
    df_ct_sorted = df_ct_corr.sort_values("pearson_r", ascending=True)

    fig, ax = plt.subplots(figsize=(10, max(8, len(common_subclasses) * 0.4)))
    colors = [CLASS_COLORS.get(SUBCLASS_TO_CLASS.get(ct, "Unknown"), "#999999")
              for ct in df_ct_sorted.index]
    bars = ax.barh(range(len(df_ct_sorted)), df_ct_sorted["pearson_r"].values,
                   color=colors, edgecolor="black", linewidth=0.5)
    ax.set_yticks(range(len(df_ct_sorted)))
    ax.set_yticklabels(df_ct_sorted.index, fontsize=14)
    ax.set_xlabel("Pearson r (across genes)", fontsize=18)
    ax.set_title("Per-Cell-Type Correlation\n(snRNAseq vs MERSCOPE pseudobulk)",
                 fontsize=20)
    ax.axvline(x=0, color="black", linewidth=0.5)

    # Add r values as text
    for i, (idx, row) in enumerate(df_ct_sorted.iterrows()):
        r_val = row["pearson_r"]
        ax.text(r_val + 0.01 if r_val >= 0 else r_val - 0.01,
                i, f"{r_val:.3f}", va="center",
                ha="left" if r_val >= 0 else "right", fontsize=12)

    # Legend
    for cls_name, cls_color in CLASS_COLORS.items():
        ax.barh([], [], color=cls_color, label=cls_name, edgecolor="black")
    ax.legend(loc="lower right", fontsize=14)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "snrnaseq_vs_merscope_celltype_corr.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("    Saved: snrnaseq_vs_merscope_celltype_corr.png")

    # ── 6c. Per-gene correlation histogram ─────────────────────────────
    print("  6c. Per-gene correlation histogram...")
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Pearson
    axes[0].hist(df_gene_corr["pearson_r"], bins=50, color="#4C72B0",
                 edgecolor="black", linewidth=0.5, alpha=0.8)
    axes[0].axvline(df_gene_corr["pearson_r"].median(), color="red",
                    linestyle="--", linewidth=2,
                    label=f"Median = {df_gene_corr['pearson_r'].median():.3f}")
    axes[0].set_xlabel("Pearson r", fontsize=18)
    axes[0].set_ylabel("Number of genes", fontsize=18)
    axes[0].set_title("Per-Gene Pearson Correlation", fontsize=20)
    axes[0].legend(fontsize=14)

    # Spearman
    axes[1].hist(df_gene_corr["spearman_rho"], bins=50, color="#DD8452",
                 edgecolor="black", linewidth=0.5, alpha=0.8)
    axes[1].axvline(df_gene_corr["spearman_rho"].median(), color="red",
                    linestyle="--", linewidth=2,
                    label=f"Median = {df_gene_corr['spearman_rho'].median():.3f}")
    axes[1].set_xlabel("Spearman rho", fontsize=18)
    axes[1].set_ylabel("Number of genes", fontsize=18)
    axes[1].set_title("Per-Gene Spearman Correlation", fontsize=20)
    axes[1].legend(fontsize=14)

    fig.suptitle(
        f"Gene-Level Correlation (snRNAseq vs MERSCOPE)\n"
        f"{len(df_gene_corr)} genes across {len(common_subclasses)} subclasses",
        fontsize=20, y=1.02,
    )
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "snrnaseq_vs_merscope_gene_corr_hist.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("    Saved: snrnaseq_vs_merscope_gene_corr_hist.png")

    # ── 6d. Heatmap comparison of top variable genes ───────────────────
    print("  6d. Heatmap comparison of top variable genes...")

    # Select top variable genes: highest variance across subclasses in snRNAseq
    sn_var = pb_sn.var(axis=0)
    top_var_genes = sn_var.nlargest(50).index.tolist()

    # Also get top correlated genes
    top_corr_genes = df_gene_corr.nlargest(50, "pearson_r").index.tolist()

    # Use top variable genes for heatmap
    fig, axes = plt.subplots(1, 2, figsize=(22, 14))

    # snRNAseq heatmap
    sn_heatmap = pb_sn[top_var_genes]
    # Z-score per gene for visualization
    sn_z = (sn_heatmap - sn_heatmap.mean(axis=0)) / (sn_heatmap.std(axis=0) + 1e-10)

    im0 = axes[0].imshow(sn_z.values, aspect="auto", cmap="RdBu_r",
                          vmin=-2, vmax=2)
    axes[0].set_yticks(range(len(sn_z)))
    axes[0].set_yticklabels(sn_z.index, fontsize=11)
    axes[0].set_xticks(range(len(top_var_genes)))
    axes[0].set_xticklabels(top_var_genes, fontsize=8, rotation=90)
    axes[0].set_title("snRNAseq (z-scored)", fontsize=20)
    plt.colorbar(im0, ax=axes[0], shrink=0.5, label="Z-score")

    # MERSCOPE heatmap (same genes, same order)
    mer_heatmap = pb_merscope[top_var_genes]
    mer_z = (mer_heatmap - mer_heatmap.mean(axis=0)) / (mer_heatmap.std(axis=0) + 1e-10)

    im1 = axes[1].imshow(mer_z.values, aspect="auto", cmap="RdBu_r",
                          vmin=-2, vmax=2)
    axes[1].set_yticks(range(len(mer_z)))
    axes[1].set_yticklabels(mer_z.index, fontsize=11)
    axes[1].set_xticks(range(len(top_var_genes)))
    axes[1].set_xticklabels(top_var_genes, fontsize=8, rotation=90)
    axes[1].set_title("MERSCOPE (z-scored)", fontsize=20)
    plt.colorbar(im1, ax=axes[1], shrink=0.5, label="Z-score")

    fig.suptitle(
        "Top 50 Variable Genes: snRNAseq vs MERSCOPE\n(z-scored pseudobulk expression)",
        fontsize=20, y=1.01,
    )
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "snrnaseq_vs_merscope_heatmap_top_genes.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("    Saved: snrnaseq_vs_merscope_heatmap_top_genes.png")

    # ── 6e. Per-cell-type scatter grid ─────────────────────────────────
    print("  6e. Per-cell-type scatter grid...")
    n_cts = len(common_subclasses)
    ncols = 5
    nrows = int(np.ceil(n_cts / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 4))
    axes = axes.flatten()

    for i, ct in enumerate(common_subclasses):
        ax = axes[i]
        cell_class = SUBCLASS_TO_CLASS.get(ct, "Unknown")
        color = CLASS_COLORS.get(cell_class, "#999999")
        sn_vals = pb_sn.loc[ct].values
        mer_vals = pb_merscope.loc[ct].values
        r_val = celltype_corrs[ct]["pearson_r"]

        ax.scatter(sn_vals, mer_vals, c=color, alpha=0.3, s=6, rasterized=True)
        lims = [0, max(sn_vals.max(), mer_vals.max()) * 1.1]
        ax.plot(lims, lims, "k--", alpha=0.3, linewidth=0.8)
        ax.set_title(f"{ct}\nr = {r_val:.3f}", fontsize=12)
        ax.set_xlabel("snRNAseq", fontsize=10)
        ax.set_ylabel("MERSCOPE", fontsize=10)
        ax.tick_params(labelsize=8)

    # Hide unused axes
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(
        "Per-Cell-Type Pseudobulk Scatter\n(snRNAseq vs MERSCOPE, each dot = 1 gene)",
        fontsize=20, y=1.01,
    )
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "snrnaseq_vs_merscope_celltype_scatter_grid.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("    Saved: snrnaseq_vs_merscope_celltype_scatter_grid.png")

    # ------------------------------------------------------------------
    # 7. Save outputs
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STEP 7: Saving outputs")
    print("=" * 70)

    # Save per-cell-type correlations
    df_ct_corr.to_csv(OUTPUT_DIR / "snrnaseq_vs_merscope_celltype_corr.csv")
    print(f"  Saved: snrnaseq_vs_merscope_celltype_corr.csv")

    # Save per-gene correlations
    df_gene_corr.to_csv(OUTPUT_DIR / "snrnaseq_vs_merscope_gene_corr.csv")
    print(f"  Saved: snrnaseq_vs_merscope_gene_corr.csv")

    # Save pseudobulk matrices
    pb_sn.to_csv(OUTPUT_DIR / "pseudobulk_snrnaseq_by_subclass.csv")
    pb_merscope.to_csv(OUTPUT_DIR / "pseudobulk_merscope_by_subclass.csv")
    print(f"  Saved: pseudobulk_snrnaseq_by_subclass.csv")
    print(f"  Saved: pseudobulk_merscope_by_subclass.csv")

    # ------------------------------------------------------------------
    # Summary statistics
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Shared genes: {len(shared_genes)}")
    print(f"  Shared subclasses: {len(common_subclasses)}")
    print(f"  Overall Pearson r: {r_overall:.4f}")
    print(f"  Overall Spearman rho: {rho_overall:.4f}")
    print(f"  Per-cell-type Pearson r:")
    print(f"    Mean: {df_ct_corr['pearson_r'].astype(float).mean():.3f}")
    print(f"    Median: {df_ct_corr['pearson_r'].astype(float).median():.3f}")
    print(f"    Min: {df_ct_corr['pearson_r'].astype(float).min():.3f} "
          f"({df_ct_corr['pearson_r'].astype(float).idxmin()})")
    print(f"    Max: {df_ct_corr['pearson_r'].astype(float).max():.3f} "
          f"({df_ct_corr['pearson_r'].astype(float).idxmax()})")
    print(f"  Per-gene Pearson r:")
    print(f"    Mean: {df_gene_corr['pearson_r'].mean():.3f}")
    print(f"    Median: {df_gene_corr['pearson_r'].median():.3f}")
    print(f"    Fraction r > 0.5: {(df_gene_corr['pearson_r'] > 0.5).mean():.1%}")
    print(f"    Fraction r > 0.8: {(df_gene_corr['pearson_r'] > 0.8).mean():.1%}")

    # Top 10 and bottom 10 genes by correlation
    print(f"\n  Top 10 genes by Pearson r:")
    for gene, row in df_gene_corr.nlargest(10, "pearson_r").iterrows():
        print(f"    {gene:15s}  r={row['pearson_r']:.3f}")

    print(f"\n  Bottom 10 genes by Pearson r:")
    for gene, row in df_gene_corr.nsmallest(10, "pearson_r").iterrows():
        print(f"    {gene:15s}  r={row['pearson_r']:.3f}")

    elapsed = time.time() - t0
    print(f"\n  Total runtime: {elapsed:.1f}s")
    print("  Done!")


if __name__ == "__main__":
    main()
