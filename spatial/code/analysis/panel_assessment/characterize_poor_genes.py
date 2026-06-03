#!/usr/bin/env python3
"""
Characterize genes with poor cross-platform correlation (snRNAseq vs MERSCOPE).

Investigates gene-level properties that predict poor reproducibility:
- Expression level (mean, variance across cell types)
- Detection rate (fraction of cells expressing each gene)
- Cell-type specificity (entropy, max/mean ratio)
- Platform bias (systematic over/under-expression in MERSCOPE vs snRNAseq)
- Gene biotype (coding, lncRNA, antisense, etc. inferred from gene names)
- Gene length proxy (number of characters is crude; we look at expression ratio)

Usage:
    python characterize_poor_genes.py
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
from scipy.stats import pearsonr, spearmanr, mannwhitneyu

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ── Shared module imports ─────────────────────────────────────────────────
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "modules"))
from constants import SUBCLASS_TO_CLASS
from gene_properties import (classify_gene_biotype, compute_detection_rate,
                              compute_specificity)

# ── Paths ──────────────────────────────────────────────────────────────
PROJECT_DIR = Path("/Users/shreejoy/Github/SCZ_Xenium")
SNRNASEQ_PATH = PROJECT_DIR / "data/reference/nicole_sea_ad_snrnaseq_reference.h5ad"
MERSCOPE_DIR = PROJECT_DIR / "output/merscope_h5ad"
OUTPUT_DIR = PROJECT_DIR / "output/presentation"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Load precomputed results
GENE_CORR_PATH = OUTPUT_DIR / "snrnaseq_vs_merscope_gene_corr.csv"
PB_SN_PATH = OUTPUT_DIR / "pseudobulk_snrnaseq_by_subclass.csv"
PB_MER_PATH = OUTPUT_DIR / "pseudobulk_merscope_by_subclass.csv"


def main():
    t0 = time.time()

    # ------------------------------------------------------------------
    # 1. Load precomputed data
    # ------------------------------------------------------------------
    print("=" * 70)
    print("STEP 1: Loading precomputed data")
    print("=" * 70)

    df_gene_corr = pd.read_csv(GENE_CORR_PATH, index_col=0)
    pb_sn = pd.read_csv(PB_SN_PATH, index_col=0)
    pb_mer = pd.read_csv(PB_MER_PATH, index_col=0)

    print(f"  Gene correlations: {len(df_gene_corr)} genes")
    print(f"  Pseudobulk snRNAseq: {pb_sn.shape}")
    print(f"  Pseudobulk MERSCOPE: {pb_mer.shape}")

    # ------------------------------------------------------------------
    # 2. Load raw data for detection rates
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STEP 2: Loading raw data for detection rates")
    print("=" * 70)

    adata_sn = ad.read_h5ad(SNRNASEQ_PATH)
    print(f"  snRNAseq shape: {adata_sn.shape}")

    merscope_files = sorted(glob(str(MERSCOPE_DIR / "*_4000_*_annotated.h5ad")))
    merscope_adatas = []
    for f in merscope_files:
        adata_m = ad.read_h5ad(f)
        if "corr_qc_pass" in adata_m.obs.columns:
            adata_m = adata_m[adata_m.obs["corr_qc_pass"] == True].copy()
        if "corr_subclass" in adata_m.obs.columns:
            adata_m = adata_m[adata_m.obs["corr_subclass"] != "Unassigned"].copy()
        merscope_adatas.append(adata_m)
    adata_merscope = ad.concat(merscope_adatas, join="inner")
    print(f"  MERSCOPE shape: {adata_merscope.shape}")

    shared_genes = list(df_gene_corr.index)

    # ------------------------------------------------------------------
    # 3. Compute gene properties
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STEP 3: Computing gene properties")
    print("=" * 70)

    # 3a. Detection rates
    print("  Computing detection rates...")
    det_sn = compute_detection_rate(adata_sn, shared_genes)
    det_mer = compute_detection_rate(adata_merscope, shared_genes)

    # 3b. Mean expression across all cell types (from pseudobulk)
    mean_sn = pb_sn[shared_genes].mean(axis=0)
    mean_mer = pb_mer[shared_genes].mean(axis=0)

    # 3c. Variance across cell types
    var_sn = pb_sn[shared_genes].var(axis=0)
    var_mer = pb_mer[shared_genes].var(axis=0)

    # 3d. Expression ratio (MERSCOPE / snRNAseq) — log2 fold change
    log2fc = np.log2((mean_mer + 0.001) / (mean_sn + 0.001))

    # 3e. Cell-type specificity
    print("  Computing cell-type specificity (snRNAseq)...")
    spec_sn = compute_specificity(pb_sn[shared_genes])
    print("  Computing cell-type specificity (MERSCOPE)...")
    spec_mer = compute_specificity(pb_mer[shared_genes])

    # 3f. Gene biotype
    biotypes = pd.Series({g: classify_gene_biotype(g) for g in shared_genes})

    # 3g. Combine into master table
    print("  Building master gene properties table...")
    df = pd.DataFrame({
        'pearson_r': df_gene_corr.loc[shared_genes, 'pearson_r'],
        'spearman_rho': df_gene_corr.loc[shared_genes, 'spearman_rho'],
        'mean_sn': mean_sn,
        'mean_mer': mean_mer,
        'mean_both': (mean_sn + mean_mer) / 2,
        'var_sn': var_sn,
        'var_mer': var_mer,
        'log2fc_mer_vs_sn': log2fc,
        'abs_log2fc': np.abs(log2fc),
        'det_rate_sn': det_sn,
        'det_rate_mer': det_mer,
        'det_ratio': det_mer / (det_sn + 0.001),
        'entropy_sn': spec_sn.loc[shared_genes, 'entropy'],
        'entropy_mer': spec_mer.loc[shared_genes, 'entropy'],
        'max_over_mean_sn': spec_sn.loc[shared_genes, 'max_over_mean'],
        'cv_sn': spec_sn.loc[shared_genes, 'cv'],
        'cv_mer': spec_mer.loc[shared_genes, 'cv'],
        'top_ct_sn': spec_sn.loc[shared_genes, 'top_celltype'],
        'top_ct_mer': spec_mer.loc[shared_genes, 'top_celltype'],
        'same_top_ct': spec_sn.loc[shared_genes, 'top_celltype'].values == spec_mer.loc[shared_genes, 'top_celltype'].values,
        'biotype': biotypes,
    })

    # Define quintiles
    df['r_quintile'] = pd.qcut(df['pearson_r'], 5, labels=['Q1 (worst)', 'Q2', 'Q3', 'Q4', 'Q5 (best)'])
    df['r_bin'] = pd.cut(df['pearson_r'], bins=[-1, 0, 0.3, 0.5, 0.7, 0.9, 1.0],
                          labels=['<0', '0-0.3', '0.3-0.5', '0.5-0.7', '0.7-0.9', '>0.9'])

    n_genes = len(df)
    print(f"\n  Master table: {n_genes} genes x {df.shape[1]} properties")

    # ------------------------------------------------------------------
    # 4. Diagnostic summaries
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STEP 4: Characterizing poorly-performing genes")
    print("=" * 70)

    bottom_20 = df[df['r_quintile'] == 'Q1 (worst)']
    top_20 = df[df['r_quintile'] == 'Q5 (best)']

    print(f"\n  Bottom quintile (n={len(bottom_20)}): median r = {bottom_20['pearson_r'].median():.3f}")
    print(f"  Top quintile (n={len(top_20)}): median r = {top_20['pearson_r'].median():.3f}")

    # Compare properties
    props_to_compare = [
        ('mean_sn', 'Mean expression (snRNAseq)'),
        ('mean_mer', 'Mean expression (MERSCOPE)'),
        ('mean_both', 'Mean expression (both)'),
        ('var_sn', 'Variance (snRNAseq)'),
        ('var_mer', 'Variance (MERSCOPE)'),
        ('abs_log2fc', '|log2FC| MERSCOPE vs snRNAseq'),
        ('det_rate_sn', 'Detection rate (snRNAseq)'),
        ('det_rate_mer', 'Detection rate (MERSCOPE)'),
        ('entropy_sn', 'Specificity entropy (snRNAseq)'),
        ('cv_sn', 'CV across cell types (snRNAseq)'),
        ('cv_mer', 'CV across cell types (MERSCOPE)'),
    ]

    print(f"\n  {'Property':<40s} {'Bottom20% median':>17s} {'Top20% median':>15s} {'MWU p-val':>12s}")
    print(f"  {'-'*40} {'-'*17} {'-'*15} {'-'*12}")
    for prop, label in props_to_compare:
        bot_vals = pd.to_numeric(bottom_20[prop], errors='coerce').dropna()
        top_vals = pd.to_numeric(top_20[prop], errors='coerce').dropna()
        bot_med = bot_vals.median()
        top_med = top_vals.median()
        stat, pval = mannwhitneyu(bot_vals.values, top_vals.values, alternative='two-sided')
        sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else ""
        print(f"  {label:<40s} {bot_med:>15.4f} {top_med:>15.4f} {pval:>10.2e} {sig}")

    # Biotype breakdown
    print(f"\n  Biotype breakdown by quintile:")
    biotype_table = pd.crosstab(df['r_quintile'], df['biotype'], normalize='index') * 100
    print(biotype_table.round(1).to_string())

    # Same top cell type agreement
    print(f"\n  Top cell type agreement (snRNAseq vs MERSCOPE) by quintile:")
    for q in df['r_quintile'].cat.categories:
        subset = df[df['r_quintile'] == q]
        agree_pct = subset['same_top_ct'].mean() * 100
        print(f"    {q}: {agree_pct:.1f}% agree on top cell type (n={len(subset)})")

    # Worst genes detail
    print(f"\n  Bottom 30 genes by Pearson r:")
    worst30 = df.nsmallest(30, 'pearson_r')
    print(f"  {'Gene':<20s} {'r':>6s} {'mean_sn':>8s} {'mean_mer':>9s} {'log2FC':>7s} {'det_sn':>7s} {'det_mer':>8s} {'biotype':>15s} {'top_ct_sn':>15s}")
    for gene, row in worst30.iterrows():
        print(f"  {gene:<20s} {row['pearson_r']:>6.3f} {row['mean_sn']:>8.4f} {row['mean_mer']:>9.4f} "
              f"{row['log2fc_mer_vs_sn']:>7.2f} {row['det_rate_sn']:>7.4f} {row['det_rate_mer']:>8.4f} "
              f"{row['biotype']:>15s} {row['top_ct_sn']:>15s}")

    # Genes with negative correlation
    neg_genes = df[df['pearson_r'] < 0]
    print(f"\n  Genes with NEGATIVE correlation: {len(neg_genes)}")
    if len(neg_genes) > 0:
        for gene, row in neg_genes.sort_values('pearson_r').iterrows():
            print(f"    {gene:<20s} r={row['pearson_r']:.3f}  mean_sn={row['mean_sn']:.4f}  mean_mer={row['mean_mer']:.4f}  "
                  f"log2FC={row['log2fc_mer_vs_sn']:.2f}  biotype={row['biotype']}")

    # ------------------------------------------------------------------
    # 5. Figures
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STEP 5: Creating diagnostic figures")
    print("=" * 70)

    plt.rcParams.update({
        "font.size": 16,
        "axes.titlesize": 20,
        "axes.labelsize": 18,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "legend.fontsize": 13,
        "figure.dpi": 150,
    })

    # ── 5a. Multi-panel: key predictors of poor correlation ──
    print("  5a. Key predictors of poor gene correlation...")
    fig, axes = plt.subplots(2, 3, figsize=(22, 14))

    # Panel A: r vs mean expression (both platforms)
    ax = axes[0, 0]
    sc_plot = ax.scatter(df['mean_both'], df['pearson_r'],
                         c=df['abs_log2fc'], cmap='RdYlBu_r', s=8, alpha=0.5,
                         vmin=0, vmax=4, rasterized=True)
    ax.set_xlabel('Mean expression (avg of both platforms)')
    ax.set_ylabel('Pearson r (across cell types)')
    ax.set_title('A. Correlation vs Expression Level')
    ax.axhline(0, color='red', linestyle='--', alpha=0.3)
    ax.axhline(0.5, color='gray', linestyle='--', alpha=0.3)
    plt.colorbar(sc_plot, ax=ax, label='|log2FC| MER vs SN', shrink=0.8)

    # Panel B: r vs detection rate in MERSCOPE
    ax = axes[0, 1]
    ax.scatter(df['det_rate_mer'], df['pearson_r'],
               c='#4C72B0', s=8, alpha=0.3, rasterized=True)
    ax.set_xlabel('Detection rate (MERSCOPE)')
    ax.set_ylabel('Pearson r')
    ax.set_title('B. Correlation vs MERSCOPE Detection Rate')
    ax.axhline(0, color='red', linestyle='--', alpha=0.3)
    ax.axhline(0.5, color='gray', linestyle='--', alpha=0.3)

    # Panel C: r vs variance in snRNAseq (cell-type variance)
    ax = axes[0, 2]
    ax.scatter(df['cv_sn'], df['pearson_r'],
               c='#55A868', s=8, alpha=0.3, rasterized=True)
    ax.set_xlabel('CV across cell types (snRNAseq)')
    ax.set_ylabel('Pearson r')
    ax.set_title('C. Correlation vs Cell-Type Variability')
    ax.axhline(0, color='red', linestyle='--', alpha=0.3)
    ax.axhline(0.5, color='gray', linestyle='--', alpha=0.3)

    # Panel D: r vs |log2FC| — platform bias
    ax = axes[1, 0]
    ax.scatter(df['abs_log2fc'], df['pearson_r'],
               c='#C44E52', s=8, alpha=0.3, rasterized=True)
    ax.set_xlabel('|log2FC| (MERSCOPE vs snRNAseq)')
    ax.set_ylabel('Pearson r')
    ax.set_title('D. Correlation vs Platform Bias')
    ax.axhline(0, color='red', linestyle='--', alpha=0.3)
    ax.axhline(0.5, color='gray', linestyle='--', alpha=0.3)

    # Panel E: Boxplot of key properties by quintile
    ax = axes[1, 1]
    quintile_data = []
    for q in df['r_quintile'].cat.categories:
        subset = df[df['r_quintile'] == q]
        quintile_data.append(subset['mean_both'].values)
    bp = ax.boxplot(quintile_data, labels=['Q1\n(worst)', 'Q2', 'Q3', 'Q4', 'Q5\n(best)'],
                    patch_artist=True, showfliers=False)
    colors_q = ['#d73027', '#fc8d59', '#fee090', '#91bfdb', '#4575b4']
    for patch, color in zip(bp['boxes'], colors_q):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax.set_xlabel('Correlation quintile')
    ax.set_ylabel('Mean expression (avg both platforms)')
    ax.set_title('E. Expression Level by Correlation Quintile')

    # Panel F: detection rate ratio by quintile
    ax = axes[1, 2]
    quintile_det = []
    for q in df['r_quintile'].cat.categories:
        subset = df[df['r_quintile'] == q]
        quintile_det.append(subset['det_rate_mer'].values)
    bp2 = ax.boxplot(quintile_det, labels=['Q1\n(worst)', 'Q2', 'Q3', 'Q4', 'Q5\n(best)'],
                     patch_artist=True, showfliers=False)
    for patch, color in zip(bp2['boxes'], colors_q):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax.set_xlabel('Correlation quintile')
    ax.set_ylabel('Detection rate (MERSCOPE)')
    ax.set_title('F. MERSCOPE Detection by Correlation Quintile')

    fig.suptitle('Predictors of Poor Gene Correlation (snRNAseq vs MERSCOPE)\n'
                 f'{n_genes} genes across 24 subclasses',
                 fontsize=22, y=1.02)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "gene_corr_predictors_multipanel.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("    Saved: gene_corr_predictors_multipanel.png")

    # ── 5b. Scatter: snRNAseq vs MERSCOPE mean expression colored by r ──
    print("  5b. Mean expression scatter colored by correlation...")
    fig, ax = plt.subplots(figsize=(10, 9))
    sc2 = ax.scatter(df['mean_sn'], df['mean_mer'],
                     c=df['pearson_r'], cmap='RdYlGn', s=10, alpha=0.6,
                     vmin=-0.3, vmax=1.0, rasterized=True)
    ax.plot([0, df['mean_sn'].max()], [0, df['mean_sn'].max()],
            'k--', alpha=0.3, linewidth=1)
    ax.set_xlabel('Mean expression (snRNAseq)')
    ax.set_ylabel('Mean expression (MERSCOPE)')
    ax.set_title('Gene Expression: snRNAseq vs MERSCOPE\n(colored by cross-platform Pearson r)')
    plt.colorbar(sc2, ax=ax, label='Pearson r', shrink=0.8)

    # Label worst genes
    worst15 = df.nsmallest(15, 'pearson_r')
    for gene, row in worst15.iterrows():
        ax.annotate(gene, (row['mean_sn'], row['mean_mer']),
                     fontsize=7, alpha=0.8,
                     xytext=(4, 4), textcoords='offset points')

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "gene_mean_expr_scatter_by_corr.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("    Saved: gene_mean_expr_scatter_by_corr.png")

    # ── 5c. Biotype breakdown by correlation bin ──
    print("  5c. Biotype breakdown by correlation bin...")
    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    # Stacked bar: biotype proportions by r quintile
    ax = axes[0]
    biotype_counts = pd.crosstab(df['r_quintile'], df['biotype'])
    biotype_pcts = biotype_counts.div(biotype_counts.sum(axis=1), axis=0) * 100
    biotype_pcts.plot(kind='bar', stacked=True, ax=ax,
                       color=['#4C72B0', '#DD8452', '#55A868', '#C44E52', '#8172B3'],
                       edgecolor='black', linewidth=0.5)
    ax.set_xlabel('Correlation quintile')
    ax.set_ylabel('Percentage of genes (%)')
    ax.set_title('Biotype Composition by Correlation Quintile')
    ax.legend(title='Biotype', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=11)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

    # Histogram of r by biotype
    ax = axes[1]
    for bt in df['biotype'].unique():
        subset = df[df['biotype'] == bt]
        if len(subset) >= 10:
            ax.hist(subset['pearson_r'], bins=30, alpha=0.5, label=f"{bt} (n={len(subset)})",
                    edgecolor='black', linewidth=0.3)
    ax.set_xlabel('Pearson r')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of Per-Gene Correlation by Biotype')
    ax.legend(fontsize=11)
    ax.axvline(0, color='red', linestyle='--', alpha=0.3)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "gene_corr_by_biotype.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("    Saved: gene_corr_by_biotype.png")

    # ── 5d. MERSCOPE-enriched vs depleted genes ──
    print("  5d. Platform bias analysis...")
    fig, axes = plt.subplots(1, 2, figsize=(18, 8))

    # Panel A: log2FC vs r
    ax = axes[0]
    colors_bt = {'protein_coding': '#4C72B0', 'antisense': '#DD8452',
                 'lncRNA/novel': '#55A868', 'unannotated': '#C44E52', 'other': '#8172B3'}
    for bt, color in colors_bt.items():
        subset = df[df['biotype'] == bt]
        ax.scatter(subset['log2fc_mer_vs_sn'], subset['pearson_r'],
                   c=color, s=10, alpha=0.4, label=f"{bt} (n={len(subset)})",
                   rasterized=True)
    ax.set_xlabel('log2(MERSCOPE / snRNAseq) mean expression')
    ax.set_ylabel('Pearson r')
    ax.set_title('Platform Bias vs Correlation')
    ax.axvline(0, color='black', linestyle='-', alpha=0.3, linewidth=0.5)
    ax.axhline(0, color='red', linestyle='--', alpha=0.3)
    ax.axhline(0.5, color='gray', linestyle='--', alpha=0.3)
    ax.legend(fontsize=11, markerscale=3)

    # Panel B: genes systematically over/under-expressed in MERSCOPE
    ax = axes[1]
    merscope_enriched = df[df['log2fc_mer_vs_sn'] > 2].sort_values('log2fc_mer_vs_sn', ascending=False)
    merscope_depleted = df[df['log2fc_mer_vs_sn'] < -2].sort_values('log2fc_mer_vs_sn')

    # Top 20 most enriched and depleted
    top_enriched = merscope_enriched.head(20)
    top_depleted = merscope_depleted.head(20)
    combined = pd.concat([top_depleted, top_enriched])

    colors_bars = ['#d73027' if fc < 0 else '#4575b4' for fc in combined['log2fc_mer_vs_sn']]
    bars = ax.barh(range(len(combined)), combined['log2fc_mer_vs_sn'].values,
                    color=colors_bars, edgecolor='black', linewidth=0.3)
    ax.set_yticks(range(len(combined)))
    ax.set_yticklabels(combined.index, fontsize=9)
    ax.set_xlabel('log2(MERSCOPE / snRNAseq)')
    ax.set_title('Most Extreme Platform Biases')
    ax.axvline(0, color='black', linewidth=0.8)

    # Add r values
    for i, (gene, row) in enumerate(combined.iterrows()):
        fc = row['log2fc_mer_vs_sn']
        r_val = row['pearson_r']
        offset = 0.1 if fc > 0 else -0.1
        ha = 'left' if fc > 0 else 'right'
        ax.text(fc + offset, i, f"r={r_val:.2f}", va='center', ha=ha, fontsize=8)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "gene_platform_bias_analysis.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("    Saved: gene_platform_bias_analysis.png")

    # ── 5e. Detailed scatter for worst genes ──
    print("  5e. Scatter plots for worst genes...")
    worst12 = df.nsmallest(12, 'pearson_r')
    fig, axes_grid = plt.subplots(3, 4, figsize=(20, 15))
    axes_flat = axes_grid.flatten()

    for i, (gene, row) in enumerate(worst12.iterrows()):
        ax = axes_flat[i]
        sn_vals = pb_sn[gene].values
        mer_vals = pb_mer[gene].values

        # Color by class
        for j, ct in enumerate(pb_sn.index):
            cls = SUBCLASS_TO_CLASS.get(ct, "Unknown")
            color = {'Glutamatergic': '#00ADF8', 'GABAergic': '#F05A28',
                     'Non-neuronal': '#808080'}.get(cls, '#999')
            ax.scatter(sn_vals[j], mer_vals[j], c=color, s=50, edgecolors='black',
                       linewidths=0.5, zorder=5)
            ax.annotate(ct, (sn_vals[j], mer_vals[j]), fontsize=5,
                         xytext=(3, 3), textcoords='offset points')

        lims = [0, max(max(sn_vals), max(mer_vals)) * 1.15]
        ax.plot(lims, lims, 'k--', alpha=0.3, linewidth=0.8)
        ax.set_title(f"{gene}\nr={row['pearson_r']:.3f}, log2FC={row['log2fc_mer_vs_sn']:.1f}",
                     fontsize=12)
        ax.set_xlabel('snRNAseq', fontsize=10)
        ax.set_ylabel('MERSCOPE', fontsize=10)
        ax.tick_params(labelsize=8)

    fig.suptitle('Worst 12 Genes: Per-Cell-Type Expression\n(each dot = 1 cell type, colored by class)',
                 fontsize=20, y=1.02)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "worst_genes_celltype_scatter.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("    Saved: worst_genes_celltype_scatter.png")

    # ── 5f. Summary: what predicts bad correlation? ──
    print("  5f. Correlation matrix of gene properties with r...")
    props_for_corr = ['mean_sn', 'mean_mer', 'mean_both', 'var_sn', 'var_mer',
                      'abs_log2fc', 'det_rate_sn', 'det_rate_mer',
                      'entropy_sn', 'cv_sn', 'cv_mer']

    corr_with_r = {}
    for prop in props_for_corr:
        valid = df[[prop, 'pearson_r']].copy()
        valid[prop] = pd.to_numeric(valid[prop], errors='coerce')
        valid = valid.dropna()
        rho, p = spearmanr(valid[prop].values, valid['pearson_r'].values)
        corr_with_r[prop] = {'spearman_rho_with_r': rho, 'p_value': p}

    df_corr_summary = pd.DataFrame(corr_with_r).T.sort_values('spearman_rho_with_r')

    fig, ax = plt.subplots(figsize=(10, 7))
    colors_corr = ['#d73027' if v < 0 else '#4575b4' for v in df_corr_summary['spearman_rho_with_r']]
    bars = ax.barh(range(len(df_corr_summary)), df_corr_summary['spearman_rho_with_r'].values,
                    color=colors_corr, edgecolor='black', linewidth=0.5)
    ax.set_yticks(range(len(df_corr_summary)))
    pretty_labels = {
        'mean_sn': 'Mean expr (snRNAseq)',
        'mean_mer': 'Mean expr (MERSCOPE)',
        'mean_both': 'Mean expr (both)',
        'var_sn': 'Variance (snRNAseq)',
        'var_mer': 'Variance (MERSCOPE)',
        'abs_log2fc': '|log2FC| platform bias',
        'det_rate_sn': 'Detection rate (snRNAseq)',
        'det_rate_mer': 'Detection rate (MERSCOPE)',
        'entropy_sn': 'Expression entropy (snRNAseq)',
        'cv_sn': 'CV across types (snRNAseq)',
        'cv_mer': 'CV across types (MERSCOPE)',
    }
    ax.set_yticklabels([pretty_labels.get(x, x) for x in df_corr_summary.index], fontsize=14)
    ax.set_xlabel('Spearman rho with per-gene Pearson r', fontsize=16)
    ax.set_title('What Predicts Good Cross-Platform Gene Correlation?\n'
                 '(Spearman rank correlation of gene property with cross-platform r)',
                 fontsize=18)
    ax.axvline(0, color='black', linewidth=0.8)

    for i, (idx, row_val) in enumerate(df_corr_summary.iterrows()):
        v = row_val['spearman_rho_with_r']
        offset = 0.005 if v > 0 else -0.005
        ha = 'left' if v > 0 else 'right'
        ax.text(v + offset, i, f"{v:.3f}", va='center', ha=ha, fontsize=12, fontweight='bold')

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "gene_property_correlation_with_r.png",
                dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("    Saved: gene_property_correlation_with_r.png")

    # ------------------------------------------------------------------
    # 6. Save master table
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STEP 6: Saving results")
    print("=" * 70)

    df.to_csv(OUTPUT_DIR / "gene_properties_vs_correlation.csv")
    print(f"  Saved: gene_properties_vs_correlation.csv ({len(df)} genes)")

    df_corr_summary.to_csv(OUTPUT_DIR / "gene_property_predictors_of_r.csv")
    print(f"  Saved: gene_property_predictors_of_r.csv")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("SUMMARY: Key findings about poorly-performing genes")
    print("=" * 70)

    print(f"\n  Total genes analyzed: {n_genes}")
    print(f"  Genes with r < 0 (negative correlation): {(df['pearson_r'] < 0).sum()}")
    print(f"  Genes with r < 0.3: {(df['pearson_r'] < 0.3).sum()}")
    print(f"  Genes with r > 0.8: {(df['pearson_r'] > 0.8).sum()}")

    print(f"\n  Strongest predictors of poor correlation (Spearman rho with r):")
    for prop, row_val in df_corr_summary.iterrows():
        rho = row_val['spearman_rho_with_r']
        if abs(rho) > 0.1:
            direction = "↑ better r" if rho > 0 else "↓ worse r"
            print(f"    {pretty_labels.get(prop, prop):40s} rho={rho:+.3f} ({direction})")

    # Median properties by quintile
    print(f"\n  Median properties by correlation quintile:")
    for prop in ['mean_both', 'det_rate_mer', 'abs_log2fc', 'cv_sn']:
        print(f"\n    {pretty_labels.get(prop, prop)}:")
        for q in df['r_quintile'].cat.categories:
            subset = df[df['r_quintile'] == q]
            print(f"      {q}: {subset[prop].median():.4f}")

    elapsed = time.time() - t0
    print(f"\n  Total runtime: {elapsed:.1f}s")
    print("  Done!")


if __name__ == "__main__":
    main()
