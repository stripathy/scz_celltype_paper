#!/usr/bin/env python3
"""
Cross-platform generalization of gene correlation findings.

Computes per-gene snRNAseq vs spatial correlations for:
1. MERSCOPE 4K (already computed — load from CSV)
2. SEA-AD MERFISH (180 genes, ~1.9M cells)
3. SCZ Xenium (300 genes, ~500K cells across 24 samples)

Then asks: do the same genes perform poorly across platforms?
Do the same gene properties predict poor correlation?

Usage:
    python cross_platform_gene_corr.py
"""

import os
import sys
import time
import warnings
from pathlib import Path
from glob import glob

import anndata as ad
import matplotlib.pyplot as plt
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
MERFISH_PATH = PROJECT_DIR / "data/reference/SEAAD_MTG_MERFISH.2024-12-11.h5ad"
MERSCOPE_DIR = PROJECT_DIR / "output/merscope_h5ad"
XENIUM_DIR = PROJECT_DIR / "output/h5ad"
OUTPUT_DIR = PROJECT_DIR / "output/presentation"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Precomputed MERSCOPE results
MERSCOPE_GENE_CORR = OUTPUT_DIR / "snrnaseq_vs_merscope_gene_corr.csv"
MERSCOPE_PROPS = OUTPUT_DIR / "gene_properties_vs_correlation.csv"

# ── Constants ──────────────────────────────────────────────────────────
PLATFORM_COLORS = {
    "MERSCOPE 4K": "#2bae89",
    "SEA-AD MERFISH": "#e94560",
    "SCZ Xenium": "#4C72B0",
}


def compute_gene_correlations(pb_ref, pb_spatial, platform_name):
    """Compute per-gene Pearson r between reference and spatial pseudobulk."""
    shared_cts = sorted(set(pb_ref.index) & set(pb_spatial.index))
    shared_genes = sorted(set(pb_ref.columns) & set(pb_spatial.columns))

    print(f"  {platform_name}: {len(shared_genes)} shared genes, {len(shared_cts)} shared cell types")

    pb_ref_aligned = pb_ref.loc[shared_cts, shared_genes]
    pb_spatial_aligned = pb_spatial.loc[shared_cts, shared_genes]

    gene_corrs = {}
    for gene in shared_genes:
        sn_vals = pb_ref_aligned[gene].values
        sp_vals = pb_spatial_aligned[gene].values
        if np.std(sn_vals) < 1e-10 or np.std(sp_vals) < 1e-10:
            continue
        r, p = pearsonr(sn_vals, sp_vals)
        rho, p_sp = spearmanr(sn_vals, sp_vals)
        gene_corrs[gene] = {
            'pearson_r': r,
            'spearman_rho': rho,
            'mean_ref': np.mean(sn_vals),
            'mean_spatial': np.mean(sp_vals),
            'var_ref': np.var(sn_vals),
            'var_spatial': np.var(sp_vals),
        }

    df = pd.DataFrame(gene_corrs).T
    df.index.name = 'gene'
    return df


def main():
    t0 = time.time()

    # ------------------------------------------------------------------
    # 1. Load precomputed MERSCOPE 4K results
    # ------------------------------------------------------------------
    print("=" * 70)
    print("STEP 1: Loading precomputed MERSCOPE 4K results")
    print("=" * 70)

    df_merscope = pd.read_csv(MERSCOPE_GENE_CORR, index_col=0)
    df_merscope_props = pd.read_csv(MERSCOPE_PROPS, index_col=0)
    print(f"  MERSCOPE 4K: {len(df_merscope)} genes with correlations")

    # ------------------------------------------------------------------
    # 2. Load snRNAseq reference
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STEP 2: Loading snRNAseq reference")
    print("=" * 70)

    adata_sn = load_and_normalize_reference(str(SNRNASEQ_PATH),
                                              normalize=True, min_cells=0)

    # Compute snRNAseq pseudobulk for shared subclasses
    sn_subclass_col = "Subclass"
    valid_subclasses = [s for s in adata_sn.obs[sn_subclass_col].unique()
                        if s in SUBCLASS_TO_CLASS]
    adata_sn_sub = adata_sn[adata_sn.obs[sn_subclass_col].isin(valid_subclasses)].copy()

    print(f"\n  Computing snRNAseq pseudobulk ({len(valid_subclasses)} subclasses)...")
    pb_sn = compute_pseudobulk_mean(adata_sn_sub, sn_subclass_col)
    print(f"  snRNAseq pseudobulk shape: {pb_sn.shape}")

    del adata_sn  # free memory

    # ------------------------------------------------------------------
    # 3. SEA-AD MERFISH
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STEP 3: SEA-AD MERFISH (180 genes)")
    print("=" * 70)

    adata_merfish = ad.read_h5ad(MERFISH_PATH)
    print(f"  Shape: {adata_merfish.shape}")

    # MERFISH is already normalized (ln(spots/10K + 1))
    # Filter to valid subclasses
    merfish_subclass_col = "Subclass"
    valid_merfish = [s for s in adata_merfish.obs[merfish_subclass_col].unique()
                     if s in SUBCLASS_TO_CLASS]
    adata_merfish_sub = adata_merfish[
        adata_merfish.obs[merfish_subclass_col].isin(valid_merfish)
    ].copy()
    print(f"  Valid subclasses: {len(valid_merfish)}")
    print(f"  Cells after filtering: {adata_merfish_sub.shape[0]}")

    print("  Computing MERFISH pseudobulk...")
    pb_merfish = compute_pseudobulk_mean(adata_merfish_sub, merfish_subclass_col)
    print(f"  MERFISH pseudobulk shape: {pb_merfish.shape}")

    # Compute per-gene correlations
    df_merfish = compute_gene_correlations(pb_sn, pb_merfish, "SEA-AD MERFISH")
    print(f"  Valid gene correlations: {len(df_merfish)}")
    print(f"  Median Pearson r: {df_merfish['pearson_r'].median():.3f}")

    del adata_merfish, adata_merfish_sub  # free memory

    # ------------------------------------------------------------------
    # 4. SCZ Xenium (300 genes)
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STEP 4: SCZ Xenium (300 genes)")
    print("=" * 70)

    xenium_files = sorted(glob(str(XENIUM_DIR / "Br*_annotated.h5ad")))
    print(f"  Found {len(xenium_files)} Xenium sample files")

    xenium_adatas = []
    for f in xenium_files:
        fname = Path(f).name
        adata_x = ad.read_h5ad(f)

        # Filter to QC-pass cells
        if "corr_qc_pass" in adata_x.obs.columns:
            adata_x = adata_x[adata_x.obs["corr_qc_pass"] == True].copy()
        elif "qc_pass" in adata_x.obs.columns:
            adata_x = adata_x[adata_x.obs["qc_pass"] == True].copy()

        # Filter out Unassigned
        subclass_col = "corr_subclass_label" if "corr_subclass_label" in adata_x.obs.columns else "subclass_label"
        if subclass_col in adata_x.obs.columns:
            adata_x = adata_x[adata_x.obs[subclass_col] != "Unassigned"].copy()

        if adata_x.shape[0] > 0:
            xenium_adatas.append(adata_x)
            print(f"  {fname}: {adata_x.shape[0]} QC-pass cells, {adata_x.shape[1]} genes")

    adata_xenium = ad.concat(xenium_adatas, join="inner")
    print(f"\n  Combined Xenium shape: {adata_xenium.shape}")

    # Normalize Xenium (raw counts)
    X_sample_x = adata_xenium.X[:100, :100]
    if sp.issparse(X_sample_x):
        X_sample_x = X_sample_x.toarray()
    if np.all(X_sample_x == X_sample_x.astype(int)):
        print("  Normalizing raw counts...")
        sc.pp.normalize_total(adata_xenium, target_sum=1e4)
        sc.pp.log1p(adata_xenium)

    # Determine subclass column
    xenium_subclass_col = "corr_subclass_label" if "corr_subclass_label" in adata_xenium.obs.columns else "subclass_label"
    valid_xenium = [s for s in adata_xenium.obs[xenium_subclass_col].unique()
                    if s in SUBCLASS_TO_CLASS]
    adata_xenium_sub = adata_xenium[
        adata_xenium.obs[xenium_subclass_col].isin(valid_xenium)
    ].copy()
    print(f"  Valid subclasses: {len(valid_xenium)}")
    print(f"  Cells after filtering: {adata_xenium_sub.shape[0]}")

    print("  Computing Xenium pseudobulk...")
    pb_xenium = compute_pseudobulk_mean(adata_xenium_sub, xenium_subclass_col)
    print(f"  Xenium pseudobulk shape: {pb_xenium.shape}")

    # Compute per-gene correlations
    df_xenium = compute_gene_correlations(pb_sn, pb_xenium, "SCZ Xenium")
    print(f"  Valid gene correlations: {len(df_xenium)}")
    print(f"  Median Pearson r: {df_xenium['pearson_r'].median():.3f}")

    del adata_xenium, adata_xenium_sub, xenium_adatas

    # ------------------------------------------------------------------
    # 5. Cross-platform comparison
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STEP 5: Cross-platform gene correlation comparison")
    print("=" * 70)

    # Find genes shared across all three platforms
    merscope_genes = set(df_merscope.index)
    merfish_genes = set(df_merfish.index)
    xenium_genes = set(df_xenium.index)

    # Pairwise overlaps
    mer_mf = sorted(merscope_genes & merfish_genes)
    mer_xen = sorted(merscope_genes & xenium_genes)
    mf_xen = sorted(merfish_genes & xenium_genes)
    all_three = sorted(merscope_genes & merfish_genes & xenium_genes)

    print(f"\n  Gene overlaps:")
    print(f"    MERSCOPE 4K genes with valid r: {len(merscope_genes)}")
    print(f"    MERFISH genes with valid r: {len(merfish_genes)}")
    print(f"    Xenium genes with valid r: {len(xenium_genes)}")
    print(f"    MERSCOPE ∩ MERFISH: {len(mer_mf)}")
    print(f"    MERSCOPE ∩ Xenium: {len(mer_xen)}")
    print(f"    MERFISH ∩ Xenium: {len(mf_xen)}")
    print(f"    All three: {len(all_three)}")

    # Build combined table for shared genes
    # MERSCOPE vs MERFISH
    print(f"\n  Correlation of gene-level r: MERSCOPE 4K vs SEA-AD MERFISH ({len(mer_mf)} shared genes):")
    r_mer_mf, p_mer_mf = pearsonr(df_merscope.loc[mer_mf, 'pearson_r'], df_merfish.loc[mer_mf, 'pearson_r'])
    rho_mer_mf, _ = spearmanr(df_merscope.loc[mer_mf, 'pearson_r'], df_merfish.loc[mer_mf, 'pearson_r'])
    print(f"    Pearson r = {r_mer_mf:.3f}, Spearman rho = {rho_mer_mf:.3f}")

    # MERSCOPE vs Xenium
    print(f"\n  Correlation of gene-level r: MERSCOPE 4K vs SCZ Xenium ({len(mer_xen)} shared genes):")
    r_mer_xen, p_mer_xen = pearsonr(df_merscope.loc[mer_xen, 'pearson_r'], df_xenium.loc[mer_xen, 'pearson_r'])
    rho_mer_xen, _ = spearmanr(df_merscope.loc[mer_xen, 'pearson_r'], df_xenium.loc[mer_xen, 'pearson_r'])
    print(f"    Pearson r = {r_mer_xen:.3f}, Spearman rho = {rho_mer_xen:.3f}")

    # MERFISH vs Xenium
    print(f"\n  Correlation of gene-level r: SEA-AD MERFISH vs SCZ Xenium ({len(mf_xen)} shared genes):")
    r_mf_xen, p_mf_xen = pearsonr(df_merfish.loc[mf_xen, 'pearson_r'], df_xenium.loc[mf_xen, 'pearson_r'])
    rho_mf_xen, _ = spearmanr(df_merfish.loc[mf_xen, 'pearson_r'], df_xenium.loc[mf_xen, 'pearson_r'])
    print(f"    Pearson r = {r_mf_xen:.3f}, Spearman rho = {rho_mf_xen:.3f}")

    # Summary stats per platform
    print(f"\n  Summary statistics per platform:")
    for name, df_corr in [("MERSCOPE 4K", df_merscope), ("SEA-AD MERFISH", df_merfish), ("SCZ Xenium", df_xenium)]:
        print(f"\n    {name} ({len(df_corr)} genes):")
        print(f"      Median Pearson r:  {df_corr['pearson_r'].median():.3f}")
        print(f"      Mean Pearson r:    {df_corr['pearson_r'].mean():.3f}")
        print(f"      Fraction r > 0.5:  {(df_corr['pearson_r'] > 0.5).mean():.1%}")
        print(f"      Fraction r > 0.8:  {(df_corr['pearson_r'] > 0.8).mean():.1%}")
        print(f"      Fraction r < 0:    {(df_corr['pearson_r'] < 0).mean():.1%}")

    # ------------------------------------------------------------------
    # 6. Do the same gene properties predict poor r across platforms?
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STEP 6: Do gene properties generalize as predictors?")
    print("=" * 70)

    # Use MERSCOPE gene properties and check if they predict r in other platforms
    print("\n  Using MERSCOPE 4K gene properties to predict r in other platforms:")
    props_to_test = ['mean_sn', 'mean_mer', 'var_sn', 'var_mer',
                     'abs_log2fc', 'det_rate_sn', 'det_rate_mer',
                     'entropy_sn', 'cv_sn', 'cv_mer']

    pretty_labels = {
        'mean_sn': 'Mean expr (snRNAseq)',
        'mean_mer': 'Mean expr (MERSCOPE)',
        'var_sn': 'Variance (snRNAseq)',
        'var_mer': 'Variance (MERSCOPE)',
        'abs_log2fc': '|log2FC| platform bias',
        'det_rate_sn': 'Detection rate (snRNAseq)',
        'det_rate_mer': 'Detection rate (MERSCOPE)',
        'entropy_sn': 'Expression entropy (snRNAseq)',
        'cv_sn': 'CV across types (snRNAseq)',
        'cv_mer': 'CV across types (MERSCOPE)',
    }

    # For MERFISH shared genes
    print(f"\n  MERSCOPE properties predicting MERFISH r ({len(mer_mf)} genes):")
    merfish_pred = {}
    for prop in props_to_test:
        shared = [g for g in mer_mf if g in df_merscope_props.index]
        prop_vals = pd.to_numeric(df_merscope_props.loc[shared, prop], errors='coerce')
        merfish_r = df_merfish.loc[shared, 'pearson_r']
        valid = prop_vals.notna() & merfish_r.notna()
        if valid.sum() > 10:
            rho, p = spearmanr(prop_vals[valid].values, merfish_r[valid].values)
            merfish_pred[prop] = rho
            print(f"    {pretty_labels.get(prop, prop):40s} rho = {rho:+.3f}")

    # For Xenium shared genes
    print(f"\n  MERSCOPE properties predicting Xenium r ({len(mer_xen)} genes):")
    xenium_pred = {}
    for prop in props_to_test:
        shared = [g for g in mer_xen if g in df_merscope_props.index]
        prop_vals = pd.to_numeric(df_merscope_props.loc[shared, prop], errors='coerce')
        xenium_r = df_xenium.loc[shared, 'pearson_r']
        valid = prop_vals.notna() & xenium_r.notna()
        if valid.sum() > 10:
            rho, p = spearmanr(prop_vals[valid].values, xenium_r[valid].values)
            xenium_pred[prop] = rho
            print(f"    {pretty_labels.get(prop, prop):40s} rho = {rho:+.3f}")

    # Also: reference-only properties (independent of spatial platform)
    # snRNAseq variance and entropy should predict r in ALL platforms
    print(f"\n  snRNAseq-only properties predicting r (platform-independent predictors):")
    ref_props = ['mean_sn', 'var_sn', 'det_rate_sn', 'entropy_sn', 'cv_sn']
    print(f"  {'Property':<35s} {'MERSCOPE':>10s} {'MERFISH':>10s} {'Xenium':>10s}")
    print(f"  {'-'*35} {'-'*10} {'-'*10} {'-'*10}")

    for prop in ref_props:
        # MERSCOPE
        vals_m = pd.to_numeric(df_merscope_props.loc[df_merscope.index, prop], errors='coerce')
        valid_m = vals_m.notna()
        rho_m, _ = spearmanr(vals_m[valid_m].values, df_merscope.loc[vals_m[valid_m].index, 'pearson_r'].values) if valid_m.sum() > 10 else (np.nan, np.nan)

        # MERFISH
        shared_mf = [g for g in mer_mf if g in df_merscope_props.index]
        vals_mf = pd.to_numeric(df_merscope_props.loc[shared_mf, prop], errors='coerce')
        valid_mf = vals_mf.notna()
        rho_mf, _ = spearmanr(vals_mf[valid_mf].values, df_merfish.loc[vals_mf[valid_mf].index, 'pearson_r'].values) if valid_mf.sum() > 10 else (np.nan, np.nan)

        # Xenium
        shared_xn = [g for g in mer_xen if g in df_merscope_props.index]
        vals_xn = pd.to_numeric(df_merscope_props.loc[shared_xn, prop], errors='coerce')
        valid_xn = vals_xn.notna()
        rho_xn, _ = spearmanr(vals_xn[valid_xn].values, df_xenium.loc[vals_xn[valid_xn].index, 'pearson_r'].values) if valid_xn.sum() > 10 else (np.nan, np.nan)

        print(f"  {pretty_labels.get(prop, prop):<35s} {rho_m:>+10.3f} {rho_mf:>+10.3f} {rho_xn:>+10.3f}")

    # ------------------------------------------------------------------
    # 7. Figures
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STEP 7: Creating figures")
    print("=" * 70)

    plt.rcParams.update({
        "font.size": 16, "axes.titlesize": 20, "axes.labelsize": 18,
        "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 13,
        "figure.dpi": 150,
    })

    # ── 7a. Pairwise scatter of gene-level r across platforms ──
    print("  7a. Pairwise scatter of gene-level r...")
    fig, axes = plt.subplots(1, 3, figsize=(21, 6.5))

    # MERSCOPE vs MERFISH
    ax = axes[0]
    ax.scatter(df_merscope.loc[mer_mf, 'pearson_r'],
               df_merfish.loc[mer_mf, 'pearson_r'],
               s=20, alpha=0.5, c='#2bae89', edgecolors='none', rasterized=True)
    ax.plot([-0.6, 1], [-0.6, 1], 'k--', alpha=0.3)
    ax.axhline(0, color='gray', linewidth=0.5, alpha=0.3)
    ax.axvline(0, color='gray', linewidth=0.5, alpha=0.3)
    ax.set_xlabel('MERSCOPE 4K: Pearson r')
    ax.set_ylabel('SEA-AD MERFISH: Pearson r')
    ax.set_title(f'MERSCOPE vs MERFISH\n(n={len(mer_mf)} genes, r={r_mer_mf:.3f}, ρ={rho_mer_mf:.3f})')
    ax.set_xlim(-0.6, 1.05)
    ax.set_ylim(-0.6, 1.05)

    # MERSCOPE vs Xenium
    ax = axes[1]
    ax.scatter(df_merscope.loc[mer_xen, 'pearson_r'],
               df_xenium.loc[mer_xen, 'pearson_r'],
               s=20, alpha=0.5, c='#4C72B0', edgecolors='none', rasterized=True)
    ax.plot([-0.6, 1], [-0.6, 1], 'k--', alpha=0.3)
    ax.axhline(0, color='gray', linewidth=0.5, alpha=0.3)
    ax.axvline(0, color='gray', linewidth=0.5, alpha=0.3)
    ax.set_xlabel('MERSCOPE 4K: Pearson r')
    ax.set_ylabel('SCZ Xenium: Pearson r')
    ax.set_title(f'MERSCOPE vs Xenium\n(n={len(mer_xen)} genes, r={r_mer_xen:.3f}, ρ={rho_mer_xen:.3f})')
    ax.set_xlim(-0.6, 1.05)
    ax.set_ylim(-0.6, 1.05)

    # MERFISH vs Xenium
    ax = axes[2]
    ax.scatter(df_merfish.loc[mf_xen, 'pearson_r'],
               df_xenium.loc[mf_xen, 'pearson_r'],
               s=20, alpha=0.5, c='#e94560', edgecolors='none', rasterized=True)
    ax.plot([-0.6, 1], [-0.6, 1], 'k--', alpha=0.3)
    ax.axhline(0, color='gray', linewidth=0.5, alpha=0.3)
    ax.axvline(0, color='gray', linewidth=0.5, alpha=0.3)
    ax.set_xlabel('SEA-AD MERFISH: Pearson r')
    ax.set_ylabel('SCZ Xenium: Pearson r')
    ax.set_title(f'MERFISH vs Xenium\n(n={len(mf_xen)} genes, r={r_mf_xen:.3f}, ρ={rho_mf_xen:.3f})')
    ax.set_xlim(-0.6, 1.05)
    ax.set_ylim(-0.6, 1.05)

    fig.suptitle('Cross-Platform Generalization of Per-Gene Correlation\n'
                 '(each point = 1 gene, comparing its snRNAseq vs spatial Pearson r across platforms)',
                 fontsize=18, y=1.05)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "cross_platform_gene_r_scatter.png", dpi=200, bbox_inches='tight')
    plt.close(fig)
    print("    Saved: cross_platform_gene_r_scatter.png")

    # ── 7b. Histogram comparison of r distributions ──
    print("  7b. Histogram comparison of r distributions...")
    fig, ax = plt.subplots(figsize=(12, 7))

    for name, df_corr, color in [("MERSCOPE 4K", df_merscope, "#2bae89"),
                                   ("SEA-AD MERFISH", df_merfish, "#e94560"),
                                   ("SCZ Xenium", df_xenium, "#4C72B0")]:
        ax.hist(df_corr['pearson_r'], bins=40, alpha=0.45, color=color,
                edgecolor='black', linewidth=0.3,
                label=f"{name} (n={len(df_corr)}, med={df_corr['pearson_r'].median():.3f})")
        ax.axvline(df_corr['pearson_r'].median(), color=color, linestyle='--', linewidth=2, alpha=0.8)

    ax.set_xlabel('Pearson r (gene vs snRNAseq, across cell types)')
    ax.set_ylabel('Number of genes')
    ax.set_title('Distribution of Per-Gene Cross-Platform Correlation\nacross 3 spatial platforms')
    ax.legend(fontsize=14, loc='upper left')
    ax.axvline(0, color='black', linewidth=0.8, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "cross_platform_gene_r_histograms.png", dpi=200, bbox_inches='tight')
    plt.close(fig)
    print("    Saved: cross_platform_gene_r_histograms.png")

    # ── 7c. Grouped bar: snRNAseq property predicts r across all platforms ──
    print("  7c. Property generalization bar chart...")

    ref_props = ['var_sn', 'cv_sn', 'entropy_sn', 'mean_sn', 'det_rate_sn']
    prop_labels = ['Variance\n(snRNAseq)', 'CV across\ntypes', 'Expression\nentropy',
                   'Mean\nexpression', 'Detection\nrate']

    platform_rhos = {pname: [] for pname in ["MERSCOPE 4K", "SEA-AD MERFISH", "SCZ Xenium"]}

    for prop in ref_props:
        # MERSCOPE
        vals = pd.to_numeric(df_merscope_props.loc[df_merscope.index, prop], errors='coerce')
        valid = vals.notna()
        rho_m, _ = spearmanr(vals[valid].values, df_merscope.loc[vals[valid].index, 'pearson_r'].values)
        platform_rhos["MERSCOPE 4K"].append(rho_m)

        # MERFISH
        shared = [g for g in mer_mf if g in df_merscope_props.index]
        vals_mf = pd.to_numeric(df_merscope_props.loc[shared, prop], errors='coerce')
        valid_mf = vals_mf.notna()
        rho_mf, _ = spearmanr(vals_mf[valid_mf].values, df_merfish.loc[vals_mf[valid_mf].index, 'pearson_r'].values)
        platform_rhos["SEA-AD MERFISH"].append(rho_mf)

        # Xenium
        shared_x = [g for g in mer_xen if g in df_merscope_props.index]
        vals_xn = pd.to_numeric(df_merscope_props.loc[shared_x, prop], errors='coerce')
        valid_xn = vals_xn.notna()
        rho_xn, _ = spearmanr(vals_xn[valid_xn].values, df_xenium.loc[vals_xn[valid_xn].index, 'pearson_r'].values)
        platform_rhos["SCZ Xenium"].append(rho_xn)

    fig, ax = plt.subplots(figsize=(14, 7))
    x = np.arange(len(ref_props))
    width = 0.25

    for i, (pname, rhos) in enumerate(platform_rhos.items()):
        color = PLATFORM_COLORS[pname]
        bars = ax.bar(x + i * width, rhos, width, label=pname, color=color,
                       edgecolor='black', linewidth=0.5, alpha=0.8)
        for j, v in enumerate(rhos):
            ax.text(x[j] + i * width, v + 0.01 * np.sign(v), f"{v:.2f}",
                    ha='center', va='bottom' if v > 0 else 'top', fontsize=10, fontweight='bold')

    ax.set_xticks(x + width)
    ax.set_xticklabels(prop_labels, fontsize=14)
    ax.set_ylabel('Spearman rho with per-gene Pearson r', fontsize=16)
    ax.set_title('snRNAseq Gene Properties Predict Cross-Platform Correlation\nAcross All 3 Spatial Platforms',
                 fontsize=20)
    ax.legend(fontsize=14)
    ax.axhline(0, color='black', linewidth=0.8)
    ax.set_ylim(-0.6, 0.7)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "cross_platform_property_generalization.png", dpi=200, bbox_inches='tight')
    plt.close(fig)
    print("    Saved: cross_platform_property_generalization.png")

    # ── 7d. Genes shared across all 3: worst and best performers ──
    if len(all_three) > 20:
        print(f"  7d. Shared genes across all 3 platforms (n={len(all_three)})...")

        df_shared = pd.DataFrame({
            'r_MERSCOPE': df_merscope.loc[all_three, 'pearson_r'],
            'r_MERFISH': df_merfish.loc[all_three, 'pearson_r'],
            'r_Xenium': df_xenium.loc[all_three, 'pearson_r'],
        })
        df_shared['mean_r'] = df_shared.mean(axis=1)
        df_shared = df_shared.sort_values('mean_r')

        # Heatmap-style figure
        n_show = min(30, len(df_shared))
        worst = df_shared.head(n_show)
        best = df_shared.tail(n_show)
        combined_show = pd.concat([worst, best])

        fig, ax = plt.subplots(figsize=(10, max(10, len(combined_show) * 0.22)))
        y_pos = range(len(combined_show))

        for i, (gene, row) in enumerate(combined_show.iterrows()):
            ax.scatter(row['r_MERSCOPE'], i, c=PLATFORM_COLORS['MERSCOPE 4K'],
                       s=60, zorder=5, edgecolors='black', linewidths=0.3)
            ax.scatter(row['r_MERFISH'], i, c=PLATFORM_COLORS['SEA-AD MERFISH'],
                       s=60, zorder=5, edgecolors='black', linewidths=0.3)
            ax.scatter(row['r_Xenium'], i, c=PLATFORM_COLORS['SCZ Xenium'],
                       s=60, zorder=5, edgecolors='black', linewidths=0.3)

        ax.set_yticks(range(len(combined_show)))
        ax.set_yticklabels(combined_show.index, fontsize=9)
        ax.set_xlabel('Pearson r (vs snRNAseq)')
        ax.axvline(0, color='red', linestyle='--', alpha=0.3)
        ax.axvline(0.5, color='gray', linestyle='--', alpha=0.2)

        # Add separator
        if n_show < len(df_shared):
            ax.axhline(n_show - 0.5, color='black', linewidth=1.5, linestyle='-')
            ax.text(0.95, n_show - 0.5, '— gap —', ha='right', va='center',
                    fontsize=10, style='italic', transform=ax.get_yaxis_transform())

        # Legend
        for pname, color in PLATFORM_COLORS.items():
            ax.scatter([], [], c=color, s=60, label=pname, edgecolors='black', linewidths=0.3)
        ax.legend(loc='lower right', fontsize=12)

        ax.set_title(f'Worst {n_show} & Best {n_show} Genes Shared Across All 3 Platforms\n'
                     f'(n={len(all_three)} genes shared)',
                     fontsize=16)
        fig.tight_layout()
        fig.savefig(OUTPUT_DIR / "cross_platform_shared_genes_worst_best.png",
                    dpi=200, bbox_inches='tight')
        plt.close(fig)
        print("    Saved: cross_platform_shared_genes_worst_best.png")

        # Print table
        print(f"\n  Worst {n_show} genes across all 3 platforms:")
        print(f"  {'Gene':<15s} {'MERSCOPE':>10s} {'MERFISH':>10s} {'Xenium':>10s} {'Mean r':>8s}")
        for gene, row in worst.iterrows():
            print(f"  {gene:<15s} {row['r_MERSCOPE']:>10.3f} {row['r_MERFISH']:>10.3f} "
                  f"{row['r_Xenium']:>10.3f} {row['mean_r']:>8.3f}")

        print(f"\n  Best {n_show} genes across all 3 platforms:")
        for gene, row in best.iterrows():
            print(f"  {gene:<15s} {row['r_MERSCOPE']:>10.3f} {row['r_MERFISH']:>10.3f} "
                  f"{row['r_Xenium']:>10.3f} {row['mean_r']:>8.3f}")

        # Correlation of mean_r across platforms
        print(f"\n  Consistency: fraction of genes that are in the same quintile across platforms:")
        for q_label, q_range in [('worst 20%', (0, 0.2)), ('best 20%', (0.8, 1.0))]:
            mer_q = df_shared['r_MERSCOPE'].rank(pct=True).between(*q_range)
            mf_q = df_shared['r_MERFISH'].rank(pct=True).between(*q_range)
            xen_q = df_shared['r_Xenium'].rank(pct=True).between(*q_range)

            # How many are in worst/best quintile in ALL platforms?
            all_q = (mer_q & mf_q & xen_q).sum()
            any_two = ((mer_q & mf_q) | (mer_q & xen_q) | (mf_q & xen_q)).sum()
            expected = len(df_shared) * 0.2 * 0.2 * 0.2  # random chance

            print(f"    {q_label}: {all_q}/{len(df_shared)} in all 3 "
                  f"(expected by chance: {expected:.1f}), "
                  f"{any_two}/{len(df_shared)} in ≥2")

    # ------------------------------------------------------------------
    # 8. Save results
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("STEP 8: Saving results")
    print("=" * 70)

    df_merfish.to_csv(OUTPUT_DIR / "snrnaseq_vs_merfish_gene_corr.csv")
    df_xenium.to_csv(OUTPUT_DIR / "snrnaseq_vs_xenium_gene_corr.csv")
    print(f"  Saved: snrnaseq_vs_merfish_gene_corr.csv ({len(df_merfish)} genes)")
    print(f"  Saved: snrnaseq_vs_xenium_gene_corr.csv ({len(df_xenium)} genes)")

    if len(all_three) > 0:
        df_shared_out = pd.DataFrame({
            'r_MERSCOPE': df_merscope.reindex(all_three)['pearson_r'],
            'r_MERFISH': df_merfish.reindex(all_three)['pearson_r'],
            'r_Xenium': df_xenium.reindex(all_three)['pearson_r'],
        })
        df_shared_out['mean_r'] = df_shared_out.mean(axis=1)
        df_shared_out.to_csv(OUTPUT_DIR / "cross_platform_shared_gene_corr.csv")
        print(f"  Saved: cross_platform_shared_gene_corr.csv ({len(df_shared_out)} genes)")

    elapsed = time.time() - t0
    print(f"\n  Total runtime: {elapsed:.1f}s")
    print("  Done!")


if __name__ == "__main__":
    main()
