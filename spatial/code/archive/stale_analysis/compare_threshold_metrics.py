#!/usr/bin/env python3
"""
Compare margin threshold variants on three cross-platform metrics.

Metrics (all Pearson):
  1. snRNAseq logFC concordance: crumblr SCZ logFC vs Nicole's meta-analysis betas
  2. Proportion similarity: median cell-type proportions, Xenium vs MERFISH
  3. Depth similarity: median predicted_norm_depth per cell type, Xenium vs MERFISH

Variants: corr (1% pctl), pctl01, pctl05, pctl10

Output:
  output/presentation/margin_threshold_metrics_comparison.png
  output/presentation/margin_depth_similarity.csv
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

BASE_DIR = os.path.expanduser("~/Github/SCZ_Xenium")
sys.path.insert(0, os.path.join(BASE_DIR, "code", "analysis"))
from config import PRESENTATION_DIR

CRUMBLR_DIR = os.path.join(BASE_DIR, "output", "crumblr")
H5AD_DIR = os.path.join(BASE_DIR, "output", "h5ad")
MERFISH_PATH = os.path.join(BASE_DIR, "output", "merfish_benchmark", "merfish_reclassified.h5ad")
NICOLE_PATH = os.path.join(BASE_DIR, "data", "nicole_scz_snrnaseq_betas", "scz_coefs.xlsx")
OUT_DIR = PRESENTATION_DIR

VARIANTS = {
    'corr': {'label': 'corr (1%)', 'suffix': '_corr'},
    'pctl01': {'label': '1st pctl', 'suffix': '_pctl01'},
    'pctl05': {'label': '5th pctl', 'suffix': '_pctl05'},
    'pctl10': {'label': '10th pctl', 'suffix': '_pctl10'},
}


# ──────────────────────────────────────────────────────────────────────
# Metric 1: snRNAseq logFC concordance
# ──────────────────────────────────────────────────────────────────────

def compute_snrnaseq_concordance():
    """Spearman r between crumblr supertype logFC and Nicole's snRNAseq betas."""
    print("\n=== Metric 1: snRNAseq logFC concordance (Spearman) ===")
    nicole = pd.read_excel(NICOLE_PATH)
    nicole = nicole.rename(columns={'CellType': 'celltype', 'estimate': 'beta_snrnaseq'})

    results = {}
    for var_name, var_info in VARIANTS.items():
        suffix = var_info['suffix']
        fpath = os.path.join(CRUMBLR_DIR, f"crumblr_results_supertype{suffix}.csv")
        if not os.path.exists(fpath):
            print(f"  {var_name}: MISSING {os.path.basename(fpath)}")
            continue

        xen = pd.read_csv(fpath)
        merged = pd.merge(
            nicole[['celltype', 'beta_snrnaseq']],
            xen[['celltype', 'logFC']].rename(columns={'logFC': 'logFC_xenium'}),
            on='celltype', how='inner'
        )
        rho, pval = pearsonr(merged['beta_snrnaseq'], merged['logFC_xenium'])
        results[var_name] = {'r': rho, 'p': pval, 'n': len(merged)}
        print(f"  {var_name}: r={rho:.3f} (p={pval:.3e}, n={len(merged)})")

    return results


# ──────────────────────────────────────────────────────────────────────
# Metric 2: Proportion similarity (Xenium vs MERFISH)
# ──────────────────────────────────────────────────────────────────────

def compute_proportion_similarity():
    """Spearman r of median subclass proportions between Xenium and MERFISH."""
    print("\n=== Metric 2: Proportion similarity (Spearman) ===")

    # MERFISH reference proportions from crumblr input (or compute from h5ad)
    # Use crumblr input files which have per-donor cell counts
    merfish_props = _get_merfish_proportions()

    results = {}
    for var_name, var_info in VARIANTS.items():
        suffix = var_info['suffix']
        fpath = os.path.join(CRUMBLR_DIR, f"crumblr_input_subclass{suffix}.csv")
        if not os.path.exists(fpath):
            print(f"  {var_name}: MISSING {os.path.basename(fpath)}")
            continue

        xen = pd.read_csv(fpath)
        xen_props = _compute_median_proportions(xen)

        merged = pd.merge(merfish_props, xen_props, on='celltype',
                          suffixes=('_merfish', '_xenium'), how='inner')
        rho, pval = pearsonr(merged['prop_merfish'], merged['prop_xenium'])
        results[var_name] = {'r': rho, 'p': pval, 'n': len(merged)}
        n_cells = xen['count'].sum()
        results[var_name]['n_cells'] = n_cells
        print(f"  {var_name}: r={rho:.3f} (p={pval:.3e}, n_types={len(merged)}, n_cells={n_cells:,})")

    return results


def _get_merfish_proportions():
    """Compute median per-donor subclass proportions from MERFISH."""
    print("  Loading MERFISH proportions...")
    adata = ad.read_h5ad(MERFISH_PATH, backed='r')
    obs = adata.obs
    # Filter to cortical neurons + glia (exclude WM-like using depth)
    # Use Subclass as the label
    df = obs[['Donor ID', 'Subclass']].copy()
    df.columns = ['donor', 'celltype']
    df['celltype'] = df['celltype'].astype(str)
    adata.file.close()

    counts = df.groupby(['donor', 'celltype']).size().reset_index(name='count')
    totals = df.groupby('donor')['celltype'].count().reset_index(name='total')
    counts = counts.merge(totals, on='donor')
    counts['prop'] = counts['count'] / counts['total']
    median_props = counts.groupby('celltype')['prop'].median().reset_index(name='prop')
    print(f"  MERFISH: {len(median_props)} subclasses from {counts['donor'].nunique()} donors")
    return median_props


def _compute_median_proportions(df):
    """Compute median per-donor proportions from crumblr input CSV."""
    if 'total' not in df.columns:
        totals = df.groupby('donor')['count'].sum().reset_index(name='total')
        df = df.merge(totals, on='donor')
    df = df.copy()
    df['prop'] = df['count'] / df['total']
    return df.groupby('celltype')['prop'].median().reset_index(name='prop')


# ──────────────────────────────────────────────────────────────────────
# Metric 3: Depth similarity (Xenium vs MERFISH)
# ──────────────────────────────────────────────────────────────────────

def compute_depth_similarity():
    """Spearman r of median depth per subclass between Xenium and MERFISH."""
    print("\n=== Metric 3: Depth similarity (Spearman) ===")

    # MERFISH reference: median predicted_norm_depth per subclass
    merfish_depth = _get_merfish_depth()

    # Load all Xenium h5ad files once, extract needed columns
    xenium_obs = _load_xenium_obs_all()

    results = {}
    for var_name, var_info in VARIANTS.items():
        # Filter Xenium cells based on the variant's QC gate
        xen_filtered = _apply_variant_filter(xenium_obs, var_name)
        xen_depth = xen_filtered.groupby('subclass_label')['predicted_norm_depth'].median().reset_index()
        xen_depth.columns = ['celltype', 'depth']

        merged = pd.merge(merfish_depth, xen_depth, on='celltype',
                          suffixes=('_merfish', '_xenium'), how='inner')
        rho, pval = pearsonr(merged['depth_merfish'], merged['depth_xenium'])
        results[var_name] = {'r': rho, 'p': pval, 'n': len(merged)}
        print(f"  {var_name}: r={rho:.3f} (p={pval:.3e}, n_types={len(merged)})")

    # Save depth comparison for the default variant
    xen_default = _apply_variant_filter(xenium_obs, 'corr')
    xen_depth_default = xen_default.groupby('subclass_label')['predicted_norm_depth'].median().reset_index()
    xen_depth_default.columns = ['celltype', 'depth']
    merged_save = pd.merge(merfish_depth, xen_depth_default, on='celltype',
                           suffixes=('_merfish', '_xenium'), how='inner')
    merged_save.to_csv(os.path.join(OUT_DIR, 'margin_depth_similarity.csv'), index=False)
    print(f"  Saved margin_depth_similarity.csv")

    return results


def _get_merfish_depth():
    """Median predicted_norm_depth per Subclass from MERFISH."""
    print("  Loading MERFISH depth data...")
    adata = ad.read_h5ad(MERFISH_PATH, backed='r')
    obs = adata.obs[['Subclass', 'predicted_norm_depth']].copy()
    obs.columns = ['celltype', 'depth']
    obs['celltype'] = obs['celltype'].astype(str)
    obs['depth'] = obs['depth'].astype(float)
    adata.file.close()
    median_depth = obs.groupby('celltype')['depth'].median().reset_index()
    print(f"  MERFISH depth: {len(median_depth)} subclasses")
    return median_depth


def _load_xenium_obs_all():
    """Load obs from all Xenium h5ad files (just the columns we need)."""
    print("  Loading Xenium obs data...")
    h5ad_files = sorted(glob.glob(os.path.join(H5AD_DIR, "*_annotated.h5ad")))
    EXCLUDE = set()  # No samples excluded
    all_obs = []
    for fpath in h5ad_files:
        sid = os.path.basename(fpath).replace("_annotated.h5ad", "")
        if sid in EXCLUDE:
            continue
        adata = ad.read_h5ad(fpath, backed='r')
        obs = adata.obs
        cols_needed = ['sample_id', 'subclass_label', 'predicted_norm_depth',
                       'corr_subclass_margin', 'qc_pass', 'corr_qc_pass',
                       'doublet_suspect', 'spatial_domain', 'layer']
        sub = obs[cols_needed].copy()
        all_obs.append(sub)
        adata.file.close()
        print(f"    {sid}: {len(sub):,} cells")

    combined = pd.concat(all_obs, ignore_index=True)
    print(f"  Total Xenium: {len(combined):,} cells")
    return combined


def _apply_variant_filter(obs, var_name):
    """Apply QC filter for a given variant, returning cortical non-WM cells."""
    qc_pass = obs['qc_pass'].astype(bool)
    not_doublet = ~obs['doublet_suspect'].astype(bool)
    cortical = (obs['spatial_domain'] == 'Cortical') & (obs['layer'].astype(str) != 'WM')
    margin = obs['corr_subclass_margin'].astype(float)

    if var_name == 'corr':
        # Use existing corr_qc_pass column (Step 01 + 02b, 1st pctl)
        mask = obs['corr_qc_pass'].astype(bool) & cortical
    else:
        # Parse percentile from name
        pctl_val = int(var_name.replace('pctl', ''))
        # Compute per-sample percentile threshold
        thresholds = obs.loc[qc_pass, :].groupby('sample_id')['corr_subclass_margin'].quantile(pctl_val / 100)
        # Map threshold per cell
        sample_thresh = obs['sample_id'].map(thresholds).astype(float)
        margin_pass = margin >= sample_thresh
        mask = qc_pass & not_doublet & margin_pass & cortical

    return obs.loc[mask].copy()


# ──────────────────────────────────────────────────────────────────────
# Combined figure
# ──────────────────────────────────────────────────────────────────────

def plot_combined_metrics(snrnaseq_results, prop_results, depth_results):
    """4-panel comparison: logFC concordance, proportion similarity, depth similarity, cells retained."""
    print("\n=== Plotting combined metrics comparison ===")

    var_order = ['corr', 'pctl01', 'pctl05', 'pctl10']
    labels = [VARIANTS[v]['label'] for v in var_order]
    colors = ['#4C72B0', '#55A868', '#C44E52', '#8172B2']

    fig, axes = plt.subplots(1, 4, figsize=(20, 5))

    # Panel 1: snRNAseq logFC Pearson r
    ax = axes[0]
    rhos = [snrnaseq_results[v]['r'] for v in var_order]
    bars = ax.bar(range(len(var_order)), rhos, color=colors, edgecolor='white', linewidth=1.5)
    for i, (bar, rho) in enumerate(zip(bars, rhos)):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                f'{rho:.3f}', ha='center', va='bottom', fontsize=12, fontweight='bold')
    ax.set_xticks(range(len(var_order)))
    ax.set_xticklabels(labels, fontsize=12)
    ax.set_ylabel('Pearson r', fontsize=14)
    ax.set_title('snRNAseq logFC\nconcordance', fontsize=14, fontweight='bold')
    ax.set_ylim(0, max(rhos) * 1.15)

    # Panel 2: Proportion similarity
    ax = axes[1]
    rhos = [prop_results[v]['r'] for v in var_order]
    bars = ax.bar(range(len(var_order)), rhos, color=colors, edgecolor='white', linewidth=1.5)
    for i, (bar, rho) in enumerate(zip(bars, rhos)):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                f'{rho:.3f}', ha='center', va='bottom', fontsize=12, fontweight='bold')
    ax.set_xticks(range(len(var_order)))
    ax.set_xticklabels(labels, fontsize=12)
    ax.set_ylabel('Pearson r', fontsize=14)
    ax.set_title('MERFISH proportion\nsimilarity', fontsize=14, fontweight='bold')
    ax.set_ylim(0, max(rhos) * 1.15)

    # Panel 3: Depth similarity
    ax = axes[2]
    rhos = [depth_results[v]['r'] for v in var_order]
    bars = ax.bar(range(len(var_order)), rhos, color=colors, edgecolor='white', linewidth=1.5)
    for i, (bar, rho) in enumerate(zip(bars, rhos)):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                f'{rho:.3f}', ha='center', va='bottom', fontsize=12, fontweight='bold')
    ax.set_xticks(range(len(var_order)))
    ax.set_xticklabels(labels, fontsize=12)
    ax.set_ylabel('Pearson r', fontsize=14)
    ax.set_title('MERFISH depth\nsimilarity', fontsize=14, fontweight='bold')
    ax.set_ylim(0, max(rhos) * 1.15)

    # Panel 4: Cells retained
    ax = axes[3]
    n_cells = [prop_results[v]['n_cells'] for v in var_order]
    bars = ax.bar(range(len(var_order)), [n/1e6 for n in n_cells],
                  color=colors, edgecolor='white', linewidth=1.5)
    for i, (bar, n) in enumerate(zip(bars, n_cells)):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f'{n/1e6:.2f}M', ha='center', va='bottom', fontsize=12, fontweight='bold')
    ax.set_xticks(range(len(var_order)))
    ax.set_xticklabels(labels, fontsize=12)
    ax.set_ylabel('Cells (millions)', fontsize=14)
    ax.set_title('Cortical cells\nretained', fontsize=14, fontweight='bold')

    plt.tight_layout()
    outpath = os.path.join(OUT_DIR, 'margin_threshold_metrics_comparison.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    snrnaseq_results = compute_snrnaseq_concordance()
    prop_results = compute_proportion_similarity()
    depth_results = compute_depth_similarity()
    plot_combined_metrics(snrnaseq_results, prop_results, depth_results)
    print("\nDone!")
