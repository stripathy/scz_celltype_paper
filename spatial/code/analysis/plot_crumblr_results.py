#!/usr/bin/env python3
"""
Generate publication figures from crumblr results.

Reads:
  output/crumblr/crumblr_results_subclass.csv
  output/crumblr/crumblr_results_supertype.csv
  output/cross_platform_comparison/median_depth_per_{subclass,supertype}.csv

Generates:
  output/crumblr/volcano_subclass.png
  output/crumblr/volcano_supertype.png
  output/crumblr/effect_bar_subclass.png
  output/crumblr/effect_vs_depth_supertype.png
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    BASE_DIR, CRUMBLR_DIR,
    GABA_PREFIXES, GLUT_PREFIXES, NN_PREFIXES, infer_class,
)

RESULTS_DIR = CRUMBLR_DIR
DEPTH_DIR = os.path.join(BASE_DIR, "output", "cross_platform_comparison")

CLASS_COLORS = {'Glut': '#00ADF8', 'GABA': '#F05A28', 'NN': '#808080', 'Other': '#999999'}


def plot_volcano(results, level, out_dir):
    """Volcano plot: logFC vs -log10(p), colored by class."""
    df = results.copy()
    df['class'] = df['celltype'].apply(infer_class)
    df['nlog10p'] = -np.log10(df['P.Value'])

    fig, ax = plt.subplots(figsize=(8, 7))

    for cls in ['Glut', 'GABA', 'NN', 'Other']:
        mask = df['class'] == cls
        if mask.sum() == 0:
            continue
        sub = df[mask]
        ax.scatter(sub['logFC'], sub['nlog10p'],
                   c=CLASS_COLORS[cls], s=60, alpha=0.7, label=cls,
                   edgecolors='white', linewidth=0.5, zorder=5)

    # Label FDR-significant types
    for _, row in df.iterrows():
        if row['FDR'] < 0.05:
            ax.annotate(row['celltype'],
                       (row['logFC'], row['nlog10p']),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=9, fontweight='bold', alpha=0.9)
        elif row['FDR'] < 0.10:
            ax.annotate(row['celltype'],
                       (row['logFC'], row['nlog10p']),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8, alpha=0.7)

    ax.axhline(-np.log10(0.05), color='grey', linestyle='--', alpha=0.5, linewidth=1)
    ax.axvline(0, color='grey', alpha=0.3, linewidth=0.5)
    ax.text(ax.get_xlim()[1] * 0.95, -np.log10(0.05) + 0.05, 'p=0.05',
            ha='right', fontsize=9, color='grey', fontstyle='italic')

    ax.set_xlabel('logFC (SCZ vs Control)', fontsize=14)
    ax.set_ylabel('-log10(p-value)', fontsize=14)
    ax.set_title(f'Xenium SCZ: {level.capitalize()} Composition\n'
                 f'Cropped (cortical only), whole composition',
                 fontsize=14, fontweight='bold')
    ax.legend(fontsize=11, loc='upper left')

    n_fdr = (df['FDR'] < 0.05).sum()
    ax.text(0.98, 0.98, f'n={len(df)} types\n{n_fdr} FDR<0.05',
            transform=ax.transAxes, ha='right', va='top', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    fname = f'volcano_{level}.png'
    plt.savefig(os.path.join(out_dir, fname), dpi=100)
    plt.close()
    print(f"  Saved {fname}")


def plot_effect_bars(results, level, out_dir):
    """Horizontal bar plot of logFC per subclass, sorted by effect size."""
    df = results.copy()
    df['class'] = df['celltype'].apply(infer_class)
    df = df.sort_values('logFC')

    fig, ax = plt.subplots(figsize=(8, max(5, len(df) * 0.3)))

    colors = [CLASS_COLORS[c] for c in df['class']]
    bars = ax.barh(range(len(df)), df['logFC'], color=colors, alpha=0.8,
                   edgecolor='white', linewidth=0.5)

    # Mark significant types
    for i, (_, row) in enumerate(df.iterrows()):
        if row['FDR'] < 0.05:
            ax.text(row['logFC'] + (0.01 if row['logFC'] >= 0 else -0.01),
                    i, '***', ha='left' if row['logFC'] >= 0 else 'right',
                    va='center', fontsize=10, fontweight='bold')
        elif row['FDR'] < 0.10:
            ax.text(row['logFC'] + (0.01 if row['logFC'] >= 0 else -0.01),
                    i, '**', ha='left' if row['logFC'] >= 0 else 'right',
                    va='center', fontsize=10)
        elif row['P.Value'] < 0.05:
            ax.text(row['logFC'] + (0.01 if row['logFC'] >= 0 else -0.01),
                    i, '*', ha='left' if row['logFC'] >= 0 else 'right',
                    va='center', fontsize=10, color='grey')

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df['celltype'], fontsize=11)
    ax.axvline(0, color='black', linewidth=0.5)
    ax.set_xlabel('logFC (SCZ vs Control)', fontsize=14)
    ax.set_title(f'SCZ Effect Sizes: {level.capitalize()}\n'
                 f'*** FDR<0.05  ** FDR<0.10  * nom p<0.05',
                 fontsize=14, fontweight='bold')

    # Legend
    for cls, color in CLASS_COLORS.items():
        if cls in df['class'].values:
            ax.barh([], [], color=color, label=cls)
    ax.legend(fontsize=10, loc='lower right')

    plt.tight_layout()
    fname = f'effect_bar_{level}.png'
    plt.savefig(os.path.join(out_dir, fname), dpi=100)
    plt.close()
    print(f"  Saved {fname}")


def plot_depth_vs_effect(results, level, depth_file, out_dir):
    """Scatter: median cortical depth vs SCZ logFC."""
    if not os.path.exists(depth_file):
        print(f"  Skipping depth plot: {depth_file} not found")
        return

    depth = pd.read_csv(depth_file, index_col=0)
    df = results.copy()
    df['class'] = df['celltype'].apply(infer_class)

    merged = df.merge(depth[['median_depth', 'n_cells']],
                      left_on='celltype', right_index=True, how='inner')

    if len(merged) < 5:
        print(f"  Skipping depth plot: only {len(merged)} matched types")
        return

    fig, ax = plt.subplots(figsize=(8, 7))

    for cls in ['Glut', 'GABA', 'NN']:
        mask = merged['class'] == cls
        if mask.sum() == 0:
            continue
        sub = merged[mask]
        ax.scatter(sub['median_depth'], sub['logFC'],
                   c=CLASS_COLORS[cls],
                   s=sub['n_cells'].clip(upper=5000) / 50 + 20,
                   alpha=0.6, label=f'{cls} (n={mask.sum()})',
                   edgecolors='white', linewidth=0.5, zorder=5)

    # Label significant and large-effect types
    for _, row in merged.iterrows():
        if row['FDR'] < 0.05 or abs(row['logFC']) > 0.4 or row['P.Value'] < 0.01:
            ax.annotate(row['celltype'],
                       (row['median_depth'], row['logFC']),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=7, alpha=0.8)

    ax.axhline(0, color='grey', alpha=0.3, linewidth=1)

    # Depth region shading
    for lname, d0, d1 in [('L1', 0, 0.10), ('L2/3', 0.10, 0.40),
                           ('L4', 0.40, 0.55), ('L5', 0.55, 0.70),
                           ('L6', 0.70, 0.90)]:
        ax.axvspan(d0, d1, alpha=0.04, color='grey')
        ax.text((d0 + d1) / 2, ax.get_ylim()[0] if ax.get_ylim()[0] != 0 else -0.5,
                lname, ha='center', va='bottom', fontsize=9, alpha=0.4, fontstyle='italic')

    # Correlation
    r_all, p_all = pearsonr(merged['median_depth'], merged['logFC'])
    neur = merged[merged['class'].isin(['Glut', 'GABA'])]
    if len(neur) > 3:
        r_neur, p_neur = pearsonr(neur['median_depth'], neur['logFC'])
    else:
        r_neur, p_neur = np.nan, np.nan

    # Regression line
    try:
        z = np.polyfit(merged['median_depth'], merged['logFC'], 1)
        x_line = np.linspace(merged['median_depth'].min(), merged['median_depth'].max(), 100)
        ax.plot(x_line, np.polyval(z, x_line), 'k-', alpha=0.3, linewidth=1.5)
    except np.linalg.LinAlgError:
        pass  # skip trendline if SVD fails

    ax.set_xlim(-0.02, 0.90)
    ax.set_xlabel('Median normalized depth from pia\n(0 = pia, 1 = WM)', fontsize=12)
    ax.set_ylabel('logFC (SCZ vs Control)', fontsize=12)
    ax.set_title(f'SCZ Effect vs Cortical Depth ({level.capitalize()})\n'
                 f'All: r={r_all:.3f} (p={p_all:.3f}) | '
                 f'Neuronal: r={r_neur:.3f} (p={p_neur:.3f})',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=9, loc='upper left')

    plt.tight_layout()
    fname = f'effect_vs_depth_{level}.png'
    plt.savefig(os.path.join(out_dir, fname), dpi=100)
    plt.close()
    print(f"  Saved {fname}")


def load_nicole_stratified(neuronal_path, nonneuronal_path):
    """Load Nicole's stratified snRNAseq meta-analysis results (neuronal + non-neuronal).

    Returns combined DataFrame with columns: celltype, beta_snrnaseq, se, pval_snrnaseq, padj_snrnaseq, analysis_stratum.
    Excludes _1-SEAAD suffixed types (alternative taxonomy definitions that don't match Xenium types).
    """
    dfs = []
    for path, stratum in [(neuronal_path, 'neuronal'), (nonneuronal_path, 'non-neuronal')]:
        if not os.path.exists(path):
            print(f"  WARNING: {path} not found, skipping {stratum}")
            continue
        df = pd.read_csv(path)
        df = df[~df['CellType'].str.contains('SEAAD', na=False)]  # exclude alt taxonomy types
        df['analysis_stratum'] = stratum
        dfs.append(df)

    if not dfs:
        return None

    combined = pd.concat(dfs, ignore_index=True)
    combined = combined.rename(columns={'CellType': 'celltype', 'estimate': 'beta_snrnaseq',
                                         'pval': 'pval_snrnaseq', 'padj': 'padj_snrnaseq'})
    print(f"  Loaded Nicole stratified data: {len(combined)} types "
          f"({sum(combined['analysis_stratum']=='neuronal')} neuronal, "
          f"{sum(combined['analysis_stratum']=='non-neuronal')} non-neuronal)")
    return combined


def plot_snrnaseq_comparison(results, nicole_neuronal_path, nicole_nonneuronal_path, out_dir,
                             level='supertype'):
    """Scatter: snRNAseq meta-analysis betas vs Xenium logFC at cluster level.

    Uses stratified Xenium results (neuronal and non-neuronal analyzed
    separately) to match Nicole's stratified snRNAseq meta-analysis.
    Falls back to whole-composition Xenium results if stratified files
    are not available.
    """
    nicole = load_nicole_stratified(nicole_neuronal_path, nicole_nonneuronal_path)
    if nicole is None:
        print(f"  Skipping snRNAseq comparison ({level}): no Nicole data found")
        return

    # Try to load stratified Xenium results
    xen_neuronal_path = os.path.join(out_dir, f'crumblr_results_{level}_neuronal.csv')
    xen_nonneuronal_path = os.path.join(out_dir, f'crumblr_results_{level}_nonneuronal.csv')

    if os.path.exists(xen_neuronal_path) and os.path.exists(xen_nonneuronal_path):
        print("  Using stratified Xenium crumblr results (neuronal + non-neuronal)")
        xen_n = pd.read_csv(xen_neuronal_path)
        xen_nn = pd.read_csv(xen_nonneuronal_path)
        xen = pd.concat([xen_n, xen_nn], ignore_index=True)
        stratified = True
    else:
        print("  Using whole-composition Xenium crumblr results (stratified not found)")
        xen = results
        stratified = False

    xen = xen.rename(columns={'logFC': 'logFC_xenium', 'P.Value': 'pval_xenium',
                                   'FDR': 'FDR_xenium'})
    se_col_xen = 'SE' if 'SE' in xen.columns else None

    merged = pd.merge(
        nicole[['celltype', 'beta_snrnaseq', 'se', 'pval_snrnaseq', 'padj_snrnaseq']],
        xen[['celltype', 'logFC_xenium', 'pval_xenium', 'FDR_xenium'] +
            ([se_col_xen] if se_col_xen else [])],
        on='celltype', how='inner'
    )
    if se_col_xen:
        merged = merged.rename(columns={se_col_xen: 'SE_xenium'})

    merged['class'] = merged['celltype'].apply(infer_class)
    print(f"  snRNAseq comparison: {len(merged)} shared cluster types")

    # Correlations
    r_all, p_all = pearsonr(merged['beta_snrnaseq'], merged['logFC_xenium'])
    neur = merged[merged['class'].isin(['Glut', 'GABA'])]
    r_neur, p_neur = pearsonr(neur['beta_snrnaseq'], neur['logFC_xenium']) if len(neur) > 3 else (np.nan, np.nan)
    nn = merged[merged['class'] == 'NN']
    r_nn, p_nn = pearsonr(nn['beta_snrnaseq'], nn['logFC_xenium']) if len(nn) > 3 else (np.nan, np.nan)

    fig, ax = plt.subplots(figsize=(9, 8))

    for cls in ['Glut', 'GABA', 'NN', 'Other']:
        mask = merged['class'] == cls
        if mask.sum() == 0:
            continue
        sub = merged[mask]
        ax.scatter(sub['beta_snrnaseq'], sub['logFC_xenium'],
                   c=CLASS_COLORS[cls], s=60, alpha=0.7, label=f'{cls} (n={mask.sum()})',
                   edgecolors='white', linewidth=0.5, zorder=5)

    # Horizontal error bars for snRNAseq SE
    for _, row in merged.iterrows():
        ax.plot([row['beta_snrnaseq'] - row['se'], row['beta_snrnaseq'] + row['se']],
                [row['logFC_xenium'], row['logFC_xenium']],
                color=CLASS_COLORS[row['class']], alpha=0.15, linewidth=1, zorder=3)

    # Vertical error bars for Xenium SE (if available)
    if 'SE_xenium' in merged.columns:
        for _, row in merged.iterrows():
            ax.plot([row['beta_snrnaseq'], row['beta_snrnaseq']],
                    [row['logFC_xenium'] - row['SE_xenium'], row['logFC_xenium'] + row['SE_xenium']],
                    color=CLASS_COLORS[row['class']], alpha=0.15, linewidth=1, zorder=3)

    # Label notable types
    for _, row in merged.iterrows():
        is_sig_either = (row['padj_snrnaseq'] < 0.05) or (row['FDR_xenium'] < 0.05)
        is_large = abs(row['beta_snrnaseq']) > 0.2 or abs(row['logFC_xenium']) > 0.5
        is_nom_both = (row['pval_snrnaseq'] < 0.05) and (row['pval_xenium'] < 0.05)
        if is_sig_either or is_large or is_nom_both:
            ax.annotate(row['celltype'],
                       (row['beta_snrnaseq'], row['logFC_xenium']),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=7, alpha=0.8)

    ax.axhline(0, color='grey', alpha=0.2, linewidth=0.5)
    ax.axvline(0, color='grey', alpha=0.2, linewidth=0.5)

    # Regression line
    z = np.polyfit(merged['beta_snrnaseq'], merged['logFC_xenium'], 1)
    lim_x = max(abs(merged['beta_snrnaseq']).max(), 0.3) * 1.3
    x_line = np.linspace(-lim_x, lim_x, 100)
    ax.plot(x_line, np.polyval(z, x_line), 'k-', alpha=0.3, linewidth=1.5)

    ylabel = ('Xenium spatial logFC (SCZ vs Control)\n'
              '(stratified: neuronal and non-neuronal analyzed separately)'
              if stratified else
              'Xenium spatial logFC (SCZ vs Control)\n(whole composition)')
    ax.set_xlabel('snRNAseq meta-analysis beta (SCZ effect)', fontsize=13)
    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_title(f'SCZ Cluster Effects: snRNAseq vs Xenium (stratified)\n'
                 f'All: r = {r_all:.3f} (p = {p_all:.1e}) | '
                 f'Neuronal: r = {r_neur:.3f} | '
                 f'Non-neuronal: r = {r_nn:.3f} | n = {len(merged)}',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=10, loc='lower right')

    plt.tight_layout()
    fname = f'snrnaseq_vs_xenium_{level}.png'
    plt.savefig(os.path.join(out_dir, fname), dpi=100)
    plt.close()
    print(f"  Saved {fname}")

    # Save comparison table
    csv_name = f'snrnaseq_vs_xenium_comparison{"_subclass" if level == "subclass" else ""}.csv'
    merged.to_csv(os.path.join(out_dir, csv_name), index=False)
    print(f"  Saved {csv_name}")

    # Print notable concordances
    sig_either = merged[(merged['padj_snrnaseq'] < 0.05) | (merged['FDR_xenium'] < 0.05)]
    if len(sig_either) > 0:
        print(f"  FDR-significant in either platform:")
        for _, row in sig_either.sort_values('celltype').iterrows():
            concordant = 'Y' if (row['beta_snrnaseq'] * row['logFC_xenium']) > 0 else 'N'
            flags = ''
            if row['padj_snrnaseq'] < 0.05: flags += ' FDR_snRNA'
            if row['FDR_xenium'] < 0.05: flags += ' FDR_Xen'
            print(f"    {concordant} {row['celltype']:25s} "
                  f"snRNA={row['beta_snrnaseq']:+.3f}  Xen={row['logFC_xenium']:+.3f}{flags}")


NICOLE_NEURONAL_PATH = os.path.join(BASE_DIR, "data", "nicole_scz_snrnaseq_betas",
                                     "final_results_crumblr_7_cohorts.csv")
NICOLE_NONNEURONAL_PATH = os.path.join(BASE_DIR, "data", "nicole_scz_snrnaseq_betas",
                                        "final_results_crumblr_7_nonN_cohorts.csv")
NICOLE_NEURONAL_SUBCLASS_PATH = os.path.join(BASE_DIR, "data", "nicole_scz_snrnaseq_betas",
                                              "final_results_crumblr_7_cohorts_subclass.csv")
NICOLE_NONNEURONAL_SUBCLASS_PATH = os.path.join(BASE_DIR, "data", "nicole_scz_snrnaseq_betas",
                                                  "final_results_crumblr_7_nonN_cohorts_subclass.csv")


def main():
    print("=" * 60)
    print("CRUMBLR RESULT FIGURES")
    print("=" * 60)

    for level in ['subclass', 'supertype']:
        fpath = os.path.join(RESULTS_DIR, f'crumblr_results_{level}.csv')
        if not os.path.exists(fpath):
            print(f"\n  Skipping {level}: {fpath} not found")
            continue

        results = pd.read_csv(fpath)
        print(f"\n  {level}: {len(results)} types")

        # Volcano plot
        plot_volcano(results, level, RESULTS_DIR)

        # Effect bar plot (subclass only — too many types for supertype)
        if level == 'subclass':
            plot_effect_bars(results, level, RESULTS_DIR)

        # Depth-vs-effect scatter
        depth_file = os.path.join(DEPTH_DIR, f'median_depth_per_{level}.csv')
        plot_depth_vs_effect(results, level, depth_file, RESULTS_DIR)

    # snRNAseq meta-analysis comparison (supertype level)
    supertype_path = os.path.join(RESULTS_DIR, 'crumblr_results_supertype.csv')
    if os.path.exists(supertype_path):
        print(f"\n  --- snRNAseq meta-analysis comparison (supertype) ---")
        supertype_results = pd.read_csv(supertype_path)
        plot_snrnaseq_comparison(supertype_results, NICOLE_NEURONAL_PATH,
                                NICOLE_NONNEURONAL_PATH, RESULTS_DIR)

    # snRNAseq meta-analysis comparison (subclass level)
    subclass_path = os.path.join(RESULTS_DIR, 'crumblr_results_subclass.csv')
    if os.path.exists(subclass_path) and os.path.exists(NICOLE_NEURONAL_SUBCLASS_PATH):
        print(f"\n  --- snRNAseq meta-analysis comparison (subclass) ---")
        subclass_results = pd.read_csv(subclass_path)
        plot_snrnaseq_comparison(subclass_results, NICOLE_NEURONAL_SUBCLASS_PATH,
                                NICOLE_NONNEURONAL_SUBCLASS_PATH, RESULTS_DIR,
                                level='subclass')

    print("\nDone!")


if __name__ == '__main__':
    main()
