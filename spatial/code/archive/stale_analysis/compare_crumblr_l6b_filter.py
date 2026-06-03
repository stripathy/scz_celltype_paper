#!/usr/bin/env python3
"""
Compare crumblr results: baseline vs L6b margin-filtered.
Also compare both to snRNAseq reference effect sizes.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import BASE_DIR, SUBCLASS_TO_CLASS

CRUMBLR_DIR = os.path.join(BASE_DIR, "output", "crumblr")
OUT_DIR = os.path.join(BASE_DIR, "output", "l6b_diagnostics")
SNRNASEQ_PATH = os.path.join(BASE_DIR, "data", "nicole_scz_snrnaseq_betas", "scz_coefs.xlsx")

plt.rcParams.update({
    'font.size': 14, 'axes.titlesize': 18, 'axes.labelsize': 16,
    'xtick.labelsize': 12, 'ytick.labelsize': 12,
    'figure.facecolor': 'white',
})


def load_results():
    """Load baseline and L6b-filtered crumblr results."""
    # Try the double-suffix filenames that run_crumblr.R produces
    base_sub = os.path.join(CRUMBLR_DIR, "crumblr_results_subclass_baseline_baseline.csv")
    filt_sub = os.path.join(CRUMBLR_DIR, "crumblr_results_subclass_l6b_margin_l6b_margin.csv")
    base_sup = os.path.join(CRUMBLR_DIR, "crumblr_results_supertype_baseline_baseline.csv")
    filt_sup = os.path.join(CRUMBLR_DIR, "crumblr_results_supertype_l6b_margin_l6b_margin.csv")

    results = {
        'subclass_baseline': pd.read_csv(base_sub),
        'subclass_filtered': pd.read_csv(filt_sub),
        'supertype_baseline': pd.read_csv(base_sup),
        'supertype_filtered': pd.read_csv(filt_sup),
    }

    for k, v in results.items():
        print(f"  {k}: {len(v)} types, {(v['FDR'] < 0.05).sum()} FDR<0.05, "
              f"{(v['FDR'] < 0.10).sum()} FDR<0.10")

    return results


def load_snrnaseq():
    """Load Nicole's snRNAseq SCZ effect sizes."""
    df = pd.read_excel(SNRNASEQ_PATH)
    print(f"  snRNAseq: {len(df)} cell types")
    # Columns: CellType, estimate, se, zval, pval, ci.lb, ci.ub, padj
    df = df.rename(columns={'CellType': 'celltype', 'estimate': 'beta_snrnaseq'})
    return df


def main():
    print("Loading results...")
    results = load_results()
    snrna = load_snrnaseq()

    base_sub = results['subclass_baseline']
    filt_sub = results['subclass_filtered']
    base_sup = results['supertype_baseline']
    filt_sup = results['supertype_filtered']

    # ══════════════════════════════════════════════════════════════
    # COMPARISON 1: Effect size changes at subclass level
    # ══════════════════════════════════════════════════════════════
    print("\n" + "=" * 70)
    print("SUBCLASS LEVEL: Baseline vs L6b-Filtered")
    print("=" * 70)

    merged_sub = base_sub[['celltype', 'logFC', 'P.Value', 'FDR', 'SE']].merge(
        filt_sub[['celltype', 'logFC', 'P.Value', 'FDR', 'SE']],
        on='celltype', suffixes=('_base', '_filt')
    )

    print(f"\n{'celltype':25s} {'logFC_base':>10s} {'logFC_filt':>10s} {'delta':>8s} "
          f"{'FDR_base':>10s} {'FDR_filt':>10s}")
    print("-" * 85)
    for _, row in merged_sub.sort_values('P.Value_base').iterrows():
        delta = row['logFC_filt'] - row['logFC_base']
        print(f"{row['celltype']:25s} {row['logFC_base']:>+10.4f} {row['logFC_filt']:>+10.4f} "
              f"{delta:>+8.4f} {row['FDR_base']:>10.4f} {row['FDR_filt']:>10.4f}")

    r_sub, p_sub = stats.pearsonr(merged_sub['logFC_base'], merged_sub['logFC_filt'])
    print(f"\nlogFC correlation (baseline vs filtered): r={r_sub:.4f}, p={p_sub:.2e}")

    # ══════════════════════════════════════════════════════════════
    # COMPARISON 2: Effect size changes at supertype level
    # ══════════════════════════════════════════════════════════════
    print("\n" + "=" * 70)
    print("SUPERTYPE LEVEL: Baseline vs L6b-Filtered")
    print("=" * 70)

    merged_sup = base_sup[['celltype', 'logFC', 'P.Value', 'FDR']].merge(
        filt_sup[['celltype', 'logFC', 'P.Value', 'FDR']],
        on='celltype', suffixes=('_base', '_filt')
    )

    r_sup, p_sup = stats.pearsonr(merged_sup['logFC_base'], merged_sup['logFC_filt'])
    print(f"logFC correlation (baseline vs filtered): r={r_sup:.4f}, p={p_sup:.2e}")
    print(f"FDR<0.05 baseline: {(merged_sup['FDR_base'] < 0.05).sum()}")
    print(f"FDR<0.05 filtered: {(merged_sup['FDR_filt'] < 0.05).sum()}")

    # Show L6b supertypes specifically
    l6b_types = merged_sup[merged_sup['celltype'].str.startswith('L6b')]
    if len(l6b_types) > 0:
        print(f"\nL6b supertypes:")
        for _, row in l6b_types.iterrows():
            print(f"  {row['celltype']:20s} logFC: {row['logFC_base']:+.4f} -> {row['logFC_filt']:+.4f}  "
                  f"FDR: {row['FDR_base']:.4f} -> {row['FDR_filt']:.4f}")

    # Show FDR<0.10 hits that changed
    sig = merged_sup[(merged_sup['FDR_base'] < 0.10) | (merged_sup['FDR_filt'] < 0.10)]
    print(f"\nAll FDR<0.10 hits (in either):")
    for _, row in sig.sort_values('FDR_base').iterrows():
        delta = row['logFC_filt'] - row['logFC_base']
        gained = "NEW" if row['FDR_base'] >= 0.10 and row['FDR_filt'] < 0.10 else ""
        lost = "LOST" if row['FDR_base'] < 0.10 and row['FDR_filt'] >= 0.10 else ""
        flag = gained or lost or ""
        print(f"  {row['celltype']:25s} logFC: {row['logFC_base']:+.4f} -> {row['logFC_filt']:+.4f} "
              f"({delta:+.4f})  FDR: {row['FDR_base']:.4f} -> {row['FDR_filt']:.4f} {flag}")

    # ══════════════════════════════════════════════════════════════
    # COMPARISON 3: Correlation with snRNAseq
    # ══════════════════════════════════════════════════════════════
    print("\n" + "=" * 70)
    print("CORRELATION WITH snRNAseq EFFECT SIZES")
    print("=" * 70)

    for level, base_df, filt_df in [
        ('subclass', base_sub, filt_sub),
        ('supertype', base_sup, filt_sup),
    ]:
        # Merge with snRNAseq
        m_base = base_df.merge(snrna[['celltype', 'beta_snrnaseq']], on='celltype', how='inner')
        m_filt = filt_df.merge(snrna[['celltype', 'beta_snrnaseq']], on='celltype', how='inner')

        if len(m_base) > 2:
            r_base, p_base = stats.pearsonr(m_base['logFC'], m_base['beta_snrnaseq'])
            r_filt, p_filt = stats.pearsonr(m_filt['logFC'], m_filt['beta_snrnaseq'])
            print(f"\n  {level} (n={len(m_base)} shared types):")
            print(f"    Baseline vs snRNAseq: r={r_base:.4f}, p={p_base:.4e}")
            print(f"    Filtered vs snRNAseq: r={r_filt:.4f}, p={p_filt:.4e}")
            print(f"    Delta r: {r_filt - r_base:+.4f}")

    # ══════════════════════════════════════════════════════════════
    # PLOTS
    # ══════════════════════════════════════════════════════════════

    fig, axes = plt.subplots(2, 3, figsize=(20, 14))

    # Row 1: Baseline vs Filtered effect sizes
    # Panel A: Subclass logFC scatter
    ax = axes[0, 0]
    for _, row in merged_sub.iterrows():
        cls = SUBCLASS_TO_CLASS.get(row['celltype'], 'Unknown')
        color = '#1f77b4' if cls == 'Glutamatergic' else '#d62728' if cls == 'GABAergic' else '#7f7f7f'
        marker = '*' if row['celltype'] == 'L6b' else 'o'
        size = 200 if row['celltype'] == 'L6b' else 50
        ax.scatter(row['logFC_base'], row['logFC_filt'], color=color, s=size,
                   marker=marker, zorder=10 if marker == '*' else 5, edgecolor='white')
        if abs(row['logFC_base']) > 0.15 or row['celltype'] == 'L6b':
            ax.annotate(row['celltype'], (row['logFC_base'], row['logFC_filt']),
                        fontsize=9, ha='left')
    lim = max(abs(merged_sub['logFC_base']).max(), abs(merged_sub['logFC_filt']).max()) * 1.2
    ax.plot([-lim, lim], [-lim, lim], 'k--', alpha=0.3)
    ax.axhline(0, color='gray', lw=0.5)
    ax.axvline(0, color='gray', lw=0.5)
    ax.set_xlabel('logFC (baseline)')
    ax.set_ylabel('logFC (L6b filtered)')
    ax.set_title(f'Subclass logFC\nr={r_sub:.4f}')

    # Panel B: Supertype logFC scatter
    ax = axes[0, 1]
    for _, row in merged_sup.iterrows():
        is_l6b = row['celltype'].startswith('L6b')
        color = '#d62728' if is_l6b else '#1f77b4'
        size = 80 if is_l6b else 15
        alpha = 1.0 if is_l6b else 0.3
        ax.scatter(row['logFC_base'], row['logFC_filt'], color=color, s=size,
                   alpha=alpha, edgecolor='white' if is_l6b else 'none')
        if is_l6b:
            ax.annotate(row['celltype'], (row['logFC_base'], row['logFC_filt']),
                        fontsize=8, ha='left', color='#d62728')
    lim = max(abs(merged_sup['logFC_base']).max(), abs(merged_sup['logFC_filt']).max()) * 1.1
    ax.plot([-lim, lim], [-lim, lim], 'k--', alpha=0.3)
    ax.axhline(0, color='gray', lw=0.5)
    ax.axvline(0, color='gray', lw=0.5)
    ax.set_xlabel('logFC (baseline)')
    ax.set_ylabel('logFC (L6b filtered)')
    ax.set_title(f'Supertype logFC\nr={r_sup:.4f} (L6b types in red)')

    # Panel C: Delta logFC distribution
    ax = axes[0, 2]
    delta_sup = merged_sup['logFC_filt'] - merged_sup['logFC_base']
    ax.hist(delta_sup, bins=40, color='#1f77b4', alpha=0.7, edgecolor='white')
    # Mark L6b types
    l6b_deltas = delta_sup[merged_sup['celltype'].str.startswith('L6b')]
    for d in l6b_deltas:
        ax.axvline(d, color='#d62728', lw=2, alpha=0.7)
    ax.axvline(0, color='black', lw=1, ls='--')
    ax.set_xlabel('Δ logFC (filtered - baseline)')
    ax.set_ylabel('Count')
    ax.set_title(f'Effect Size Changes\n(median Δ={delta_sup.median():+.4f})')

    # Row 2: Correlation with snRNAseq
    for idx, (level, base_df, filt_df) in enumerate([
        ('subclass', base_sub, filt_sub),
        ('supertype', base_sup, filt_sup),
    ]):
        m_base = base_df.merge(snrna[['celltype', 'beta_snrnaseq']], on='celltype', how='inner')
        m_filt = filt_df.merge(snrna[['celltype', 'beta_snrnaseq']], on='celltype', how='inner')

        if len(m_base) < 3:
            continue

        r_b, p_b = stats.pearsonr(m_base['logFC'], m_base['beta_snrnaseq'])
        r_f, p_f = stats.pearsonr(m_filt['logFC'], m_filt['beta_snrnaseq'])

        # Baseline vs snRNAseq
        ax = axes[1, idx]
        for _, row in m_base.iterrows():
            is_l6b = row['celltype'].startswith('L6b') or row['celltype'] == 'L6b'
            color = '#d62728' if is_l6b else '#1f77b4'
            size = 80 if is_l6b else 30
            ax.scatter(row['beta_snrnaseq'], row['logFC'], color=color, s=size,
                       alpha=0.6, edgecolor='white' if is_l6b else 'none',
                       label='_' if not is_l6b else None)
            if is_l6b:
                ax.annotate(row['celltype'], (row['beta_snrnaseq'], row['logFC']),
                            fontsize=8, color='#d62728')

        # Overlay filtered
        for _, row in m_filt.iterrows():
            is_l6b = row['celltype'].startswith('L6b') or row['celltype'] == 'L6b'
            if is_l6b:
                ax.scatter(row['beta_snrnaseq'], row['logFC'], color='#ff7f0e', s=80,
                           marker='D', edgecolor='white', zorder=10)
                ax.annotate(f"  (filt)", (row['beta_snrnaseq'], row['logFC']),
                            fontsize=7, color='#ff7f0e')

        lim = max(abs(m_base['beta_snrnaseq']).max(), abs(m_base['logFC']).max()) * 1.2
        ax.plot([-lim, lim], [-lim, lim], 'k--', alpha=0.2)
        ax.axhline(0, color='gray', lw=0.5)
        ax.axvline(0, color='gray', lw=0.5)
        ax.set_xlabel('snRNAseq beta')
        ax.set_ylabel('Xenium logFC')
        ax.set_title(f'{level.title()} vs snRNAseq\n'
                     f'base r={r_b:.3f} → filt r={r_f:.3f}')

    # Panel F: FDR comparison
    ax = axes[1, 2]
    sig_all = merged_sup[(merged_sup['FDR_base'] < 0.15) | (merged_sup['FDR_filt'] < 0.15)]
    if len(sig_all) > 0:
        sig_all = sig_all.sort_values('FDR_base')
        y = np.arange(len(sig_all))
        ax.barh(y - 0.15, -np.log10(sig_all['FDR_base']), 0.3, color='#1f77b4',
                alpha=0.7, label='Baseline')
        ax.barh(y + 0.15, -np.log10(sig_all['FDR_filt']), 0.3, color='#ff7f0e',
                alpha=0.7, label='L6b filtered')
        ax.axvline(-np.log10(0.05), color='red', ls='--', lw=1.5, alpha=0.7, label='FDR=0.05')
        ax.axvline(-np.log10(0.10), color='orange', ls='--', lw=1.5, alpha=0.7, label='FDR=0.10')
        ax.set_yticks(y)
        ax.set_yticklabels(sig_all['celltype'], fontsize=9)
        ax.set_xlabel('-log10(FDR)')
        ax.set_title('Top Hits: FDR Comparison')
        ax.legend(fontsize=10)

    fig.suptitle('Crumblr Results: Baseline vs L6b Margin Filter', fontsize=22, y=1.02)
    fig.tight_layout()
    path = os.path.join(OUT_DIR, 'crumblr_baseline_vs_l6b_filtered.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {path}")


if __name__ == '__main__':
    main()
