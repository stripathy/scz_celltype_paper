#!/usr/bin/env python3
"""
Plot crumblr results for SEA-AD pseudoprogression analysis.

Generates figures for both "all cell types" and "neurons only" analyses:
  1. Volcano plots per source/level
  2. Effect size comparison: snRNAseq vs MERFISH (orig + reclassified)
  3. Summary bar chart of significant hits per source
  4. Forest plot of top shared effects
  5. Direction-of-effect agreement

Usage:
    python3 -u plot_seaad_crumblr_results.py
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['font.size'] = 14

BASE_DIR = os.path.expanduser("~/Github/SCZ_Xenium")
IO_DIR = os.path.join(BASE_DIR, "output", "crumblr_seaad")

# Nice labels for sources
SOURCE_LABELS = {
    'snrnaseq': 'snRNAseq (n=84)',
    'merfish_orig': 'MERFISH original (n=27)',
    'merfish_recl_qc': 'MERFISH reclassified (n=27)',
    'merfish_depth': 'MERFISH depth-annot (n=27)',
}

# Colors per source
SOURCE_COLORS = {
    'snrnaseq': '#2166ac',
    'merfish_orig': '#b2182b',
    'merfish_recl_qc': '#d6604d',
    'merfish_depth': '#4daf4a',
}


def get_base_source(source):
    """Strip level suffix to get base source name."""
    # Order matters: check longer prefixes first to avoid partial matches
    for base in ['merfish_recl_qc', 'merfish_depth', 'merfish_orig', 'snrnaseq']:
        if source.startswith(base):
            return base
    return source


def load_all_results():
    """Load the combined crumblr results."""
    fpath = os.path.join(IO_DIR, "crumblr_results_all.csv")
    df = pd.read_csv(fpath)
    print(f"Loaded {len(df)} rows from crumblr_results_all.csv")
    print(f"  Sources: {df['source'].unique()}")
    print(f"  Levels: {df['level'].unique()}")
    return df


def plot_volcanos(df, level, suffix, title_extra=''):
    """Volcano plots: one per source for a given level."""
    sub = df[df['level'] == level].copy()
    sources = sorted(sub['source'].unique())

    fig, axes = plt.subplots(1, len(sources), figsize=(7 * len(sources), 6),
                              squeeze=False)

    for i, source in enumerate(sources):
        ax = axes[0, i]
        d = sub[sub['source'] == source].copy()
        d['-log10p'] = -np.log10(d['P.Value'])

        colors = []
        for _, row in d.iterrows():
            if row['FDR'] < 0.05:
                colors.append('red')
            elif row['FDR'] < 0.10:
                colors.append('orange')
            elif row['P.Value'] < 0.05:
                colors.append('#4393c3')
            else:
                colors.append('grey')

        s_size = 50 if 'subclass' in level else 30
        ax.scatter(d['logFC'], d['-log10p'], c=colors, s=s_size, alpha=0.8,
                   edgecolors='k', linewidth=0.5)

        # Label significant or top hits
        for _, row in d.iterrows():
            if row['FDR'] < 0.10 or row['P.Value'] < 0.01:
                fs = 9 if 'subclass' in level else 7
                ax.annotate(row['celltype'],
                           (row['logFC'], row['-log10p']),
                           fontsize=fs, ha='center', va='bottom',
                           xytext=(0, 5), textcoords='offset points')

        if d['P.Value'].min() < 0.05:
            ax.axhline(-np.log10(0.05), color='grey', linestyle='--',
                       alpha=0.5, linewidth=0.8)
        ax.axvline(0, color='grey', linestyle='-', alpha=0.3)

        base_src = get_base_source(source)
        ax.set_xlabel('logFC (CPS effect)', fontsize=14)
        ax.set_ylabel('-log10(p-value)', fontsize=14)
        ax.set_title(SOURCE_LABELS.get(base_src, source), fontsize=16,
                     fontweight='bold')

    level_label = level.replace('_neurons', ' — Neurons Only').capitalize()
    plt.suptitle(f'Cell Type Proportion Changes with CPS ({level_label}){title_extra}',
                 fontsize=18, fontweight='bold', y=1.02)
    plt.tight_layout()
    out = os.path.join(IO_DIR, f"volcano_{suffix}.png")
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out}")


def plot_effect_comparison(df, level, suffix, title_extra=''):
    """Scatter: snRNAseq effect vs MERFISH effect for a given level."""
    sub = df[df['level'] == level].copy()

    # Find the snRNAseq source for this level
    snr_sources = [s for s in sub['source'].unique() if 'snrnaseq' in s]
    if not snr_sources:
        print(f"  Skipping effect comparison for {level}: no snRNAseq source found")
        return
    snr_source = snr_sources[0]

    merfish_pairs = []
    for s in sorted(sub['source'].unique()):
        base = get_base_source(s)
        if base == 'merfish_depth':
            merfish_pairs.append((s, 'MERFISH depth-annotated'))
        elif base == 'merfish_orig':
            merfish_pairs.append((s, 'MERFISH original'))
        elif base == 'merfish_recl_qc':
            merfish_pairs.append((s, 'MERFISH reclassified'))

    if not merfish_pairs:
        return

    snr = sub[sub['source'] == snr_source][['celltype', 'logFC', 'FDR']].copy()
    snr = snr.rename(columns={'logFC': 'logFC_snrnaseq', 'FDR': 'FDR_snrnaseq'})

    fig, axes = plt.subplots(1, len(merfish_pairs), figsize=(7 * len(merfish_pairs), 6))
    if len(merfish_pairs) == 1:
        axes = [axes]

    for i, (ms, ml) in enumerate(merfish_pairs):
        ax = axes[i]
        mer = sub[sub['source'] == ms][['celltype', 'logFC', 'FDR']].copy()
        mer = mer.rename(columns={'logFC': 'logFC_merfish', 'FDR': 'FDR_merfish'})

        merged = snr.merge(mer, on='celltype', how='inner')

        # Color by joint significance
        colors = []
        for _, row in merged.iterrows():
            snr_sig = row['FDR_snrnaseq'] < 0.10
            mer_sig = row['FDR_merfish'] < 0.10
            if snr_sig and mer_sig:
                colors.append('red')
            elif snr_sig:
                colors.append('#2166ac')
            elif mer_sig:
                colors.append('#b2182b')
            else:
                colors.append('grey')

        s_size = 80 if 'subclass' in level else 20
        ax.scatter(merged['logFC_snrnaseq'], merged['logFC_merfish'],
                  c=colors, s=s_size, alpha=0.8, edgecolors='k', linewidth=0.5)

        # Label points
        if 'subclass' in level:
            for _, row in merged.iterrows():
                ax.annotate(row['celltype'],
                           (row['logFC_snrnaseq'], row['logFC_merfish']),
                           fontsize=8, ha='center', va='bottom',
                           xytext=(0, 5), textcoords='offset points')
        else:
            top_snr = merged.nsmallest(10, 'FDR_snrnaseq')
            for _, row in top_snr.iterrows():
                ax.annotate(row['celltype'],
                           (row['logFC_snrnaseq'], row['logFC_merfish']),
                           fontsize=6, ha='center', va='bottom',
                           xytext=(0, 3), textcoords='offset points')

        r = merged['logFC_snrnaseq'].corr(merged['logFC_merfish'])
        n_shared = len(merged)
        ax.set_title(f'{ml}\nr = {r:.3f} (n={n_shared} types)',
                     fontsize=16, fontweight='bold')

        lim = max(abs(merged['logFC_snrnaseq']).max(),
                  abs(merged['logFC_merfish']).max()) * 1.2
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.plot([-lim, lim], [-lim, lim], 'k--', alpha=0.3, linewidth=0.8)
        ax.axhline(0, color='grey', linestyle='-', alpha=0.2)
        ax.axvline(0, color='grey', linestyle='-', alpha=0.2)

        ax.set_xlabel('logFC (snRNAseq)', fontsize=14)
        ax.set_ylabel(f'logFC ({ml})', fontsize=14)

    level_label = level.replace('_neurons', ' — Neurons Only').capitalize()
    plt.suptitle(f'CPS Effect Size: snRNAseq vs MERFISH ({level_label}){title_extra}',
                 fontsize=18, fontweight='bold', y=1.02)
    plt.tight_layout()
    out = os.path.join(IO_DIR, f"effect_comparison_{suffix}.png")
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out}")


def plot_significance_summary(df, levels, suffix, title_extra=''):
    """Bar chart: number of significant hits per source/level."""
    summary = []
    for (source, level), g in df.groupby(['source', 'level']):
        if level not in levels:
            continue
        n_total = len(g)
        n_fdr05 = (g['FDR'] < 0.05).sum()
        n_fdr10 = (g['FDR'] < 0.10).sum()
        n_nom05 = (g['P.Value'] < 0.05).sum()
        base_src = get_base_source(source)
        summary.append({
            'source': source, 'base_source': base_src, 'level': level,
            'n_total': n_total,
            'FDR < 0.05': n_fdr05,
            'FDR < 0.10': n_fdr10,
            'nom p < 0.05': n_nom05,
        })
    summary = pd.DataFrame(summary)

    n_panels = len(levels)
    fig, axes = plt.subplots(1, n_panels, figsize=(7 * n_panels, 5))
    if n_panels == 1:
        axes = [axes]

    for i, level in enumerate(levels):
        ax = axes[i]
        s = summary[summary['level'] == level].copy()
        s = s.sort_values('source')

        x = np.arange(len(s))
        w = 0.25
        bars1 = ax.bar(x - w, s['nom p < 0.05'], w, label='nom p < 0.05',
                       color='#92c5de', edgecolor='k', linewidth=0.5)
        bars2 = ax.bar(x, s['FDR < 0.10'], w, label='FDR < 0.10',
                       color='#f4a582', edgecolor='k', linewidth=0.5)
        bars3 = ax.bar(x + w, s['FDR < 0.05'], w, label='FDR < 0.05',
                       color='#d73027', edgecolor='k', linewidth=0.5)

        for bars in [bars1, bars2, bars3]:
            for bar in bars:
                h = bar.get_height()
                if h > 0:
                    ax.text(bar.get_x() + bar.get_width()/2., h + 0.3,
                           f'{int(h)}', ha='center', va='bottom', fontsize=11)

        ax.set_xticks(x)
        xlabels = [SOURCE_LABELS.get(src, src) for src in s['base_source']]
        ax.set_xticklabels(xlabels, rotation=15, ha='right', fontsize=11)
        ax.set_ylabel('Number of significant cell types', fontsize=14)
        level_label = level.replace('_neurons', ' neurons').capitalize()
        ax.set_title(f'{level_label}\n(n = {s["n_total"].iloc[0]} types)',
                     fontsize=16, fontweight='bold')
        ax.legend(fontsize=10, loc='upper left')

    plt.suptitle(f'Significant CPS-associated Cell Types{title_extra}',
                 fontsize=18, fontweight='bold', y=1.02)
    plt.tight_layout()
    out = os.path.join(IO_DIR, f"significance_summary_{suffix}.png")
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out}")


def plot_forest_top_snrnaseq(df, level, suffix, title_extra=''):
    """Forest plot: top snRNAseq hits with MERFISH comparison."""
    sub = df[df['level'] == level].copy()

    # Find sources
    snr_source = [s for s in sub['source'].unique() if 'snrnaseq' in s][0]
    snr = sub[sub['source'] == snr_source].sort_values('P.Value')

    n_show = min(18, len(snr))
    top_types = snr.head(n_show)['celltype'].tolist()

    fig, ax = plt.subplots(figsize=(10, max(6, n_show * 0.5)))

    sources_in_data = sorted(sub['source'].unique())
    # Order: snrnaseq first, then merfish
    sources_order = [snr_source] + [s for s in sources_in_data if s != snr_source]

    n_src = len(sources_order)
    offsets = np.linspace(0.2, -0.2, n_src)
    markers_list = ['o', 's', 'D', 'v']

    for j, ct in enumerate(reversed(top_types)):
        for k, src in enumerate(sources_order):
            row = sub[(sub['source'] == src) & (sub['celltype'] == ct)]
            if len(row) == 0:
                continue
            row = row.iloc[0]
            y = j + offsets[k]
            se = row.get('SE', row['logFC'] / row['t'] if row['t'] != 0 else 0)

            base_src = get_base_source(src)
            color = SOURCE_COLORS.get(base_src, 'grey')
            marker = markers_list[k % len(markers_list)]
            alpha = 1.0 if row['FDR'] < 0.10 else 0.4

            ax.errorbar(row['logFC'], y, xerr=abs(se) * 1.96,
                       fmt=marker, color=color, markersize=8,
                       capsize=3, alpha=alpha, linewidth=1.5)

    # Legend
    for k, src in enumerate(sources_order):
        base_src = get_base_source(src)
        ax.plot([], [], markers_list[k % len(markers_list)],
               color=SOURCE_COLORS.get(base_src, 'grey'),
               label=SOURCE_LABELS.get(base_src, src), markersize=8)

    ax.set_yticks(range(len(top_types)))
    ax.set_yticklabels(list(reversed(top_types)), fontsize=12)
    ax.axvline(0, color='grey', linestyle='-', alpha=0.3)
    ax.set_xlabel('logFC (CPS effect ± 95% CI)', fontsize=14)
    level_label = level.replace('_neurons', ' — Neurons Only').capitalize()
    ax.set_title(f'Top snRNAseq CPS Effects: Cross-Platform ({level_label}){title_extra}',
                 fontsize=16, fontweight='bold')
    ax.legend(loc='lower right', fontsize=11)

    plt.tight_layout()
    out = os.path.join(IO_DIR, f"forest_top_{suffix}.png")
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out}")


def compute_direction_agreement(df, levels):
    """Compute direction-of-effect agreement between snRNAseq and MERFISH."""
    results = []
    for level in levels:
        sub = df[df['level'] == level]
        snr_source = [s for s in sub['source'].unique() if 'snrnaseq' in s]
        if not snr_source:
            continue
        snr_source = snr_source[0]

        snr = sub[sub['source'] == snr_source][['celltype', 'logFC', 'P.Value', 'FDR']].copy()
        snr = snr.rename(columns={'logFC': 'logFC_snr', 'P.Value': 'p_snr', 'FDR': 'fdr_snr'})

        for ms in sub['source'].unique():
            if 'merfish' not in ms:
                continue
            mer = sub[sub['source'] == ms][['celltype', 'logFC', 'P.Value', 'FDR']].copy()
            mer = mer.rename(columns={'logFC': 'logFC_mer', 'P.Value': 'p_mer', 'FDR': 'fdr_mer'})

            merged = snr.merge(mer, on='celltype', how='inner')

            same_dir_all = (np.sign(merged['logFC_snr']) == np.sign(merged['logFC_mer'])).mean()

            sig = merged[merged['p_snr'] < 0.05]
            same_dir_sig = (np.sign(sig['logFC_snr']) == np.sign(sig['logFC_mer'])).mean() if len(sig) > 0 else np.nan

            fdr = merged[merged['fdr_snr'] < 0.10]
            same_dir_fdr = (np.sign(fdr['logFC_snr']) == np.sign(fdr['logFC_mer'])).mean() if len(fdr) > 0 else np.nan

            base_ms = get_base_source(ms)
            results.append({
                'level': level,
                'merfish_source': base_ms,
                'n_shared': len(merged),
                'direction_agree_all': same_dir_all,
                'direction_agree_nom': same_dir_sig,
                'n_nom': len(sig),
                'direction_agree_fdr': same_dir_fdr,
                'n_fdr': len(fdr),
            })

    return pd.DataFrame(results)


def plot_direction_agreement(df, levels, suffix, title_extra=''):
    """Plot direction-of-effect agreement."""
    results = compute_direction_agreement(df, levels)

    print(f"\n  Direction-of-effect agreement ({suffix}):")
    for _, row in results.iterrows():
        print(f"    {row['level']} / {row['merfish_source']}:")
        print(f"      All: {row['direction_agree_all']:.1%} ({row['n_shared']} types)")
        print(f"      snRNAseq nom p<0.05: {row['direction_agree_nom']:.1%} ({row['n_nom']} types)")
        if row['n_fdr'] > 0:
            print(f"      snRNAseq FDR<0.10: {row['direction_agree_fdr']:.1%} ({row['n_fdr']} types)")

    fig, ax = plt.subplots(figsize=(max(8, len(results) * 2.5), 5))

    labels = []
    agree_all = []
    agree_sig = []
    for _, row in results.iterrows():
        level_label = row['level'].replace('_neurons', '\nneurons').capitalize()
        ms_label = SOURCE_LABELS.get(row['merfish_source'], row['merfish_source']).split('(')[0].strip()
        labels.append(f"{level_label}\n{ms_label}")
        agree_all.append(row['direction_agree_all'] * 100)
        agree_sig.append(row['direction_agree_nom'] * 100 if not np.isnan(row['direction_agree_nom']) else 0)

    x = np.arange(len(labels))
    w = 0.35
    ax.bar(x - w/2, agree_all, w, label='All cell types', color='#92c5de',
           edgecolor='k', linewidth=0.5)
    ax.bar(x + w/2, agree_sig, w, label='snRNAseq nom p < 0.05', color='#f4a582',
           edgecolor='k', linewidth=0.5)

    for xi, (a, s) in enumerate(zip(agree_all, agree_sig)):
        ax.text(xi - w/2, a + 1, f'{a:.0f}%', ha='center', va='bottom', fontsize=10)
        if s > 0:
            ax.text(xi + w/2, s + 1, f'{s:.0f}%', ha='center', va='bottom', fontsize=10)

    ax.axhline(50, color='grey', linestyle='--', alpha=0.5, label='Chance (50%)')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel('Direction agreement (%)', fontsize=14)
    ax.set_ylim(0, 110)
    ax.set_title(f'Direction-of-Effect Agreement: snRNAseq vs MERFISH{title_extra}',
                 fontsize=16, fontweight='bold')
    ax.legend(fontsize=11, loc='upper right')

    plt.tight_layout()
    out = os.path.join(IO_DIR, f"direction_agreement_{suffix}.png")
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out}")

    results.to_csv(os.path.join(IO_DIR, f"direction_agreement_{suffix}.csv"), index=False)
    return results


def run_analysis_set(df, levels, suffix_base, title_extra=''):
    """Run full set of plots for a given set of levels."""
    for level in levels:
        level_suffix = level.replace(' ', '_')
        print(f"\n--- Volcano: {level} ---")
        plot_volcanos(df, level, level_suffix, title_extra)

        print(f"\n--- Effect comparison: {level} ---")
        plot_effect_comparison(df, level, level_suffix, title_extra)

        print(f"\n--- Forest plot: {level} ---")
        plot_forest_top_snrnaseq(df, level, level_suffix, title_extra)

    print(f"\n--- Significance summary: {suffix_base} ---")
    plot_significance_summary(df, levels, suffix_base, title_extra)

    print(f"\n--- Direction agreement: {suffix_base} ---")
    plot_direction_agreement(df, levels, suffix_base, title_extra)


def main():
    print("=" * 70)
    print("Plot crumblr results: SEA-AD pseudoprogression")
    print("=" * 70)

    df = load_all_results()

    # Identify what levels we have
    all_levels = sorted(df['level'].unique())
    neuron_levels = [l for l in all_levels if 'neurons' in l]
    all_type_levels = [l for l in all_levels if 'neurons' not in l]

    print(f"\n  All-type levels: {all_type_levels}")
    print(f"  Neuron-only levels: {neuron_levels}")

    # Run all-types analysis
    if all_type_levels:
        print(f"\n{'='*70}")
        print("ALL CELL TYPES")
        print(f"{'='*70}")
        run_analysis_set(df, all_type_levels, 'all_types', '')

    # Run neuron-only analysis
    if neuron_levels:
        print(f"\n{'='*70}")
        print("NEURONS ONLY")
        print(f"{'='*70}")
        run_analysis_set(df, neuron_levels, 'neurons', '\n(Neurons Only)')

    # Print overall summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    for (source, level), g in df.groupby(['source', 'level']):
        n_fdr05 = (g['FDR'] < 0.05).sum()
        n_fdr10 = (g['FDR'] < 0.10).sum()
        n_total = len(g)
        base_src = get_base_source(source)
        label = SOURCE_LABELS.get(base_src, source)
        print(f"  {label} / {level}: {n_fdr05} FDR<0.05, {n_fdr10} FDR<0.10 of {n_total}")

    print("\nDone!")


if __name__ == "__main__":
    main()
