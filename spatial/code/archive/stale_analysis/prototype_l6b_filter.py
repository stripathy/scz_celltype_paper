#!/usr/bin/env python3
"""
Prototype L6b filtering strategies and compare their impact.

Strategy A: Depth-gated filter — reassign L6b cells with depth < 0.5 to 2nd-best type
Strategy B: Margin threshold — reassign L6b cells where margin < threshold to 2nd-best type
Strategy C: Combined — reassign L6b where margin < threshold OR (depth < 0.5 AND margin < relaxed_threshold)

For each strategy, we evaluate:
  1. How many cells are reassigned
  2. What they become
  3. Does the L6b depth distribution match MERFISH better?
  4. Does it affect diagnosis comparisons?

Uses the pre-computed l6b_top3_matches.csv from the rank analysis.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    SAMPLE_TO_DX, EXCLUDE_SAMPLES, SUBCLASS_TO_CLASS,
    load_merfish_cortical,
)

OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), "output", "l6b_diagnostics")

plt.rcParams.update({
    'font.size': 14, 'axes.titlesize': 18, 'axes.labelsize': 16,
    'xtick.labelsize': 12, 'ytick.labelsize': 12,
    'figure.facecolor': 'white',
})


def load_data():
    """Load pre-computed rank data and MERFISH reference."""
    csv_path = os.path.join(OUT_DIR, 'l6b_top3_matches.csv')
    df = pd.read_csv(csv_path, index_col=0)
    print(f"Loaded {len(df):,} L6b cells with top-3 matches")

    # Add diagnosis
    df['dx'] = df['sample_id'].map(SAMPLE_TO_DX)

    # Load MERFISH for comparison
    print("Loading MERFISH reference...")
    merfish = load_merfish_cortical()
    mer_l6b = merfish[merfish['subclass'] == 'L6b']['depth']
    print(f"MERFISH L6b: {len(mer_l6b):,} cells")

    return df, mer_l6b


def apply_strategy(df, name, mask):
    """Apply a filtering strategy. Returns modified df with reassignments."""
    result = df.copy()
    result['reassigned'] = False
    result['final_type'] = 'L6b'

    # For cells flagged by this strategy, reassign to rank2 if rank1 is L6b,
    # otherwise rank1 (since rank1 is already not L6b for some cells)
    flagged = df[mask].copy()
    n_flagged = mask.sum()

    # Reassign: use rank2 if rank1 is L6b, else use rank1
    for idx in flagged.index:
        r1 = df.loc[idx, 'rank1_type']
        r2 = df.loc[idx, 'rank2_type']
        if r1 == 'L6b':
            result.loc[idx, 'final_type'] = r2
        else:
            result.loc[idx, 'final_type'] = r1

    result.loc[mask, 'reassigned'] = True

    print(f"\n{'─'*60}")
    print(f"Strategy: {name}")
    print(f"{'─'*60}")
    print(f"Flagged: {n_flagged:,} / {len(df):,} ({100*n_flagged/len(df):.1f}%)")
    print(f"Remaining L6b: {(result['final_type'] == 'L6b').sum():,}")

    # What do reassigned cells become?
    reassigned = result[result['reassigned']]
    if len(reassigned) > 0:
        new_types = reassigned['final_type'].value_counts()
        print(f"\nReassigned cells become:")
        for t, n in new_types.head(10).items():
            print(f"  {t:20s}: {n:>5,} ({100*n/len(reassigned):.1f}%)")

    return result


def compute_depth_kl(xen_depths, mer_depths, bins):
    """Compute KL divergence between Xenium and MERFISH depth distributions."""
    h_xen, _ = np.histogram(xen_depths, bins=bins, density=True)
    h_mer, _ = np.histogram(mer_depths, bins=bins, density=True)
    # Add small epsilon to avoid log(0)
    eps = 1e-10
    h_xen = h_xen / h_xen.sum() + eps
    h_mer = h_mer / h_mer.sum() + eps
    kl = np.sum(h_xen * np.log(h_xen / h_mer))
    return kl


def main():
    df, mer_l6b_depths = load_data()

    # ══════════════════════════════════════════════════════════════
    # Define strategies
    # ══════════════════════════════════════════════════════════════

    # Scan a range of margin thresholds
    margin_thresholds = [0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05]

    strategies = {}

    # Strategy A: Depth-only (depth < 0.5)
    strategies['A: Depth<0.5'] = df['depth'] < 0.5

    # Strategy B: Margin-only at various thresholds
    for mt in margin_thresholds:
        strategies[f'B: Margin<{mt}'] = df['margin_1v2'] < mt

    # Strategy C: Combined (margin OR shallow+relaxed_margin)
    strategies['C: Margin<0.02 OR (depth<0.5 & margin<0.04)'] = (
        (df['margin_1v2'] < 0.02) |
        ((df['depth'] < 0.5) & (df['margin_1v2'] < 0.04))
    )

    # ══════════════════════════════════════════════════════════════
    # Apply all strategies and collect metrics
    # ══════════════════════════════════════════════════════════════
    bins = np.linspace(-0.1, 1.2, 60)
    results = {}
    metrics = []

    for name, mask in strategies.items():
        res = apply_strategy(df, name, mask)
        results[name] = res

        # Remaining L6b depth distribution
        remaining_l6b = res[res['final_type'] == 'L6b']
        remaining_depths = remaining_l6b['depth'].dropna()

        # Metrics
        n_remaining = len(remaining_l6b)
        n_reassigned = mask.sum()
        frac_reassigned = n_reassigned / len(df)

        # Fraction of remaining L6b that are misplaced
        if len(remaining_depths) > 0:
            frac_misplaced_after = (remaining_depths < 0.5).mean()
            kl_div = compute_depth_kl(remaining_depths, mer_l6b_depths, bins)
            median_depth = remaining_depths.median()
        else:
            frac_misplaced_after = 0
            kl_div = np.inf
            median_depth = np.nan

        # Original misplacement rate
        frac_misplaced_before = (df['depth'] < 0.5).mean()

        metrics.append({
            'strategy': name,
            'n_reassigned': n_reassigned,
            'pct_reassigned': 100 * frac_reassigned,
            'n_remaining_l6b': n_remaining,
            'frac_misplaced_before': frac_misplaced_before,
            'frac_misplaced_after': frac_misplaced_after,
            'misplacement_reduction': 1 - frac_misplaced_after / frac_misplaced_before if frac_misplaced_before > 0 else 0,
            'kl_vs_merfish': kl_div,
            'median_depth_remaining': median_depth,
        })

    metrics_df = pd.DataFrame(metrics)
    print(f"\n{'='*80}")
    print("STRATEGY COMPARISON")
    print(f"{'='*80}")
    print(metrics_df.to_string(index=False))

    csv_path = os.path.join(OUT_DIR, 'l6b_filter_strategy_comparison.csv')
    metrics_df.to_csv(csv_path, index=False)
    print(f"\nSaved: {csv_path}")

    # ══════════════════════════════════════════════════════════════
    # PLOT 1: Strategy comparison — key metrics
    # ══════════════════════════════════════════════════════════════
    fig, axes = plt.subplots(1, 3, figsize=(20, 7))

    # Extract margin-only strategies for x-axis
    margin_rows = metrics_df[metrics_df['strategy'].str.startswith('B:')]
    mt_vals = margin_thresholds[:len(margin_rows)]

    # Panel A: % reassigned vs margin threshold
    ax = axes[0]
    ax.plot(mt_vals, margin_rows['pct_reassigned'], 'o-', color='#1f77b4',
            lw=2, markersize=8, label='Margin-only (B)')
    # Add depth-only as horizontal line
    depth_row = metrics_df[metrics_df['strategy'] == 'A: Depth<0.5'].iloc[0]
    ax.axhline(depth_row['pct_reassigned'], color='#d62728', ls='--', lw=2,
               label=f"Depth<0.5 only (A): {depth_row['pct_reassigned']:.1f}%")
    ax.set_xlabel('Margin threshold')
    ax.set_ylabel('% L6b cells reassigned')
    ax.set_title('Cells Reassigned')
    ax.legend(fontsize=11)

    # Panel B: Misplacement rate after filtering
    ax = axes[1]
    ax.plot(mt_vals, 100 * margin_rows['frac_misplaced_after'], 'o-',
            color='#1f77b4', lw=2, markersize=8, label='Margin-only (B)')
    ax.axhline(100 * depth_row['frac_misplaced_after'], color='#d62728', ls='--', lw=2,
               label=f"Depth<0.5 only (A)")
    ax.axhline(100 * metrics_df.iloc[0]['frac_misplaced_before'], color='gray', ls=':',
               lw=2, label='Before filtering')
    ax.set_xlabel('Margin threshold')
    ax.set_ylabel('% remaining L6b misplaced (depth<0.5)')
    ax.set_title('Misplacement Rate After Filtering')
    ax.legend(fontsize=11)

    # Panel C: KL divergence vs MERFISH
    ax = axes[2]
    ax.plot(mt_vals, margin_rows['kl_vs_merfish'], 'o-', color='#1f77b4',
            lw=2, markersize=8, label='Margin-only (B)')
    ax.axhline(depth_row['kl_vs_merfish'], color='#d62728', ls='--', lw=2,
               label='Depth<0.5 only (A)')
    # Add unfiltered KL
    unfiltered_kl = compute_depth_kl(df['depth'].dropna(), mer_l6b_depths,
                                      np.linspace(-0.1, 1.2, 60))
    ax.axhline(unfiltered_kl, color='gray', ls=':', lw=2, label='Unfiltered')
    ax.set_xlabel('Margin threshold')
    ax.set_ylabel('KL divergence (Xenium vs MERFISH)')
    ax.set_title('Similarity to MERFISH L6b Distribution')
    ax.legend(fontsize=11)

    fig.suptitle('L6b Filter Strategy Comparison', fontsize=22, y=1.02)
    fig.tight_layout()
    path = os.path.join(OUT_DIR, 'l6b_filter_strategies_metrics.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {path}")

    # ══════════════════════════════════════════════════════════════
    # PLOT 2: Depth distributions — before/after for key strategies
    # ══════════════════════════════════════════════════════════════
    key_strategies = [
        'A: Depth<0.5',
        'B: Margin<0.02',
        'B: Margin<0.03',
        'B: Margin<0.04',
    ]

    fig, axes = plt.subplots(1, len(key_strategies) + 1, figsize=(5 * (len(key_strategies) + 1), 5))
    bins = np.linspace(-0.1, 1.2, 50)

    # First panel: unfiltered + MERFISH
    ax = axes[0]
    ax.hist(mer_l6b_depths, bins=bins, density=True, alpha=0.5,
            color='#2ca02c', label='MERFISH', edgecolor='white')
    ax.hist(df['depth'].dropna(), bins=bins, density=True, alpha=0.5,
            color='#1f77b4', label='Xenium (unfiltered)', edgecolor='white')
    ax.axvline(0.5, color='orange', ls='--', lw=1.5, alpha=0.7)
    ax.set_title(f'Unfiltered\n(n={len(df):,}, KL={unfiltered_kl:.3f})')
    ax.set_xlabel('Depth')
    ax.set_ylabel('Density')
    ax.legend(fontsize=9)

    for i, sname in enumerate(key_strategies):
        ax = axes[i + 1]
        res = results[sname]
        remaining = res[res['final_type'] == 'L6b']['depth'].dropna()
        row = metrics_df[metrics_df['strategy'] == sname].iloc[0]

        ax.hist(mer_l6b_depths, bins=bins, density=True, alpha=0.5,
                color='#2ca02c', label='MERFISH', edgecolor='white')
        ax.hist(remaining, bins=bins, density=True, alpha=0.5,
                color='#1f77b4', label=f'Xenium filtered', edgecolor='white')
        ax.axvline(0.5, color='orange', ls='--', lw=1.5, alpha=0.7)

        short_name = sname.split(': ')[1] if ': ' in sname else sname
        ax.set_title(f'{short_name}\n(n={row["n_remaining_l6b"]:,.0f}, '
                     f'KL={row["kl_vs_merfish"]:.3f})')
        ax.set_xlabel('Depth')
        ax.legend(fontsize=9)

    fig.suptitle('L6b Depth Distribution After Filtering (vs MERFISH)', fontsize=20, y=1.03)
    fig.tight_layout()
    path = os.path.join(OUT_DIR, 'l6b_filter_depth_comparison.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {path}")

    # ══════════════════════════════════════════════════════════════
    # PLOT 3: What do reassigned cells become? (for margin=0.02 and 0.03)
    # ══════════════════════════════════════════════════════════════
    fig, axes = plt.subplots(1, 3, figsize=(20, 7))

    for idx, sname in enumerate(['B: Margin<0.02', 'B: Margin<0.03', 'B: Margin<0.04']):
        ax = axes[idx]
        res = results[sname]
        reassigned = res[res['reassigned']]

        if len(reassigned) == 0:
            continue

        # Group by new type, split by depth
        new_types = reassigned['final_type'].value_counts()
        top_types = new_types.head(8).index.tolist()

        # Split into shallow vs deep reassignments
        shallow = reassigned[reassigned['depth'] < 0.5]
        deep = reassigned[reassigned['depth'] >= 0.5]

        x = np.arange(len(top_types))
        w = 0.35

        shallow_counts = [len(shallow[shallow['final_type'] == t]) for t in top_types]
        deep_counts = [len(deep[deep['final_type'] == t]) for t in top_types]

        ax.bar(x - w/2, shallow_counts, w, color='#d62728', alpha=0.7,
               label=f'Shallow (<0.5, n={len(shallow)})')
        ax.bar(x + w/2, deep_counts, w, color='#2ca02c', alpha=0.7,
               label=f'Deep (≥0.5, n={len(deep)})')

        ax.set_xticks(x)
        ax.set_xticklabels(top_types, rotation=45, ha='right')
        ax.set_ylabel('Count')

        short = sname.split(': ')[1]
        ax.set_title(f'{short}\n({len(reassigned)} reassigned)')
        ax.legend(fontsize=10)

    fig.suptitle('What Do Reassigned L6b Cells Become?', fontsize=20, y=1.03)
    fig.tight_layout()
    path = os.path.join(OUT_DIR, 'l6b_reassignment_breakdown.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {path}")

    # ══════════════════════════════════════════════════════════════
    # PLOT 4: Diagnosis impact — does filtering change SCZ vs Ctrl?
    # ══════════════════════════════════════════════════════════════
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for idx, sname in enumerate(['B: Margin<0.02', 'B: Margin<0.03']):
        ax = axes[idx]
        res = results[sname]

        # Per-sample L6b counts before and after
        samples = sorted(df['sample_id'].unique())
        before_counts = df.groupby('sample_id').size()
        after_counts = res[res['final_type'] == 'L6b'].groupby('sample_id').size()

        dx_colors = {'Control': '#4A90D9', 'SCZ': '#D94A4A'}
        for s in samples:
            dx = SAMPLE_TO_DX.get(s, 'Unknown')
            b = before_counts.get(s, 0)
            a = after_counts.get(s, 0)
            pct_lost = 100 * (1 - a / b) if b > 0 else 0
            ax.scatter(b, a, color=dx_colors.get(dx, 'gray'), s=60, zorder=5,
                       edgecolor='white')
            ax.annotate(f'{s}\n(-{pct_lost:.0f}%)', (b, a), fontsize=7,
                        ha='center', va='bottom')

        # Identity line
        lim = max(before_counts.max(), after_counts.max()) * 1.1
        ax.plot([0, lim], [0, lim], 'k--', alpha=0.3)
        ax.set_xlabel('L6b count before')
        ax.set_ylabel('L6b count after')
        short = sname.split(': ')[1]
        ax.set_title(f'{short}')
        ax.legend(handles=[
            Line2D([0], [0], marker='o', color='w', markerfacecolor='#4A90D9',
                   markersize=10, label='Control'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='#D94A4A',
                   markersize=10, label='SCZ'),
        ], fontsize=11)

    fig.suptitle('Per-Sample L6b Count: Before vs After Filtering', fontsize=20, y=1.02)
    fig.tight_layout()
    path = os.path.join(OUT_DIR, 'l6b_filter_dx_impact.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {path}")

    print(f"\nAll outputs in: {OUT_DIR}")


if __name__ == '__main__':
    main()
