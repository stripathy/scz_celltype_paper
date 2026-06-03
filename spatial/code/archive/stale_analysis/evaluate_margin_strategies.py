#!/usr/bin/env python3
"""
Evaluate margin filtering strategies for the full dataset.

Compares:
  1. Current: 5th percentile per-sample margin filter (corr_qc_pass)
  2. Absolute margin < 0.02 replacing the 5th percentile filter
  3. Both combined (5th pctl + absolute margin < 0.02)
  4. 5th percentile + absolute margin < 0.02 for L6b ONLY
  5. Unfiltered (no margin filter, just base QC)

For each, computes:
  - Total cells lost per sample, per cell type
  - Overlap between filters
  - Depth distributions vs MERFISH
  - Cell type proportions vs MERFISH
  - Builds crumblr input CSVs for downstream R analysis
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from scipy import stats, sparse
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    H5AD_DIR, SAMPLE_TO_DX, EXCLUDE_SAMPLES,
    SUBCLASS_TO_CLASS, load_sample_adata, load_merfish_cortical,
)

OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), "output", "l6b_diagnostics")
os.makedirs(OUT_DIR, exist_ok=True)

SAMPLES = [s for s in SAMPLE_TO_DX if s not in EXCLUDE_SAMPLES]
ABS_MARGIN_THRESH = 0.02

plt.rcParams.update({
    'font.size': 14, 'axes.titlesize': 18, 'axes.labelsize': 16,
    'xtick.labelsize': 11, 'ytick.labelsize': 11,
    'figure.facecolor': 'white',
})


def load_all_samples():
    """Load obs from all samples with both QC columns + margins."""
    print("Loading all samples...")
    rows = []
    for i, sample in enumerate(SAMPLES):
        print(f"  [{i+1}/{len(SAMPLES)}] {sample}...", end=" ", flush=True)
        adata = load_sample_adata(sample)
        obs = adata.obs.copy()

        # Base QC (spatial + UMI filters, before margin filter)
        obs['base_qc'] = obs['qc_pass'].values.astype(bool)
        obs['corr_qc'] = obs['corr_qc_pass'].values.astype(bool)
        obs['sample_id'] = sample
        obs['dx'] = SAMPLE_TO_DX[sample]

        # Only keep base-QC-pass cells (spatial + UMI filters)
        obs = obs[obs['base_qc']].copy()

        n_total = len(obs)
        n_corr_qc = obs['corr_qc'].sum()
        n_margin_fail = n_total - n_corr_qc
        print(f"{n_total} base QC, {n_corr_qc} corr QC "
              f"({n_margin_fail} margin-filtered = {100*n_margin_fail/n_total:.1f}%)")
        rows.append(obs)

    df = pd.concat(rows, axis=0)
    print(f"\nTotal base-QC cells: {len(df):,}")
    print(f"Total corr-QC cells: {df['corr_qc'].sum():,}")
    return df


def define_strategies(df):
    """Define filtering strategies as boolean masks (True = KEEP)."""
    strategies = {}

    # 0. Unfiltered (base QC only — no margin filter at all)
    strategies['Unfiltered\n(base QC only)'] = pd.Series(True, index=df.index)

    # 1. Current: 5th percentile per-sample margin (corr_qc_pass)
    strategies['Current\n(5th pctl)'] = df['corr_qc'].astype(bool)

    # 2. Absolute margin >= 0.02 only (replacing 5th pctl)
    strategies[f'Absolute\n(margin≥{ABS_MARGIN_THRESH})'] = (
        df['corr_subclass_margin'] >= ABS_MARGIN_THRESH
    )

    # 3. Both: 5th pctl AND absolute margin >= 0.02
    strategies[f'Both\n(5th pctl + margin≥{ABS_MARGIN_THRESH})'] = (
        df['corr_qc'].astype(bool) &
        (df['corr_subclass_margin'] >= ABS_MARGIN_THRESH)
    )

    # 4. 5th pctl + absolute margin for L6b only
    is_l6b = df['corr_subclass'] == 'L6b'
    l6b_margin_fail = is_l6b & (df['corr_subclass_margin'] < ABS_MARGIN_THRESH)
    strategies[f'5th pctl +\nL6b margin≥{ABS_MARGIN_THRESH}'] = (
        df['corr_qc'].astype(bool) & ~l6b_margin_fail
    )

    return strategies


def analyze_overlap(df):
    """Analyze overlap between existing 5th pctl filter and absolute margin filter."""
    print("\n" + "=" * 70)
    print("FILTER OVERLAP ANALYSIS")
    print("=" * 70)

    pctl_fail = ~df['corr_qc'].astype(bool)
    abs_fail = df['corr_subclass_margin'] < ABS_MARGIN_THRESH

    both_fail = pctl_fail & abs_fail
    only_pctl = pctl_fail & ~abs_fail
    only_abs = ~pctl_fail & abs_fail
    neither = ~pctl_fail & ~abs_fail

    total = len(df)
    print(f"Total base-QC cells: {total:,}")
    print(f"\n  5th pctl fails:     {pctl_fail.sum():>7,} ({100*pctl_fail.mean():.1f}%)")
    print(f"  Absolute<0.02 fails: {abs_fail.sum():>7,} ({100*abs_fail.mean():.1f}%)")
    print(f"\n  Both fail:          {both_fail.sum():>7,} ({100*both_fail.mean():.1f}%)")
    print(f"  Only 5th pctl fail: {only_pctl.sum():>7,} ({100*only_pctl.mean():.1f}%)")
    print(f"  Only absolute fail: {only_abs.sum():>7,} ({100*only_abs.mean():.1f}%)")
    print(f"  Neither fail:       {neither.sum():>7,} ({100*neither.mean():.1f}%)")

    # Overlap: of cells caught by 5th pctl, how many also caught by absolute?
    if pctl_fail.sum() > 0:
        overlap = both_fail.sum() / pctl_fail.sum()
        print(f"\n  Of 5th-pctl failures, {100*overlap:.1f}% also caught by absolute<0.02")
    if abs_fail.sum() > 0:
        overlap2 = both_fail.sum() / abs_fail.sum()
        print(f"  Of absolute<0.02 failures, {100*overlap2:.1f}% also caught by 5th-pctl")

    # Per cell type breakdown
    print(f"\n── Per-subclass: absolute<0.02 failure rate ──")
    for st in sorted(df['corr_subclass'].unique()):
        mask = df['corr_subclass'] == st
        n = mask.sum()
        n_fail = (mask & abs_fail).sum()
        n_pctl_fail = (mask & pctl_fail).sum()
        if n > 0:
            print(f"  {st:20s}: {n_fail:>5,}/{n:>6,} abs fail ({100*n_fail/n:5.1f}%)  |  "
                  f"{n_pctl_fail:>5,} pctl fail ({100*n_pctl_fail/n:5.1f}%)")

    return pctl_fail, abs_fail


def compute_strategy_metrics(df, strategies, mer_depths_by_subclass, mer_proportions):
    """Compute key metrics for each strategy."""
    print("\n" + "=" * 70)
    print("STRATEGY METRICS")
    print("=" * 70)

    all_metrics = []
    strategy_cell_counts = {}

    for sname, keep_mask in strategies.items():
        filtered = df[keep_mask].copy()
        n_kept = len(filtered)
        n_dropped = len(df) - n_kept
        pct_dropped = 100 * n_dropped / len(df)

        # Per-sample counts
        per_sample = filtered.groupby('sample_id').size()

        # Per-subclass counts
        per_subclass = filtered[filtered['spatial_domain'] == 'Cortical'].groupby(
            'corr_subclass').size()
        total_cortical = per_subclass.sum()

        # Proportions
        xen_proportions = per_subclass / total_cortical

        # Compare proportions to MERFISH
        shared = sorted(set(xen_proportions.index) & set(mer_proportions.index))
        if len(shared) > 1:
            r_prop, p_prop = stats.pearsonr(
                [xen_proportions.get(s, 0) for s in shared],
                [mer_proportions.get(s, 0) for s in shared]
            )
        else:
            r_prop, p_prop = np.nan, np.nan

        # L6b specific: depth distribution vs MERFISH
        l6b = filtered[(filtered['corr_subclass'] == 'L6b') &
                        (filtered['spatial_domain'] == 'Cortical')]
        l6b_depths = l6b['predicted_norm_depth'].dropna()
        mer_l6b_depths = mer_depths_by_subclass.get('L6b', pd.Series())

        if len(l6b_depths) > 0 and len(mer_l6b_depths) > 0:
            frac_misplaced = (l6b_depths < 0.5).mean()
            # KL divergence
            bins = np.linspace(-0.1, 1.2, 50)
            h_x, _ = np.histogram(l6b_depths, bins=bins, density=True)
            h_m, _ = np.histogram(mer_l6b_depths, bins=bins, density=True)
            eps = 1e-10
            h_x = h_x / h_x.sum() + eps
            h_m = h_m / h_m.sum() + eps
            kl = np.sum(h_x * np.log(h_x / h_m))
        else:
            frac_misplaced = np.nan
            kl = np.nan

        # Overall depth correlation with MERFISH (all subclasses)
        depth_corrs = []
        for st in shared:
            xen_d = filtered[(filtered['corr_subclass'] == st) &
                             (filtered['spatial_domain'] == 'Cortical')]['predicted_norm_depth'].dropna()
            mer_d = mer_depths_by_subclass.get(st, pd.Series())
            if len(xen_d) > 10 and len(mer_d) > 10:
                depth_corrs.append({
                    'subclass': st,
                    'xen_median_depth': xen_d.median(),
                    'mer_median_depth': mer_d.median(),
                })
        if len(depth_corrs) > 2:
            dc = pd.DataFrame(depth_corrs)
            r_depth, p_depth = stats.pearsonr(dc['xen_median_depth'], dc['mer_median_depth'])
        else:
            r_depth, p_depth = np.nan, np.nan

        metrics = {
            'strategy': sname.replace('\n', ' '),
            'n_kept': n_kept,
            'n_dropped': n_dropped,
            'pct_dropped': pct_dropped,
            'n_l6b': len(l6b),
            'l6b_frac_misplaced': frac_misplaced,
            'l6b_kl_vs_merfish': kl,
            'proportion_r_vs_merfish': r_prop,
            'depth_r_vs_merfish': r_depth,
        }
        all_metrics.append(metrics)
        strategy_cell_counts[sname] = per_subclass

        label = sname.replace('\n', ' ')
        print(f"\n{'─'*60}")
        print(f"  {label}")
        print(f"{'─'*60}")
        print(f"  Kept: {n_kept:,} / {len(df):,} ({100-pct_dropped:.1f}%)")
        print(f"  Dropped: {n_dropped:,} ({pct_dropped:.1f}%)")
        print(f"  L6b cells: {len(l6b):,}, misplaced: {100*frac_misplaced:.1f}%")
        print(f"  L6b KL vs MERFISH: {kl:.4f}")
        print(f"  Proportion r vs MERFISH: {r_prop:.4f}")
        print(f"  Median depth r vs MERFISH: {r_depth:.4f}")

    return pd.DataFrame(all_metrics), strategy_cell_counts


def plot_strategy_comparison(metrics_df, strategy_counts, df, strategies,
                              mer_proportions, mer_depths_by_subclass):
    """Generate comprehensive comparison plots."""

    # ══════════════════════════════════════════════════════════════
    # PLOT 1: Summary metrics bar chart
    # ══════════════════════════════════════════════════════════════
    fig, axes = plt.subplots(1, 4, figsize=(22, 6))

    labels = metrics_df['strategy'].tolist()
    x = np.arange(len(labels))
    colors = ['#aaaaaa', '#4A90D9', '#2ca02c', '#d62728', '#ff7f0e']

    ax = axes[0]
    ax.bar(x, metrics_df['pct_dropped'], color=colors)
    ax.set_ylabel('% cells dropped')
    ax.set_title('Cells Dropped')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha='right', fontsize=9)

    ax = axes[1]
    ax.bar(x, 100 * metrics_df['l6b_frac_misplaced'], color=colors)
    ax.set_ylabel('% L6b misplaced (depth<0.5)')
    ax.set_title('L6b Misplacement Rate')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha='right', fontsize=9)

    ax = axes[2]
    ax.bar(x, metrics_df['proportion_r_vs_merfish'], color=colors)
    ax.set_ylabel('Pearson r')
    ax.set_title('Proportion Correlation\nwith MERFISH')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha='right', fontsize=9)
    ax.set_ylim(0.8, 1.0)

    ax = axes[3]
    ax.bar(x, metrics_df['l6b_kl_vs_merfish'], color=colors)
    ax.set_ylabel('KL divergence')
    ax.set_title('L6b Depth KL\nvs MERFISH (lower=better)')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha='right', fontsize=9)

    fig.suptitle('Margin Filter Strategy Comparison', fontsize=22, y=1.03)
    fig.tight_layout()
    path = os.path.join(OUT_DIR, 'margin_strategy_summary.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {path}")

    # ══════════════════════════════════════════════════════════════
    # PLOT 2: Per-sample cell loss across strategies
    # ══════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(18, 7))

    strat_names = list(strategies.keys())
    # Skip unfiltered
    strat_names_plot = strat_names[1:]
    dx_colors = {'Control': '#4A90D9', 'SCZ': '#D94A4A'}

    bar_width = 0.18
    samples_sorted = sorted(SAMPLES, key=lambda s: SAMPLE_TO_DX[s])

    for si, sname in enumerate(strat_names_plot):
        keep = strategies[sname]
        per_sample_kept = df[keep].groupby('sample_id').size()
        per_sample_total = df.groupby('sample_id').size()
        pct_kept = 100 * per_sample_kept / per_sample_total

        positions = np.arange(len(samples_sorted)) + si * bar_width
        vals = [pct_kept.get(s, 0) for s in samples_sorted]
        label = sname.replace('\n', ' ')
        ax.bar(positions, vals, bar_width, label=label, alpha=0.8)

    ax.set_xticks(np.arange(len(samples_sorted)) + bar_width * 1.5)
    ax.set_xticklabels(samples_sorted, rotation=45, ha='right', fontsize=10)
    for i, s in enumerate(samples_sorted):
        ax.get_xticklabels()[i].set_color(dx_colors[SAMPLE_TO_DX[s]])
        ax.get_xticklabels()[i].set_fontweight('bold')

    ax.set_ylabel('% cells kept')
    ax.set_title('Per-Sample Cell Retention by Strategy', fontsize=20)
    ax.legend(fontsize=10, loc='lower left')
    ax.set_ylim(85, 101)

    fig.tight_layout()
    path = os.path.join(OUT_DIR, 'margin_strategy_per_sample.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {path}")

    # ══════════════════════════════════════════════════════════════
    # PLOT 3: Per-subclass cell loss
    # ══════════════════════════════════════════════════════════════
    fig, ax = plt.subplots(figsize=(18, 7))

    # All cortical subclasses
    cortical = df[df['spatial_domain'] == 'Cortical']
    all_subclasses = sorted(cortical['corr_subclass'].unique())

    bar_width = 0.18
    for si, sname in enumerate(strat_names_plot):
        keep = strategies[sname]
        filtered = df[keep & (df['spatial_domain'] == 'Cortical')]
        per_type = filtered.groupby('corr_subclass').size()
        per_type_total = cortical.groupby('corr_subclass').size()
        pct_kept = 100 * per_type / per_type_total

        positions = np.arange(len(all_subclasses)) + si * bar_width
        vals = [pct_kept.get(s, 0) for s in all_subclasses]
        label = sname.replace('\n', ' ')
        ax.bar(positions, vals, bar_width, label=label, alpha=0.8)

    ax.set_xticks(np.arange(len(all_subclasses)) + bar_width * 1.5)
    ax.set_xticklabels(all_subclasses, rotation=45, ha='right', fontsize=10)
    ax.set_ylabel('% cells kept')
    ax.set_title('Per-Subclass Cell Retention by Strategy', fontsize=20)
    ax.legend(fontsize=10, loc='lower left')
    ax.set_ylim(50, 101)

    fig.tight_layout()
    path = os.path.join(OUT_DIR, 'margin_strategy_per_subclass.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {path}")

    # ══════════════════════════════════════════════════════════════
    # PLOT 4: L6b depth distributions for each strategy vs MERFISH
    # ══════════════════════════════════════════════════════════════
    mer_l6b = mer_depths_by_subclass.get('L6b', pd.Series())
    bins = np.linspace(-0.1, 1.2, 50)

    fig, axes = plt.subplots(1, len(strat_names), figsize=(5 * len(strat_names), 5))

    for i, sname in enumerate(strat_names):
        ax = axes[i]
        keep = strategies[sname]
        l6b = df[keep & (df['corr_subclass'] == 'L6b') &
                  (df['spatial_domain'] == 'Cortical')]
        depths = l6b['predicted_norm_depth'].dropna()

        ax.hist(mer_l6b, bins=bins, density=True, alpha=0.5,
                color='#2ca02c', label='MERFISH', edgecolor='white')
        ax.hist(depths, bins=bins, density=True, alpha=0.5,
                color='#1f77b4', label='Xenium', edgecolor='white')
        ax.axvline(0.5, color='orange', ls='--', lw=1.5, alpha=0.7)

        row = metrics_df[metrics_df['strategy'] == sname.replace('\n', ' ')].iloc[0]
        short = sname.replace('\n', ' ')
        ax.set_title(f'{short}\n(n={len(l6b):,}, KL={row["l6b_kl_vs_merfish"]:.3f})',
                     fontsize=13)
        ax.set_xlabel('Depth')
        if i == 0:
            ax.set_ylabel('Density')
        ax.legend(fontsize=9)

    fig.suptitle('L6b Depth Distribution by Strategy (vs MERFISH)', fontsize=20, y=1.03)
    fig.tight_layout()
    path = os.path.join(OUT_DIR, 'margin_strategy_l6b_depth.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {path}")

    # ══════════════════════════════════════════════════════════════
    # PLOT 5: Proportion scatter vs MERFISH for key strategies
    # ══════════════════════════════════════════════════════════════
    key_strats = [strat_names[1], strat_names[2], strat_names[4]]  # current, absolute, L6b-only

    fig, axes = plt.subplots(1, len(key_strats), figsize=(6 * len(key_strats), 6))

    for i, sname in enumerate(key_strats):
        ax = axes[i]
        keep = strategies[sname]
        filtered = df[keep & (df['spatial_domain'] == 'Cortical')]
        per_type = filtered.groupby('corr_subclass').size()
        xen_prop = per_type / per_type.sum()

        shared = sorted(set(xen_prop.index) & set(mer_proportions.index))

        for st in shared:
            cls = SUBCLASS_TO_CLASS.get(st, 'Unknown')
            color = '#1f77b4' if cls == 'Glutamatergic' else '#d62728' if cls == 'GABAergic' else '#7f7f7f'
            ax.scatter(mer_proportions[st], xen_prop.get(st, 0),
                       color=color, s=40, zorder=5)
            if mer_proportions[st] > 0.05 or xen_prop.get(st, 0) > 0.05:
                ax.annotate(st, (mer_proportions[st], xen_prop.get(st, 0)),
                            fontsize=7, ha='left')

        # Diagonal
        lim = max(mer_proportions[shared].max(), xen_prop[shared].max()) * 1.1
        ax.plot([0, lim], [0, lim], 'k--', alpha=0.3)

        r, p = stats.pearsonr(
            [mer_proportions[s] for s in shared],
            [xen_prop.get(s, 0) for s in shared]
        )
        short = sname.replace('\n', ' ')
        ax.set_title(f'{short}\nr={r:.4f}', fontsize=14)
        ax.set_xlabel('MERFISH proportion')
        ax.set_ylabel('Xenium proportion')

    fig.suptitle('Subclass Proportions: Xenium vs MERFISH', fontsize=20, y=1.03)
    fig.tight_layout()
    path = os.path.join(OUT_DIR, 'margin_strategy_proportions.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {path}")


def build_crumblr_inputs(df, strategies):
    """Build crumblr input CSVs for each strategy (for R analysis)."""
    print("\n" + "=" * 70)
    print("BUILDING CRUMBLR INPUTS")
    print("=" * 70)

    crumblr_dir = os.path.join(OUT_DIR, 'crumblr_inputs')
    os.makedirs(crumblr_dir, exist_ok=True)

    # Load sample metadata (age, sex)
    sample_meta = {}
    for sample in SAMPLES:
        sample_meta[sample] = {
            'dx': SAMPLE_TO_DX[sample],
        }

    for sname, keep_mask in strategies.items():
        filtered = df[keep_mask & (df['spatial_domain'] == 'Cortical')].copy()

        # Build count table: donor × supertype
        counts = filtered.groupby(['sample_id', 'corr_subclass']).size().reset_index(name='count')
        totals = filtered.groupby('sample_id').size().reset_index(name='total')
        counts = counts.merge(totals, on='sample_id')
        counts['diagnosis'] = counts['sample_id'].map(SAMPLE_TO_DX)

        safe_name = sname.replace('\n', '_').replace(' ', '_').replace('≥', 'ge').replace('+', '_')
        safe_name = ''.join(c for c in safe_name if c.isalnum() or c == '_')
        path = os.path.join(crumblr_dir, f'crumblr_input_{safe_name}.csv')
        counts.to_csv(path, index=False)
        print(f"  Saved: {path} ({len(counts)} rows, {counts['sample_id'].nunique()} samples)")


def main():
    t0 = time.time()

    # Load data
    df = load_all_samples()

    # Load MERFISH reference
    print("\nLoading MERFISH reference...")
    merfish = load_merfish_cortical()
    mer_proportions = merfish['subclass'].value_counts(normalize=True)
    mer_depths_by_subclass = {
        st: merfish.loc[merfish['subclass'] == st, 'depth']
        for st in merfish['subclass'].unique()
    }
    print(f"MERFISH: {len(merfish):,} cells, {len(mer_proportions)} subclasses")

    # Analyze overlap
    analyze_overlap(df)

    # Define strategies
    strategies = define_strategies(df)

    # Compute metrics
    metrics_df, strategy_counts = compute_strategy_metrics(
        df, strategies, mer_depths_by_subclass, mer_proportions
    )

    # Save metrics
    csv_path = os.path.join(OUT_DIR, 'margin_strategy_metrics.csv')
    metrics_df.to_csv(csv_path, index=False)
    print(f"\nSaved: {csv_path}")

    # Plots
    plot_strategy_comparison(metrics_df, strategy_counts, df, strategies,
                              mer_proportions, mer_depths_by_subclass)

    # Build crumblr inputs
    build_crumblr_inputs(df, strategies)

    print(f"\nDone in {time.time()-t0:.0f}s")
    print(f"All outputs: {OUT_DIR}")


if __name__ == '__main__':
    main()
