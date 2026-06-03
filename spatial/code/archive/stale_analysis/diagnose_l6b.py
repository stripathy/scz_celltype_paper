#!/usr/bin/env python3
"""
Diagnostic analysis of L6b cell type: depth distribution, transcriptome quality,
and sample-level variation.

Hypothesis: L6b acts as a "garbage" cell type — neurons with poor transcriptomes
get annotated as L6b. "Misplaced" L6b cells (in upper layers) should have worse
QC metrics than correctly-placed L6b cells (deep layers).

Outputs (saved to output/l6b_diagnostics/):
  1. l6b_depth_xenium_vs_merfish.png — depth distributions compared
  2. l6b_qc_by_depth.png — QC metrics stratified by depth bin
  3. l6b_per_sample.png — per-sample L6b depth distributions + fraction misplaced
  4. l6b_dx_comparison.png — Control vs SCZ misplacement rates
  5. l6b_gene_expression_by_depth.png — marker gene expression by depth
  6. l6b_classifier_confidence.png — classifier confidence for misplaced vs correct
  7. l6b_diagnostics_summary.csv — tabular summary
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
    H5AD_DIR, SAMPLE_TO_DX, EXCLUDE_SAMPLES,
    SUBCLASS_TO_CLASS, load_cells, load_merfish_cortical,
    load_sample_adata,
)

# ── Config ──
OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), "output", "l6b_diagnostics")
os.makedirs(OUT_DIR, exist_ok=True)

# Depth threshold: cells above this are "misplaced" (upper cortex)
# L6 starts at 0.65 in the depth model, so anything < 0.5 is clearly wrong
MISPLACED_THRESHOLD = 0.5

# All samples
SAMPLES = [s for s in SAMPLE_TO_DX if s not in EXCLUDE_SAMPLES]

plt.rcParams.update({
    'font.size': 14,
    'axes.titlesize': 18,
    'axes.labelsize': 16,
    'xtick.labelsize': 13,
    'ytick.labelsize': 13,
    'figure.facecolor': 'white',
})


def load_all_l6b_with_qc():
    """Load L6b cells from all samples with full QC metrics."""
    print("Loading L6b cells from all samples...")
    rows = []
    for i, sample in enumerate(SAMPLES):
        print(f"  [{i+1}/{len(SAMPLES)}] {sample}...", end=" ")
        try:
            adata = load_sample_adata(sample)
            # Get QC-pass cells
            mask = adata.obs['corr_qc_pass'].values.astype(bool)
            obs = adata.obs.loc[mask].copy()

            # Extract L6b cells
            l6b_mask = obs['corr_subclass'] == 'L6b'
            l6b = obs.loc[l6b_mask].copy()

            # Also get all glutamatergic neurons for comparison
            glut_mask = obs['corr_class'] == 'Glutamatergic'

            l6b['sample_id'] = sample
            l6b['dx'] = SAMPLE_TO_DX[sample]

            # Get total cell counts for this sample
            n_total = mask.sum()
            n_l6b = l6b_mask.sum()
            n_glut = glut_mask.sum()

            print(f"{n_l6b} L6b / {n_glut} Glut / {n_total} total")
            rows.append(l6b)
        except Exception as e:
            print(f"ERROR: {e}")

    df = pd.concat(rows, axis=0)
    print(f"\nTotal L6b cells: {len(df)}")
    return df


def load_all_excitatory():
    """Load all excitatory neurons for baseline comparison."""
    print("\nLoading all excitatory neurons for comparison...")
    glut_subclasses = [s for s, c in SUBCLASS_TO_CLASS.items() if c == 'Glutamatergic']
    print(f"  Glutamatergic subclasses: {glut_subclasses}")
    rows = []
    for i, sample in enumerate(SAMPLES):
        print(f"  [{i+1}/{len(SAMPLES)}] {sample}...", end=" ")
        try:
            adata = load_sample_adata(sample)
            mask = adata.obs['corr_qc_pass'].values.astype(bool)
            obs = adata.obs.loc[mask].copy()
            # Filter to cortical
            if 'spatial_domain' in obs.columns:
                obs = obs[obs['spatial_domain'] == 'Cortical']
            glut = obs[obs['corr_subclass'].isin(glut_subclasses)].copy()
            glut['sample_id'] = sample
            glut['dx'] = SAMPLE_TO_DX[sample]
            print(f"{len(glut)} excitatory neurons")
            rows.append(glut)
        except Exception as e:
            print(f"ERROR: {e}")

    df = pd.concat(rows, axis=0)
    print(f"Total excitatory neurons: {len(df)}")
    return df


# ══════════════════════════════════════════════════════════════════════
# PLOT 1: L6b depth distribution — Xenium vs MERFISH
# ══════════════════════════════════════════════════════════════════════
def plot_depth_comparison(l6b_df, merfish_df):
    """Compare L6b depth distributions between Xenium and MERFISH."""
    print("\n── Plot 1: Depth distribution comparison ──")

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    bins = np.linspace(-0.1, 1.2, 60)

    # Panel A: MERFISH L6b depth
    mer_l6b = merfish_df[merfish_df['subclass'] == 'L6b']
    axes[0].hist(mer_l6b['depth'], bins=bins, color='#2ca02c',
                 alpha=0.7, edgecolor='white', linewidth=0.5)
    axes[0].axvline(0.65, color='red', ls='--', lw=2, label='L6 boundary (0.65)')
    axes[0].axvline(MISPLACED_THRESHOLD, color='orange', ls='--', lw=2,
                    label=f'Misplaced threshold ({MISPLACED_THRESHOLD})')
    axes[0].set_title(f'MERFISH L6b (n={len(mer_l6b):,})')
    axes[0].set_xlabel('Normalized depth (0=pia, 1=WM)')
    axes[0].set_ylabel('Cell count')
    axes[0].legend(fontsize=11)

    # Panel B: Xenium L6b depth
    xen_depths = l6b_df['predicted_norm_depth'].dropna()
    axes[1].hist(xen_depths, bins=bins, color='#1f77b4', alpha=0.7,
                 edgecolor='white', linewidth=0.5)
    axes[1].axvline(0.65, color='red', ls='--', lw=2)
    axes[1].axvline(MISPLACED_THRESHOLD, color='orange', ls='--', lw=2)
    frac_misplaced = (xen_depths < MISPLACED_THRESHOLD).mean()
    axes[1].set_title(f'Xenium L6b (n={len(xen_depths):,})\n'
                      f'{frac_misplaced:.1%} above depth {MISPLACED_THRESHOLD}')
    axes[1].set_xlabel('Normalized depth (0=pia, 1=WM)')

    # Panel C: Overlay (normalized density)
    if len(mer_l6b) > 0:
        axes[2].hist(mer_l6b['depth'], bins=bins, density=True,
                     color='#2ca02c', alpha=0.5, label='MERFISH', edgecolor='white')
    axes[2].hist(xen_depths, bins=bins, density=True,
                 color='#1f77b4', alpha=0.5, label='Xenium', edgecolor='white')
    axes[2].axvline(0.65, color='red', ls='--', lw=2)
    axes[2].axvline(MISPLACED_THRESHOLD, color='orange', ls='--', lw=2)
    axes[2].set_title('Density overlay')
    axes[2].set_xlabel('Normalized depth (0=pia, 1=WM)')
    axes[2].set_ylabel('Density')
    axes[2].legend(fontsize=13)

    fig.suptitle('L6b Depth Distribution: Xenium vs MERFISH', fontsize=20, y=1.02)
    fig.tight_layout()
    path = os.path.join(OUT_DIR, 'l6b_depth_xenium_vs_merfish.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


# ══════════════════════════════════════════════════════════════════════
# PLOT 2: QC metrics by depth bin (misplaced vs correctly placed)
# ══════════════════════════════════════════════════════════════════════
def plot_qc_by_depth(l6b_df):
    """Compare QC metrics between misplaced and correctly-placed L6b cells."""
    print("\n── Plot 2: QC metrics by depth ──")

    df = l6b_df.copy()
    df['depth'] = df['predicted_norm_depth']
    df = df.dropna(subset=['depth'])

    # Categorize
    df['location'] = pd.cut(df['depth'],
                            bins=[-np.inf, 0.3, 0.5, 0.65, 0.85, np.inf],
                            labels=['Very shallow\n(<0.3)', 'Shallow\n(0.3-0.5)',
                                    'Transition\n(0.5-0.65)', 'L6\n(0.65-0.85)',
                                    'Deep/WM\n(>0.85)'])

    metrics = {
        'total_counts': 'Total UMI counts',
        'n_genes': 'Genes detected',
        'corr_subclass_corr': 'Classifier correlation',
        'corr_subclass_margin': 'Classifier margin',
    }

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flat

    for idx, (col, label) in enumerate(metrics.items()):
        ax = axes[idx]
        if col not in df.columns:
            ax.text(0.5, 0.5, f'{col}\nnot available', transform=ax.transAxes,
                    ha='center', va='center', fontsize=14)
            continue

        data_by_loc = [df.loc[df['location'] == loc, col].dropna()
                       for loc in df['location'].cat.categories]

        bp = ax.boxplot(data_by_loc, labels=df['location'].cat.categories,
                        patch_artist=True, showfliers=False, widths=0.6)

        colors = ['#d62728', '#ff7f0e', '#ffbb78', '#2ca02c', '#1f77b4']
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)

        # Add sample sizes
        for i, d in enumerate(data_by_loc):
            ax.text(i + 1, ax.get_ylim()[0], f'n={len(d):,}',
                    ha='center', va='bottom', fontsize=10, color='gray')

        ax.set_ylabel(label)
        ax.set_title(label, fontsize=16)

        # Statistical test: very shallow vs L6
        if len(data_by_loc[0]) > 5 and len(data_by_loc[3]) > 5:
            stat, pval = stats.mannwhitneyu(data_by_loc[0], data_by_loc[3],
                                             alternative='two-sided')
            stars = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else 'ns'
            ax.text(0.98, 0.98, f'Shallow vs L6: p={pval:.2e} {stars}',
                    transform=ax.transAxes, ha='right', va='top', fontsize=11,
                    bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    fig.suptitle('L6b Transcriptome Quality by Cortical Depth', fontsize=20, y=1.01)
    fig.tight_layout()
    path = os.path.join(OUT_DIR, 'l6b_qc_by_depth.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


# ══════════════════════════════════════════════════════════════════════
# PLOT 3: Per-sample L6b depth + fraction misplaced
# ══════════════════════════════════════════════════════════════════════
def plot_per_sample(l6b_df):
    """Per-sample L6b depth distributions and misplacement rates."""
    print("\n── Plot 3: Per-sample analysis ──")

    df = l6b_df.dropna(subset=['predicted_norm_depth']).copy()

    # Compute per-sample stats
    sample_stats = []
    for sample in SAMPLES:
        s = df[df['sample_id'] == sample]
        if len(s) == 0:
            continue
        depths = s['predicted_norm_depth']
        sample_stats.append({
            'sample': sample,
            'dx': SAMPLE_TO_DX[sample],
            'n_l6b': len(s),
            'frac_misplaced': (depths < MISPLACED_THRESHOLD).mean(),
            'median_depth': depths.median(),
            'mean_depth': depths.mean(),
            'frac_very_shallow': (depths < 0.3).mean(),
            'median_corr': s['corr_subclass_corr'].median() if 'corr_subclass_corr' in s.columns else np.nan,
            'median_umi': s['total_counts'].median() if 'total_counts' in s.columns else np.nan,
            'median_ngenes': s['n_genes'].median() if 'n_genes' in s.columns else np.nan,
        })

    stats_df = pd.DataFrame(sample_stats).sort_values('frac_misplaced', ascending=False)

    fig, axes = plt.subplots(2, 1, figsize=(18, 12))

    # Panel A: Violin/strip of depth per sample
    ax = axes[0]
    samples_ordered = stats_df['sample'].tolist()
    dx_colors = {'Control': '#4A90D9', 'SCZ': '#D94A4A'}

    positions = range(len(samples_ordered))
    for i, sample in enumerate(samples_ordered):
        s = df[df['sample_id'] == sample]['predicted_norm_depth']
        dx = SAMPLE_TO_DX[sample]

        vp = ax.violinplot([s.values], positions=[i], showmedians=True,
                           showextrema=False, widths=0.7)
        for pc in vp['bodies']:
            pc.set_facecolor(dx_colors[dx])
            pc.set_alpha(0.6)
        vp['cmedians'].set_color('black')

    ax.axhline(0.65, color='red', ls='--', lw=2, alpha=0.7, label='L6 boundary')
    ax.axhline(MISPLACED_THRESHOLD, color='orange', ls='--', lw=2, alpha=0.7,
               label=f'Misplaced threshold ({MISPLACED_THRESHOLD})')
    ax.set_xticks(positions)
    ax.set_xticklabels(samples_ordered, rotation=45, ha='right', fontsize=11)
    ax.set_ylabel('Normalized depth')
    ax.set_title('L6b Depth Distribution per Sample (sorted by misplacement rate)')
    ax.legend(fontsize=12)

    # Color x-tick labels by diagnosis
    for i, label in enumerate(ax.get_xticklabels()):
        dx = SAMPLE_TO_DX[samples_ordered[i]]
        label.set_color(dx_colors[dx])
        label.set_fontweight('bold')

    # Panel B: Bar chart of fraction misplaced
    ax = axes[1]
    colors = [dx_colors[stats_df.iloc[i]['dx']] for i in range(len(stats_df))]
    bars = ax.bar(positions, stats_df['frac_misplaced'], color=colors, alpha=0.7,
                  edgecolor='white', linewidth=0.5)

    # Add count labels
    for i, (_, row) in enumerate(stats_df.iterrows()):
        ax.text(i, row['frac_misplaced'] + 0.01, f"n={row['n_l6b']}",
                ha='center', va='bottom', fontsize=9, rotation=45)

    ax.set_xticks(positions)
    ax.set_xticklabels(samples_ordered, rotation=45, ha='right', fontsize=11)
    for i, label in enumerate(ax.get_xticklabels()):
        dx = SAMPLE_TO_DX[samples_ordered[i]]
        label.set_color(dx_colors[dx])
        label.set_fontweight('bold')

    ax.set_ylabel(f'Fraction L6b with depth < {MISPLACED_THRESHOLD}')
    ax.set_title('L6b Misplacement Rate per Sample')
    ax.legend(handles=[
        Line2D([0], [0], marker='s', color='w', markerfacecolor='#4A90D9',
               markersize=12, label='Control'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='#D94A4A',
               markersize=12, label='SCZ'),
    ], fontsize=13)

    fig.tight_layout()
    path = os.path.join(OUT_DIR, 'l6b_per_sample.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")

    # Save summary table
    csv_path = os.path.join(OUT_DIR, 'l6b_per_sample_stats.csv')
    stats_df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")

    return stats_df


# ══════════════════════════════════════════════════════════════════════
# PLOT 4: Diagnosis comparison
# ══════════════════════════════════════════════════════════════════════
def plot_dx_comparison(l6b_df, sample_stats_df):
    """Compare L6b misplacement between Control and SCZ."""
    print("\n── Plot 4: Diagnosis comparison ──")

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    dx_colors = {'Control': '#4A90D9', 'SCZ': '#D94A4A'}

    # Panel A: Depth density by diagnosis
    ax = axes[0]
    bins = np.linspace(-0.1, 1.2, 50)
    for dx in ['Control', 'SCZ']:
        depths = l6b_df.loc[l6b_df['dx'] == dx, 'predicted_norm_depth'].dropna()
        ax.hist(depths, bins=bins, density=True, alpha=0.5, color=dx_colors[dx],
                label=f'{dx} (n={len(depths):,})', edgecolor='white')
    ax.axvline(0.65, color='red', ls='--', lw=2, alpha=0.7)
    ax.axvline(MISPLACED_THRESHOLD, color='orange', ls='--', lw=2, alpha=0.7)
    ax.set_xlabel('Normalized depth')
    ax.set_ylabel('Density')
    ax.set_title('L6b Depth by Diagnosis')
    ax.legend(fontsize=12)

    # Panel B: Misplacement rate by diagnosis (sample-level)
    ax = axes[1]
    ctrl = sample_stats_df[sample_stats_df['dx'] == 'Control']['frac_misplaced']
    scz = sample_stats_df[sample_stats_df['dx'] == 'SCZ']['frac_misplaced']

    bp = ax.boxplot([ctrl, scz], labels=['Control', 'SCZ'],
                    patch_artist=True, showfliers=True, widths=0.5)
    bp['boxes'][0].set_facecolor(dx_colors['Control'])
    bp['boxes'][1].set_facecolor(dx_colors['SCZ'])
    for b in bp['boxes']:
        b.set_alpha(0.6)

    # Individual points
    jitter = 0.05
    ax.scatter(np.ones(len(ctrl)) + np.random.uniform(-jitter, jitter, len(ctrl)),
               ctrl, color=dx_colors['Control'], s=40, zorder=5, edgecolor='white')
    ax.scatter(np.ones(len(scz)) * 2 + np.random.uniform(-jitter, jitter, len(scz)),
               scz, color=dx_colors['SCZ'], s=40, zorder=5, edgecolor='white')

    stat, pval = stats.mannwhitneyu(ctrl, scz, alternative='two-sided')
    ax.set_title(f'Misplacement Rate\n(Mann-Whitney p={pval:.3f})')
    ax.set_ylabel(f'Fraction L6b with depth < {MISPLACED_THRESHOLD}')

    # Panel C: Classifier confidence by diagnosis × placement
    ax = axes[2]
    if 'corr_subclass_corr' in l6b_df.columns:
        df = l6b_df.dropna(subset=['predicted_norm_depth', 'corr_subclass_corr']).copy()
        df['misplaced'] = df['predicted_norm_depth'] < MISPLACED_THRESHOLD

        groups = {
            'Ctrl\ncorrect': df[(df['dx'] == 'Control') & ~df['misplaced']]['corr_subclass_corr'],
            'Ctrl\nmisplaced': df[(df['dx'] == 'Control') & df['misplaced']]['corr_subclass_corr'],
            'SCZ\ncorrect': df[(df['dx'] == 'SCZ') & ~df['misplaced']]['corr_subclass_corr'],
            'SCZ\nmisplaced': df[(df['dx'] == 'SCZ') & df['misplaced']]['corr_subclass_corr'],
        }

        bp = ax.boxplot(list(groups.values()), labels=list(groups.keys()),
                        patch_artist=True, showfliers=False, widths=0.6)
        box_colors = [dx_colors['Control'], dx_colors['Control'],
                      dx_colors['SCZ'], dx_colors['SCZ']]
        alphas = [0.7, 0.3, 0.7, 0.3]
        for patch, color, alpha in zip(bp['boxes'], box_colors, alphas):
            patch.set_facecolor(color)
            patch.set_alpha(alpha)

        for i, (name, vals) in enumerate(groups.items()):
            ax.text(i + 1, ax.get_ylim()[0], f'n={len(vals):,}',
                    ha='center', va='bottom', fontsize=10, color='gray')

        ax.set_ylabel('Classifier correlation')
        ax.set_title('Classifier Confidence\n(Correct vs Misplaced)')

    fig.suptitle('L6b: Control vs SCZ Comparison', fontsize=20, y=1.02)
    fig.tight_layout()
    path = os.path.join(OUT_DIR, 'l6b_dx_comparison.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


# ══════════════════════════════════════════════════════════════════════
# PLOT 5: Classifier confidence — L6b vs other excitatory types
# ══════════════════════════════════════════════════════════════════════
def plot_classifier_confidence(l6b_df, excit_df):
    """Compare L6b classifier confidence to other excitatory types."""
    print("\n── Plot 5: Classifier confidence comparison ──")

    if 'corr_subclass_corr' not in excit_df.columns:
        print("  SKIP: corr_subclass_corr not available")
        return

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # Panel A: Boxplot of correlation by subclass
    ax = axes[0]
    excit_types = ['L2/3 IT', 'L4 IT', 'L5 IT', 'L5 ET', 'L5/6 NP',
                   'L6 IT', 'L6 IT Car3', 'L6 CT', 'L6b']
    data = []
    labels = []
    for st in excit_types:
        vals = excit_df.loc[excit_df['corr_subclass'] == st, 'corr_subclass_corr'].dropna()
        if len(vals) > 0:
            data.append(vals.values)
            labels.append(f'{st}\n(n={len(vals):,})')

    bp = ax.boxplot(data, labels=labels, patch_artist=True, showfliers=False, widths=0.6)
    for i, patch in enumerate(bp['boxes']):
        color = '#d62728' if 'L6b' in labels[i] else '#1f77b4'
        patch.set_facecolor(color)
        patch.set_alpha(0.6)

    ax.set_ylabel('Classifier correlation')
    ax.set_title('Classifier Confidence by Excitatory Subclass')
    ax.tick_params(axis='x', rotation=45)

    # Panel B: L6b confidence distribution — misplaced vs correct, with comparison to
    # what the misplaced cells' SECOND-best match looks like
    ax = axes[1]
    df = l6b_df.dropna(subset=['predicted_norm_depth', 'corr_subclass_corr']).copy()

    correct = df[df['predicted_norm_depth'] >= MISPLACED_THRESHOLD]['corr_subclass_corr']
    misplaced = df[df['predicted_norm_depth'] < MISPLACED_THRESHOLD]['corr_subclass_corr']
    very_shallow = df[df['predicted_norm_depth'] < 0.3]['corr_subclass_corr']

    bins = np.linspace(0, 1, 50)
    ax.hist(correct, bins=bins, density=True, alpha=0.5, color='#2ca02c',
            label=f'Correct (≥{MISPLACED_THRESHOLD}, n={len(correct):,})')
    ax.hist(misplaced, bins=bins, density=True, alpha=0.5, color='#d62728',
            label=f'Misplaced (<{MISPLACED_THRESHOLD}, n={len(misplaced):,})')
    ax.hist(very_shallow, bins=bins, density=True, alpha=0.3, color='#ff7f0e',
            label=f'Very shallow (<0.3, n={len(very_shallow):,})')

    ax.set_xlabel('Classifier correlation')
    ax.set_ylabel('Density')
    ax.set_title('L6b Classifier Correlation:\nCorrectly Placed vs Misplaced')
    ax.legend(fontsize=11)

    fig.tight_layout()
    path = os.path.join(OUT_DIR, 'l6b_classifier_confidence.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


# ══════════════════════════════════════════════════════════════════════
# PLOT 6: What subclass would misplaced L6b cells be if not L6b?
# ══════════════════════════════════════════════════════════════════════
def plot_l6b_vs_other_excitatory_depth(excit_df):
    """Show depth distributions of ALL excitatory types on one panel,
    highlighting L6b's anomalous spread."""
    print("\n── Plot 6: All excitatory depth distributions ──")

    excit_types = ['L2/3 IT', 'L4 IT', 'L5 IT', 'L5 ET', 'L5/6 NP',
                   'L6 IT', 'L6 IT Car3', 'L6 CT', 'L6b']

    fig, ax = plt.subplots(figsize=(14, 8))

    bins = np.linspace(-0.1, 1.2, 60)

    for st in excit_types:
        depths = excit_df.loc[excit_df['corr_subclass'] == st,
                              'predicted_norm_depth'].dropna()
        if len(depths) == 0:
            continue
        lw = 3 if st == 'L6b' else 1.5
        alpha = 1.0 if st == 'L6b' else 0.6
        color = '#d62728' if st == 'L6b' else None

        counts, edges = np.histogram(depths, bins=bins, density=True)
        centers = (edges[:-1] + edges[1:]) / 2
        kwargs = {'lw': lw, 'alpha': alpha, 'label': f'{st} (n={len(depths):,})'}
        if color:
            kwargs['color'] = color
        ax.plot(centers, counts, **kwargs)

    ax.axvline(0.65, color='gray', ls='--', lw=1.5, alpha=0.5)
    ax.axvline(MISPLACED_THRESHOLD, color='orange', ls='--', lw=1.5, alpha=0.5)
    ax.set_xlabel('Normalized depth (0=pia, 1=WM)')
    ax.set_ylabel('Density')
    ax.set_title('Depth Distributions of All Excitatory Subclasses\n(L6b highlighted in red)',
                 fontsize=18)
    ax.legend(fontsize=11, loc='upper left')

    fig.tight_layout()
    path = os.path.join(OUT_DIR, 'l6b_vs_all_excitatory_depth.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


# ══════════════════════════════════════════════════════════════════════
# PLOT 7: UMI/gene scatter colored by depth — are shallow L6b low-quality?
# ══════════════════════════════════════════════════════════════════════
def plot_qc_scatter(l6b_df):
    """Scatter of UMI vs genes, colored by depth."""
    print("\n── Plot 7: QC scatter colored by depth ──")

    df = l6b_df.dropna(subset=['predicted_norm_depth', 'total_counts', 'n_genes']).copy()

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # Panel A: UMI vs genes colored by depth
    ax = axes[0]
    sc = ax.scatter(df['n_genes'], df['total_counts'],
                    c=df['predicted_norm_depth'], cmap='RdYlBu_r',
                    s=3, alpha=0.3, vmin=0, vmax=1)
    plt.colorbar(sc, ax=ax, label='Normalized depth', shrink=0.8)
    ax.set_xlabel('Genes detected')
    ax.set_ylabel('Total UMI counts')
    ax.set_title('L6b: UMI vs Genes (colored by depth)')

    # Panel B: Classifier correlation vs depth
    ax = axes[1]
    if 'corr_subclass_corr' in df.columns:
        ax.scatter(df['predicted_norm_depth'], df['corr_subclass_corr'],
                   c=df['total_counts'], cmap='viridis', s=3, alpha=0.3,
                   norm=plt.Normalize(vmin=0, vmax=df['total_counts'].quantile(0.95)))
        plt.colorbar(ax.collections[0], ax=ax, label='Total UMI', shrink=0.8)
        ax.axvline(MISPLACED_THRESHOLD, color='orange', ls='--', lw=2)
        ax.axvline(0.65, color='red', ls='--', lw=2, alpha=0.5)
        ax.set_xlabel('Normalized depth')
        ax.set_ylabel('Classifier correlation')
        ax.set_title('Depth vs Classifier Confidence\n(colored by UMI)')

    fig.suptitle('L6b Transcriptome Quality Indicators', fontsize=20, y=1.01)
    fig.tight_layout()
    path = os.path.join(OUT_DIR, 'l6b_qc_scatter.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


# ══════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    import time
    t0 = time.time()

    # Load data
    l6b_df = load_all_l6b_with_qc()
    excit_df = load_all_excitatory()

    # Load MERFISH reference
    print("\nLoading MERFISH reference...")
    try:
        merfish = load_merfish_cortical()
        print(f"  MERFISH: {len(merfish)} cells")
        has_merfish = True
    except Exception as e:
        print(f"  Could not load MERFISH: {e}")
        merfish = pd.DataFrame()
        has_merfish = False

    # Generate plots
    if has_merfish:
        plot_depth_comparison(l6b_df, merfish)

    plot_qc_by_depth(l6b_df)
    sample_stats = plot_per_sample(l6b_df)
    plot_dx_comparison(l6b_df, sample_stats)
    plot_classifier_confidence(l6b_df, excit_df)
    plot_l6b_vs_other_excitatory_depth(excit_df)
    plot_qc_scatter(l6b_df)

    # Print summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    depths = l6b_df['predicted_norm_depth'].dropna()
    print(f"Total L6b cells: {len(l6b_df):,}")
    print(f"Fraction misplaced (depth < {MISPLACED_THRESHOLD}): {(depths < MISPLACED_THRESHOLD).mean():.1%}")
    print(f"Fraction very shallow (depth < 0.3): {(depths < 0.3).mean():.1%}")
    print(f"Median depth: {depths.median():.3f}")

    if 'corr_subclass_corr' in l6b_df.columns:
        correct = l6b_df[l6b_df['predicted_norm_depth'] >= MISPLACED_THRESHOLD]['corr_subclass_corr']
        misplaced = l6b_df[l6b_df['predicted_norm_depth'] < MISPLACED_THRESHOLD]['corr_subclass_corr']
        print(f"\nClassifier correlation (median):")
        print(f"  Correctly placed: {correct.median():.3f}")
        print(f"  Misplaced:        {misplaced.median():.3f}")

    print(f"\nDiagnosis breakdown:")
    for dx in ['Control', 'SCZ']:
        d = l6b_df[l6b_df['dx'] == dx]
        mis = (d['predicted_norm_depth'] < MISPLACED_THRESHOLD).mean()
        print(f"  {dx}: {len(d):,} L6b cells, {mis:.1%} misplaced")

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.0f}s. Output: {OUT_DIR}")
