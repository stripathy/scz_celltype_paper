#!/usr/bin/env python3
"""
Compare cell type proportions across SEA-AD datasets.

Per-donor matched comparison of subclass proportions between:
  1. snRNAseq (1.38M cells, 89 donors)
  2. MERFISH original labels (1.89M cells, 27 donors)
  3. MERFISH reclassified by our pipeline (same cells, our labels)

For the 27 donors with matched data from both technologies, computes:
  - Per-donor, per-subclass proportions
  - Correlation (Pearson, Spearman) between technologies per subclass
  - Scatter plots: snRNAseq vs MERFISH proportions per donor
  - Overall R² and bias for each label set

Also examines class-level proportions and key known biases between
spatial (MERFISH) and dissociated (snRNAseq) technologies.

Usage:
    python3 -u compare_seaad_proportions.py

Output: output/merfish_benchmark/
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

BASE_DIR = os.path.expanduser("~/Github/SCZ_Xenium")
REF_DIR = os.path.join(BASE_DIR, "data", "reference")
OUT_DIR = os.path.join(BASE_DIR, "output", "merfish_benchmark")
os.makedirs(OUT_DIR, exist_ok=True)

# Add analysis config
sys.path.insert(0, os.path.join(BASE_DIR, "code", "analysis"))
from config import SUBCLASS_TO_CLASS, CLASS_COLORS


def load_proportion_data():
    """Load all cell count CSVs and compute matched per-donor proportions."""

    # snRNAseq
    snr_sub = pd.read_csv(os.path.join(REF_DIR, "SEAAD_snrnaseq_counts_by_subclass.csv"))
    snr_sub = snr_sub.rename(columns={"Subclass": "subclass"})

    # MERFISH original
    mer_sub = pd.read_csv(os.path.join(REF_DIR, "SEAAD_merfish_original_counts_by_subclass.csv"))
    mer_sub = mer_sub.rename(columns={"Subclass": "subclass"})

    # MERFISH reclassified (all cells)
    recl_sub = pd.read_csv(os.path.join(REF_DIR, "SEAAD_merfish_reclassified_counts_by_subclass.csv"))
    recl_sub = recl_sub.rename(columns={"Subclass": "subclass"})

    # MERFISH reclassified (QC-pass only)
    recl_qc_sub = pd.read_csv(os.path.join(REF_DIR, "SEAAD_merfish_reclassified_qcpass_counts_by_subclass.csv"))
    recl_qc_sub = recl_qc_sub.rename(columns={"Subclass": "subclass"})

    # MERFISH depth-annotated only (19.5% of cells — used in original SEA-AD paper)
    mer_depth_sub = pd.read_csv(os.path.join(REF_DIR, "SEAAD_merfish_depth_annotated_counts_by_subclass.csv"))
    mer_depth_sub = mer_depth_sub.rename(columns={"Subclass": "subclass"})

    # Find overlapping donors (all 27 MERFISH donors should be in snRNAseq)
    overlap_donors = sorted(set(snr_sub['Donor ID']) & set(mer_sub['Donor ID']))
    print(f"Overlapping donors: {len(overlap_donors)}")

    return {
        'snrnaseq': snr_sub,
        'merfish_orig': mer_sub,
        'merfish_depth': mer_depth_sub,
        'merfish_recl': recl_sub,
        'merfish_recl_qc': recl_qc_sub,
        'overlap_donors': overlap_donors,
    }


def load_proportion_data_supertype():
    """Load supertype-level count CSVs."""

    snr_sup = pd.read_csv(os.path.join(REF_DIR, "SEAAD_snrnaseq_counts_by_supertype.csv"))
    snr_sup = snr_sup.rename(columns={"Supertype": "supertype"})

    mer_sup = pd.read_csv(os.path.join(REF_DIR, "SEAAD_merfish_original_counts_by_supertype.csv"))
    mer_sup = mer_sup.rename(columns={"Supertype": "supertype"})

    recl_sup = pd.read_csv(os.path.join(REF_DIR, "SEAAD_merfish_reclassified_counts_by_supertype.csv"))
    recl_sup = recl_sup.rename(columns={"Supertype": "supertype"})

    recl_qc_sup = pd.read_csv(os.path.join(REF_DIR, "SEAAD_merfish_reclassified_qcpass_counts_by_supertype.csv"))
    recl_qc_sup = recl_qc_sup.rename(columns={"Supertype": "supertype"})

    mer_depth_sup = pd.read_csv(os.path.join(REF_DIR, "SEAAD_merfish_depth_annotated_counts_by_supertype.csv"))
    mer_depth_sup = mer_depth_sup.rename(columns={"Supertype": "supertype"})

    overlap_donors = sorted(set(snr_sup['Donor ID']) & set(mer_sup['Donor ID']))

    return {
        'snrnaseq': snr_sup,
        'merfish_orig': mer_sup,
        'merfish_depth': mer_depth_sup,
        'merfish_recl': recl_sup,
        'merfish_recl_qc': recl_qc_sup,
        'overlap_donors': overlap_donors,
    }


def build_matched_proportions(data, level='subclass'):
    """Build matched proportion tables for overlapping donors.

    Returns a dict of DataFrames: {source: DataFrame(donor x subclass proportions)}
    """
    donors = data['overlap_donors']
    all_types = sorted(set(
        data['snrnaseq'][level].unique().tolist() +
        data['merfish_orig'][level].unique().tolist()
    ))

    result = {}
    for source_name, df in [('snrnaseq', data['snrnaseq']),
                             ('merfish_orig', data['merfish_orig']),
                             ('merfish_depth', data['merfish_depth']),
                             ('merfish_recl', data['merfish_recl']),
                             ('merfish_recl_qc', data['merfish_recl_qc'])]:
        # Filter to overlapping donors
        df_filt = df[df['Donor ID'].isin(donors)].copy()

        # Pivot to donor x subclass proportion matrix
        pivot = df_filt.pivot_table(
            index='Donor ID', columns=level, values='proportion',
            fill_value=0
        )

        # Ensure all types present
        for t in all_types:
            if t not in pivot.columns:
                pivot[t] = 0.0
        pivot = pivot[all_types]

        result[source_name] = pivot

    return result, all_types


def compute_per_subclass_correlations(props, all_types):
    """Compute per-subclass correlation between snRNAseq and each MERFISH label set."""

    records = []
    for merfish_source in ['merfish_orig', 'merfish_depth', 'merfish_recl', 'merfish_recl_qc']:
        for sub in all_types:
            snr_vals = props['snrnaseq'][sub].values
            mer_vals = props[merfish_source][sub].values

            # Skip types not present in either
            if snr_vals.sum() == 0 and mer_vals.sum() == 0:
                continue

            r_pearson, p_pearson = pearsonr(snr_vals, mer_vals) if len(snr_vals) > 2 else (np.nan, np.nan)
            r_spearman, p_spearman = spearmanr(snr_vals, mer_vals) if len(snr_vals) > 2 else (np.nan, np.nan)

            # Mean proportions
            snr_mean = snr_vals.mean()
            mer_mean = mer_vals.mean()
            bias = mer_mean - snr_mean  # positive = MERFISH overrepresents

            cls = SUBCLASS_TO_CLASS.get(sub, 'Unknown')

            records.append({
                'subclass': sub,
                'class': cls,
                'merfish_labels': merfish_source.replace('merfish_', ''),
                'pearson_r': r_pearson,
                'pearson_p': p_pearson,
                'spearman_r': r_spearman,
                'spearman_p': p_spearman,
                'snrnaseq_mean_prop': snr_mean,
                'merfish_mean_prop': mer_mean,
                'bias': bias,
                'fold_change': mer_mean / snr_mean if snr_mean > 0 else np.inf,
            })

    return pd.DataFrame(records)


def plot_scatter_all_subclasses(props, all_types, corr_df,
                                suffix='', ylabel_level='Subclass'):
    """Grand scatter: all proportions across all donors, snRNAseq vs MERFISH."""

    fig, axes = plt.subplots(1, 4, figsize=(32, 7))
    sources = [
        ('merfish_orig', 'MERFISH Original (all cells)'),
        ('merfish_depth', 'MERFISH Original (depth-annotated)'),
        ('merfish_recl', 'MERFISH Our Pipeline (all cells)'),
        ('merfish_recl_qc', 'MERFISH Our Pipeline (QC-pass)'),
    ]

    for ax, (src, title) in zip(axes, sources):
        # Collect all (snRNAseq proportion, MERFISH proportion) pairs
        x_all, y_all, colors, labels = [], [], [], []
        for sub in all_types:
            snr = props['snrnaseq'][sub].values
            mer = props[src][sub].values
            cls = SUBCLASS_TO_CLASS.get(sub, 'Unknown')
            color = CLASS_COLORS.get(cls, '#888888')
            for s, m in zip(snr, mer):
                x_all.append(s)
                y_all.append(m)
                colors.append(color)
                labels.append(sub)

        x_all = np.array(x_all)
        y_all = np.array(y_all)

        ax.scatter(x_all, y_all, c=colors, s=15 if not suffix else 8,
                   alpha=0.5, edgecolors='none')

        # Diagonal
        lim = max(x_all.max(), y_all.max()) * 1.05
        ax.plot([0, lim], [0, lim], 'k--', alpha=0.3, lw=1)
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)

        # Overall correlation
        r, p = pearsonr(x_all, y_all)
        ax.set_xlabel("snRNAseq Proportion (per donor)", fontsize=14)
        ax.set_ylabel("MERFISH Proportion (per donor)", fontsize=14)
        ax.set_title(f"{title}\nPearson r = {r:.3f}", fontsize=14)
        ax.set_aspect('equal')

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=CLASS_COLORS[c], label=c)
                       for c in ['Glutamatergic', 'GABAergic', 'Non-neuronal']]
    axes[-1].legend(handles=legend_elements, loc='upper left', fontsize=11)

    plt.suptitle(f"Per-Donor {ylabel_level} Proportions: snRNAseq vs MERFISH (27 matched donors)",
                 fontsize=16, y=1.02)
    plt.tight_layout()
    path = os.path.join(OUT_DIR, f"proportion_scatter_snrnaseq_vs_merfish{suffix}.png")
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {path}")


def plot_per_subclass_scatter(props, subclasses_to_plot):
    """Individual scatter plots for selected subclasses, one panel each."""

    n_sub = len(subclasses_to_plot)
    n_cols = min(6, n_sub)
    n_rows = (n_sub + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    axes_flat = axes.flatten()

    for i, sub in enumerate(subclasses_to_plot):
        ax = axes_flat[i]
        cls = SUBCLASS_TO_CLASS.get(sub, 'Unknown')
        color = CLASS_COLORS.get(cls, '#888888')

        snr = props['snrnaseq'][sub].values
        mer_orig = props['merfish_orig'][sub].values
        mer_recl = props['merfish_recl_qc'][sub].values

        # snRNAseq vs MERFISH original
        ax.scatter(snr, mer_orig, c=color, s=30, alpha=0.6, edgecolors='none',
                   label='Original', marker='o')
        # snRNAseq vs MERFISH reclassified
        ax.scatter(snr, mer_recl, c=color, s=30, alpha=0.6, edgecolors='black',
                   linewidths=0.5, label='Our pipeline', marker='s')

        # Diagonal
        lim = max(snr.max(), mer_orig.max(), mer_recl.max()) * 1.1 + 0.001
        ax.plot([0, lim], [0, lim], 'k--', alpha=0.3, lw=1)
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)

        # Correlations
        r_orig, _ = pearsonr(snr, mer_orig) if len(snr) > 2 else (0, 1)
        r_recl, _ = pearsonr(snr, mer_recl) if len(snr) > 2 else (0, 1)

        ax.set_title(f"{sub}\nr_orig={r_orig:.2f}, r_ours={r_recl:.2f}",
                     fontsize=11)
        ax.set_xlabel("snRNAseq", fontsize=9)
        ax.set_ylabel("MERFISH", fontsize=9)
        ax.set_aspect('equal')

        if i == 0:
            ax.legend(fontsize=8, loc='upper left')

    # Hide unused axes
    for j in range(i + 1, len(axes_flat)):
        axes_flat[j].set_visible(False)

    plt.suptitle("Per-Donor Subclass Proportions: snRNAseq vs MERFISH\n"
                 "Circles = original labels, Squares = our pipeline (QC-pass)",
                 fontsize=14, y=1.02)
    plt.tight_layout()
    path = os.path.join(OUT_DIR, "proportion_scatter_per_subclass.png")
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {path}")


def plot_bias_comparison(corr_df):
    """Compare bias (MERFISH - snRNAseq mean proportion) across label sets."""

    # Focus on original vs reclassified_qc
    orig = corr_df[corr_df['merfish_labels'] == 'orig'].set_index('subclass')
    recl = corr_df[corr_df['merfish_labels'] == 'recl_qc'].set_index('subclass')

    common = sorted(set(orig.index) & set(recl.index))
    orig = orig.loc[common]
    recl = recl.loc[common]

    fig, axes = plt.subplots(1, 2, figsize=(20, 8))

    # Panel 1: Bias (fold change)
    ax = axes[0]
    x = np.arange(len(common))
    width = 0.35

    fc_orig = np.clip(orig['fold_change'].values, 0.01, 100)
    fc_recl = np.clip(recl['fold_change'].values, 0.01, 100)

    colors_orig = [CLASS_COLORS.get(SUBCLASS_TO_CLASS.get(s, 'Unknown'), '#888')
                   for s in common]
    colors_recl = colors_orig  # same colors

    bars1 = ax.barh(x - width/2, np.log2(fc_orig), width, color=colors_orig,
                     alpha=0.6, label='Original labels')
    bars2 = ax.barh(x + width/2, np.log2(fc_recl), width, color=colors_recl,
                     alpha=1.0, edgecolor='black', linewidth=0.5,
                     label='Our pipeline (QC-pass)')

    ax.set_yticks(x)
    ax.set_yticklabels(common, fontsize=10)
    ax.axvline(0, color='black', lw=0.5)
    ax.set_xlabel("log₂(MERFISH / snRNAseq) mean proportion", fontsize=12)
    ax.set_title("Proportion Bias: MERFISH vs snRNAseq", fontsize=14)
    ax.legend(fontsize=11)
    ax.invert_yaxis()

    # Panel 2: Per-subclass Pearson r
    ax = axes[1]
    r_orig = orig['pearson_r'].values
    r_recl = recl['pearson_r'].values

    bars1 = ax.barh(x - width/2, r_orig, width, color=colors_orig,
                     alpha=0.6, label='Original labels')
    bars2 = ax.barh(x + width/2, r_recl, width, color=colors_recl,
                     alpha=1.0, edgecolor='black', linewidth=0.5,
                     label='Our pipeline (QC-pass)')

    ax.set_yticks(x)
    ax.set_yticklabels(common, fontsize=10)
    ax.set_xlabel("Pearson r (per-donor proportion, n=27)", fontsize=12)
    ax.set_title("Per-Subclass Correlation with snRNAseq", fontsize=14)
    ax.legend(fontsize=11)
    ax.set_xlim(-0.3, 1.05)
    ax.axvline(0, color='gray', lw=0.5, ls='--')
    ax.invert_yaxis()

    plt.tight_layout()
    path = os.path.join(OUT_DIR, "proportion_bias_and_correlation.png")
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {path}")


def plot_class_proportions(data):
    """Bar chart comparing class-level proportions per donor (matched)."""

    donors = data['overlap_donors']

    # Map MERFISH reclassified class names to match original naming
    class_map_recl = {
        'Glutamatergic': 'Neuronal: Glutamatergic',
        'GABAergic': 'Neuronal: GABAergic',
        'Non-neuronal': 'Non-neuronal and Non-neural',
    }

    sources = {
        'snRNAseq': pd.read_csv(os.path.join(REF_DIR, "SEAAD_snrnaseq_counts_by_class.csv")),
        'MERFISH Original': pd.read_csv(os.path.join(REF_DIR, "SEAAD_merfish_original_counts_by_class.csv")),
        'MERFISH Our Pipeline': pd.read_csv(os.path.join(REF_DIR, "SEAAD_merfish_reclassified_qcpass_counts_by_class.csv")),
    }

    # Harmonize class names
    sources['MERFISH Our Pipeline']['Class'] = sources['MERFISH Our Pipeline']['Class'].map(class_map_recl)

    classes = ['Neuronal: Glutamatergic', 'Neuronal: GABAergic', 'Non-neuronal and Non-neural']
    short_names = ['Glutamatergic', 'GABAergic', 'Non-neuronal']

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    for ax, cls, short in zip(axes, classes, short_names):
        color = CLASS_COLORS.get(short, '#888888')
        positions = []
        for i, (src_name, df) in enumerate(sources.items()):
            df_filt = df[df['Donor ID'].isin(donors) & (df['Class'] == cls)]
            donor_props = df_filt.set_index('Donor ID')['proportion']

            # Ensure all donors present
            vals = [donor_props.get(d, 0.0) for d in donors]
            x = np.arange(len(donors)) + i * 0.25
            ax.bar(x, vals, width=0.23, label=src_name,
                   alpha=0.7 if 'snRNA' in src_name else 1.0,
                   color=color if 'snRNA' not in src_name else '#aaaaaa',
                   edgecolor='black' if 'Our' in src_name else 'none',
                   linewidth=0.5)

        ax.set_xticks(np.arange(len(donors)) + 0.25)
        ax.set_xticklabels([d.split('.')[-1] for d in donors], rotation=90, fontsize=7)
        ax.set_ylabel("Proportion", fontsize=12)
        ax.set_title(short, fontsize=14, fontweight='bold')
        if ax == axes[0]:
            ax.legend(fontsize=9)

    plt.suptitle("Class-Level Proportions per Matched Donor (27 donors)",
                 fontsize=16, y=1.02)
    plt.tight_layout()
    path = os.path.join(OUT_DIR, "class_proportions_per_donor.png")
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {path}")


def run_comparison(data, level='subclass', level_col='subclass'):
    """Run full comparison at a given taxonomy level."""

    print(f"\n{'='*70}")
    print(f"  COMPARISON AT {level.upper()} LEVEL")
    print(f"{'='*70}")

    # Build matched proportion tables
    print(f"\nBuilding matched proportion tables at {level} level...")
    props, all_types = build_matched_proportions(data, level=level_col)
    print(f"  {len(all_types)} types, {len(data['overlap_donors'])} matched donors")

    # Per-type correlations
    print(f"\nComputing per-{level} correlations...")
    corr_df = compute_per_subclass_correlations(props, all_types)
    corr_path = os.path.join(OUT_DIR, f"proportion_correlations_{level}.csv")
    corr_df.to_csv(corr_path, index=False)
    print(f"  Saved: {corr_path}")

    # Print summary table
    print(f"\n  Per-{level} Pearson r (snRNAseq vs MERFISH, 27 donors):")
    print(f"  {level.capitalize():30s} {'Class':15s} {'r_orig':>8s} {'r_depth':>8s} {'r_ours':>8s} "
          f"{'snr_prop':>9s} {'mer_orig':>9s} {'mer_depth':>9s} {'mer_ours':>9s}")
    print("  " + "-" * 130)

    for sub in sorted(all_types):
        cls = SUBCLASS_TO_CLASS.get(sub, 'Unknown')
        row_orig = corr_df[(corr_df['subclass'] == sub) & (corr_df['merfish_labels'] == 'orig')]
        row_depth = corr_df[(corr_df['subclass'] == sub) & (corr_df['merfish_labels'] == 'depth')]
        row_recl = corr_df[(corr_df['subclass'] == sub) & (corr_df['merfish_labels'] == 'recl_qc')]

        if len(row_orig) == 0:
            continue

        ro = row_orig.iloc[0]
        rd = row_depth.iloc[0] if len(row_depth) > 0 else pd.Series(
            {'pearson_r': np.nan, 'merfish_mean_prop': 0})
        rr = row_recl.iloc[0] if len(row_recl) > 0 else pd.Series(
            {'pearson_r': np.nan, 'merfish_mean_prop': 0})

        # Only print types with meaningful abundance
        if ro['snrnaseq_mean_prop'] < 0.001 and ro['merfish_mean_prop'] < 0.001:
            continue

        print(f"  {sub:30s} {cls:15s} {ro['pearson_r']:8.3f} {rd['pearson_r']:8.3f} {rr['pearson_r']:8.3f} "
              f"{ro['snrnaseq_mean_prop']:9.4f} {ro['merfish_mean_prop']:9.4f} {rd['merfish_mean_prop']:9.4f} {rr['merfish_mean_prop']:9.4f}")

    # Overall correlation summary
    print(f"\n  Overall Pearson r across all (donor, {level}) pairs:")
    for src in ['merfish_orig', 'merfish_depth', 'merfish_recl', 'merfish_recl_qc']:
        x = props['snrnaseq'].values.flatten()
        y = props[src].values.flatten()
        r, p = pearsonr(x, y)
        label = src.replace('merfish_', '')
        print(f"    {label:20s}: r = {r:.4f}, p = {p:.2e}")

    return props, all_types, corr_df


def main():
    print("=" * 70)
    print("Compare SEA-AD Proportions: snRNAseq vs MERFISH")
    print("=" * 70)

    # Load data
    data = load_proportion_data()

    # Also load supertype-level data
    data_sup = load_proportion_data_supertype()

    # ── Subclass-level comparison ──
    props_sub, types_sub, corr_sub = run_comparison(data, level='subclass', level_col='subclass')

    # ── Supertype-level comparison ──
    props_sup, types_sup, corr_sup = run_comparison(data_sup, level='supertype', level_col='supertype')

    # ── Generate plots ──
    print("\n" + "=" * 70)
    print("Generating plots...")
    print("=" * 70)

    # Subclass plots
    plot_scatter_all_subclasses(props_sub, types_sub, corr_sub)
    plot_per_subclass_scatter(props_sub, types_sub)
    plot_bias_comparison(corr_sub)
    plot_class_proportions(data)

    # Supertype-level scatter
    plot_scatter_all_subclasses(props_sup, types_sup, corr_sup,
                                suffix='_supertype', ylabel_level='Supertype')

    print("\nDone!")


if __name__ == "__main__":
    main()
