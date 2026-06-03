#!/usr/bin/env python3
"""
Compare cell type classification using nuclear-only vs whole-cell counts.

Hypothesis: Spatial doublets arise from cytoplasmic mRNA spillover. Nuclear-only
counts should reduce doublet artifacts and yield cleaner cell type assignments.

Steps:
1. Load centroids from all 24 samples (same as pipeline 02b)
2. Run correlation classifier on Br8667 nuclear counts
3. Compare with existing whole-cell classifications
4. Generate 5 diagnostic figures + summary CSV

Requires: build_nuclear_counts.py must have been run first.

Output figures:
  nuclear_fraction_histogram.png
  nuclear_fraction_by_subclass.png
  nuclear_vs_wholecell_concordance.png
  nuclear_vs_wholecell_doublet_rates.png
  nuclear_classification_quality.png
  nuclear_vs_wholecell_comparison.csv
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
from scipy import sparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    H5AD_DIR, PRESENTATION_DIR, BG_COLOR,
    SAMPLE_TO_DX, EXCLUDE_SAMPLES, SUBCLASS_TO_CLASS,
    CLASS_COLORS,
)

# Classifier module
MODULES_DIR = os.path.expanduser('~/Github/SCZ_Xenium/code/modules')
sys.path.insert(0, MODULES_DIR)
from correlation_classifier import (
    build_subclass_centroids,
    build_supertype_centroids,
    run_two_stage_classifier,
    flag_doublet_cells,
)

# ── Config ──
SAMPLE = 'Br8667'
H5AD_WC = os.path.join(H5AD_DIR, f'{SAMPLE}_annotated.h5ad')
H5AD_NUC = os.path.join(H5AD_DIR, f'{SAMPLE}_nuclear_counts.h5ad')
OUT_DIR = PRESENTATION_DIR
TOP_N = 100  # exemplars per type for centroid building

# GABA/Glut markers for doublet scoring
GABA_MARKERS = ['GAD1', 'GAD2', 'SLC32A1', 'SST', 'PVALB', 'VIP', 'LAMP5']
GLUT_MARKERS = ['CUX2', 'RORB', 'GRIN2A', 'THEMIS']


def load_centroids():
    """Load all 24 samples and build centroids (replicating 02b)."""
    print("=" * 60)
    print("Loading samples for centroid building...")
    print("=" * 60)

    sample_ids = sorted(set(SAMPLE_TO_DX.keys()) - EXCLUDE_SAMPLES)
    adatas = []
    for i, sid in enumerate(sample_ids):
        fpath = os.path.join(H5AD_DIR, f"{sid}_annotated.h5ad")
        if not os.path.exists(fpath):
            print(f"  WARNING: {fpath} not found, skipping")
            continue
        adata = ad.read_h5ad(fpath)
        if 'qc_pass' in adata.obs.columns:
            adata = adata[adata.obs['qc_pass'].values.astype(bool)].copy()
        print(f"  [{i+1:2d}/{len(sample_ids)}] {sid}: {adata.n_obs:,} cells")
        adatas.append(adata)

    print("  Concatenating...")
    combined = ad.concat(adatas, join='outer')
    print(f"  Total: {combined.n_obs:,} cells × {combined.n_vars} genes")

    print(f"\nBuilding subclass centroids (top-{TOP_N})...")
    sub_centroids, sub_counts, gene_names = build_subclass_centroids(
        combined, top_n=TOP_N)

    print(f"Building supertype centroids (top-{TOP_N})...")
    sup_centroids, sup_counts = build_supertype_centroids(
        combined, top_n=TOP_N)

    return sub_centroids, sup_centroids, gene_names


def classify_nuclear(adata_nuc, sub_centroids, sup_centroids, gene_names):
    """Run two-stage classifier on nuclear count matrix."""
    print("\n" + "=" * 60)
    print("Classifying nuclear counts...")
    print("=" * 60)

    # Filter to QC-pass cells (use whole-cell QC status)
    qc_mask = adata_nuc.obs['qc_pass'].values.astype(bool)
    adata_qc = adata_nuc[qc_mask].copy()
    print(f"  {adata_qc.n_obs:,} QC-pass cells")

    # Run classifier
    results = run_two_stage_classifier(
        adata_qc, sub_centroids, sup_centroids, gene_names)

    # Derive class
    results['corr_class'] = results['corr_subclass'].map(
        lambda x: SUBCLASS_TO_CLASS.get(x, 'Unknown'))

    # Detect doublets
    class_labels = results['corr_class'].values
    doublet_suspect, doublet_type, doublet_stats = flag_doublet_cells(
        adata_qc, class_labels, SUBCLASS_TO_CLASS)

    results['doublet_suspect'] = doublet_suspect
    results['doublet_type'] = doublet_type

    n_gg = (doublet_type == 'Glut+GABA').sum()
    n_ga = (doublet_type == 'GABA+GABA').sum()
    print(f"  Nuclear doublets: {doublet_suspect.sum():,} "
          f"({n_gg} Glut+GABA, {n_ga} GABA+GABA)")

    return results, qc_mask


def build_comparison(adata_wc, adata_nuc, results_nuc, qc_mask):
    """Merge whole-cell and nuclear classifications into a comparison DataFrame."""
    obs_wc = adata_wc.obs.copy()
    obs_wc = obs_wc.loc[qc_mask]

    df = pd.DataFrame(index=obs_wc.index)
    df['wc_subclass'] = obs_wc['corr_subclass'].astype(str).values
    df['wc_class'] = obs_wc['corr_class'].astype(str).values
    df['wc_doublet'] = obs_wc['doublet_suspect'].values.astype(bool)
    df['wc_doublet_type'] = obs_wc['doublet_type'].astype(str).values
    df['wc_corr'] = obs_wc['corr_subclass_corr'].values.astype(float)
    df['wc_margin'] = obs_wc['corr_subclass_margin'].values.astype(float)

    df['nuc_subclass'] = results_nuc['corr_subclass'].values
    df['nuc_class'] = results_nuc['corr_class'].values
    df['nuc_doublet'] = results_nuc['doublet_suspect'].values.astype(bool)
    df['nuc_doublet_type'] = results_nuc['doublet_type'].values
    df['nuc_corr'] = results_nuc['corr_subclass_corr'].values.astype(float)
    df['nuc_margin'] = results_nuc['corr_subclass_margin'].values.astype(float)

    # Nuclear fraction and counts
    nuc_obs = adata_nuc.obs.loc[qc_mask]
    df['total_counts'] = obs_wc['total_counts'].values.astype(int)
    df['nuclear_total_counts'] = nuc_obs['nuclear_total_counts'].values.astype(int)
    df['nuclear_fraction'] = nuc_obs['nuclear_fraction'].values.astype(float)

    # Concordance flags
    df['concordant_class'] = df['wc_class'] == df['nuc_class']
    df['concordant_subclass'] = df['wc_subclass'] == df['nuc_subclass']

    return df


# ──────────────────────────────────────────────────────────────────────
# Figures
# ──────────────────────────────────────────────────────────────────────

def fig1_nuclear_fraction_histogram(df):
    """Nuclear fraction distribution by class."""
    fig, ax = plt.subplots(figsize=(10, 6), facecolor='white')

    classes = ['Glutamatergic', 'GABAergic', 'Non-neuronal']
    for cls in classes:
        mask = df['wc_class'] == cls
        vals = df.loc[mask, 'nuclear_fraction']
        ax.hist(vals, bins=80, alpha=0.55, label=f'{cls} (n={mask.sum():,})',
                color=CLASS_COLORS.get(cls, 'gray'), density=True)
        med = vals.median()
        ax.axvline(med, color=CLASS_COLORS.get(cls, 'gray'),
                   linestyle='--', linewidth=2, alpha=0.8)
        ax.text(med + 0.01, ax.get_ylim()[1] * 0.85 - classes.index(cls) * 0.8,
                f'med={med:.2f}', fontsize=12,
                color=CLASS_COLORS.get(cls, 'gray'), fontweight='bold')

    ax.set_xlabel('Nuclear Fraction (nuclear UMI / total UMI)', fontsize=16)
    ax.set_ylabel('Density', fontsize=16)
    ax.set_title('Distribution of Nuclear Transcript Fraction by Cell Class',
                 fontsize=20, fontweight='bold', pad=12)
    ax.legend(fontsize=13, loc='upper right')
    ax.tick_params(labelsize=13)
    ax.set_xlim(0, 1)

    # Overall stats annotation
    med_all = df['nuclear_fraction'].median()
    ax.text(0.02, 0.95, f'Overall median: {med_all:.3f}\n'
            f'N cells: {len(df):,}',
            transform=ax.transAxes, fontsize=12, va='top',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    plt.tight_layout()
    path = os.path.join(OUT_DIR, 'nuclear_fraction_histogram.png')
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def fig2_nuclear_fraction_by_subclass(df):
    """Boxplot of nuclear fraction per subclass."""
    # Compute medians for ordering
    sub_medians = df.groupby('wc_subclass')['nuclear_fraction'].median().sort_values()
    order = sub_medians.index.tolist()

    fig, ax = plt.subplots(figsize=(14, 8), facecolor='white')

    positions = range(len(order))
    colors = [CLASS_COLORS.get(SUBCLASS_TO_CLASS.get(s, ''), 'gray') for s in order]
    counts = [df[df['wc_subclass'] == s].shape[0] for s in order]

    bp = ax.boxplot(
        [df.loc[df['wc_subclass'] == s, 'nuclear_fraction'].values for s in order],
        positions=positions, vert=False, widths=0.7,
        patch_artist=True, showfliers=False,
        medianprops=dict(color='black', linewidth=2),
    )
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)

    ax.set_yticks(positions)
    ylabels = [f'{s}  (n={c:,})' for s, c in zip(order, counts)]
    ax.set_yticklabels(ylabels, fontsize=11)
    ax.set_xlabel('Nuclear Fraction', fontsize=16)
    ax.set_title('Nuclear Fraction by Cell Subclass', fontsize=20,
                 fontweight='bold', pad=12)
    ax.tick_params(labelsize=12)

    # Add median annotations
    for i, s in enumerate(order):
        med = sub_medians[s]
        ax.text(med + 0.005, i, f'{med:.2f}', va='center', fontsize=9,
                fontweight='bold')

    # Legend for class colors
    legend_items = [Line2D([0], [0], marker='s', color='w',
                           markerfacecolor=CLASS_COLORS[c], markersize=12,
                           label=c) for c in classes if c in CLASS_COLORS]
    ax.legend(handles=legend_items, fontsize=12, loc='lower right')

    plt.tight_layout()
    path = os.path.join(OUT_DIR, 'nuclear_fraction_by_subclass.png')
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


classes = ['Glutamatergic', 'GABAergic', 'Non-neuronal']


def fig3_concordance(df):
    """Class-level confusion matrix + subclass concordance bars."""
    fig = plt.figure(figsize=(18, 8), facecolor='white')
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1.5], wspace=0.35)

    # ── Panel A: Class confusion matrix ──
    ax_a = fig.add_subplot(gs[0])
    class_order = ['Glutamatergic', 'GABAergic', 'Non-neuronal']
    n_cls = len(class_order)
    conf = np.zeros((n_cls, n_cls), dtype=int)
    for i, wc in enumerate(class_order):
        for j, nc in enumerate(class_order):
            conf[i, j] = ((df['wc_class'] == wc) & (df['nuc_class'] == nc)).sum()

    im = ax_a.imshow(conf, cmap='Blues', aspect='auto')
    for i in range(n_cls):
        row_total = conf[i].sum()
        for j in range(n_cls):
            pct = 100 * conf[i, j] / row_total if row_total > 0 else 0
            txt = f'{conf[i,j]:,}\n({pct:.1f}%)'
            color = 'white' if conf[i, j] > conf.max() * 0.5 else 'black'
            ax_a.text(j, i, txt, ha='center', va='center', fontsize=11,
                     fontweight='bold', color=color)

    ax_a.set_xticks(range(n_cls))
    ax_a.set_xticklabels([c[:4] for c in class_order], fontsize=13)
    ax_a.set_yticks(range(n_cls))
    ax_a.set_yticklabels([c[:4] for c in class_order], fontsize=13)
    ax_a.set_xlabel('Nuclear Classification', fontsize=14, fontweight='bold')
    ax_a.set_ylabel('Whole-Cell Classification', fontsize=14, fontweight='bold')
    ax_a.set_title('A) Class-Level Concordance', fontsize=16, fontweight='bold')

    # Overall concordance
    total_agree = sum(conf[i, i] for i in range(n_cls))
    total = conf.sum()
    ax_a.text(0.5, -0.18, f'Overall: {total_agree:,}/{total:,} = {100*total_agree/total:.1f}%',
              transform=ax_a.transAxes, ha='center', fontsize=13, fontweight='bold')

    # ── Panel B: Per-subclass concordance ──
    ax_b = fig.add_subplot(gs[1])
    sub_conc = df.groupby('wc_subclass')['concordant_subclass'].agg(['mean', 'count'])
    sub_conc = sub_conc.sort_values('mean', ascending=True)

    y_pos = range(len(sub_conc))
    bar_colors = [CLASS_COLORS.get(SUBCLASS_TO_CLASS.get(s, ''), 'gray')
                  for s in sub_conc.index]

    ax_b.barh(y_pos, sub_conc['mean'] * 100, color=bar_colors, alpha=0.75,
              edgecolor='white', linewidth=0.5)
    ax_b.set_yticks(y_pos)
    ax_b.set_yticklabels(
        [f'{s} (n={int(sub_conc.loc[s, "count"]):,})' for s in sub_conc.index],
        fontsize=10)
    ax_b.set_xlabel('Subclass Concordance (%)', fontsize=14, fontweight='bold')
    ax_b.set_title('B) Per-Subclass Concordance Rate', fontsize=16,
                   fontweight='bold')
    ax_b.set_xlim(0, 105)

    # Annotate percentages
    for i, (s, row) in enumerate(sub_conc.iterrows()):
        ax_b.text(row['mean'] * 100 + 1, i, f'{row["mean"]*100:.1f}%',
                 va='center', fontsize=9, fontweight='bold')

    # Overall subclass concordance
    overall = df['concordant_subclass'].mean()
    ax_b.axvline(overall * 100, color='red', linestyle='--', linewidth=2, alpha=0.7)
    ax_b.text(overall * 100 + 1, len(sub_conc) - 1,
              f'Overall: {overall*100:.1f}%', fontsize=11,
              fontweight='bold', color='red')

    plt.tight_layout()
    path = os.path.join(OUT_DIR, 'nuclear_vs_wholecell_concordance.png')
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def fig4_doublet_rates(df, adata_wc, adata_nuc, qc_mask):
    """Doublet rate comparison: whole-cell vs nuclear."""
    fig = plt.figure(figsize=(16, 7), facecolor='white')
    gs = gridspec.GridSpec(1, 2, wspace=0.3)

    # ── Panel A: Grouped bar chart ──
    ax_a = fig.add_subplot(gs[0])

    wc_gg = (df['wc_doublet_type'] == 'Glut+GABA').sum()
    wc_ga = (df['wc_doublet_type'] == 'GABA+GABA').sum()
    wc_tot = df['wc_doublet'].sum()
    nuc_gg = (df['nuc_doublet_type'] == 'Glut+GABA').sum()
    nuc_ga = (df['nuc_doublet_type'] == 'GABA+GABA').sum()
    nuc_tot = df['nuc_doublet'].sum()

    categories = ['Glut+GABA', 'GABA+GABA', 'Total']
    wc_vals = [wc_gg, wc_ga, wc_tot]
    nuc_vals = [nuc_gg, nuc_ga, nuc_tot]

    x = np.arange(len(categories))
    w = 0.35
    bars1 = ax_a.bar(x - w/2, wc_vals, w, label='Whole-Cell', color='#FF6B6B',
                     alpha=0.8, edgecolor='white')
    bars2 = ax_a.bar(x + w/2, nuc_vals, w, label='Nuclear', color='#44AAFF',
                     alpha=0.8, edgecolor='white')

    # Annotate counts + % reduction
    for i in range(len(categories)):
        ax_a.text(x[i] - w/2, wc_vals[i] + 15, f'{wc_vals[i]:,}',
                 ha='center', va='bottom', fontsize=12, fontweight='bold')
        ax_a.text(x[i] + w/2, nuc_vals[i] + 15, f'{nuc_vals[i]:,}',
                 ha='center', va='bottom', fontsize=12, fontweight='bold',
                 color='#0066CC')
        if wc_vals[i] > 0:
            reduction = 100 * (1 - nuc_vals[i] / wc_vals[i])
            ax_a.text(x[i], max(wc_vals[i], nuc_vals[i]) + 50,
                     f'{reduction:+.0f}%', ha='center', fontsize=11,
                     color='green' if reduction > 0 else 'red',
                     fontweight='bold')

    ax_a.set_xticks(x)
    ax_a.set_xticklabels(categories, fontsize=14)
    ax_a.set_ylabel('Number of Doublets', fontsize=14)
    ax_a.set_title('A) Doublet Detection: Whole-Cell vs Nuclear',
                   fontsize=16, fontweight='bold')
    ax_a.legend(fontsize=13)
    ax_a.tick_params(labelsize=12)

    # ── Panel B: GABA marker score scatter ──
    ax_b = fig.add_subplot(gs[1])

    # Compute GABA marker score for whole-cell and nuclear
    X_wc = adata_wc.X
    X_nuc = adata_nuc.X
    if sparse.issparse(X_wc):
        X_wc = X_wc.toarray()
    if sparse.issparse(X_nuc):
        X_nuc = X_nuc.toarray()

    gene_names = list(adata_wc.var_names)
    gaba_idx = [gene_names.index(g) for g in GABA_MARKERS if g in gene_names]

    wc_gaba_score = (X_wc[:, gaba_idx] > 0).sum(axis=1)
    nuc_gaba_score = (X_nuc[:, gaba_idx] > 0).sum(axis=1)

    # Filter to QC-pass and Glut cells (where Glut+GABA doublets are detected)
    qc_idx = np.where(qc_mask)[0]
    glut_mask = df['wc_class'].values == 'Glutamatergic'
    plot_idx = qc_idx[glut_mask]

    wc_scores = wc_gaba_score[plot_idx]
    nuc_scores = nuc_gaba_score[plot_idx]
    is_wc_doublet = df.loc[glut_mask, 'wc_doublet'].values
    is_nuc_doublet = df.loc[glut_mask, 'nuc_doublet'].values

    # Add jitter for visibility
    jitter = np.random.RandomState(42).uniform(-0.25, 0.25, len(wc_scores))

    # Plot non-doublets (gray, small)
    neither = ~is_wc_doublet & ~is_nuc_doublet
    ax_b.scatter(wc_scores[neither] + jitter[neither],
                 nuc_scores[neither] + jitter[neither],
                 c='#AAAAAA', s=3, alpha=0.15, rasterized=True)

    # WC-only doublets (red circles)
    wc_only = is_wc_doublet & ~is_nuc_doublet
    ax_b.scatter(wc_scores[wc_only] + jitter[wc_only],
                 nuc_scores[wc_only] + jitter[wc_only],
                 c='#FF4444', s=30, alpha=0.7, edgecolors='white',
                 linewidths=0.5, label=f'WC-only doublet (n={wc_only.sum()})',
                 zorder=5)

    # Both doublets (purple)
    both = is_wc_doublet & is_nuc_doublet
    ax_b.scatter(wc_scores[both] + jitter[both],
                 nuc_scores[both] + jitter[both],
                 c='#9944FF', s=30, alpha=0.7, edgecolors='white',
                 linewidths=0.5, label=f'Both doublet (n={both.sum()})',
                 zorder=5)

    # Nuclear-only doublets (blue)
    nuc_only = ~is_wc_doublet & is_nuc_doublet
    if nuc_only.sum() > 0:
        ax_b.scatter(wc_scores[nuc_only] + jitter[nuc_only],
                     nuc_scores[nuc_only] + jitter[nuc_only],
                     c='#44AAFF', s=30, alpha=0.7, edgecolors='white',
                     linewidths=0.5, label=f'Nuclear-only doublet (n={nuc_only.sum()})',
                     zorder=5)

    # Threshold line
    ax_b.axhline(4, color='#44AAFF', linestyle='--', alpha=0.5, linewidth=1)
    ax_b.axvline(4, color='#FF4444', linestyle='--', alpha=0.5, linewidth=1)
    ax_b.text(6.5, 3.5, 'WC threshold', fontsize=9, color='#FF4444', alpha=0.7)
    ax_b.text(3.5, 6.5, 'Nuclear\nthreshold', fontsize=9, color='#44AAFF',
              alpha=0.7, ha='right')

    ax_b.set_xlabel('Whole-Cell GABA Marker Score (# positive)', fontsize=14)
    ax_b.set_ylabel('Nuclear GABA Marker Score (# positive)', fontsize=14)
    ax_b.set_title('B) GABA Marker Scores in Glutamatergic Cells',
                   fontsize=16, fontweight='bold')
    ax_b.legend(fontsize=10, loc='upper left')
    ax_b.tick_params(labelsize=12)
    ax_b.plot([0, 7], [0, 7], 'k--', alpha=0.3, linewidth=1)

    plt.tight_layout()
    path = os.path.join(OUT_DIR, 'nuclear_vs_wholecell_doublet_rates.png')
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def fig5_classification_quality(df):
    """Correlation score comparison: whole-cell vs nuclear."""
    fig = plt.figure(figsize=(16, 7), facecolor='white')
    gs = gridspec.GridSpec(1, 2, wspace=0.3)

    # ── Panel A: Paired scatter ──
    ax_a = fig.add_subplot(gs[0])

    concordant = df['concordant_subclass'].values
    ax_a.scatter(df.loc[concordant, 'wc_corr'],
                 df.loc[concordant, 'nuc_corr'],
                 c='#44AAFF', s=2, alpha=0.08, rasterized=True,
                 label=f'Concordant (n={concordant.sum():,})')
    ax_a.scatter(df.loc[~concordant, 'wc_corr'],
                 df.loc[~concordant, 'nuc_corr'],
                 c='#FF4444', s=5, alpha=0.3, rasterized=True,
                 label=f'Discordant (n={(~concordant).sum():,})')

    ax_a.plot([0, 1], [0, 1], 'k--', alpha=0.4, linewidth=1)
    ax_a.set_xlabel('Whole-Cell Subclass Correlation', fontsize=14)
    ax_a.set_ylabel('Nuclear Subclass Correlation', fontsize=14)
    ax_a.set_title('A) Classification Confidence: WC vs Nuclear',
                   fontsize=16, fontweight='bold')
    ax_a.legend(fontsize=12, loc='lower right')
    ax_a.tick_params(labelsize=12)

    # Median drop
    med_wc = df['wc_corr'].median()
    med_nuc = df['nuc_corr'].median()
    ax_a.text(0.02, 0.95,
              f'Median corr: WC={med_wc:.3f}, Nuc={med_nuc:.3f}\n'
              f'Drop: {med_nuc - med_wc:+.3f}',
              transform=ax_a.transAxes, fontsize=11, va='top',
              bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    # ── Panel B: Correlation distribution by concordance ──
    ax_b = fig.add_subplot(gs[1])

    ax_b.hist(df.loc[concordant, 'nuc_corr'], bins=60, alpha=0.6,
              color='#44AAFF', density=True,
              label=f'Concordant (n={concordant.sum():,})')
    ax_b.hist(df.loc[~concordant, 'nuc_corr'], bins=40, alpha=0.6,
              color='#FF4444', density=True,
              label=f'Discordant (n={(~concordant).sum():,})')

    ax_b.axvline(df.loc[concordant, 'nuc_corr'].median(),
                 color='#0066CC', linestyle='--', linewidth=2)
    ax_b.axvline(df.loc[~concordant, 'nuc_corr'].median(),
                 color='#CC0000', linestyle='--', linewidth=2)

    ax_b.set_xlabel('Nuclear Subclass Correlation', fontsize=14)
    ax_b.set_ylabel('Density', fontsize=16)
    ax_b.set_title('B) Nuclear Correlation by Concordance',
                   fontsize=16, fontweight='bold')
    ax_b.legend(fontsize=12)
    ax_b.tick_params(labelsize=12)

    plt.tight_layout()
    path = os.path.join(OUT_DIR, 'nuclear_classification_quality.png')
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def save_summary_csv(df):
    """Save summary metrics table."""
    metrics = []

    n = len(df)
    metrics.append(('Total cells (QC-pass)', n, n, 0))
    metrics.append(('Mean UMI per cell',
                    df['total_counts'].mean(),
                    df['nuclear_total_counts'].mean(),
                    df['nuclear_total_counts'].mean() - df['total_counts'].mean()))
    metrics.append(('Median nuclear fraction', '-',
                    f"{df['nuclear_fraction'].median():.3f}", '-'))

    cls_conc = df['concordant_class'].mean()
    sub_conc = df['concordant_subclass'].mean()
    metrics.append(('Class concordance (%)', '-', '-', f'{100*cls_conc:.1f}'))
    metrics.append(('Subclass concordance (%)', '-', '-', f'{100*sub_conc:.1f}'))

    wc_gg = (df['wc_doublet_type'] == 'Glut+GABA').sum()
    nuc_gg = (df['nuc_doublet_type'] == 'Glut+GABA').sum()
    metrics.append(('Glut+GABA doublets', wc_gg, nuc_gg, nuc_gg - wc_gg))

    wc_ga = (df['wc_doublet_type'] == 'GABA+GABA').sum()
    nuc_ga = (df['nuc_doublet_type'] == 'GABA+GABA').sum()
    metrics.append(('GABA+GABA doublets', wc_ga, nuc_ga, nuc_ga - wc_ga))

    wc_tot = df['wc_doublet'].sum()
    nuc_tot = df['nuc_doublet'].sum()
    metrics.append(('Total doublets', wc_tot, nuc_tot, nuc_tot - wc_tot))

    metrics.append(('Median subclass correlation',
                    f"{df['wc_corr'].median():.4f}",
                    f"{df['nuc_corr'].median():.4f}",
                    f"{df['nuc_corr'].median() - df['wc_corr'].median():+.4f}"))

    metrics_df = pd.DataFrame(metrics,
                              columns=['Metric', 'Whole-Cell', 'Nuclear', 'Delta'])
    path = os.path.join(OUT_DIR, 'nuclear_vs_wholecell_comparison.csv')
    metrics_df.to_csv(path, index=False)
    print(f"  Saved: {path}")

    # Print summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(metrics_df.to_string(index=False))


def main():
    t_start = time.time()

    # 1. Check nuclear h5ad exists
    if not os.path.exists(H5AD_NUC):
        print(f"ERROR: Nuclear h5ad not found: {H5AD_NUC}")
        print("Run build_nuclear_counts.py first.")
        sys.exit(1)

    # 2. Load centroids
    sub_centroids, sup_centroids, gene_names = load_centroids()

    # 3. Load data
    print("\nLoading Br8667 whole-cell and nuclear h5ads...")
    adata_wc = ad.read_h5ad(H5AD_WC)
    adata_nuc = ad.read_h5ad(H5AD_NUC)
    print(f"  Whole-cell: {adata_wc.n_obs:,} × {adata_wc.n_vars}")
    print(f"  Nuclear:    {adata_nuc.n_obs:,} × {adata_nuc.n_vars}")

    # 4. Classify nuclear counts
    results_nuc, qc_mask = classify_nuclear(
        adata_nuc, sub_centroids, sup_centroids, gene_names)

    # 5. Build comparison
    print("\nBuilding comparison DataFrame...")
    df = build_comparison(adata_wc, adata_nuc, results_nuc, qc_mask)
    print(f"  {len(df):,} cells compared")
    print(f"  Class concordance: {100*df['concordant_class'].mean():.1f}%")
    print(f"  Subclass concordance: {100*df['concordant_subclass'].mean():.1f}%")

    # 6. Generate figures
    print("\nGenerating figures...")
    fig1_nuclear_fraction_histogram(df)
    fig2_nuclear_fraction_by_subclass(df)
    fig3_concordance(df)
    fig4_doublet_rates(df, adata_wc, adata_nuc, qc_mask)
    fig5_classification_quality(df)

    # 7. Summary CSV
    save_summary_csv(df)

    elapsed = time.time() - t_start
    print(f"\nTotal elapsed: {elapsed:.0f}s ({elapsed/60:.1f} min)")


if __name__ == '__main__':
    main()
