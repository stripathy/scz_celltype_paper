#!/usr/bin/env python3
"""
Regenerate the three key presentation plots using the standard pipeline QC
filters (corr_qc_pass: 5th-percentile margin + L6b threshold + doublet exclusion).

Plots:
  1. snRNAseq meta-analysis beta vs Xenium logFC (supertype level)
  2. MERFISH vs Xenium subclass proportions (controls only)
  3. MERFISH vs Xenium median depth by cell type (subclass + supertype)

Output: output/presentation/slide_*.png
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
from matplotlib.lines import Line2D
from scipy.stats import pearsonr
from adjustText import adjust_text

BASE_DIR = os.path.expanduser("~/Github/SCZ_Xenium")
sys.path.insert(0, os.path.join(BASE_DIR, "code", "analysis"))
from config import (
    BG_COLOR, H5AD_DIR, MERFISH_PATH, PRESENTATION_DIR, CRUMBLR_DIR,
    SUBCLASS_TO_CLASS, CLASS_COLORS, EXCLUDE_SAMPLES, infer_class,
    load_cells,
    load_merfish_cortical as _load_merfish_cortical,
)

OUT_DIR = PRESENTATION_DIR
BG = BG_COLOR
SCATTER_COLORS = {'Glut': '#00ADF8', 'GABA': '#F05A28', 'NN': '#808080', 'Other': '#999999'}


# ══════════════════════════════════════════════════════════════════════
# Helpers
# ══════════════════════════════════════════════════════════════════════

def load_xenium_cortical(extra_cols=None):
    """Load all Xenium cortical cells using standard pipeline QC (corr_qc_pass)."""
    h5ad_files = sorted(glob.glob(os.path.join(H5AD_DIR, "*_annotated.h5ad")))
    all_obs = []
    for fpath in h5ad_files:
        sid = os.path.basename(fpath).replace("_annotated.h5ad", "")
        if sid in EXCLUDE_SAMPLES:
            continue
        obs = load_cells(sid, cortical_only=True, extra_obs_columns=extra_cols)
        all_obs.append(obs)
        print(f"  {sid}: {len(obs):,} cells")

    combined = pd.concat(all_obs, ignore_index=True)
    print(f"  Total: {len(combined):,} cells")
    return combined


# ══════════════════════════════════════════════════════════════════════
# Plot 1: snRNAseq vs Xenium logFC
# ══════════════════════════════════════════════════════════════════════

def plot_snrnaseq_scatter():
    """Dark-background scatter: snRNAseq beta vs Xenium logFC at supertype level."""
    print("\n=== Plot 1: snRNAseq vs Xenium logFC ===")
    csv_path = os.path.join(CRUMBLR_DIR, "snrnaseq_vs_xenium_comparison.csv")
    df = pd.read_csv(csv_path)
    df['class'] = df['celltype'].apply(infer_class)
    print(f"  {len(df)} shared supertypes")

    snrna_sig = set(df[df['padj_snrnaseq'] < 0.1]['celltype'].values)
    snrna_nom = set(df[(df['pval_snrnaseq'] < 0.05) & (df['padj_snrnaseq'] >= 0.1)]['celltype'].values)

    r_all, p_all = pearsonr(df['beta_snrnaseq'], df['logFC_xenium'])
    neur = df[df['class'].isin(['Glut', 'GABA'])]
    r_neur, p_neur = pearsonr(neur['beta_snrnaseq'], neur['logFC_xenium'])
    print(f"  All: r={r_all:.3f} (p={p_all:.1e}), Neuronal: r={r_neur:.3f} (p={p_neur:.1e})")

    fig, ax = plt.subplots(figsize=(11, 9), facecolor=BG)
    ax.set_facecolor(BG)

    # Error bars
    for _, row in df.iterrows():
        c = SCATTER_COLORS[row['class']]
        ax.plot([row['beta_snrnaseq'] - row['se'], row['beta_snrnaseq'] + row['se']],
                [row['logFC_xenium'], row['logFC_xenium']],
                color=c, alpha=0.15, linewidth=1.2, zorder=2, solid_capstyle='round')
        if pd.notna(row.get('SE_xenium')):
            ax.plot([row['beta_snrnaseq'], row['beta_snrnaseq']],
                    [row['logFC_xenium'] - row['SE_xenium'],
                     row['logFC_xenium'] + row['SE_xenium']],
                    color=c, alpha=0.15, linewidth=1.2, zorder=2, solid_capstyle='round')

    # Scatter
    for cls in ['Glut', 'GABA', 'NN', 'Other']:
        mask = df['class'] == cls
        if mask.sum() == 0:
            continue
        sub = df[mask]
        ax.scatter(sub['beta_snrnaseq'], sub['logFC_xenium'],
                   c=SCATTER_COLORS[cls], s=70, alpha=0.8,
                   edgecolors='white', linewidth=0.5, zorder=5,
                   label=f'{cls} (n={mask.sum()})')

    # Regression line
    z = np.polyfit(df['beta_snrnaseq'], df['logFC_xenium'], 1)
    lim_x = max(abs(df['beta_snrnaseq']).max(), 0.3) * 1.3
    x_line = np.linspace(-lim_x, lim_x, 100)
    ax.plot(x_line, np.polyval(z, x_line), color='#888888', alpha=0.5,
            linewidth=1.5, linestyle='-', zorder=1)

    ax.axhline(0, color='#555555', alpha=0.4, linewidth=0.8, zorder=1)
    ax.axvline(0, color='#555555', alpha=0.4, linewidth=0.8, zorder=1)

    # Labels
    texts = []
    for _, row in df.iterrows():
        ct = row['celltype']
        if ct in snrna_sig:
            txt = ax.text(row['beta_snrnaseq'], row['logFC_xenium'], f"  {ct}",
                          fontsize=13, fontweight='bold', color='white', alpha=0.95, zorder=10)
            texts.append(txt)
        elif ct in snrna_nom:
            txt = ax.text(row['beta_snrnaseq'], row['logFC_xenium'], f"  {ct}",
                          fontsize=11, fontweight='normal', color='#bbbbbb', alpha=0.85, zorder=10)
            texts.append(txt)
    if texts:
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle='-', color='#aaaaaa', alpha=0.5, linewidth=0.7),
                    force_text=(0.9, 0.9), force_points=(0.6, 0.6),
                    expand_text=(1.4, 1.4), expand_points=(1.6, 1.6))

    ax.text(0.97, 0.04,
            f"All: r = {r_all:.2f} (p = {p_all:.1e})\n"
            f"Neuronal: r = {r_neur:.2f} (p = {p_neur:.1e})\n"
            f"n = {len(df)} shared supertypes",
            transform=ax.transAxes, ha='right', va='bottom',
            fontsize=14, color='#dddddd',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='#222222',
                      edgecolor='#555555', alpha=0.9))

    ax.set_xlabel('snRNAseq meta-analysis beta (SCZ effect)', fontsize=18, color='white')
    ax.set_ylabel('Xenium spatial logFC (SCZ vs Control)', fontsize=18, color='white')
    ax.set_title('SCZ compositional effects: snRNAseq vs Xenium spatial',
                 fontsize=22, fontweight='bold', color='white', pad=12)
    ax.tick_params(colors='white', labelsize=14)
    for spine in ax.spines.values():
        spine.set_color('#555555')
    ax.grid(True, alpha=0.15, color='#555555')

    legend_elements = [
        Line2D([0], [0], marker='o', color=BG, markerfacecolor=SCATTER_COLORS['Glut'],
               markersize=10, label=f"Glutamatergic (n={(df['class']=='Glut').sum()})", linewidth=0),
        Line2D([0], [0], marker='o', color=BG, markerfacecolor=SCATTER_COLORS['GABA'],
               markersize=10, label=f"GABAergic (n={(df['class']=='GABA').sum()})", linewidth=0),
        Line2D([0], [0], marker='o', color=BG, markerfacecolor=SCATTER_COLORS['NN'],
               markersize=10, label=f"Non-neuronal (n={(df['class']=='NN').sum()})", linewidth=0),
    ]
    ax.legend(handles=legend_elements, fontsize=13, loc='upper left',
              facecolor='#1a1a1a', edgecolor='#555555', labelcolor='white')

    plt.tight_layout()
    outpath = os.path.join(OUT_DIR, 'slide_snrnaseq_vs_xenium.png')
    plt.savefig(outpath, dpi=150, facecolor=BG)
    plt.close()
    print(f"  Saved: {outpath}")


# ══════════════════════════════════════════════════════════════════════
# Plot 2: Proportion scatter (MERFISH vs Xenium)
# ══════════════════════════════════════════════════════════════════════

def plot_proportion_scatter():
    """Subclass proportion scatter: MERFISH vs Xenium (controls only)."""
    print("\n=== Plot 2: Proportion scatter ===")

    # MERFISH reference proportions (all donors)
    print("  Loading MERFISH proportions...")
    merfish = ad.read_h5ad(MERFISH_PATH, backed='r')
    obs_m = merfish.obs[['Donor ID', 'Subclass']].copy()
    obs_m.columns = ['donor', 'celltype']
    obs_m['celltype'] = obs_m['celltype'].astype(str)
    merfish.file.close()

    m_counts = obs_m.groupby(['donor', 'celltype']).size().reset_index(name='count')
    m_totals = obs_m.groupby('donor').size().reset_index(name='total')
    m_counts = m_counts.merge(m_totals, on='donor')
    m_counts['prop'] = m_counts['count'] / m_counts['total']
    merfish_props = m_counts.groupby('celltype')['prop'].median().reset_index(name='merfish_prop')

    # Xenium proportions from crumblr input (default = standard pipeline QC)
    xen_path = os.path.join(CRUMBLR_DIR, 'crumblr_input_subclass.csv')
    xen = pd.read_csv(xen_path)
    xen['prop'] = xen['count'] / xen['total']
    xenium_props = xen.groupby('celltype')['prop'].median().reset_index(name='xenium_prop')

    merged = pd.merge(merfish_props, xenium_props, on='celltype', how='inner')
    merged['class'] = merged['celltype'].map(SUBCLASS_TO_CLASS).fillna('Non-neuronal')
    print(f"  {len(merged)} shared subclasses")

    r, p = pearsonr(merged['merfish_prop'], merged['xenium_prop'])
    print(f"  Pearson r={r:.3f} (p={p:.1e})")

    fig, ax = plt.subplots(figsize=(10, 9), facecolor=BG)
    ax.set_facecolor(BG)

    for cls in ['Glutamatergic', 'GABAergic', 'Non-neuronal']:
        mask = merged['class'] == cls
        if mask.sum() == 0:
            continue
        sub = merged[mask]
        short = cls[:4] if cls != 'Non-neuronal' else 'NN'
        c = CLASS_COLORS.get(cls, '#888888')
        ax.scatter(sub['merfish_prop'], sub['xenium_prop'], c=c, s=80, alpha=0.8,
                   edgecolors='white', linewidth=0.5, zorder=5,
                   label=f'{cls} (n={mask.sum()})')

    # Identity line
    lims = [0, max(merged['merfish_prop'].max(), merged['xenium_prop'].max()) * 1.1]
    ax.plot(lims, lims, '--', color='#666666', linewidth=1, zorder=1)

    # Labels
    texts = []
    for _, row in merged.iterrows():
        texts.append(ax.text(row['merfish_prop'], row['xenium_prop'],
                             row['celltype'], fontsize=11, color='#dddddd'))
    if texts:
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle='-', color='#555555', lw=0.5),
                    force_text=(0.3, 0.3))

    ax.text(0.03, 0.97, f"r = {r:.3f}\nn = {len(merged)}",
            transform=ax.transAxes, ha='left', va='top', fontsize=14, color='white',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#1a1a1a',
                      edgecolor='#444444', alpha=0.85))

    ax.set_xlabel('MERFISH median proportion', fontsize=18, color='white')
    ax.set_ylabel('Xenium median proportion', fontsize=18, color='white')
    ax.set_title('Subclass proportions: MERFISH vs Xenium',
                 fontsize=20, fontweight='bold', color='white', pad=12)
    ax.tick_params(colors='white', labelsize=14)
    for spine in ax.spines.values():
        spine.set_color('#555555')
    ax.legend(fontsize=12, loc='lower right', facecolor='#1a1a1a',
              edgecolor='#555555', labelcolor='white')

    plt.tight_layout()
    outpath = os.path.join(OUT_DIR, 'slide_proportion_scatter.png')
    plt.savefig(outpath, dpi=150, facecolor=BG)
    plt.close()
    print(f"  Saved: {outpath}")


# ══════════════════════════════════════════════════════════════════════
# Plot 3: Depth scatter (MERFISH vs Xenium)
# ══════════════════════════════════════════════════════════════════════

def plot_depth_scatter():
    """Median depth scatter: MERFISH vs Xenium (subclass + supertype)."""
    print("\n=== Plot 3: Depth scatter ===")

    # MERFISH reference depth (manual annotations, cortical only)
    print("  Loading MERFISH cortical depth...")
    merfish_df = _load_merfish_cortical()
    merfish_df = merfish_df.rename(columns={'predicted_norm_depth': 'depth'})
    # Use manual depth where available
    if 'manual_depth' in merfish_df.columns:
        has_manual = ~np.isnan(merfish_df['manual_depth'])
        merfish_df.loc[has_manual, 'depth'] = merfish_df.loc[has_manual, 'manual_depth']

    merfish_sub = merfish_df.groupby('subclass')['depth'].agg(['median', 'count']).reset_index()
    merfish_sub.columns = ['celltype', 'median_depth', 'n_cells']
    merfish_sub = merfish_sub[merfish_sub['n_cells'] >= 20]

    merfish_sup = merfish_df.groupby('supertype')['depth'].agg(['median', 'count']).reset_index()
    merfish_sup.columns = ['celltype', 'median_depth', 'n_cells']
    merfish_sup = merfish_sup[merfish_sup['n_cells'] >= 20]

    # Xenium (standard pipeline QC)
    print("  Loading Xenium cortical cells...")
    xen_df = load_xenium_cortical(extra_cols=['predicted_norm_depth'])
    xen_df['depth'] = xen_df['predicted_norm_depth'].astype(float)

    xen_sub = xen_df.groupby('subclass_label')['depth'].agg(['median', 'count']).reset_index()
    xen_sub.columns = ['celltype', 'median_depth', 'n_cells']
    xen_sub = xen_sub[xen_sub['n_cells'] >= 20]

    xen_sup = xen_df.groupby('supertype_label')['depth'].agg(['median', 'count']).reset_index()
    xen_sup.columns = ['celltype', 'median_depth', 'n_cells']
    xen_sup = xen_sup[xen_sup['n_cells'] >= 20]

    print(f"  MERFISH subclasses: {len(merfish_sub)}, supertypes: {len(merfish_sup)}")
    print(f"  Xenium subclasses: {len(xen_sub)}, supertypes: {len(xen_sup)}")

    fig, axes = plt.subplots(1, 2, figsize=(18, 8), facecolor=BG)

    def _scatter(ax, mdf, xdf, level, label_all=True, label_thresh=0.1):
        merged = mdf.merge(xdf, on='celltype', suffixes=('_merfish', '_xenium'))
        if len(merged) == 0:
            ax.text(0.5, 0.5, 'No shared cell types', transform=ax.transAxes,
                    ha='center', color='white', fontsize=14)
            return

        x = merged['median_depth_merfish'].values
        y = merged['median_depth_xenium'].values
        names = merged['celltype'].values

        if level == 'subclass':
            colors = [CLASS_COLORS.get(SUBCLASS_TO_CLASS.get(n, 'Non-neuronal'), '#888888')
                      for n in names]
        else:
            # Infer from prefix
            colors = []
            for n in names:
                matched = False
                for sub, cls in SUBCLASS_TO_CLASS.items():
                    if n.startswith(sub.split()[0]):
                        colors.append(CLASS_COLORS.get(cls, '#888888'))
                        matched = True
                        break
                if not matched:
                    colors.append('#888888')

        min_n = np.minimum(merged['n_cells_merfish'].values, merged['n_cells_xenium'].values)
        sizes = np.clip(np.log10(min_n) * 25, 15, 120)

        ax.scatter(x, y, c=colors, s=sizes, alpha=0.8, edgecolors='white',
                   linewidths=0.5, zorder=3)

        lims = [min(x.min(), y.min()) - 0.02, max(x.max(), y.max()) + 0.02]
        ax.plot(lims, lims, '--', color='#666666', linewidth=1, zorder=1)

        r, p = pearsonr(x, y)
        ax.text(0.03, 0.97, f"r = {r:.3f}\nn = {len(merged)}",
                transform=ax.transAxes, ha='left', va='top', fontsize=13, color='white',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#1a1a1a',
                          edgecolor='#444444', alpha=0.85))

        texts = []
        for i, name in enumerate(names):
            if level == 'supertype' and not label_all:
                if abs(x[i] - y[i]) < label_thresh:
                    continue
            fs = 8 if level == 'supertype' else 11
            texts.append(ax.text(x[i], y[i], name, fontsize=fs, color='#dddddd'))
        if texts:
            adjust_text(texts, ax=ax,
                        arrowprops=dict(arrowstyle='-', color='#555555', lw=0.5),
                        force_text=(0.3, 0.3))

        ax.set_xlabel('MERFISH median depth from pia', fontsize=14, color='white')
        ax.set_ylabel('Xenium median depth from pia', fontsize=14, color='white')
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.tick_params(colors='white', labelsize=12)
        ax.set_facecolor(BG)
        for spine in ax.spines.values():
            spine.set_color('#555555')

    _scatter(axes[0], merfish_sub, xen_sub, 'subclass')
    axes[0].set_title('Subclass level', fontsize=18, fontweight='bold', color='white', pad=10)

    _scatter(axes[1], merfish_sup, xen_sup, 'supertype', label_all=False, label_thresh=0.08)
    axes[1].set_title('Supertype level', fontsize=18, fontweight='bold', color='white', pad=10)

    # Class legend
    legend_elements = [
        Line2D([0], [0], marker='o', color=BG, markerfacecolor=CLASS_COLORS.get('Glutamatergic', '#00ADF8'),
               markersize=10, label='Glutamatergic', linewidth=0),
        Line2D([0], [0], marker='o', color=BG, markerfacecolor=CLASS_COLORS.get('GABAergic', '#F05A28'),
               markersize=10, label='GABAergic', linewidth=0),
        Line2D([0], [0], marker='o', color=BG, markerfacecolor=CLASS_COLORS.get('Non-neuronal', '#808080'),
               markersize=10, label='Non-neuronal', linewidth=0),
    ]
    axes[1].legend(handles=legend_elements, fontsize=12, loc='lower right',
                   facecolor='#1a1a1a', edgecolor='#555555', labelcolor='white')

    fig.suptitle('Median depth: MERFISH vs Xenium',
                 fontsize=22, fontweight='bold', color='white', y=1.02)
    plt.tight_layout()
    outpath = os.path.join(OUT_DIR, 'slide_median_depth_by_celltype.png')
    plt.savefig(outpath, dpi=150, facecolor=BG, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outpath}")


# ══════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    plot_snrnaseq_scatter()
    plot_proportion_scatter()
    plot_depth_scatter()
    print("\nAll three presentation plots regenerated.")
