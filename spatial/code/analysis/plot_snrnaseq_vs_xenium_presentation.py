#!/usr/bin/env python3
"""
Presentation-quality scatter: snRNAseq meta-analysis beta vs Xenium spatial logFC.

Two-panel layout: Neuronal (left) | Non-neuronal (right).
Generates both supertype-level and subclass-level figures.

Output:
  output/presentation/slide_snrnaseq_vs_xenium.png          (supertype)
  output/presentation/slide_snrnaseq_vs_xenium_subclass.png (subclass)
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import pearsonr
from adjustText import adjust_text

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    BG_COLOR, CRUMBLR_DIR, PRESENTATION_DIR, BASE_DIR,
    GABA_PREFIXES, GLUT_PREFIXES, NN_PREFIXES, infer_class,
)

RESULTS_DIR = CRUMBLR_DIR
OUT_DIR = PRESENTATION_DIR

BG = BG_COLOR

CLASS_COLORS = {'Glut': '#00ADF8', 'GABA': '#F05A28', 'NN': '#808080', 'Other': '#999999'}

# Nicole subclass paths
NICOLE_NEURONAL_SUBCLASS = os.path.join(BASE_DIR, "data", "nicole_scz_snrnaseq_betas",
                                         "final_results_crumblr_7_cohorts_subclass.csv")
NICOLE_NONNEURONAL_SUBCLASS = os.path.join(BASE_DIR, "data", "nicole_scz_snrnaseq_betas",
                                            "final_results_crumblr_7_nonN_cohorts_subclass.csv")


def plot_panel(ax, panel_df, snrna_sig, snrna_nom, title, show_ylabel=True,
               xcol='beta_snrnaseq', ycol='logFC_xenium',
               se_x_col='se', se_y_col='SE_xenium',
               ylabel='Xenium spatial logFC (SCZ vs Control)',
               label_all=False):
    """Plot a single scatter panel with error bars, labels, and regression line."""
    ax.set_facecolor(BG)

    # Error bars first (behind dots)
    for _, row in panel_df.iterrows():
        c = CLASS_COLORS.get(row['class'], '#999999')
        if pd.notna(row.get(se_x_col)):
            ax.plot([row[xcol] - row[se_x_col], row[xcol] + row[se_x_col]],
                    [row[ycol], row[ycol]],
                    color=c, alpha=0.15, linewidth=1.2, zorder=2, solid_capstyle='round')
        if pd.notna(row.get(se_y_col)):
            ax.plot([row[xcol], row[xcol]],
                    [row[ycol] - row[se_y_col], row[ycol] + row[se_y_col]],
                    color=c, alpha=0.15, linewidth=1.2, zorder=2, solid_capstyle='round')

    # Scatter dots
    for cls in ['Glut', 'GABA', 'NN', 'Other']:
        mask = panel_df['class'] == cls
        if mask.sum() == 0:
            continue
        sub = panel_df[mask]
        ax.scatter(sub[xcol], sub[ycol],
                   c=CLASS_COLORS[cls], s=80, alpha=0.8,
                   edgecolors='white', linewidth=0.5, zorder=5,
                   label=f'{cls} (n={mask.sum()})')

    # Regression line
    valid = panel_df[xcol].notna() & panel_df[ycol].notna()
    m = panel_df[valid]
    if len(m) > 2:
        z = np.polyfit(m[xcol], m[ycol], 1)
        lim_x = max(abs(m[xcol]).max(), 0.3) * 1.3
        x_line = np.linspace(-lim_x, lim_x, 100)
        ax.plot(x_line, np.polyval(z, x_line), color='#888888', alpha=0.5,
                linewidth=1.5, linestyle='-', zorder=1)

    # Reference lines
    ax.axhline(0, color='#555555', alpha=0.4, linewidth=0.8, zorder=1)
    ax.axvline(0, color='#555555', alpha=0.4, linewidth=0.8, zorder=1)

    # Labels: two tiers (or label all if few points)
    texts = []
    for _, row in panel_df.iterrows():
        ct = row['celltype']
        if label_all or ct in snrna_sig:
            txt = ax.text(row[xcol], row[ycol],
                          f"  {ct}",
                          fontsize=14, fontweight='bold',
                          color='white', alpha=0.95, zorder=10)
            texts.append(txt)
        elif ct in snrna_nom:
            txt = ax.text(row[xcol], row[ycol],
                          f"  {ct}",
                          fontsize=12, fontweight='normal',
                          color='#bbbbbb', alpha=0.85, zorder=10)
            texts.append(txt)

    if texts:
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle='-', color='#aaaaaa',
                                   alpha=0.5, linewidth=0.7),
                    force_text=(0.9, 0.9),
                    force_points=(0.6, 0.6),
                    expand_text=(1.3, 1.3),
                    expand_points=(1.5, 1.5))

    # Correlation annotation
    if len(m) > 2:
        r_val, p_val = pearsonr(m[xcol], m[ycol])
        ax.text(0.97, 0.04,
                f"r = {r_val:.2f} (p = {p_val:.1e})\nn = {len(m)}",
                transform=ax.transAxes, ha='right', va='bottom',
                fontsize=16, color='#dddddd',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='#222222',
                          edgecolor='#555555', alpha=0.9))

    ax.set_xlabel('snRNAseq meta-analysis beta', fontsize=20, color='white')
    if show_ylabel:
        ax.set_ylabel(ylabel, fontsize=20, color='white')
    ax.set_title(title, fontsize=22, fontweight='bold', color='white', pad=10)

    ax.tick_params(colors='white', labelsize=16)
    for spine in ax.spines.values():
        spine.set_color('#555555')
    ax.grid(True, alpha=0.15, color='#555555')

    # Note about labeled types
    ax.text(0.03, 0.04,
            'Bold = FDR < 0.1  |  Light = nom. p < 0.05',
            transform=ax.transAxes, ha='left', va='bottom',
            fontsize=13, color='#aaaaaa', fontstyle='italic')


def make_two_panel_figure(df, out_path, suptitle):
    """Create a 2-panel (neuronal | non-neuronal) scatter figure."""
    df['class'] = df['celltype'].apply(infer_class)

    # Split into neuronal and non-neuronal
    neuronal = df[df['class'].isin(['Glut', 'GABA'])].copy()
    nonneuronal = df[df['class'] == 'NN'].copy()
    print(f"  Neuronal: {len(neuronal)}, Non-neuronal: {len(nonneuronal)}")

    # Identify significance tiers
    snrna_sig = set(df[df['padj_snrnaseq'] < 0.1]['celltype'].values)
    snrna_nom = set(df[(df['pval_snrnaseq'] < 0.05) &
                       (df['padj_snrnaseq'] >= 0.1)]['celltype'].values)
    print(f"  snRNAseq FDR < 0.1: {len(snrna_sig)} types")

    # Overall correlations
    r_all, p_all = pearsonr(df['beta_snrnaseq'], df['logFC_xenium'])
    if len(neuronal) > 2:
        r_neur, p_neur = pearsonr(neuronal['beta_snrnaseq'], neuronal['logFC_xenium'])
    else:
        r_neur, p_neur = np.nan, np.nan
    if len(nonneuronal) > 2:
        r_nn, p_nn = pearsonr(nonneuronal['beta_snrnaseq'], nonneuronal['logFC_xenium'])
    else:
        r_nn, p_nn = np.nan, np.nan
    print(f"  All: r={r_all:.3f} (p={p_all:.1e})")
    print(f"  Neuronal: r={r_neur:.3f} (p={p_neur:.1e})")
    print(f"  Non-neuronal: r={r_nn:.3f} (p={p_nn:.1e})")

    # Two-panel figure
    fig, (ax_n, ax_nn) = plt.subplots(1, 2, figsize=(22, 10), facecolor=BG)

    # For non-neuronal or if few points, label all
    label_all_nn = len(nonneuronal) <= 20
    label_all_neur = len(neuronal) <= 25

    plot_panel(ax_n, neuronal, snrna_sig, snrna_nom,
              title=f'Neuronal (n={len(neuronal)})',
              show_ylabel=True, label_all=label_all_neur)

    plot_panel(ax_nn, nonneuronal, snrna_sig, snrna_nom,
              title=f'Non-neuronal (n={len(nonneuronal)})',
              show_ylabel=False, label_all=label_all_nn)

    # Legend on left panel
    legend_elements = [
        Line2D([0], [0], marker='o', color=BG, markerfacecolor=CLASS_COLORS['Glut'],
               markersize=12, label=f"Glutamatergic (n={(neuronal['class']=='Glut').sum()})", linewidth=0),
        Line2D([0], [0], marker='o', color=BG, markerfacecolor=CLASS_COLORS['GABA'],
               markersize=12, label=f"GABAergic (n={(neuronal['class']=='GABA').sum()})", linewidth=0),
    ]
    leg = ax_n.legend(handles=legend_elements, loc='upper left', fontsize=15,
                      frameon=True, fancybox=True, framealpha=0.85,
                      edgecolor='#555555', labelcolor='white')
    leg.get_frame().set_facecolor('#222222')

    fig.suptitle(suptitle, fontsize=26, fontweight='bold', color='white', y=1.02)

    plt.tight_layout(pad=2.0)
    plt.savefig(out_path, dpi=200, facecolor=BG, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")


def main():
    # === SUPERTYPE LEVEL ===
    csv_path = os.path.join(RESULTS_DIR, "snrnaseq_vs_xenium_comparison.csv")
    if os.path.exists(csv_path):
        print("=== Supertype-level composition scatter ===")
        df = pd.read_csv(csv_path)
        print(f"  Loaded {len(df)} cell types")
        outpath = os.path.join(OUT_DIR, "slide_snrnaseq_vs_xenium.png")
        make_two_panel_figure(df, outpath,
                              'SCZ compositional effects: snRNAseq vs Xenium (supertype)')
    else:
        print(f"  Supertype comparison CSV not found: {csv_path}")

    # === SUBCLASS LEVEL ===
    # Build subclass comparison from stratified Xenium + Nicole subclass data
    subclass_csv = os.path.join(RESULTS_DIR, "snrnaseq_vs_xenium_comparison_subclass.csv")
    if os.path.exists(subclass_csv):
        print("\n=== Subclass-level composition scatter ===")
        df_sub = pd.read_csv(subclass_csv)
        print(f"  Loaded {len(df_sub)} subclasses")
        outpath_sub = os.path.join(OUT_DIR, "slide_snrnaseq_vs_xenium_subclass.png")
        make_two_panel_figure(df_sub, outpath_sub,
                              'SCZ compositional effects: snRNAseq vs Xenium (subclass)')
    else:
        print(f"\n  Subclass comparison CSV not found: {subclass_csv}")
        print("  Run plot_crumblr_results.py first to generate it.")


if __name__ == '__main__':
    main()
