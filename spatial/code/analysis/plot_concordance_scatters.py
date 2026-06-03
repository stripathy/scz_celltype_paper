#!/usr/bin/env python3
"""
Generate all cross-platform concordance scatter plots (snRNAseq vs Xenium).

Produces 4 two-panel figures (neuronal | non-neuronal):
  1. Subclass composition    → output/presentation/concordance_subclass_composition.png
  2. Subclass density        → output/presentation/concordance_subclass_density.png
  3. Supertype composition   → output/presentation/concordance_supertype_composition.png
  4. Supertype density       → output/presentation/concordance_supertype_density.png

Each figure shows neuronal types (left) and non-neuronal types (right) as
separate panels with independent correlation statistics.

Inputs:
  - Xenium stratified crumblr results (composition)
  - Xenium density results (aggregated to subclass for subclass-level)
  - Nicole's snRNAseq meta-analysis results (supertype + subclass, neuronal + non-neuronal)

Usage:
    cd ~/Github/SCZ_Xenium && python3 -u code/analysis/plot_concordance_scatters.py
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import pearsonr, ttest_ind
from adjustText import adjust_text

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    BASE_DIR, CRUMBLR_DIR, PRESENTATION_DIR, EXCLUDE_SAMPLES,
    BG_COLOR, infer_class,
)

BG = BG_COLOR
CLASS_COLORS = {'Glut': '#00ADF8', 'GABA': '#F05A28', 'NN': '#808080', 'Other': '#999999'}

# ── Data paths ──
NICOLE_DIR = os.path.join(BASE_DIR, "data", "nicole_scz_snrnaseq_betas")
DENSITY_DIR = os.path.join(BASE_DIR, "output", "density_analysis")

# Xenium density uses abbreviated subclass names; Nicole uses full names.
SUBCLASS_NAME_MAP = {
    "Astro": "Astrocyte",
    "Endo": "Endothelial",
    "Micro-PVM": "Microglia-PVM",
    "Oligo": "Oligodendrocyte",
}


# ──────────────────────────────────────────────────────────────────────
# Data loading
# ──────────────────────────────────────────────────────────────────────

def load_nicole(level='supertype'):
    """Load Nicole's stratified snRNAseq meta-analysis results.

    Returns DataFrame with: celltype, beta_snrnaseq, se, pval_snrnaseq, padj_snrnaseq
    """
    suffix = '_subclass' if level == 'subclass' else ''
    neuronal_path = os.path.join(NICOLE_DIR, f"final_results_crumblr_7_cohorts{suffix}.csv")
    nonneuronal_path = os.path.join(NICOLE_DIR, f"final_results_crumblr_7_nonN_cohorts{suffix}.csv")

    dfs = []
    for path in [neuronal_path, nonneuronal_path]:
        if not os.path.exists(path):
            print(f"  WARNING: {path} not found")
            continue
        df = pd.read_csv(path)
        df = df[~df['CellType'].str.contains('SEAAD', na=False)]
        dfs.append(df)

    if not dfs:
        return None

    combined = pd.concat(dfs, ignore_index=True)
    combined = combined.rename(columns={
        'CellType': 'celltype', 'estimate': 'beta_snrnaseq',
        'pval': 'pval_snrnaseq', 'padj': 'padj_snrnaseq',
    })
    return combined


def load_xenium_composition(level='supertype'):
    """Load stratified Xenium crumblr results (neuronal + non-neuronal)."""
    n_path = os.path.join(CRUMBLR_DIR, f'crumblr_results_{level}_neuronal.csv')
    nn_path = os.path.join(CRUMBLR_DIR, f'crumblr_results_{level}_nonneuronal.csv')

    if not (os.path.exists(n_path) and os.path.exists(nn_path)):
        print(f"  Stratified crumblr results not found for {level}")
        return None

    xen = pd.concat([pd.read_csv(n_path), pd.read_csv(nn_path)], ignore_index=True)
    xen = xen.rename(columns={
        'logFC': 'logFC_xenium', 'P.Value': 'pval_xenium', 'FDR': 'FDR_xenium',
    })
    return xen


def supertype_to_subclass(st):
    """Map supertype name to subclass by stripping trailing _N."""
    parts = st.rsplit("_", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return st


def load_xenium_density_subclass():
    """Load Xenium density results aggregated to subclass level."""
    raw_path = os.path.join(DENSITY_DIR, "density_per_sample_supertype.csv")
    if not os.path.exists(raw_path):
        return None

    raw = pd.read_csv(raw_path)
    raw = raw[raw["region"] == "cortical"]
    raw = raw[~raw["sample_id"].isin(EXCLUDE_SAMPLES)]
    raw["subclass"] = raw["supertype"].apply(supertype_to_subclass)

    sub_density = (raw.groupby(["sample_id", "diagnosis", "subclass"])
                   .agg(density_per_mm2=("density_per_mm2", "sum"),
                        count=("count", "sum"),
                        area_mm2=("area_mm2", "first"))
                   .reset_index())

    results = []
    for sc in sorted(sub_density["subclass"].unique()):
        sc_data = sub_density[sub_density["subclass"] == sc]
        ctrl = sc_data[sc_data["diagnosis"] == "Control"]["density_per_mm2"]
        scz = sc_data[sc_data["diagnosis"] == "SCZ"]["density_per_mm2"]
        if len(ctrl) < 3 or len(scz) < 3:
            continue
        ctrl_mean, scz_mean = ctrl.mean(), scz.mean()
        logFC = np.log(scz_mean / ctrl_mean) if ctrl_mean > 0 else np.nan
        _, pval = ttest_ind(np.log1p(ctrl), np.log1p(scz), equal_var=False)
        se = np.log1p(sc_data["density_per_mm2"]).std() / np.sqrt(len(sc_data))
        results.append({
            "celltype": sc, "logFC_density": logFC, "pval_density": pval, "se_density": se,
        })

    df = pd.DataFrame(results)
    df["celltype"] = df["celltype"].replace(SUBCLASS_NAME_MAP)
    df = df[~df["celltype"].str.contains("SEAAD", na=False)]
    return df


def merge_composition(level='supertype'):
    """Merge Nicole snRNAseq + Xenium composition for a given level."""
    nicole = load_nicole(level)
    xen = load_xenium_composition(level)
    if nicole is None or xen is None:
        return None

    se_col = 'SE' if 'SE' in xen.columns else None
    merge_cols = ['celltype', 'logFC_xenium', 'pval_xenium', 'FDR_xenium']
    if se_col:
        merge_cols.append(se_col)

    merged = pd.merge(
        nicole[['celltype', 'beta_snrnaseq', 'se', 'pval_snrnaseq', 'padj_snrnaseq']],
        xen[merge_cols], on='celltype', how='inner'
    )
    if se_col:
        merged = merged.rename(columns={se_col: 'SE_xenium'})
    merged['class'] = merged['celltype'].apply(infer_class)
    return merged


def merge_density(level='supertype'):
    """Merge Nicole snRNAseq + Xenium density for a given level."""
    nicole = load_nicole(level)
    if nicole is None:
        return None

    if level == 'subclass':
        density = load_xenium_density_subclass()
        nicole_renamed = nicole.rename(columns={'se': 'se_snrnaseq'})
    else:
        density_path = os.path.join(DENSITY_DIR, "density_results_supertype_cortical.csv")
        if not os.path.exists(density_path):
            return None
        density = pd.read_csv(density_path)
        density = density.rename(columns={
            "supertype": "celltype", "logFC": "logFC_density",
            "pval": "pval_density", "fdr": "fdr_density", "se": "se_density",
        })
        nicole_renamed = nicole.rename(columns={'se': 'se_snrnaseq'})

    if density is None:
        return None

    merged = pd.merge(
        nicole_renamed[['celltype', 'beta_snrnaseq', 'se_snrnaseq', 'pval_snrnaseq', 'padj_snrnaseq']],
        density, on='celltype', how='inner'
    )
    merged['class'] = merged['celltype'].apply(infer_class)
    return merged


# ──────────────────────────────────────────────────────────────────────
# Plotting
# ──────────────────────────────────────────────────────────────────────

def plot_panel(ax, panel_df, xcol, ycol, se_x_col, se_y_col,
               snrna_sig, snrna_nom, title, ylabel, show_ylabel=True,
               label_all=False):
    """Plot a single scatter panel with error bars, labels, and regression line."""
    ax.set_facecolor(BG)

    # Error bars
    for _, row in panel_df.iterrows():
        c = CLASS_COLORS.get(row['class'], '#999999')
        if se_x_col and pd.notna(row.get(se_x_col)):
            ax.plot([row[xcol] - row[se_x_col], row[xcol] + row[se_x_col]],
                    [row[ycol], row[ycol]],
                    color=c, alpha=0.15, linewidth=1.2, zorder=2, solid_capstyle='round')
        if se_y_col and pd.notna(row.get(se_y_col)):
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
                   edgecolors='white', linewidth=0.5, zorder=5)

    # Regression line
    valid = panel_df[xcol].notna() & panel_df[ycol].notna()
    m = panel_df[valid]
    if len(m) > 2:
        z = np.polyfit(m[xcol], m[ycol], 1)
        lim_x = max(abs(m[xcol]).max(), 0.3) * 1.3
        x_line = np.linspace(-lim_x, lim_x, 100)
        ax.plot(x_line, np.polyval(z, x_line), color='#888888', alpha=0.5,
                linewidth=1.5, linestyle='-', zorder=1)

    ax.axhline(0, color='#555555', alpha=0.4, linewidth=0.8, zorder=1)
    ax.axvline(0, color='#555555', alpha=0.4, linewidth=0.8, zorder=1)

    # Labels
    texts = []
    for _, row in panel_df.iterrows():
        ct = row['celltype']
        if label_all or ct in snrna_sig:
            txt = ax.text(row[xcol], row[ycol], f"  {ct}",
                          fontsize=14, fontweight='bold',
                          color='white', alpha=0.95, zorder=10)
            texts.append(txt)
        elif ct in snrna_nom:
            txt = ax.text(row[xcol], row[ycol], f"  {ct}",
                          fontsize=12, fontweight='normal',
                          color='#bbbbbb', alpha=0.85, zorder=10)
            texts.append(txt)

    if texts:
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle='-', color='#aaaaaa',
                                   alpha=0.5, linewidth=0.7),
                    force_text=(0.9, 0.9), force_points=(0.6, 0.6),
                    expand_text=(1.3, 1.3), expand_points=(1.5, 1.5))

    # Correlation box
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

    # Significance note
    ax.text(0.03, 0.04,
            'Bold = FDR < 0.1  |  Light = nom. p < 0.05',
            transform=ax.transAxes, ha='left', va='bottom',
            fontsize=13, color='#aaaaaa', fontstyle='italic')


def make_concordance_figure(df, xcol, ycol, se_x_col, se_y_col,
                             ylabel, suptitle, out_path, sig_col='padj_snrnaseq'):
    """Create a 2-panel (neuronal | non-neuronal) concordance scatter."""
    neuronal = df[df['class'].isin(['Glut', 'GABA'])].copy()
    nonneuronal = df[df['class'] == 'NN'].copy()

    # Significance tiers
    snrna_sig = set(df[df[sig_col] < 0.1]['celltype'].values)
    snrna_nom = set(df[(df['pval_snrnaseq'] < 0.05) & (df[sig_col] >= 0.1)]['celltype'].values)

    # Correlations
    def safe_corr(sub, col_x, col_y):
        valid = sub[col_x].notna() & sub[col_y].notna()
        m = sub[valid]
        if len(m) > 2:
            return pearsonr(m[col_x], m[col_y])
        return (np.nan, np.nan)

    r_n, p_n = safe_corr(neuronal, xcol, ycol)
    r_nn, p_nn = safe_corr(nonneuronal, xcol, ycol)
    r_all, p_all = safe_corr(df, xcol, ycol)
    print(f"  All: r={r_all:.3f} (p={p_all:.1e}, n={len(df)})")
    print(f"  Neuronal: r={r_n:.3f} (p={p_n:.1e}, n={len(neuronal)})")
    print(f"  Non-neuronal: r={r_nn:.3f} (p={p_nn:.1e}, n={len(nonneuronal)})")

    # Figure
    fig, (ax_n, ax_nn) = plt.subplots(1, 2, figsize=(22, 10), facecolor=BG)

    label_all_n = len(neuronal) <= 25
    label_all_nn = len(nonneuronal) <= 20

    plot_panel(ax_n, neuronal, xcol, ycol, se_x_col, se_y_col,
              snrna_sig, snrna_nom,
              title=f'Neuronal (n={len(neuronal)})', ylabel=ylabel,
              show_ylabel=True, label_all=label_all_n)

    plot_panel(ax_nn, nonneuronal, xcol, ycol, se_x_col, se_y_col,
              snrna_sig, snrna_nom,
              title=f'Non-neuronal (n={len(nonneuronal)})', ylabel=ylabel,
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


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────

def main():
    os.makedirs(PRESENTATION_DIR, exist_ok=True)

    # ── 1. Subclass composition ──
    print("\n=== Subclass composition ===")
    sub_comp = merge_composition('subclass')
    if sub_comp is not None:
        print(f"  {len(sub_comp)} shared subclasses")
        make_concordance_figure(
            sub_comp, 'beta_snrnaseq', 'logFC_xenium', 'se', 'SE_xenium',
            ylabel='Xenium compositional logFC',
            suptitle='SCZ compositional effects: snRNAseq vs Xenium (subclass)',
            out_path=os.path.join(PRESENTATION_DIR, 'concordance_subclass_composition.png'),
        )
        sub_comp.to_csv(os.path.join(CRUMBLR_DIR, 'snrnaseq_vs_xenium_comparison_subclass.csv'), index=False)

    # ── 2. Subclass density ──
    print("\n=== Subclass density ===")
    sub_dens = merge_density('subclass')
    if sub_dens is not None:
        print(f"  {len(sub_dens)} shared subclasses")
        make_concordance_figure(
            sub_dens, 'beta_snrnaseq', 'logFC_density', 'se_snrnaseq', 'se_density',
            ylabel='Xenium density logFC (cells/mm²)',
            suptitle='SCZ density effects: snRNAseq vs Xenium (subclass)',
            out_path=os.path.join(PRESENTATION_DIR, 'concordance_subclass_density.png'),
        )
        sub_dens.to_csv(os.path.join(DENSITY_DIR, 'snrnaseq_vs_density_subclass.csv'), index=False)

    # ── 3. Supertype composition ──
    print("\n=== Supertype composition ===")
    sup_comp = merge_composition('supertype')
    if sup_comp is not None:
        print(f"  {len(sup_comp)} shared supertypes")
        make_concordance_figure(
            sup_comp, 'beta_snrnaseq', 'logFC_xenium', 'se', 'SE_xenium',
            ylabel='Xenium compositional logFC',
            suptitle='SCZ compositional effects: snRNAseq vs Xenium (supertype)',
            out_path=os.path.join(PRESENTATION_DIR, 'concordance_supertype_composition.png'),
        )
        sup_comp.to_csv(os.path.join(CRUMBLR_DIR, 'snrnaseq_vs_xenium_comparison.csv'), index=False)

    # ── 4. Supertype density ──
    print("\n=== Supertype density ===")
    sup_dens = merge_density('supertype')
    if sup_dens is not None:
        print(f"  {len(sup_dens)} shared supertypes")
        make_concordance_figure(
            sup_dens, 'beta_snrnaseq', 'logFC_density', 'se_snrnaseq', 'se_density',
            ylabel='Xenium density logFC (cells/mm²)',
            suptitle='SCZ density effects: snRNAseq vs Xenium (supertype)',
            out_path=os.path.join(PRESENTATION_DIR, 'concordance_supertype_density.png'),
        )
        sup_dens.to_csv(os.path.join(DENSITY_DIR, 'snrnaseq_vs_density_comparison.csv'), index=False)

    print("\nDone!")


if __name__ == '__main__':
    main()
