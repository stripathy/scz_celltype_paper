#!/usr/bin/env python3
"""
Compare crumblr results: current default vs proposed QC (5th pctl margin filter).
Also compare both to snRNAseq betas.

Output: output/presentation/crumblr_current_vs_proposed.png
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import BASE_DIR

CRUMBLR_DIR = os.path.join(BASE_DIR, "output", "crumblr")
PRES_DIR = os.path.join(BASE_DIR, "output", "presentation")
SNRNASEQ_PATH = os.path.join(BASE_DIR, "data", "nicole_scz_snrnaseq_betas", "scz_coefs.xlsx")

# Cell class mapping for coloring
SUBCLASS_TO_CLASS = {
    'L2/3 IT': 'Glut', 'L4 IT': 'Glut', 'L5 ET': 'Glut', 'L5 IT': 'Glut',
    'L5/6 NP': 'Glut', 'L6 CT': 'Glut', 'L6 IT': 'Glut', 'L6 IT Car3': 'Glut',
    'L6b': 'Glut',
    'Lamp5': 'GABA', 'Lamp5 Lhx6': 'GABA', 'Sncg': 'GABA', 'Sst': 'GABA',
    'Sst Chodl': 'GABA', 'Pvalb': 'GABA', 'Pax6': 'GABA', 'Chandelier': 'GABA',
    'Vip': 'GABA',
    'Astrocyte': 'NN', 'Oligodendrocyte': 'NN', 'OPC': 'NN',
    'Microglia-PVM': 'NN', 'Endothelial': 'NN', 'VLMC': 'NN',
}

CLASS_COLORS = {'Glut': '#4C8BF5', 'GABA': '#E74C3C', 'NN': '#7F8C8D'}


def load_results(suffix):
    """Load subclass + supertype results for a given suffix."""
    dfs = []
    for level in ['subclass', 'supertype']:
        # Handle double-suffix naming from run_crumblr.R
        fname = f"crumblr_results_{level}_{suffix}_{suffix}.csv"
        fpath = os.path.join(CRUMBLR_DIR, fname)
        if not os.path.exists(fpath):
            fname = f"crumblr_results_{level}_{suffix}.csv"
            fpath = os.path.join(CRUMBLR_DIR, fname)
        if os.path.exists(fpath):
            df = pd.read_csv(fpath)
            df['level'] = level
            dfs.append(df)
            print(f"  Loaded {fname}: {len(df)} rows")
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def load_snrnaseq():
    """Load snRNAseq SCZ betas."""
    df = pd.read_excel(SNRNASEQ_PATH)
    df = df.rename(columns={'CellType': 'celltype', 'estimate': 'snrnaseq_beta'})
    return df[['celltype', 'snrnaseq_beta', 'pval', 'padj']]


def main():
    os.makedirs(PRES_DIR, exist_ok=True)

    # Load results
    print("Loading current baseline results...")
    current = load_results("current_baseline")
    print("Loading proposed results...")
    proposed = load_results("proposed")
    print("Loading snRNAseq betas...")
    snrna = load_snrnaseq()

    # ── Figure: 6-panel comparison ──
    fig, axes = plt.subplots(2, 3, figsize=(22, 14))
    fig.patch.set_facecolor('white')

    for row, level in enumerate(['subclass', 'supertype']):
        cur = current[current['level'] == level].copy()
        pro = proposed[proposed['level'] == level].copy()

        # Panel 1: Current vs Proposed logFC
        ax = axes[row, 0]
        m = cur.merge(pro, on='celltype', suffixes=('_cur', '_pro'))
        r = np.corrcoef(m['logFC_cur'], m['logFC_pro'])[0, 1]

        for _, row_data in m.iterrows():
            ct = row_data['celltype']
            cls = SUBCLASS_TO_CLASS.get(ct, 'NN')
            color = CLASS_COLORS.get(cls, '#7F8C8D')
            ax.scatter(row_data['logFC_cur'], row_data['logFC_pro'],
                      c=color, s=40, alpha=0.7, edgecolors='none')

        lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
                max(ax.get_xlim()[1], ax.get_ylim()[1])]
        ax.plot(lims, lims, 'k--', alpha=0.3, lw=1)
        ax.set_xlabel('Current logFC (no margin filter)', fontsize=14)
        ax.set_ylabel('Proposed logFC (5th pctl margin)', fontsize=14)
        ax.set_title(f'{level.capitalize()}: Current vs Proposed\nr={r:.4f}', fontsize=16)
        ax.tick_params(labelsize=12)

        # Annotate cells that changed FDR significance
        for _, row_data in m.iterrows():
            cur_sig = row_data['FDR_cur'] < 0.1
            pro_sig = row_data['FDR_pro'] < 0.1
            if cur_sig != pro_sig:
                ax.annotate(row_data['celltype'],
                           (row_data['logFC_cur'], row_data['logFC_pro']),
                           fontsize=8, alpha=0.8,
                           xytext=(5, 5), textcoords='offset points')

        # Panel 2: snRNAseq correlation — current
        ax = axes[row, 1]
        m_sn_cur = cur.merge(snrna, on='celltype')
        if len(m_sn_cur) > 0:
            r_cur = np.corrcoef(m_sn_cur['snrnaseq_beta'], m_sn_cur['logFC'])[0, 1]
            # Neuron-only correlation
            m_sn_cur['class'] = m_sn_cur['celltype'].map(SUBCLASS_TO_CLASS)
            neur = m_sn_cur[m_sn_cur['class'].isin(['Glut', 'GABA'])]
            r_neur_cur = np.corrcoef(neur['snrnaseq_beta'], neur['logFC'])[0, 1] if len(neur) > 2 else np.nan

            for _, row_data in m_sn_cur.iterrows():
                cls = SUBCLASS_TO_CLASS.get(row_data['celltype'], 'NN')
                color = CLASS_COLORS.get(cls, '#7F8C8D')
                ax.scatter(row_data['snrnaseq_beta'], row_data['logFC'],
                          c=color, s=40, alpha=0.7, edgecolors='none')

            ax.axhline(0, color='gray', lw=0.5, alpha=0.3)
            ax.axvline(0, color='gray', lw=0.5, alpha=0.3)
            # Regression line
            z = np.polyfit(m_sn_cur['snrnaseq_beta'], m_sn_cur['logFC'], 1)
            x_fit = np.linspace(m_sn_cur['snrnaseq_beta'].min(), m_sn_cur['snrnaseq_beta'].max(), 100)
            ax.plot(x_fit, np.polyval(z, x_fit), 'k-', alpha=0.3, lw=1)

            ax.set_xlabel('snRNAseq beta', fontsize=14)
            ax.set_ylabel('Xenium logFC (current)', fontsize=14)
            n_glut = (m_sn_cur['class'] == 'Glut').sum()
            n_gaba = (m_sn_cur['class'] == 'GABA').sum()
            n_nn = (m_sn_cur['class'] == 'NN').sum()
            ax.set_title(f'{level.capitalize()}: Current vs snRNAseq\n'
                        f'r={r_cur:.3f} | neur r={r_neur_cur:.3f}\n'
                        f'Glut(n={n_glut}) GABA(n={n_gaba}) NN(n={n_nn})',
                        fontsize=14)
            ax.tick_params(labelsize=12)

        # Panel 3: snRNAseq correlation — proposed
        ax = axes[row, 2]
        m_sn_pro = pro.merge(snrna, on='celltype')
        if len(m_sn_pro) > 0:
            r_pro = np.corrcoef(m_sn_pro['snrnaseq_beta'], m_sn_pro['logFC'])[0, 1]
            m_sn_pro['class'] = m_sn_pro['celltype'].map(SUBCLASS_TO_CLASS)
            neur = m_sn_pro[m_sn_pro['class'].isin(['Glut', 'GABA'])]
            r_neur_pro = np.corrcoef(neur['snrnaseq_beta'], neur['logFC'])[0, 1] if len(neur) > 2 else np.nan

            for _, row_data in m_sn_pro.iterrows():
                cls = SUBCLASS_TO_CLASS.get(row_data['celltype'], 'NN')
                color = CLASS_COLORS.get(cls, '#7F8C8D')
                ax.scatter(row_data['snrnaseq_beta'], row_data['logFC'],
                          c=color, s=40, alpha=0.7, edgecolors='none')

            ax.axhline(0, color='gray', lw=0.5, alpha=0.3)
            ax.axvline(0, color='gray', lw=0.5, alpha=0.3)
            z = np.polyfit(m_sn_pro['snrnaseq_beta'], m_sn_pro['logFC'], 1)
            x_fit = np.linspace(m_sn_pro['snrnaseq_beta'].min(), m_sn_pro['snrnaseq_beta'].max(), 100)
            ax.plot(x_fit, np.polyval(z, x_fit), 'k-', alpha=0.3, lw=1)

            ax.set_xlabel('snRNAseq beta', fontsize=14)
            ax.set_ylabel('Xenium logFC (proposed)', fontsize=14)
            n_glut = (m_sn_pro['class'] == 'Glut').sum()
            n_gaba = (m_sn_pro['class'] == 'GABA').sum()
            n_nn = (m_sn_pro['class'] == 'NN').sum()
            ax.set_title(f'{level.capitalize()}: Proposed vs snRNAseq\n'
                        f'r={r_pro:.3f} | neur r={r_neur_pro:.3f}\n'
                        f'Glut(n={n_glut}) GABA(n={n_gaba}) NN(n={n_nn})',
                        fontsize=14)
            ax.tick_params(labelsize=12)

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=CLASS_COLORS['Glut'],
               markersize=10, label='Glutamatergic'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=CLASS_COLORS['GABA'],
               markersize=10, label='GABAergic'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=CLASS_COLORS['NN'],
               markersize=10, label='Non-neuronal'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=3, fontsize=14,
              bbox_to_anchor=(0.5, -0.02))

    fig.suptitle('Crumblr Compositional Analysis: Current vs Proposed QC\n'
                 'Proposed = spatial QC + 5th percentile margin filter (drops 5% cells)',
                 fontsize=18, fontweight='bold', y=1.02)
    plt.tight_layout()
    outpath = os.path.join(PRES_DIR, "crumblr_current_vs_proposed.png")
    fig.savefig(outpath, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"\nSaved: {outpath}")
    plt.close()

    # ── Summary table ──
    print("\n" + "="*80)
    print("SUMMARY: Crumblr results comparison")
    print("="*80)

    for level in ['subclass', 'supertype']:
        cur = current[current['level'] == level]
        pro = proposed[proposed['level'] == level]

        m = cur.merge(pro, on='celltype', suffixes=('_cur', '_pro'))
        r_logfc = np.corrcoef(m['logFC_cur'], m['logFC_pro'])[0, 1]

        m_sn_cur = cur.merge(snrna, on='celltype')
        m_sn_pro = pro.merge(snrna, on='celltype')
        r_sn_cur = np.corrcoef(m_sn_cur['snrnaseq_beta'], m_sn_cur['logFC'])[0, 1] if len(m_sn_cur) > 2 else np.nan
        r_sn_pro = np.corrcoef(m_sn_pro['snrnaseq_beta'], m_sn_pro['logFC'])[0, 1] if len(m_sn_pro) > 2 else np.nan

        # Neuron only
        m_sn_cur['class'] = m_sn_cur['celltype'].map(SUBCLASS_TO_CLASS)
        m_sn_pro['class'] = m_sn_pro['celltype'].map(SUBCLASS_TO_CLASS)
        neur_cur = m_sn_cur[m_sn_cur['class'].isin(['Glut', 'GABA'])]
        neur_pro = m_sn_pro[m_sn_pro['class'].isin(['Glut', 'GABA'])]
        r_neur_cur = np.corrcoef(neur_cur['snrnaseq_beta'], neur_cur['logFC'])[0, 1] if len(neur_cur) > 2 else np.nan
        r_neur_pro = np.corrcoef(neur_pro['snrnaseq_beta'], neur_pro['logFC'])[0, 1] if len(neur_pro) > 2 else np.nan

        n_fdr05_cur = (cur['FDR'] < 0.05).sum()
        n_fdr05_pro = (pro['FDR'] < 0.05).sum()
        n_fdr10_cur = (cur['FDR'] < 0.10).sum()
        n_fdr10_pro = (pro['FDR'] < 0.10).sum()

        print(f"\n{level.upper()}:")
        print(f"  Current vs Proposed logFC r = {r_logfc:.4f}")
        print(f"  snRNAseq correlation:  current r={r_sn_cur:.3f}  proposed r={r_sn_pro:.3f}  (Δ={r_sn_pro-r_sn_cur:+.3f})")
        print(f"  snRNAseq neur only:    current r={r_neur_cur:.3f}  proposed r={r_neur_pro:.3f}  (Δ={r_neur_pro-r_neur_cur:+.3f})")
        print(f"  FDR < 0.05:  current={n_fdr05_cur}  proposed={n_fdr05_pro}")
        print(f"  FDR < 0.10:  current={n_fdr10_cur}  proposed={n_fdr10_pro}")

        # Show FDR changes for hits
        for _, row_data in m.iterrows():
            cur_sig = row_data['FDR_cur'] < 0.1
            pro_sig = row_data['FDR_pro'] < 0.1
            if cur_sig or pro_sig:
                direction = "↑SCZ" if row_data['logFC_cur'] > 0 else "↓SCZ"
                print(f"    {row_data['celltype']:30s} {direction} "
                      f"logFC: {row_data['logFC_cur']:+.4f}->{row_data['logFC_pro']:+.4f}  "
                      f"FDR: {row_data['FDR_cur']:.4f}->{row_data['FDR_pro']:.4f}")

    # ── Also output the pctl05 comparison if available ──
    pctl05_sub = os.path.join(CRUMBLR_DIR, "crumblr_results_subclass_pctl05.csv")
    if os.path.exists(pctl05_sub):
        print("\n" + "="*80)
        print("COMPARISON WITH EXISTING pctl05 RESULTS")
        print("="*80)
        p05 = pd.read_csv(pctl05_sub)
        pro_sub = proposed[proposed['level'] == 'subclass']
        m = p05.merge(pro_sub, on='celltype', suffixes=('_p05', '_pro'))
        r = np.corrcoef(m['logFC_p05'], m['logFC_pro'])[0, 1]
        print(f"  pctl05 vs proposed subclass logFC r = {r:.4f}")
        print(f"  Max |diff|: {(m['logFC_p05'] - m['logFC_pro']).abs().max():.4f}")
        print(f"  pctl05 FDR<0.05: {(p05['FDR']<0.05).sum()}, proposed FDR<0.05: {(pro_sub['FDR']<0.05).sum()}")

    print("\nDone!")


if __name__ == "__main__":
    main()
