#!/usr/bin/env python3
"""
Validation figures for the hybrid nuclear doublet resolution pipeline (optional, unnumbered).

Loads all processed samples, computes aggregate statistics, and generates
6 diagnostic figures plus a summary CSV.

Requires: Nuclear doublet resolution (04_run_nuclear_doublet_resolution.py) must have been run.

Output figures (saved to output/presentation/):
  1. nuclear_doublet_resolution_summary.png
  2. nuclear_fraction_distributions.png
  3. nuclear_doublet_marker_evidence.png
  4. nuclear_doublet_aggregate_table.png
  5. nuclear_doublet_scz_vs_control.png
  6. nuclear_hybrid_qc_impact.png

Usage:
    python3 -u plot_nuclear_doublet_validation.py
"""

import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from scipy import sparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    H5AD_DIR, PRESENTATION_DIR,
    SAMPLE_TO_DX, EXCLUDE_SAMPLES, SUBCLASS_TO_CLASS,
    CLASS_COLORS, DX_COLORS,
)

# ── Constants ──
STATUS_COLORS = {
    'resolved': '#4CAF50',       # green
    'persistent': '#F44336',     # red
    'insufficient': '#9E9E9E',   # grey
    'nuclear_only': '#FF9800',   # orange
}

STATUS_ORDER = ['resolved', 'persistent', 'insufficient', 'nuclear_only']

# GABA/Glut marker genes (same as correlation_classifier.py)
GABA_MARKERS = ['GAD1', 'GAD2', 'SLC32A1', 'SST', 'PVALB', 'VIP', 'LAMP5']
GLUT_MARKERS = ['CUX2', 'RORB', 'GRIN2A', 'THEMIS']


# ══════════════════════════════════════════════════════════════════════
# Data loading
# ══════════════════════════════════════════════════════════════════════

def load_all_samples():
    """Load obs metadata from all processed samples (those with hybrid_qc_pass)."""
    sample_ids = sorted(set(SAMPLE_TO_DX.keys()) - EXCLUDE_SAMPLES)
    all_obs = []
    loaded = []

    for sid in sample_ids:
        fpath = os.path.join(H5AD_DIR, f"{sid}_annotated.h5ad")
        if not os.path.exists(fpath):
            continue
        adata = ad.read_h5ad(fpath, backed='r')

        if 'hybrid_qc_pass' not in adata.obs.columns:
            continue

        # Load relevant columns
        cols = ['sample_id', 'qc_pass', 'corr_class', 'corr_subclass',
                'corr_qc_pass', 'doublet_suspect', 'doublet_type',
                'nuclear_total_counts', 'nuclear_n_genes', 'nuclear_fraction',
                'nuclear_doublet_suspect', 'nuclear_doublet_type',
                'nuclear_doublet_status', 'hybrid_qc_pass', 'total_counts',
                'fail_total_counts_high']
        avail_cols = [c for c in cols if c in adata.obs.columns]
        obs = adata.obs[avail_cols].copy()
        obs['diagnosis'] = SAMPLE_TO_DX.get(sid, 'Unknown')
        all_obs.append(obs)
        loaded.append(sid)

    if not all_obs:
        print("ERROR: No samples with hybrid_qc_pass found. Run 04_run_nuclear_doublet_resolution.py first.")
        return None, []

    df = pd.concat(all_obs, axis=0)
    print(f"Loaded {len(loaded)} samples, {len(df):,} total cells")
    return df, loaded


def load_expression_for_doublets(sample_ids):
    """Load raw expression for doublet cells (WC and nuclear) for marker analysis."""
    records = []

    for sid in sample_ids:
        fpath = os.path.join(H5AD_DIR, f"{sid}_annotated.h5ad")
        adata = ad.read_h5ad(fpath)

        # Get doublet cells (either WC or nuclear)
        is_doublet = (adata.obs['doublet_suspect'].values.astype(bool) |
                      adata.obs['nuclear_doublet_suspect'].values.astype(bool))
        if not is_doublet.any():
            continue

        qc = adata.obs['qc_pass'].values.astype(bool)
        mask = is_doublet & qc
        if not mask.any():
            continue

        # Whole-cell expression
        X_wc = adata[mask].X
        if sparse.issparse(X_wc):
            X_wc = X_wc.toarray()

        gene_names = list(adata.var_names)
        status = adata.obs.loc[mask, 'nuclear_doublet_status'].values
        wc_class = adata.obs.loc[mask, 'corr_class'].astype(str).values

        for i in range(mask.sum()):
            for marker_list, marker_type in [(GABA_MARKERS, 'GABA'), (GLUT_MARKERS, 'Glut')]:
                for gene in marker_list:
                    if gene in gene_names:
                        gi = gene_names.index(gene)
                        records.append({
                            'sample_id': sid,
                            'status': status[i],
                            'wc_class': wc_class[i],
                            'gene': gene,
                            'marker_type': marker_type,
                            'wc_expr': float(X_wc[i, gi]),
                        })

    return pd.DataFrame(records)


# ══════════════════════════════════════════════════════════════════════
# Figure 1: Per-sample resolution summary (stacked bar)
# ══════════════════════════════════════════════════════════════════════

def plot_resolution_summary(df, loaded_samples):
    """Stacked bar chart of resolution status per sample."""
    fig, ax = plt.subplots(figsize=(16, 7), facecolor='white')

    # Count by sample and status
    qc = df['qc_pass'].astype(bool)
    doublet_df = df[qc & (df['nuclear_doublet_status'] != 'clean')]

    sample_order = sorted(loaded_samples, key=lambda s: SAMPLE_TO_DX.get(s, ''))
    n_samples = len(sample_order)

    x = np.arange(n_samples)
    width = 0.7
    bottoms = np.zeros(n_samples)

    for status in STATUS_ORDER:
        counts = []
        for sid in sample_order:
            mask = (doublet_df['sample_id'] == sid) & (doublet_df['nuclear_doublet_status'] == status)
            counts.append(mask.sum())
        counts = np.array(counts)
        ax.bar(x, counts, width, bottom=bottoms, label=status.capitalize(),
               color=STATUS_COLORS[status], edgecolor='white', linewidth=0.5)
        bottoms += counts

    ax.set_xticks(x)
    labels = [f"{s}\n({'SCZ' if SAMPLE_TO_DX.get(s)=='SCZ' else 'Ctrl'})"
              for s in sample_order]
    ax.set_xticklabels(labels, fontsize=11, rotation=45, ha='right')
    ax.set_ylabel('Number of Doublet Cells', fontsize=16)
    ax.set_title('Nuclear Doublet Resolution by Sample', fontsize=20, fontweight='bold')
    ax.legend(fontsize=13, loc='upper right')
    ax.tick_params(labelsize=12)

    plt.tight_layout()
    out = os.path.join(PRESENTATION_DIR, 'nuclear_doublet_resolution_summary.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ══════════════════════════════════════════════════════════════════════
# Figure 2: Nuclear fraction distributions by class
# ══════════════════════════════════════════════════════════════════════

def plot_nuclear_fraction_distributions(df):
    """Violin plots of nuclear fraction per class (all samples)."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 7), facecolor='white')

    qc = df['qc_pass'].astype(bool) & df['corr_qc_pass'].astype(bool)
    df_qc = df[qc].copy()
    nf = df_qc['nuclear_fraction'].astype(float)

    # Panel A: By class
    ax = axes[0]
    classes = ['Glutamatergic', 'GABAergic', 'Non-neuronal']
    positions = np.arange(len(classes))
    data_by_class = [nf[df_qc['corr_class'] == c].values for c in classes]

    parts = ax.violinplot(data_by_class, positions=positions, showmedians=True,
                          showextrema=False)
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(CLASS_COLORS.get(classes[i], '#888'))
        pc.set_alpha(0.7)
    parts['cmedians'].set_color('black')

    ax.set_xticks(positions)
    ax.set_xticklabels(classes, fontsize=14)
    ax.set_ylabel('Nuclear Fraction', fontsize=16)
    ax.set_title('Nuclear Fraction by Class', fontsize=18, fontweight='bold')
    ax.tick_params(labelsize=12)

    # Add median labels
    for i, data in enumerate(data_by_class):
        med = np.median(data)
        ax.text(i, med + 0.03, f'{med:.2f}', ha='center', fontsize=12,
                fontweight='bold')

    # Panel B: Histogram
    ax = axes[1]
    for cls in classes:
        vals = nf[df_qc['corr_class'] == cls].values
        ax.hist(vals, bins=60, range=(0, 1.5), alpha=0.6,
                color=CLASS_COLORS.get(cls, '#888'), label=cls, density=True)

    ax.set_xlabel('Nuclear Fraction', fontsize=16)
    ax.set_ylabel('Density', fontsize=16)
    ax.set_title('Nuclear Fraction Distribution', fontsize=18, fontweight='bold')
    ax.legend(fontsize=13)
    ax.tick_params(labelsize=12)
    ax.axvline(0.5, color='grey', linestyle='--', alpha=0.5)

    plt.tight_layout()
    out = os.path.join(PRESENTATION_DIR, 'nuclear_fraction_distributions.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ══════════════════════════════════════════════════════════════════════
# Figure 3: Marker evidence heatmap (resolved vs persistent)
# ══════════════════════════════════════════════════════════════════════

def plot_marker_evidence(df, loaded_samples):
    """Detection rate of GABA markers in WC for resolved vs persistent doublets."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 7), facecolor='white')

    qc = df['qc_pass'].astype(bool)

    for ax, dtype, title in [(axes[0], 'Glut+GABA', 'Glut+GABA Doublets'),
                              (axes[1], 'GABA+GABA', 'GABA+GABA Doublets')]:
        # Get resolved vs persistent for this doublet type
        type_mask = qc & (df['doublet_type'] == dtype)
        resolved = type_mask & (df['nuclear_doublet_status'] == 'resolved')
        persistent = type_mask & (df['nuclear_doublet_status'] == 'persistent')

        n_resolved = resolved.sum()
        n_persistent = persistent.sum()

        # Summary text
        total = n_resolved + n_persistent
        if total > 0:
            pct_resolved = 100 * n_resolved / total
            ax.text(0.5, 0.95,
                    f'Resolved: {n_resolved:,} ({pct_resolved:.0f}%)  |  '
                    f'Persistent: {n_persistent:,} ({100-pct_resolved:.0f}%)',
                    transform=ax.transAxes, ha='center', va='top',
                    fontsize=14, fontweight='bold',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        # Nuclear fraction distributions
        if n_resolved > 5 and n_persistent > 5:
            nf_res = df.loc[resolved, 'nuclear_fraction'].astype(float).values
            nf_per = df.loc[persistent, 'nuclear_fraction'].astype(float).values
            ax.hist(nf_res, bins=40, range=(0, 1.2), alpha=0.6,
                    color=STATUS_COLORS['resolved'], label='Resolved', density=True)
            ax.hist(nf_per, bins=40, range=(0, 1.2), alpha=0.6,
                    color=STATUS_COLORS['persistent'], label='Persistent', density=True)
            ax.set_xlabel('Nuclear Fraction', fontsize=14)
            ax.set_ylabel('Density', fontsize=14)
            ax.legend(fontsize=12)
        else:
            ax.text(0.5, 0.5, f'n={total} (too few for histogram)',
                    transform=ax.transAxes, ha='center', fontsize=14)

        ax.set_title(title, fontsize=18, fontweight='bold')
        ax.tick_params(labelsize=12)

    plt.tight_layout()
    out = os.path.join(PRESENTATION_DIR, 'nuclear_doublet_marker_evidence.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ══════════════════════════════════════════════════════════════════════
# Figure 4: Aggregate summary table
# ══════════════════════════════════════════════════════════════════════

def plot_aggregate_table(df, loaded_samples):
    """Render summary table as a figure."""
    fig, ax = plt.subplots(figsize=(18, max(5, len(loaded_samples) * 0.4 + 2)),
                           facecolor='white')
    ax.axis('off')

    qc = df['qc_pass'].astype(bool)

    # Build table data
    rows = []
    sample_order = sorted(loaded_samples, key=lambda s: SAMPLE_TO_DX.get(s, ''))

    for sid in sample_order:
        sm = (df['sample_id'] == sid) & qc
        n_cells = sm.sum()
        n_wc_dbl = (sm & df['doublet_suspect'].astype(bool)).sum()
        n_res = (sm & (df['nuclear_doublet_status'] == 'resolved')).sum()
        n_per = (sm & (df['nuclear_doublet_status'] == 'persistent')).sum()
        n_ins = (sm & (df['nuclear_doublet_status'] == 'insufficient')).sum()
        n_nuc = (sm & (df['nuclear_doublet_status'] == 'nuclear_only')).sum()
        rate = f"{100*n_res/n_wc_dbl:.0f}%" if n_wc_dbl > 0 else "—"
        n_corr = (sm & df['corr_qc_pass'].astype(bool)).sum()
        n_hyb = (sm & df['hybrid_qc_pass'].astype(bool)).sum()
        net = n_hyb - n_corr
        dx = SAMPLE_TO_DX.get(sid, '?')
        rows.append([sid, dx, f"{n_cells:,}", f"{n_wc_dbl:,}",
                      f"{n_res:,}", f"{n_per:,}", f"{n_ins:,}", f"{n_nuc:,}",
                      rate, f"{net:+,}"])

    col_labels = ['Sample', 'Dx', 'QC Cells', 'WC Dbl',
                   'Resolved', 'Persist', 'Insuff', 'Nuc Only',
                   'Res Rate', 'Net ±']

    table = ax.table(cellText=rows, colLabels=col_labels,
                      loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1.0, 1.4)

    # Style header
    for j in range(len(col_labels)):
        table[(0, j)].set_facecolor('#333333')
        table[(0, j)].set_text_props(color='white', fontweight='bold')

    # Color code rows by diagnosis
    for i, row in enumerate(rows):
        dx = row[1]
        color = '#E3F2FD' if dx == 'Control' else '#FFEBEE'
        for j in range(len(col_labels)):
            table[(i + 1, j)].set_facecolor(color)

    ax.set_title('Nuclear Doublet Resolution: Per-Sample Summary',
                  fontsize=20, fontweight='bold', pad=20)

    plt.tight_layout()
    out = os.path.join(PRESENTATION_DIR, 'nuclear_doublet_aggregate_table.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ══════════════════════════════════════════════════════════════════════
# Figure 5: SCZ vs Control comparison
# ══════════════════════════════════════════════════════════════════════

def plot_scz_vs_control(df, loaded_samples):
    """Grouped bar: doublet and resolution rates by diagnosis."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 7), facecolor='white')

    qc = df['qc_pass'].astype(bool)

    # Compute per-sample metrics
    per_sample = []
    for sid in loaded_samples:
        sm = (df['sample_id'] == sid) & qc
        n_cells = sm.sum()
        if n_cells == 0:
            continue
        n_wc_dbl = (sm & df['doublet_suspect'].astype(bool)).sum()
        n_res = (sm & (df['nuclear_doublet_status'] == 'resolved')).sum()
        n_per = (sm & (df['nuclear_doublet_status'] == 'persistent')).sum()
        dx = SAMPLE_TO_DX.get(sid, '?')
        per_sample.append({
            'sample_id': sid,
            'diagnosis': dx,
            'wc_doublet_rate': 100 * n_wc_dbl / n_cells,
            'resolution_rate': 100 * n_res / n_wc_dbl if n_wc_dbl > 0 else 0,
            'persistent_rate': 100 * n_per / n_cells,
            'n_wc_dbl': n_wc_dbl,
            'n_resolved': n_res,
            'n_persistent': n_per,
        })

    ps = pd.DataFrame(per_sample)

    # Panel A: WC doublet rate by diagnosis
    ax = axes[0]
    for i, dx in enumerate(['Control', 'SCZ']):
        vals = ps[ps['diagnosis'] == dx]['wc_doublet_rate'].values
        positions = np.random.normal(i, 0.08, len(vals))
        ax.scatter(positions, vals, c=DX_COLORS[dx], s=60, alpha=0.7, zorder=3)
        ax.bar(i, np.mean(vals), width=0.5, alpha=0.3, color=DX_COLORS[dx])
        ax.errorbar(i, np.mean(vals), yerr=np.std(vals), color='black',
                     capsize=8, capthick=2, linewidth=2, zorder=4)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Control', 'SCZ'], fontsize=16)
    ax.set_ylabel('WC Doublet Rate (%)', fontsize=16)
    ax.set_title('Whole-Cell Doublet Rate', fontsize=18, fontweight='bold')
    ax.tick_params(labelsize=12)

    # Panel B: Resolution rate by diagnosis
    ax = axes[1]
    for i, dx in enumerate(['Control', 'SCZ']):
        vals = ps[ps['diagnosis'] == dx]['resolution_rate'].values
        positions = np.random.normal(i, 0.08, len(vals))
        ax.scatter(positions, vals, c=DX_COLORS[dx], s=60, alpha=0.7, zorder=3)
        ax.bar(i, np.mean(vals), width=0.5, alpha=0.3, color=DX_COLORS[dx])
        ax.errorbar(i, np.mean(vals), yerr=np.std(vals), color='black',
                     capsize=8, capthick=2, linewidth=2, zorder=4)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Control', 'SCZ'], fontsize=16)
    ax.set_ylabel('Resolution Rate (%)', fontsize=16)
    ax.set_title('Nuclear Resolution Rate', fontsize=18, fontweight='bold')
    ax.tick_params(labelsize=12)

    # Panel C: Persistent doublet rate by diagnosis
    ax = axes[2]
    for i, dx in enumerate(['Control', 'SCZ']):
        vals = ps[ps['diagnosis'] == dx]['persistent_rate'].values
        positions = np.random.normal(i, 0.08, len(vals))
        ax.scatter(positions, vals, c=DX_COLORS[dx], s=60, alpha=0.7, zorder=3)
        ax.bar(i, np.mean(vals), width=0.5, alpha=0.3, color=DX_COLORS[dx])
        ax.errorbar(i, np.mean(vals), yerr=np.std(vals), color='black',
                     capsize=8, capthick=2, linewidth=2, zorder=4)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Control', 'SCZ'], fontsize=16)
    ax.set_ylabel('Persistent Doublet Rate (%)', fontsize=16)
    ax.set_title('Persistent Doublet Rate', fontsize=18, fontweight='bold')
    ax.tick_params(labelsize=12)

    plt.tight_layout()
    out = os.path.join(PRESENTATION_DIR, 'nuclear_doublet_scz_vs_control.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ══════════════════════════════════════════════════════════════════════
# Figure 6: Hybrid QC compositional impact
# ══════════════════════════════════════════════════════════════════════

def plot_hybrid_qc_impact(df):
    """Subclass proportions: corr_qc_pass vs hybrid_qc_pass."""
    fig, axes = plt.subplots(1, 2, figsize=(18, 8), facecolor='white')

    qc = df['qc_pass'].astype(bool)

    # Get subclass labels
    corr_pass = qc & df['corr_qc_pass'].astype(bool)
    hybrid_pass = qc & df['hybrid_qc_pass'].astype(bool)

    corr_sub = df.loc[corr_pass, 'corr_subclass'].astype(str)
    hybrid_sub = df.loc[hybrid_pass, 'corr_subclass'].astype(str)

    # Get all subclasses, ordered by class then count
    all_subs = sorted(set(corr_sub) | set(hybrid_sub))
    # Order by class
    class_order = {'Glutamatergic': 0, 'GABAergic': 1, 'Non-neuronal': 2}
    all_subs = sorted(all_subs,
                       key=lambda s: (class_order.get(SUBCLASS_TO_CLASS.get(s, 'Unknown'), 3), s))

    # Panel A: Paired bar chart of proportions
    ax = axes[0]
    corr_counts = corr_sub.value_counts()
    hybrid_counts = hybrid_sub.value_counts()
    corr_props = {s: 100 * corr_counts.get(s, 0) / len(corr_sub) for s in all_subs}
    hybrid_props = {s: 100 * hybrid_counts.get(s, 0) / len(hybrid_sub) for s in all_subs}

    y_pos = np.arange(len(all_subs))
    width = 0.35

    bars1 = ax.barh(y_pos - width/2, [corr_props[s] for s in all_subs],
                     width, label='corr_qc_pass', color='#1976D2', alpha=0.7)
    bars2 = ax.barh(y_pos + width/2, [hybrid_props[s] for s in all_subs],
                     width, label='hybrid_qc_pass', color='#4CAF50', alpha=0.7)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(all_subs, fontsize=11)
    ax.set_xlabel('Proportion (%)', fontsize=14)
    ax.set_title('Subclass Composition', fontsize=18, fontweight='bold')
    ax.legend(fontsize=12)
    ax.tick_params(labelsize=11)
    ax.invert_yaxis()

    # Panel B: Difference (hybrid - corr) in absolute proportion
    ax = axes[1]
    diffs = [hybrid_props[s] - corr_props[s] for s in all_subs]
    colors = ['#4CAF50' if d >= 0 else '#F44336' for d in diffs]
    ax.barh(y_pos, diffs, 0.6, color=colors, alpha=0.7)
    ax.axvline(0, color='black', linewidth=1)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(all_subs, fontsize=11)
    ax.set_xlabel('Δ Proportion (hybrid − corr) [pp]', fontsize=14)
    ax.set_title('Compositional Change from Hybrid QC', fontsize=18, fontweight='bold')
    ax.tick_params(labelsize=11)
    ax.invert_yaxis()

    # Add count annotations
    for i, s in enumerate(all_subs):
        corr_n = corr_counts.get(s, 0)
        hyb_n = hybrid_counts.get(s, 0)
        delta = hyb_n - corr_n
        if delta != 0:
            ax.text(diffs[i] + (0.002 if diffs[i] >= 0 else -0.002), i,
                    f'{delta:+,}', va='center',
                    ha='left' if diffs[i] >= 0 else 'right',
                    fontsize=9, color='#333')

    plt.tight_layout()
    out = os.path.join(PRESENTATION_DIR, 'nuclear_hybrid_qc_impact.png')
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ══════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════

def main():
    print("=" * 70)
    print("Nuclear Doublet Resolution: Validation Figures")
    print("=" * 70)

    # Load data
    print("\nLoading samples...")
    df, loaded_samples = load_all_samples()
    if df is None:
        return

    os.makedirs(PRESENTATION_DIR, exist_ok=True)

    # Generate figures
    print("\nGenerating figures...")

    print("\n  Fig 1: Per-sample resolution summary...")
    plot_resolution_summary(df, loaded_samples)

    print("\n  Fig 2: Nuclear fraction distributions...")
    plot_nuclear_fraction_distributions(df)

    print("\n  Fig 3: Marker evidence (resolved vs persistent)...")
    plot_marker_evidence(df, loaded_samples)

    print("\n  Fig 4: Aggregate summary table...")
    plot_aggregate_table(df, loaded_samples)

    print("\n  Fig 5: SCZ vs Control comparison...")
    plot_scz_vs_control(df, loaded_samples)

    print("\n  Fig 6: Hybrid QC compositional impact...")
    plot_hybrid_qc_impact(df)

    print("\nDone!")


if __name__ == '__main__':
    main()
