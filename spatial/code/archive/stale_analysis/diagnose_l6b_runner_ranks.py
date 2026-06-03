#!/usr/bin/env python3
"""
For misplaced L6b cells, recompute full correlation against all subclass
centroids and report 2nd/3rd best type assignments.

Question: what would these cells be called if not L6b?
"""

import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import sparse
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    H5AD_DIR, SAMPLE_TO_DX, EXCLUDE_SAMPLES,
    SUBCLASS_TO_CLASS, load_sample_adata,
)

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "modules"))
from correlation_classifier import (
    build_subclass_centroids, correlate, _normalize_adata,
)

OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), "output", "l6b_diagnostics")
os.makedirs(OUT_DIR, exist_ok=True)

SAMPLES = [s for s in SAMPLE_TO_DX if s not in EXCLUDE_SAMPLES]
MISPLACED_THRESHOLD = 0.5

plt.rcParams.update({
    'font.size': 14, 'axes.titlesize': 18, 'axes.labelsize': 16,
    'xtick.labelsize': 12, 'ytick.labelsize': 12,
    'figure.facecolor': 'white',
})


def main():
    import time
    t0 = time.time()

    # ── Step 1: Load all samples and combine for centroid building ──
    print("Loading all samples to build centroids and extract L6b cells...")
    all_adata = []
    l6b_indices = {}  # sample -> list of cell barcodes

    for i, sample in enumerate(SAMPLES):
        print(f"  [{i+1}/{len(SAMPLES)}] {sample}...", end=" ", flush=True)
        adata = load_sample_adata(sample)
        mask = adata.obs['corr_qc_pass'].values.astype(bool)
        adata = adata[mask].copy()

        # Track L6b cells
        l6b_mask = adata.obs['corr_subclass'] == 'L6b'
        n_l6b = l6b_mask.sum()
        print(f"{n_l6b} L6b / {len(adata)} total")

        all_adata.append(adata)

    # Concatenate
    print("\nConcatenating all samples...")
    combined = ad.concat(all_adata, merge='same')
    print(f"Combined: {combined.shape[0]:,} cells × {combined.shape[1]} genes")

    # ── Step 2: Build subclass centroids (same as pipeline) ──
    print("\nBuilding subclass centroids...")
    centroids, cell_counts, gene_names = build_subclass_centroids(
        combined, top_n=100,
        subclass_col='subclass_label',
        confidence_col='subclass_label_confidence',
    )
    print(f"Centroids shape: {centroids.shape}")

    # ── Step 3: Extract L6b cells and compute full correlation matrix ──
    l6b_mask = combined.obs['corr_subclass'] == 'L6b'
    l6b_adata = combined[l6b_mask].copy()
    print(f"\nL6b cells to classify: {l6b_adata.shape[0]:,}")

    # Normalize L6b expression
    print("Normalizing L6b expression...")
    sc.pp.normalize_total(l6b_adata, target_sum=1e4)
    sc.pp.log1p(l6b_adata)

    X_l6b = l6b_adata.X
    if sparse.issparse(X_l6b):
        X_l6b = X_l6b.toarray()
    X_l6b = X_l6b.astype(np.float32)
    X_l6b = np.nan_to_num(X_l6b, nan=0.0)

    print("Computing full correlation matrix...")
    corr_matrix, type_names = correlate(X_l6b, centroids)
    print(f"Correlation matrix: {corr_matrix.shape} (cells × {len(type_names)} subclasses)")

    # ── Step 4: Extract top-3 matches for each cell ──
    ranks = np.argsort(-corr_matrix, axis=1)  # descending sort

    results = pd.DataFrame(index=l6b_adata.obs.index)
    results['sample_id'] = l6b_adata.obs['sample_id'].values
    results['depth'] = l6b_adata.obs['predicted_norm_depth'].values
    results['total_counts'] = l6b_adata.obs['total_counts'].values
    results['n_genes'] = l6b_adata.obs['n_genes'].values
    results['misplaced'] = results['depth'] < MISPLACED_THRESHOLD

    for rank in range(3):
        idx = ranks[:, rank]
        results[f'rank{rank+1}_type'] = [type_names[i] for i in idx]
        results[f'rank{rank+1}_corr'] = corr_matrix[np.arange(len(idx)), idx]

    results['margin_1v2'] = results['rank1_corr'] - results['rank2_corr']
    results['margin_1v3'] = results['rank1_corr'] - results['rank3_corr']

    # ── Step 5: Analyze misplaced vs correctly placed ──
    misplaced = results[results['misplaced']]
    correct = results[~results['misplaced']]

    print(f"\n{'='*70}")
    print(f"RESULTS: Top-3 subclass matches for L6b cells")
    print(f"{'='*70}")
    print(f"Total L6b: {len(results):,}")
    print(f"  Misplaced (depth < {MISPLACED_THRESHOLD}): {len(misplaced):,}")
    print(f"  Correctly placed (depth >= {MISPLACED_THRESHOLD}): {len(correct):,}")

    print(f"\n── Rank 1 (current assignment) ──")
    print(f"  Should be 100% L6b: {(results['rank1_type'] == 'L6b').mean():.1%}")

    print(f"\n── Rank 2 (next-best type) for MISPLACED L6b ──")
    r2_mis = misplaced['rank2_type'].value_counts()
    print(r2_mis.head(10).to_string())

    print(f"\n── Rank 2 for CORRECTLY PLACED L6b ──")
    r2_cor = correct['rank2_type'].value_counts()
    print(r2_cor.head(10).to_string())

    print(f"\n── Rank 3 for MISPLACED L6b ──")
    r3_mis = misplaced['rank3_type'].value_counts()
    print(r3_mis.head(10).to_string())

    print(f"\n── Margin (L6b corr - 2nd best corr) ──")
    print(f"  Misplaced:  median={misplaced['margin_1v2'].median():.4f}, "
          f"mean={misplaced['margin_1v2'].mean():.4f}")
    print(f"  Correct:    median={correct['margin_1v2'].median():.4f}, "
          f"mean={correct['margin_1v2'].mean():.4f}")

    # ── Step 6: Plots ──
    fig, axes = plt.subplots(2, 2, figsize=(18, 14))

    # Panel A: 2nd-best type distribution for misplaced vs correct
    ax = axes[0, 0]
    all_types = sorted(set(results['rank2_type']))
    # Compute fractions
    mis_frac = misplaced['rank2_type'].value_counts(normalize=True)
    cor_frac = correct['rank2_type'].value_counts(normalize=True)
    # Top types by misplaced frequency
    top_types = mis_frac.head(10).index.tolist()
    x = np.arange(len(top_types))
    w = 0.35
    ax.bar(x - w/2, [mis_frac.get(t, 0) for t in top_types], w,
           color='#d62728', alpha=0.7, label=f'Misplaced (n={len(misplaced):,})')
    ax.bar(x + w/2, [cor_frac.get(t, 0) for t in top_types], w,
           color='#2ca02c', alpha=0.7, label=f'Correct (n={len(correct):,})')
    ax.set_xticks(x)
    ax.set_xticklabels(top_types, rotation=45, ha='right')
    ax.set_ylabel('Fraction of cells')
    ax.set_title('2nd-Best Subclass Match')
    ax.legend(fontsize=12)

    # Panel B: Margin distribution by placement
    ax = axes[0, 1]
    bins = np.linspace(0, 0.15, 50)
    ax.hist(correct['margin_1v2'], bins=bins, density=True, alpha=0.5,
            color='#2ca02c', label='Correct')
    ax.hist(misplaced['margin_1v2'], bins=bins, density=True, alpha=0.5,
            color='#d62728', label='Misplaced')
    ax.set_xlabel('Margin (L6b corr − 2nd best)')
    ax.set_ylabel('Density')
    ax.set_title('Classification Margin: L6b vs Next-Best')
    ax.legend(fontsize=12)

    # Panel C: Scatter: depth vs margin, colored by 2nd-best type
    ax = axes[1, 0]
    # Color by broad 2nd-best class
    color_map = {}
    for t in all_types:
        cls = SUBCLASS_TO_CLASS.get(t, 'Unknown')
        if cls == 'Glutamatergic':
            color_map[t] = '#1f77b4'
        elif cls == 'GABAergic':
            color_map[t] = '#d62728'
        else:
            color_map[t] = '#7f7f7f'

    colors = [color_map.get(t, '#7f7f7f') for t in results['rank2_type']]
    ax.scatter(results['depth'], results['margin_1v2'],
               c=colors, s=3, alpha=0.3)
    ax.axvline(MISPLACED_THRESHOLD, color='orange', ls='--', lw=2)
    ax.set_xlabel('Normalized depth')
    ax.set_ylabel('Margin (L6b − 2nd best)')
    ax.set_title('Depth vs Margin (colored by 2nd-best class)')
    ax.legend(handles=[
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#1f77b4',
               markersize=8, label='2nd = Glutamatergic'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#d62728',
               markersize=8, label='2nd = GABAergic'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#7f7f7f',
               markersize=8, label='2nd = Non-neuronal'),
    ], fontsize=11)

    # Panel D: For very shallow cells, what's the 2nd best type vs depth
    ax = axes[1, 1]
    shallow = results[results['depth'] < 0.65].copy()
    top5_r2 = shallow['rank2_type'].value_counts().head(6).index.tolist()
    type_colors = plt.cm.Set1(np.linspace(0, 1, len(top5_r2)))
    for i, t in enumerate(top5_r2):
        subset = shallow[shallow['rank2_type'] == t]
        ax.scatter(subset['depth'], subset['rank2_corr'],
                   s=8, alpha=0.4, color=type_colors[i], label=t)
    ax.axvline(MISPLACED_THRESHOLD, color='orange', ls='--', lw=2)
    ax.set_xlabel('Normalized depth')
    ax.set_ylabel('2nd-best correlation')
    ax.set_title('Shallow L6b: 2nd-Best Type vs Depth')
    ax.legend(fontsize=10, markerscale=2)

    fig.suptitle('L6b Next-Best Type Analysis', fontsize=22, y=1.01)
    fig.tight_layout()
    path = os.path.join(OUT_DIR, 'l6b_rank2_rank3_analysis.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {path}")

    # ── Save full results table ──
    csv_path = os.path.join(OUT_DIR, 'l6b_top3_matches.csv')
    results.to_csv(csv_path)
    print(f"Saved: {csv_path}")

    # ── Summary table: misplaced 2nd-best breakdown ──
    print(f"\n{'='*70}")
    print("DETAILED BREAKDOWN: Misplaced L6b 2nd-best type")
    print(f"{'='*70}")
    summary = misplaced.groupby('rank2_type').agg(
        n_cells=('depth', 'count'),
        median_depth=('depth', 'median'),
        median_margin=('margin_1v2', 'median'),
        median_r2_corr=('rank2_corr', 'median'),
        median_umi=('total_counts', 'median'),
    ).sort_values('n_cells', ascending=False)
    summary['pct'] = 100 * summary['n_cells'] / len(misplaced)
    print(summary.head(15).to_string())

    summary_path = os.path.join(OUT_DIR, 'l6b_misplaced_2nd_best_summary.csv')
    summary.to_csv(summary_path)
    print(f"\nSaved: {summary_path}")

    print(f"\nDone in {time.time()-t0:.0f}s")


if __name__ == '__main__':
    main()
