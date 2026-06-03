#!/usr/bin/env python3
"""
Capture-rate-weighted correlation classifier for Xenium cell types.

Builds reference expression centroids from high-confidence HANN-labeled
Xenium cells, then reclassifies all cells using weighted Pearson correlation
where each gene's contribution is weighted by its transcript capture rate
(fraction of transcripts inside cell boundaries).

Motivation: genes with low capture rates (high ambient RNA contamination)
contribute noise to classification. Downweighting these should reduce
misclassifications like L6b being assigned to upper cortical layers.
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# ── Path setup ──
BASE_DIR = os.path.expanduser("~/Github/SCZ_Xenium")
sys.path.insert(0, os.path.join(BASE_DIR, "code", "analysis"))
from config import (SAMPLE_TO_DX, EXCLUDE_SAMPLES, SUBCLASS_TO_CLASS,
                     H5AD_DIR, LAYER_COLORS, CORTICAL_LAYERS, LAYER_ORDER)

CAPTURE_RATE_CSV = os.path.join(BASE_DIR, "output", "plots",
                                 "transcript_capture_by_gene.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "plots")
os.makedirs(OUTPUT_DIR, exist_ok=True)

UPPER_LAYERS = {"L1", "L2/3", "L4"}
CONF_THRESHOLD = 0.7   # for building reference centroids
MIN_CELLS = 100         # minimum cells per subclass for centroid


# ═══════════════════════════════════════════════════════════════════════
# 1. DATA LOADING
# ═══════════════════════════════════════════════════════════════════════

def load_capture_rates():
    """Load per-gene transcript capture rates."""
    df = pd.read_csv(CAPTURE_RATE_CSV)
    rates = dict(zip(df['gene'], df['pct_inside']))
    print(f"  Loaded capture rates for {len(rates)} genes")
    print(f"  Median capture: {df['pct_inside'].median():.1f}%")
    print(f"  Range: {df['pct_inside'].min():.1f}% — {df['pct_inside'].max():.1f}%")

    # Flag MT genes
    mt_genes = [g for g in rates if g.startswith('MT-') or g.startswith('MTRNR')]
    print(f"  MT/MTRNR genes ({len(mt_genes)}): "
          f"{', '.join(f'{g}={rates[g]:.0f}%' for g in sorted(mt_genes)[:5])}...")
    return rates


def load_all_samples():
    """Load expression + metadata for all samples."""
    sample_ids = sorted(set(SAMPLE_TO_DX.keys()) - EXCLUDE_SAMPLES)
    print(f"\nLoading {len(sample_ids)} samples...")

    adatas = []
    for i, sid in enumerate(sample_ids):
        fpath = os.path.join(H5AD_DIR, f"{sid}_annotated.h5ad")
        adata = ad.read_h5ad(fpath)

        # Filter to QC-pass
        if "qc_pass" in adata.obs.columns:
            adata = adata[adata.obs["qc_pass"].values.astype(bool)].copy()

        print(f"  [{i+1:2d}/{len(sample_ids)}] {sid}: {adata.n_obs:,} cells", flush=True)
        adatas.append(adata)

    # Concatenate
    print("  Concatenating...", flush=True)
    combined = ad.concat(adatas, join='outer')
    print(f"  Total: {combined.n_obs:,} cells x {combined.n_vars} genes")
    return combined


# ═══════════════════════════════════════════════════════════════════════
# 2. BUILD REFERENCE CENTROIDS
# ═══════════════════════════════════════════════════════════════════════

def build_centroids(adata, conf_threshold=CONF_THRESHOLD, min_cells=MIN_CELLS):
    """Build per-subclass mean expression centroids from high-confidence cells.

    Parameters
    ----------
    adata : AnnData
        Combined expression data with raw counts in .X
    conf_threshold : float
        Minimum subclass_label_confidence for inclusion
    min_cells : int
        Minimum cells per subclass

    Returns
    -------
    centroids : pd.DataFrame
        (n_subclasses, n_genes) mean log-normalized expression
    counts : dict
        {subclass: n_cells} for cells used
    gene_names : list
        Gene names in order
    """
    print(f"\nBuilding reference centroids (confidence >= {conf_threshold})...")

    # Filter to high-confidence cells
    conf = adata.obs['subclass_label_confidence'].astype(float).values
    mask = conf >= conf_threshold
    n_pass = mask.sum()
    print(f"  {n_pass:,} / {adata.n_obs:,} cells pass confidence filter "
          f"({100*n_pass/adata.n_obs:.1f}%)")

    adata_hc = adata[mask].copy()

    # Normalize: counts per 10k + log1p
    sc.pp.normalize_total(adata_hc, target_sum=1e4)
    sc.pp.log1p(adata_hc)

    # Get dense matrix
    X = adata_hc.X
    if sparse.issparse(X):
        X = X.toarray()
    X = X.astype(np.float32)

    gene_names = list(adata_hc.var_names)
    labels = adata_hc.obs['subclass_label'].astype(str).values

    # Compute per-subclass means
    unique_labels = sorted(set(labels))
    centroids_dict = {}
    counts = {}

    for lab in unique_labels:
        lab_mask = labels == lab
        n = lab_mask.sum()
        if n < min_cells:
            print(f"  WARNING: {lab} has only {n} cells (< {min_cells}), skipping")
            continue
        centroids_dict[lab] = X[lab_mask].mean(axis=0)
        counts[lab] = n

    centroids = pd.DataFrame(centroids_dict, index=gene_names).T
    print(f"  Built centroids for {len(centroids)} subclasses")
    for lab in sorted(counts):
        print(f"    {lab:20s}: {counts[lab]:>8,} cells")

    return centroids, counts, gene_names


# ═══════════════════════════════════════════════════════════════════════
# 3. WEIGHTED CORRELATION ENGINE
# ═══════════════════════════════════════════════════════════════════════

def weighted_correlate(query_expr, centroids, weights, chunk_size=5000):
    """Compute weighted Pearson correlation between cells and centroids.

    Parameters
    ----------
    query_expr : np.ndarray
        (n_cells, n_genes) log-normalized expression
    centroids : pd.DataFrame
        (n_types, n_genes) centroid expression
    weights : np.ndarray
        (n_genes,) gene weights in [0, 1]
    chunk_size : int
        Process cells in chunks

    Returns
    -------
    corr_matrix : np.ndarray
        (n_cells, n_types)
    type_names : list
    """
    type_names = list(centroids.index)
    centroid_arr = centroids.values.astype(np.float64)  # (n_types, n_genes)
    centroid_arr = np.nan_to_num(centroid_arr, nan=0.0)

    n_cells = query_expr.shape[0]
    n_genes = query_expr.shape[1]
    n_types = len(type_names)

    sqrt_w = np.sqrt(weights).astype(np.float64)  # (n_genes,)

    # Weight centroids
    c_weighted = centroid_arr * sqrt_w[np.newaxis, :]  # (n_types, n_genes)

    # Standardize weighted centroids
    c_mean = c_weighted.mean(axis=1, keepdims=True)
    c_std = c_weighted.std(axis=1, keepdims=True, ddof=0)
    c_std[c_std == 0] = 1.0
    c_norm = (c_weighted - c_mean) / c_std  # (n_types, n_genes)

    corr_matrix = np.zeros((n_cells, n_types), dtype=np.float32)

    for start in range(0, n_cells, chunk_size):
        end = min(start + chunk_size, n_cells)
        chunk = query_expr[start:end].astype(np.float64)

        # Weight query
        q_weighted = chunk * sqrt_w[np.newaxis, :]

        q_mean = q_weighted.mean(axis=1, keepdims=True)
        q_std = q_weighted.std(axis=1, keepdims=True, ddof=0)
        q_std[q_std == 0] = 1.0
        q_norm = (q_weighted - q_mean) / q_std

        q_norm = np.nan_to_num(q_norm, nan=0.0)

        corr_matrix[start:end] = ((q_norm @ c_norm.T) / n_genes).astype(np.float32)

    return corr_matrix, type_names


def assign_labels(corr_matrix, type_names):
    """Assign best-match labels from correlation matrix.

    Returns
    -------
    labels : np.ndarray of str
    best_corr : np.ndarray of float
    margin : np.ndarray of float
        Difference between best and second-best correlation
    """
    sorted_corr = np.sort(corr_matrix, axis=1)
    best_idx = np.argmax(corr_matrix, axis=1)
    labels = np.array([type_names[i] for i in best_idx])
    best_corr = sorted_corr[:, -1]
    second_corr = sorted_corr[:, -2]
    margin = best_corr - second_corr
    return labels, best_corr, margin


# ═══════════════════════════════════════════════════════════════════════
# 4. MAIN CLASSIFICATION PIPELINE
# ═══════════════════════════════════════════════════════════════════════

def run_classification(adata, centroids, gene_names, capture_rates):
    """Run unweighted, linear-weighted, and quadratic-weighted classifiers.

    Returns DataFrame with per-cell results.
    """
    # Normalize query expression
    print("\nNormalizing query expression...")
    adata_q = adata.copy()
    sc.pp.normalize_total(adata_q, target_sum=1e4)
    sc.pp.log1p(adata_q)

    X_query = adata_q.X
    if sparse.issparse(X_query):
        X_query = X_query.toarray()
    X_query = X_query.astype(np.float32)
    X_query = np.nan_to_num(X_query, nan=0.0)

    # Build gene weight arrays
    n_genes = len(gene_names)
    w_uniform = np.ones(n_genes, dtype=np.float64)

    w_linear = np.array([capture_rates.get(g, 50.0) / 100.0 for g in gene_names],
                         dtype=np.float64)

    w_quadratic = w_linear ** 2

    # Print weight summary
    print(f"\nGene weight summary:")
    print(f"  Linear:    mean={w_linear.mean():.3f}, min={w_linear.min():.3f}, "
          f"max={w_linear.max():.3f}")
    print(f"  Quadratic: mean={w_quadratic.mean():.3f}, min={w_quadratic.min():.3f}, "
          f"max={w_quadratic.max():.3f}")

    # Run each classifier
    results = {}
    for name, weights in [("unweighted", w_uniform),
                           ("linear", w_linear),
                           ("quadratic", w_quadratic)]:
        print(f"\n  Running {name} classifier...", flush=True)
        t1 = time.time()
        corr_mat, type_names = weighted_correlate(X_query, centroids, weights)
        labels, best_corr, margin = assign_labels(corr_mat, type_names)
        elapsed = time.time() - t1
        print(f"    Done in {elapsed:.1f}s | Median corr: {np.median(best_corr):.3f}, "
              f"Median margin: {np.median(margin):.4f}")
        results[name] = {
            'labels': labels,
            'best_corr': best_corr,
            'margin': margin,
        }

    # Build output DataFrame
    print("\nBuilding results DataFrame...")
    df = pd.DataFrame({
        'sample_id': adata.obs['sample_id'].values.astype(str),
        'layer': adata.obs['layer'].values.astype(str),
        'depth': adata.obs['predicted_norm_depth'].values.astype(float),
        'hann_subclass': adata.obs['subclass_label'].values.astype(str),
        'hann_confidence': adata.obs['subclass_label_confidence'].values.astype(float),
        'unweighted_subclass': results['unweighted']['labels'],
        'unweighted_corr': results['unweighted']['best_corr'],
        'unweighted_margin': results['unweighted']['margin'],
        'linear_subclass': results['linear']['labels'],
        'linear_corr': results['linear']['best_corr'],
        'linear_margin': results['linear']['margin'],
        'quadratic_subclass': results['quadratic']['labels'],
        'quadratic_corr': results['quadratic']['best_corr'],
        'quadratic_margin': results['quadratic']['margin'],
    })

    return df


# ═══════════════════════════════════════════════════════════════════════
# 5. L6b EVALUATION
# ═══════════════════════════════════════════════════════════════════════

def evaluate_l6b(df):
    """Evaluate L6b misclassification across methods."""
    print("\n" + "=" * 80)
    print("L6b EVALUATION")
    print("=" * 80)

    methods = [
        ('HANN', 'hann_subclass'),
        ('Unweighted corr', 'unweighted_subclass'),
        ('Linear-weighted', 'linear_subclass'),
        ('Quadratic-weighted', 'quadratic_subclass'),
    ]

    summary = []
    for method_name, col in methods:
        l6b_mask = df[col] == 'L6b'
        n_l6b = l6b_mask.sum()
        if n_l6b == 0:
            print(f"\n  {method_name}: 0 L6b cells")
            summary.append({'method': method_name, 'n_l6b': 0,
                            'n_upper': 0, 'pct_upper': 0,
                            'pct_of_total': 0})
            continue

        l6b_upper = l6b_mask & df['layer'].isin(UPPER_LAYERS)
        n_upper = l6b_upper.sum()
        pct_upper = 100 * n_upper / n_l6b
        pct_of_total = 100 * n_l6b / len(df)

        print(f"\n  {method_name}:")
        print(f"    Total L6b: {n_l6b:,} ({pct_of_total:.2f}% of all cells)")
        print(f"    L6b in upper layers: {n_upper:,} ({pct_upper:.1f}%)")

        # Depth stats for L6b cells
        l6b_depths = df.loc[l6b_mask, 'depth']
        print(f"    L6b depth: mean={l6b_depths.mean():.3f}, "
              f"median={l6b_depths.median():.3f}")

        summary.append({
            'method': method_name,
            'n_l6b': n_l6b,
            'n_upper': n_upper,
            'pct_upper': pct_upper,
            'pct_of_total': pct_of_total,
            'mean_depth': l6b_depths.mean(),
        })

    summary_df = pd.DataFrame(summary)
    print(f"\n{summary_df.to_string(index=False)}")

    # What do reclassified HANN-L6b cells become?
    print("\n" + "-" * 60)
    print("What do HANN-L6b cells become under each method?")
    print("-" * 60)

    hann_l6b = df['hann_subclass'] == 'L6b'
    for method_name, col in methods[1:]:  # skip HANN itself
        new_labels = df.loc[hann_l6b, col].value_counts()
        stayed_l6b = new_labels.get('L6b', 0)
        changed = hann_l6b.sum() - stayed_l6b
        print(f"\n  {method_name}: {stayed_l6b:,} stayed L6b, "
              f"{changed:,} changed ({100*changed/hann_l6b.sum():.1f}%)")
        for lab, n in new_labels.head(8).items():
            print(f"    {lab:20s}: {n:>7,} ({100*n/hann_l6b.sum():.1f}%)")

    # Upper-layer HANN-L6b cells specifically
    print("\n" + "-" * 60)
    print("Upper-layer HANN-L6b cells — what do they become?")
    print("-" * 60)
    hann_l6b_upper = hann_l6b & df['layer'].isin(UPPER_LAYERS)
    n_upper = hann_l6b_upper.sum()
    print(f"  {n_upper:,} HANN-L6b cells in upper layers")

    for method_name, col in methods[1:]:
        new_labels = df.loc[hann_l6b_upper, col].value_counts()
        stayed_l6b = new_labels.get('L6b', 0)
        print(f"\n  {method_name}: {stayed_l6b:,}/{n_upper:,} "
              f"({100*stayed_l6b/n_upper:.1f}%) still L6b")
        for lab, n in new_labels.head(8).items():
            print(f"    {lab:20s}: {n:>7,} ({100*n/n_upper:.1f}%)")

    return summary_df


# ═══════════════════════════════════════════════════════════════════════
# 6. VISUALIZATION
# ═══════════════════════════════════════════════════════════════════════

def make_comparison_figure(df, summary_df, capture_rates, gene_names):
    """Generate multi-panel comparison figure."""

    fig, axes = plt.subplots(2, 3, figsize=(24, 14))
    fig.suptitle("Capture-Rate-Weighted Classifier: L6b Evaluation",
                 fontsize=22, fontweight='bold', y=0.98)

    # ── Panel A: % L6b in upper layers by method ──
    ax = axes[0, 0]
    methods = summary_df['method'].values
    pct_upper = summary_df['pct_upper'].values
    colors = ['#888888', '#4a90d9', '#e6a040', '#e05050']
    bars = ax.bar(range(len(methods)), pct_upper, color=colors[:len(methods)],
                  edgecolor='white', width=0.65)
    ax.set_xticks(range(len(methods)))
    ax.set_xticklabels(methods, fontsize=12, rotation=15, ha='right')
    ax.set_ylabel("% L6b in Upper Layers", fontsize=16)
    ax.set_title("L6b in Upper Layers (L1/L2-3/L4)", fontsize=18)
    for bar, val in zip(bars, pct_upper):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                f"{val:.1f}%", ha='center', va='bottom', fontsize=14,
                fontweight='bold')
    ax.tick_params(labelsize=13)
    ax.set_ylim(0, max(pct_upper) * 1.2)

    # ── Panel B: L6b depth distributions ──
    ax = axes[0, 1]
    method_cols = [('HANN', 'hann_subclass', '#888888'),
                   ('Unweighted', 'unweighted_subclass', '#4a90d9'),
                   ('Linear', 'linear_subclass', '#e6a040'),
                   ('Quadratic', 'quadratic_subclass', '#e05050')]
    for mname, col, color in method_cols:
        l6b_depths = df.loc[df[col] == 'L6b', 'depth'].values
        if len(l6b_depths) > 0:
            ax.hist(l6b_depths, bins=50, range=(0, 1.2), alpha=0.5,
                    label=f"{mname} (n={len(l6b_depths):,})", color=color,
                    density=True)
    ax.axvline(0.85, color='white', ls=':', lw=1, alpha=0.5)
    ax.text(0.86, ax.get_ylim()[1]*0.9, 'WM', fontsize=10, color='white', alpha=0.5)
    ax.set_xlabel("Predicted Depth", fontsize=16)
    ax.set_ylabel("Density", fontsize=16)
    ax.set_title("L6b Depth Distribution", fontsize=18)
    ax.legend(fontsize=11)
    ax.tick_params(labelsize=13)

    # ── Panel C: Total L6b proportion ──
    ax = axes[0, 2]
    n_total = len(df)
    proportions = summary_df['pct_of_total'].values
    bars = ax.bar(range(len(methods)), proportions, color=colors[:len(methods)],
                  edgecolor='white', width=0.65)
    # Add MERFISH reference line (1.77% from previous analysis)
    ax.axhline(1.77, color='#2ecc71', ls='--', lw=2, label='MERFISH ref (1.77%)')
    ax.set_xticks(range(len(methods)))
    ax.set_xticklabels(methods, fontsize=12, rotation=15, ha='right')
    ax.set_ylabel("% of All Cells", fontsize=16)
    ax.set_title("Total L6b Proportion", fontsize=18)
    ax.legend(fontsize=12)
    for bar, val in zip(bars, proportions):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
                f"{val:.2f}%", ha='center', va='bottom', fontsize=12)
    ax.tick_params(labelsize=13)

    # ── Panel D: Gene weight distribution ──
    ax = axes[1, 0]
    w_vals = np.array([capture_rates.get(g, 50.0) / 100.0 for g in gene_names])

    # Color by gene category
    colors_genes = []
    for g in gene_names:
        if g.startswith('MT-') or g.startswith('MTRNR'):
            colors_genes.append('#e74c3c')
        else:
            colors_genes.append('#4a90d9')

    ax.hist(w_vals, bins=40, color='#4a90d9', edgecolor='white', alpha=0.8)
    mt_w = [capture_rates.get(g, 50.0)/100.0 for g in gene_names
            if g.startswith('MT-') or g.startswith('MTRNR')]
    if mt_w:
        ax.hist(mt_w, bins=40, color='#e74c3c', edgecolor='white', alpha=0.8,
                label=f'MT genes (n={len(mt_w)})')
    ax.set_xlabel("Gene Weight (linear = pct_inside/100)", fontsize=16)
    ax.set_ylabel("Number of Genes", fontsize=16)
    ax.set_title("Gene Weight Distribution", fontsize=18)
    ax.legend(fontsize=12)
    ax.tick_params(labelsize=13)

    # ── Panel E: Reclassification of HANN-L6b cells (linear) ──
    ax = axes[1, 1]
    hann_l6b = df['hann_subclass'] == 'L6b'
    new_labels = df.loc[hann_l6b, 'linear_subclass'].value_counts().head(10)
    y_pos = np.arange(len(new_labels))
    bar_colors = ['#2ecc71' if lab == 'L6b' else '#e74c3c' for lab in new_labels.index]
    ax.barh(y_pos, new_labels.values, color=bar_colors, edgecolor='white')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(new_labels.index, fontsize=12)
    ax.set_xlabel("Number of Cells", fontsize=16)
    ax.set_title("HANN-L6b → Linear-Weighted Labels", fontsize=18)
    ax.invert_yaxis()
    ax.tick_params(labelsize=13)
    for i, (lab, n) in enumerate(new_labels.items()):
        ax.text(n + max(new_labels)*0.01, i,
                f"{100*n/hann_l6b.sum():.1f}%", va='center', fontsize=11)

    # ── Panel F: Per-sample L6b upper-layer comparison ──
    ax = axes[1, 2]
    samples = sorted(df['sample_id'].unique())
    hann_rates = []
    linear_rates = []
    for sid in samples:
        smask = df['sample_id'] == sid
        hann_l6b_s = smask & (df['hann_subclass'] == 'L6b')
        lin_l6b_s = smask & (df['linear_subclass'] == 'L6b')

        if hann_l6b_s.sum() > 0:
            h_upper = (hann_l6b_s & df['layer'].isin(UPPER_LAYERS)).sum()
            hann_rates.append(100 * h_upper / hann_l6b_s.sum())
        else:
            hann_rates.append(0)

        if lin_l6b_s.sum() > 0:
            l_upper = (lin_l6b_s & df['layer'].isin(UPPER_LAYERS)).sum()
            linear_rates.append(100 * l_upper / lin_l6b_s.sum())
        else:
            linear_rates.append(0)

    ax.scatter(hann_rates, linear_rates, s=60, c='#4a90d9', edgecolors='white',
               zorder=3)
    max_val = max(max(hann_rates), max(linear_rates)) * 1.1
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.3, lw=1)
    ax.set_xlabel("HANN: % L6b in Upper Layers", fontsize=16)
    ax.set_ylabel("Linear-Weighted: % L6b in Upper Layers", fontsize=16)
    ax.set_title("Per-Sample Comparison", fontsize=18)
    ax.tick_params(labelsize=13)
    # Count improvements
    n_improved = sum(1 for h, l in zip(hann_rates, linear_rates) if l < h)
    ax.text(0.05, 0.95, f"{n_improved}/{len(samples)} samples improved",
            transform=ax.transAxes, fontsize=13, va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out_path = os.path.join(OUTPUT_DIR, "weighted_classifier_comparison.png")
    plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"\nFigure saved: {out_path}")


# ═══════════════════════════════════════════════════════════════════════
# 7. MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    t0 = time.time()

    # Load data
    capture_rates = load_capture_rates()
    adata = load_all_samples()
    gene_names = list(adata.var_names)

    # Build centroids
    centroids, counts, gene_names = build_centroids(adata)

    # Run classifiers
    df = run_classification(adata, centroids, gene_names, capture_rates)

    # Evaluate L6b
    summary_df = evaluate_l6b(df)

    # Visualize
    make_comparison_figure(df, summary_df, capture_rates, gene_names)

    # Save results
    csv_path = os.path.join(OUTPUT_DIR, "weighted_classifier_results.csv")
    # Save a subset to keep file size manageable
    df_save = df[['sample_id', 'layer', 'depth', 'hann_subclass', 'hann_confidence',
                   'linear_subclass', 'linear_corr', 'quadratic_subclass',
                   'quadratic_corr']].copy()
    df_save.to_csv(csv_path, index=False)
    print(f"Results saved: {csv_path}")

    summary_path = os.path.join(OUTPUT_DIR, "weighted_classifier_l6b_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    print(f"L6b summary saved: {summary_path}")

    elapsed = time.time() - t0
    print(f"\nTotal elapsed: {elapsed:.0f}s ({elapsed/60:.1f} min)")


if __name__ == "__main__":
    main()
