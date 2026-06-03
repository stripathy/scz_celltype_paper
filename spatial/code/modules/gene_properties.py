"""
Gene property classification and quality filtering utilities.

Used by: characterize_poor_genes.py, hierarchical_probe_selection.py
"""

import re
import numpy as np
import pandas as pd
import scipy.sparse as sp


def classify_gene_biotype(gene_name):
    """Infer gene biotype from naming conventions.

    Parameters
    ----------
    gene_name : str
        Gene symbol.

    Returns
    -------
    str : one of 'protein_coding', 'antisense', 'lncRNA/novel',
          'unannotated', 'other'
    """
    if re.match(r'^AC\d+\.\d+$', gene_name) or re.match(r'^AL\d+\.\d+$', gene_name):
        return "lncRNA/novel"
    if gene_name.endswith('-AS1') or gene_name.endswith('-AS2'):
        return "antisense"
    if gene_name.startswith('LINC') or gene_name.startswith('MIR'):
        return "lncRNA/novel"
    if gene_name.startswith('LOC') or gene_name.startswith('ENSG'):
        return "unannotated"
    if re.match(r'^[A-Z][A-Z0-9]+$', gene_name):
        return "protein_coding"
    if re.match(r'^[A-Za-z]', gene_name):
        return "protein_coding"
    return "other"


def filter_eligible_genes(gene_quality_df, all_gene_names,
                          min_pearson_r_quantile=0.2,
                          exclude_biotypes=("antisense", "lncRNA/novel"),
                          min_detection_rate=0.005):
    """Build set of eligible genes for probe selection.

    Genes are excluded if they:
    - Fall in the bottom quintile of cross-platform correlation
      (but only genes PRESENT in gene_quality_df are filtered this way;
      genes not in the 4K panel are eligible by default)
    - Are antisense/lncRNA biotypes
    - Have detection rate below threshold

    Parameters
    ----------
    gene_quality_df : DataFrame
        Gene properties with pearson_r, biotype columns (from characterize_poor_genes.py).
    all_gene_names : array-like
        Full set of gene names to filter from.
    min_pearson_r_quantile : float
        Bottom quantile to exclude (default 0.2 = bottom 20%).
    exclude_biotypes : tuple of str
        Biotypes to exclude.
    min_detection_rate : float
        Minimum snRNAseq detection rate.

    Returns
    -------
    set : eligible gene names
    set : excluded gene names with reasons
    """
    eligible = set(all_gene_names)
    excluded = {}

    # 1. Filter by cross-platform correlation (only for genes in the table)
    if gene_quality_df is not None and 'pearson_r' in gene_quality_df.columns:
        r_threshold = gene_quality_df['pearson_r'].quantile(min_pearson_r_quantile)
        poor_genes = set(gene_quality_df[
            gene_quality_df['pearson_r'] < r_threshold
        ].index)
        for g in poor_genes:
            if g in eligible:
                eligible.discard(g)
                excluded[g] = f"poor_cross_platform_r (< {r_threshold:.3f})"

    # 2. Filter by biotype
    for g in list(eligible):
        bt = classify_gene_biotype(g)
        if bt in exclude_biotypes:
            eligible.discard(g)
            excluded[g] = f"biotype_{bt}"

    # 3. Filter by detection rate (if available in quality table)
    if gene_quality_df is not None and 'det_rate_sn' in gene_quality_df.columns:
        for g in list(eligible):
            if g in gene_quality_df.index:
                det = gene_quality_df.loc[g, 'det_rate_sn']
                if pd.notna(det) and det < min_detection_rate:
                    eligible.discard(g)
                    excluded[g] = f"low_detection_rate ({det:.4f})"

    n_total = len(all_gene_names)
    n_excluded = n_total - len(eligible)
    print(f"  Gene quality filter: {n_total} -> {len(eligible)} eligible "
          f"({n_excluded} excluded)")

    return eligible, excluded


def compute_detection_rate(adata, gene_names):
    """Compute detection rate (fraction of cells with expression > 0) per gene.

    Parameters
    ----------
    adata : AnnData
        Expression data (raw or normalized, either works for > 0 check).
    gene_names : list of str
        Genes to compute detection rate for.

    Returns
    -------
    pd.Series : gene_name -> fraction detected
    """
    gene_idx_map = {g: i for i, g in enumerate(adata.var_names)}
    valid_genes = [g for g in gene_names if g in gene_idx_map]
    idx = [gene_idx_map[g] for g in valid_genes]

    X = adata.X[:, idx]
    if sp.issparse(X):
        det_rate = np.array((X > 0).mean(axis=0)).flatten()
    else:
        det_rate = np.mean(X > 0, axis=0)

    return pd.Series(det_rate, index=valid_genes)


def compute_specificity(pb_df):
    """Compute cell-type specificity metrics for each gene.

    Parameters
    ----------
    pb_df : DataFrame
        Pseudobulk expression (rows = cell types, columns = genes).

    Returns
    -------
    DataFrame with columns: entropy, max_over_mean, gini, top_celltype, cv
    """
    results = {}
    for gene in pb_df.columns:
        vals = pb_df[gene].values
        vals_pos = np.maximum(vals, 0)
        total = vals_pos.sum()

        if total < 1e-10:
            results[gene] = {
                'entropy': 0, 'max_over_mean': 0, 'gini': 0,
                'top_celltype': 'none', 'cv': 0
            }
            continue

        # Shannon entropy
        probs = vals_pos / total
        probs = probs[probs > 0]
        entropy = -np.sum(probs * np.log2(probs))

        # Max/mean ratio
        mean_val = np.mean(vals_pos)
        max_val = np.max(vals_pos)
        max_over_mean = max_val / mean_val if mean_val > 0 else 0

        # Gini coefficient
        sorted_vals = np.sort(vals_pos)
        n = len(sorted_vals)
        cumvals = np.cumsum(sorted_vals)
        gini = ((n + 1 - 2 * np.sum(cumvals) / cumvals[-1]) / n
                if cumvals[-1] > 0 else 0)

        # CV
        cv = np.std(vals_pos) / mean_val if mean_val > 0 else 0

        # Top cell type
        top_ct = pb_df.index[np.argmax(vals)]

        results[gene] = {
            'entropy': entropy, 'max_over_mean': max_over_mean,
            'gini': gini, 'top_celltype': top_ct, 'cv': cv
        }

    return pd.DataFrame(results).T
