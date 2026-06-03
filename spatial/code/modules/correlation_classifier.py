"""
Two-stage hierarchical correlation classifier for Xenium cell types.

Reclassifies cells using Pearson correlation against self-built centroids
from high-confidence HANN-labeled exemplar cells. Dramatically reduces
misclassifications such as L6b appearing in upper cortical layers.

Algorithm:
  Stage 1: Classify into 24 subclasses using top-100 HANN exemplars per subclass
  Stage 2: Within each Stage-1 subclass, classify into supertypes using
           top-100 HANN exemplars per supertype

The centroids are built from the Xenium data itself (not from reference data),
using the cells with highest HANN confidence as reference exemplars.

Also includes spatial doublet detection (Section 5), validated against
Nicole's Sea-AD snRNAseq reference (137K cells x 36.6K genes). See
data/reference/nicole_sea_ad_snrnaseq_reference.h5ad.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


# ═══════════════════════════════════════════════════════════════════════
# 1. CENTROID BUILDING
# ═══════════════════════════════════════════════════════════════════════

def build_subclass_centroids(adata, top_n=100,
                              subclass_col='subclass_label',
                              confidence_col='subclass_label_confidence'):
    """Build per-subclass mean expression centroids from top-N exemplar cells.

    For each subclass, selects the top_n cells with highest HANN confidence,
    normalizes (counts per 10k + log1p), and computes mean expression.

    Parameters
    ----------
    adata : AnnData
        Combined expression data with raw counts in .X
    top_n : int
        Number of exemplar cells per subclass (uses all if < top_n available)
    subclass_col : str
        Column name for subclass labels in adata.obs
    confidence_col : str
        Column name for subclass confidence in adata.obs

    Returns
    -------
    centroids : pd.DataFrame
        (n_subclasses, n_genes) mean log-normalized expression
    cell_counts : dict
        {subclass: n_cells_used}
    gene_names : list
        Gene names in column order
    """
    labels = adata.obs[subclass_col].astype(str).values
    confidences = adata.obs[confidence_col].astype(float).values
    unique_labels = sorted(set(labels))

    # Select top-N exemplar cells for each subclass
    exemplar_indices = []
    cell_counts = {}

    for lab in unique_labels:
        lab_mask = np.where(labels == lab)[0]
        n_available = len(lab_mask)
        if n_available == 0:
            continue

        n_use = min(top_n, n_available)
        lab_conf = confidences[lab_mask]
        top_idx = lab_mask[np.argsort(lab_conf)[-n_use:]]
        exemplar_indices.append(top_idx)
        cell_counts[lab] = n_use

    all_idx = np.concatenate(exemplar_indices)
    adata_ex = adata[all_idx].copy()

    # Normalize: counts per 10k + log1p
    sc.pp.normalize_total(adata_ex, target_sum=1e4)
    sc.pp.log1p(adata_ex)

    # Get dense matrix
    X = adata_ex.X
    if sparse.issparse(X):
        X = X.toarray()
    X = X.astype(np.float32)

    gene_names = list(adata_ex.var_names)
    ex_labels = adata_ex.obs[subclass_col].astype(str).values

    # Compute per-subclass means
    centroids_dict = {}
    for lab in sorted(cell_counts.keys()):
        lab_mask = ex_labels == lab
        centroids_dict[lab] = X[lab_mask].mean(axis=0)

    centroids = pd.DataFrame(centroids_dict, index=gene_names).T

    print(f"  Built {len(centroids)} subclass centroids (top-{top_n} exemplars each)")
    for lab in sorted(cell_counts):
        print(f"    {lab:20s}: {cell_counts[lab]:>5,} cells")

    return centroids, cell_counts, gene_names


def build_supertype_centroids(adata, top_n=100,
                               subclass_col='subclass_label',
                               supertype_col='supertype_label',
                               confidence_col='supertype_label_confidence'):
    """Build per-supertype centroids grouped by parent subclass.

    For each supertype, selects the top_n cells with highest HANN supertype
    confidence, normalizes, and computes mean expression.

    Parameters
    ----------
    adata : AnnData
        Combined expression data with raw counts in .X
    top_n : int
        Number of exemplar cells per supertype
    subclass_col, supertype_col, confidence_col : str
        Column names in adata.obs

    Returns
    -------
    subclass_to_supertype_centroids : dict
        {subclass_name: pd.DataFrame(n_supertypes, n_genes)}
    supertype_to_subclass : dict
        {supertype_name: subclass_name} mapping
    """
    labels_sub = adata.obs[subclass_col].astype(str).values
    labels_sup = adata.obs[supertype_col].astype(str).values
    confidences = adata.obs[confidence_col].astype(float).values
    gene_names = list(adata.var_names)

    # Discover supertype -> subclass mapping from data
    supertype_to_subclass = {}
    for sup, sub in zip(labels_sup, labels_sub):
        if sup not in supertype_to_subclass:
            supertype_to_subclass[sup] = sub

    # Group supertypes by subclass
    from collections import defaultdict
    subclass_to_supertypes = defaultdict(list)
    for sup, sub in sorted(supertype_to_subclass.items()):
        subclass_to_supertypes[sub].append(sup)

    # Build centroids per subclass
    result = {}
    total_types = 0

    for sub in sorted(subclass_to_supertypes.keys()):
        sups = sorted(subclass_to_supertypes[sub])

        # Select exemplar cells for each supertype in this subclass
        exemplar_indices = []
        sup_counts = {}

        for sup in sups:
            sup_mask = np.where(labels_sup == sup)[0]
            if len(sup_mask) == 0:
                continue
            n_use = min(top_n, len(sup_mask))
            sup_conf = confidences[sup_mask]
            top_idx = sup_mask[np.argsort(sup_conf)[-n_use:]]
            exemplar_indices.append(top_idx)
            sup_counts[sup] = n_use

        if not exemplar_indices:
            continue

        all_idx = np.concatenate(exemplar_indices)
        adata_ex = adata[all_idx].copy()

        # Normalize
        sc.pp.normalize_total(adata_ex, target_sum=1e4)
        sc.pp.log1p(adata_ex)

        X = adata_ex.X
        if sparse.issparse(X):
            X = X.toarray()
        X = X.astype(np.float32)

        ex_sup_labels = adata_ex.obs[supertype_col].astype(str).values

        # Compute per-supertype means
        centroids_dict = {}
        for sup in sorted(sup_counts.keys()):
            sup_mask = ex_sup_labels == sup
            centroids_dict[sup] = X[sup_mask].mean(axis=0)

        result[sub] = pd.DataFrame(centroids_dict, index=gene_names).T
        total_types += len(centroids_dict)

    print(f"  Built {total_types} supertype centroids across "
          f"{len(result)} subclasses (top-{top_n} exemplars each)")

    return result, supertype_to_subclass


def _normalize_adata(adata_sub, norm_method='log_cp10k'):
    """Normalize an AnnData object in place using the specified method.

    Parameters
    ----------
    adata_sub : AnnData
        Data to normalize (modified in place).
    norm_method : str
        'log_cp10k' — counts per 10k + log1p (our default, scanpy standard)
        'log2_cpm'  — counts per million + log2(x+1) (SEA-AD / scrattch.mapping)
    """
    if norm_method == 'log2_cpm':
        sc.pp.normalize_total(adata_sub, target_sum=1e6)
        X = adata_sub.X
        if sparse.issparse(X):
            X = X.toarray()
        adata_sub.X = np.log2(X + 1).astype(np.float32)
    else:  # log_cp10k (default)
        sc.pp.normalize_total(adata_sub, target_sum=1e4)
        sc.pp.log1p(adata_sub)


def _filter_blank_genes(adata):
    """Return adata with Blank control probes removed."""
    keep = ~adata.var_names.str.startswith('Blank')
    n_removed = (~keep).sum()
    if n_removed > 0:
        print(f"    Excluding {n_removed} Blank probes ({keep.sum()} genes retained)")
    return adata[:, keep].copy()


def build_flat_centroids(adata, label_col, confidence_col, subclass_col,
                          top_n=None, agg_func='median',
                          norm_method='log_cp10k', exclude_blank=False):
    """Build flat centroids for ALL supertypes in a single DataFrame.

    Unlike the hierarchical approach, this builds one large centroid matrix
    containing all supertypes at once (matching SEA-AD's scrattch.mapping::corrMap).

    Parameters
    ----------
    adata : AnnData
        Expression data with raw counts in .X
    label_col : str
        Column for supertype labels (e.g., 'hann_supertype')
    confidence_col : str
        Column for supertype confidence (e.g., 'hann_supertype_confidence')
    subclass_col : str
        Column for subclass labels (e.g., 'hann_subclass') — used to build lookup
    top_n : int or None
        Number of exemplar cells per supertype (None = use all cells)
    agg_func : str
        'median' (SEA-AD default) or 'mean'
    norm_method : str
        'log_cp10k' or 'log2_cpm' (SEA-AD default)
    exclude_blank : bool
        If True, remove genes starting with 'Blank'

    Returns
    -------
    centroids : pd.DataFrame
        (n_supertypes, n_genes) centroid expression
    type_to_subclass : dict
        {supertype_name: subclass_name} lookup
    gene_names : list
        Gene names in column order
    cell_counts : dict
        {supertype: n_cells_used}
    """
    # Optionally filter Blank probes
    if exclude_blank:
        adata = _filter_blank_genes(adata)

    labels = adata.obs[label_col].astype(str).values
    confidences = adata.obs[confidence_col].astype(float).values
    subclass_labels = adata.obs[subclass_col].astype(str).values
    unique_labels = sorted(set(labels))

    # Build supertype -> subclass mapping
    type_to_subclass = {}
    for sup, sub in zip(labels, subclass_labels):
        if sup not in type_to_subclass:
            type_to_subclass[sup] = sub

    # Select exemplar cells for each supertype
    exemplar_indices = []
    cell_counts = {}

    for lab in unique_labels:
        lab_mask = np.where(labels == lab)[0]
        n_available = len(lab_mask)
        if n_available == 0:
            continue

        if top_n is None:
            n_use = n_available
        else:
            n_use = min(top_n, n_available)
        lab_conf = confidences[lab_mask]
        top_idx = lab_mask[np.argsort(lab_conf)[-n_use:]]
        exemplar_indices.append(top_idx)
        cell_counts[lab] = n_use

    all_idx = np.concatenate(exemplar_indices)
    adata_ex = adata[all_idx].copy()

    # Normalize using specified method
    _normalize_adata(adata_ex, norm_method)

    # Get dense matrix
    X = adata_ex.X
    if sparse.issparse(X):
        X = X.toarray()
    X = X.astype(np.float32)

    gene_names = list(adata_ex.var_names)
    ex_labels = adata_ex.obs[label_col].astype(str).values

    # Compute per-supertype centroids
    agg_fn = np.median if agg_func == 'median' else np.mean
    centroids_dict = {}
    for lab in sorted(cell_counts.keys()):
        lab_mask = ex_labels == lab
        centroids_dict[lab] = agg_fn(X[lab_mask], axis=0)

    centroids = pd.DataFrame(centroids_dict, index=gene_names).T

    top_label = "all" if top_n is None else f"top-{top_n}"
    print(f"  Built {len(centroids)} flat supertype centroids "
          f"({top_label} exemplars, {agg_func}, {norm_method}, "
          f"{'excl blanks' if exclude_blank else 'incl blanks'})")

    return centroids, type_to_subclass, gene_names, cell_counts


def run_flat_classifier(adata, supertype_centroids, type_to_subclass, gene_names,
                         norm_method='log_cp10k'):
    """Run flat (non-hierarchical) correlation classifier.

    Maps all cells directly to all supertypes at once (like SEA-AD corrMap),
    then derives subclass from the type_to_subclass lookup.

    Parameters
    ----------
    adata : AnnData
        Expression data with raw counts in .X
    supertype_centroids : pd.DataFrame
        (n_supertypes, n_genes) from build_flat_centroids()
    type_to_subclass : dict
        {supertype_name: subclass_name} mapping
    gene_names : list
        Ordered gene names matching centroid columns
    norm_method : str
        'log_cp10k' or 'log2_cpm' — must match centroid normalization

    Returns
    -------
    pd.DataFrame with columns matching run_two_stage_classifier output:
        corr_subclass, corr_subclass_corr, corr_subclass_margin,
        corr_supertype, corr_supertype_corr
    """
    import time

    # Filter adata to match gene set
    common_genes = [g for g in gene_names if g in adata.var_names]
    if len(common_genes) < len(gene_names):
        print(f"  Warning: {len(gene_names) - len(common_genes)} genes not in query")
    adata_q = adata[:, common_genes].copy()
    supertype_centroids = supertype_centroids[common_genes]

    # Normalize query expression
    print("  Normalizing query expression...", flush=True)
    _normalize_adata(adata_q, norm_method)

    X_query = adata_q.X
    if sparse.issparse(X_query):
        X_query = X_query.toarray()
    X_query = X_query.astype(np.float32)
    X_query = np.nan_to_num(X_query, nan=0.0)

    n_cells = X_query.shape[0]

    # ── Correlate against all supertypes at once ──
    print(f"  Flat classification: {n_cells:,} cells × "
          f"{len(supertype_centroids)} supertypes...", flush=True)
    t0 = time.time()
    corr_matrix, type_names = correlate(X_query, supertype_centroids)

    # Best supertype assignment
    sup_labels, sup_corr, sup_margin = assign_labels(corr_matrix, type_names)
    print(f"    Done in {time.time()-t0:.1f}s | "
          f"Median corr: {np.median(sup_corr):.3f}")

    # Derive subclass from supertype
    sub_labels = np.array([type_to_subclass.get(s, 'Unknown') for s in sup_labels])

    # ── Compute subclass-level margin ──
    # For each cell: group supertype correlations by parent subclass,
    # take max within each subclass, then margin = best - second_best
    unique_subclasses = sorted(set(type_to_subclass.values()))
    subclass_idx_map = {}  # {subclass: list of supertype column indices}
    for i, stype in enumerate(type_names):
        sc_name = type_to_subclass.get(stype, 'Unknown')
        if sc_name not in subclass_idx_map:
            subclass_idx_map[sc_name] = []
        subclass_idx_map[sc_name].append(i)

    n_subclasses = len(subclass_idx_map)
    subclass_names = sorted(subclass_idx_map.keys())
    # Max supertype correlation per subclass
    subclass_max_corr = np.full((n_cells, n_subclasses), -np.inf, dtype=np.float32)
    for j, sc_name in enumerate(subclass_names):
        sc_cols = subclass_idx_map[sc_name]
        subclass_max_corr[:, j] = corr_matrix[:, sc_cols].max(axis=1)

    # Subclass-level: best correlation and margin
    sorted_sc_corr = np.sort(subclass_max_corr, axis=1)
    best_sc_idx = np.argmax(subclass_max_corr, axis=1)
    sub_corr_vals = sorted_sc_corr[:, -1]
    second_sc_corr = sorted_sc_corr[:, -2] if n_subclasses > 1 else np.zeros_like(sub_corr_vals)
    sub_margin = sub_corr_vals - second_sc_corr

    # Build results DataFrame
    results = pd.DataFrame({
        'corr_subclass': sub_labels,
        'corr_subclass_corr': sub_corr_vals,
        'corr_subclass_margin': sub_margin,
        'corr_supertype': sup_labels.astype(str),
        'corr_supertype_corr': sup_corr,
    }, index=adata.obs.index)

    return results


# ═══════════════════════════════════════════════════════════════════════
# 2. CORRELATION ENGINE
# ═══════════════════════════════════════════════════════════════════════

def correlate(query_expr, centroids, chunk_size=10000):
    """Compute Pearson correlation between cells and centroids.

    Standard (unweighted) Pearson correlation. Z-scores each row,
    then dot product / n_genes.

    Parameters
    ----------
    query_expr : np.ndarray
        (n_cells, n_genes) log-normalized expression
    centroids : pd.DataFrame
        (n_types, n_genes) centroid expression
    chunk_size : int
        Process cells in chunks for memory efficiency

    Returns
    -------
    corr_matrix : np.ndarray
        (n_cells, n_types) Pearson correlation
    type_names : list
        Type names matching columns of corr_matrix
    """
    type_names = list(centroids.index)
    centroid_arr = centroids.values.astype(np.float64)
    centroid_arr = np.nan_to_num(centroid_arr, nan=0.0)

    n_cells = query_expr.shape[0]
    n_genes = query_expr.shape[1]
    n_types = len(type_names)

    # Standardize centroids
    c_mean = centroid_arr.mean(axis=1, keepdims=True)
    c_std = centroid_arr.std(axis=1, keepdims=True, ddof=0)
    c_std[c_std == 0] = 1.0
    c_norm = (centroid_arr - c_mean) / c_std  # (n_types, n_genes)

    corr_matrix = np.zeros((n_cells, n_types), dtype=np.float32)

    for start in range(0, n_cells, chunk_size):
        end = min(start + chunk_size, n_cells)
        chunk = query_expr[start:end].astype(np.float64)

        # Standardize query chunk
        q_mean = chunk.mean(axis=1, keepdims=True)
        q_std = chunk.std(axis=1, keepdims=True, ddof=0)
        q_std[q_std == 0] = 1.0
        q_norm = (chunk - q_mean) / q_std

        q_norm = np.nan_to_num(q_norm, nan=0.0)

        corr_matrix[start:end] = ((q_norm @ c_norm.T) / n_genes).astype(np.float32)

    return corr_matrix, type_names


def assign_labels(corr_matrix, type_names):
    """Assign best-match labels from correlation matrix.

    Returns
    -------
    labels : np.ndarray of str
    best_corr : np.ndarray of float32
    margin : np.ndarray of float32
        Difference between best and second-best correlation
    """
    sorted_corr = np.sort(corr_matrix, axis=1)
    best_idx = np.argmax(corr_matrix, axis=1)
    labels = np.array([type_names[i] for i in best_idx])
    best_corr = sorted_corr[:, -1]
    second_corr = sorted_corr[:, -2] if corr_matrix.shape[1] > 1 else np.zeros_like(best_corr)
    margin = best_corr - second_corr
    return labels, best_corr.astype(np.float32), margin.astype(np.float32)


# ═══════════════════════════════════════════════════════════════════════
# 3. TWO-STAGE CLASSIFIER
# ═══════════════════════════════════════════════════════════════════════

def run_two_stage_classifier(adata, subclass_centroids,
                              supertype_centroids_by_subclass, gene_names):
    """Run the complete two-stage hierarchical correlation classifier.

    Stage 1: Classify all cells against 24 subclass centroids.
    Stage 2: For each cell, classify against ONLY the supertype centroids
             within its Stage 1 subclass.

    Parameters
    ----------
    adata : AnnData
        Expression data with raw counts in .X (will be normalized internally)
    subclass_centroids : pd.DataFrame
        From build_subclass_centroids()
    supertype_centroids_by_subclass : dict
        From build_supertype_centroids()
    gene_names : list
        Ordered gene names matching centroid columns

    Returns
    -------
    pd.DataFrame with columns:
        corr_subclass, corr_subclass_corr, corr_subclass_margin,
        corr_supertype, corr_supertype_corr
    """
    import time

    # Normalize query expression (once)
    print("  Normalizing query expression...")
    adata_q = adata.copy()
    sc.pp.normalize_total(adata_q, target_sum=1e4)
    sc.pp.log1p(adata_q)

    X_query = adata_q.X
    if sparse.issparse(X_query):
        X_query = X_query.toarray()
    X_query = X_query.astype(np.float32)
    X_query = np.nan_to_num(X_query, nan=0.0)

    n_cells = X_query.shape[0]

    # ── Stage 1: Subclass classification ──
    print("  Stage 1: Subclass classification...", flush=True)
    t1 = time.time()
    corr_mat_sub, type_names_sub = correlate(X_query, subclass_centroids)
    sub_labels, sub_corr, sub_margin = assign_labels(corr_mat_sub, type_names_sub)
    print(f"    Done in {time.time()-t1:.1f}s | "
          f"Median corr: {np.median(sub_corr):.3f}, "
          f"Median margin: {np.median(sub_margin):.4f}")

    # ── Stage 2: Supertype classification within each subclass ──
    print("  Stage 2: Supertype classification (within subclass)...", flush=True)
    t2 = time.time()

    sup_labels = np.full(n_cells, '', dtype=object)
    sup_corr = np.zeros(n_cells, dtype=np.float32)

    unique_subs = sorted(set(sub_labels))
    for sub in unique_subs:
        sub_mask = sub_labels == sub
        n_sub_cells = sub_mask.sum()

        if sub not in supertype_centroids_by_subclass:
            # No supertype centroids for this subclass — assign subclass as supertype
            sup_labels[sub_mask] = sub
            sup_corr[sub_mask] = sub_corr[sub_mask]
            continue

        sup_centroids = supertype_centroids_by_subclass[sub]

        if len(sup_centroids) == 1:
            # Only one supertype — assign directly
            sup_labels[sub_mask] = sup_centroids.index[0]
            sup_corr[sub_mask] = sub_corr[sub_mask]
            continue

        # Correlate subset against supertype centroids
        sub_indices = np.where(sub_mask)[0]
        X_sub = X_query[sub_indices]
        corr_mat_sup, type_names_sup = correlate(X_sub, sup_centroids)
        s_labels, s_corr, _ = assign_labels(corr_mat_sup, type_names_sup)

        sup_labels[sub_indices] = s_labels
        sup_corr[sub_indices] = s_corr

    print(f"    Done in {time.time()-t2:.1f}s | "
          f"{len(unique_subs)} subclasses processed")

    # Build results DataFrame
    results = pd.DataFrame({
        'corr_subclass': sub_labels,
        'corr_subclass_corr': sub_corr,
        'corr_subclass_margin': sub_margin,
        'corr_supertype': sup_labels.astype(str),
        'corr_supertype_corr': sup_corr,
    }, index=adata.obs.index)

    return results


# ═══════════════════════════════════════════════════════════════════════
# 4. QC FLAGGING
# ═══════════════════════════════════════════════════════════════════════

def flag_low_margin_cells(margins, sample_ids, percentile=1.0,
                          subclass_labels=None, l6b_margin_threshold=None):
    """Flag bottom percentile of cells by Stage 1 subclass margin, per sample.

    Parameters
    ----------
    margins : np.ndarray
        Per-cell correlation margin (from Stage 1)
    sample_ids : np.ndarray
        Per-cell sample ID
    percentile : float
        Bottom percentage to flag (e.g., 5.0 = bottom 5%)
    subclass_labels : np.ndarray, optional
        Per-cell subclass labels (e.g., corr_subclass). Required for
        cell-type-specific thresholds.
    l6b_margin_threshold : float, optional
        Absolute margin threshold for L6b cells. L6b cells with
        margin < this value are flagged regardless of percentile.

    Returns
    -------
    corr_qc_pass : np.ndarray of bool
        True = passes QC, False = flagged as low-confidence
    thresholds : dict
        {sample_id: margin_threshold} for reporting
    """
    corr_qc_pass = np.ones(len(margins), dtype=bool)
    thresholds = {}

    for sid in sorted(set(sample_ids)):
        sid_mask = sample_ids == sid
        sid_margins = margins[sid_mask]
        threshold = np.percentile(sid_margins, percentile)
        fail_mask = sid_mask & (margins < threshold)
        corr_qc_pass[fail_mask] = False
        thresholds[sid] = threshold

    n_pctl_fail = (~corr_qc_pass).sum()
    print(f"  Percentile QC: {n_pctl_fail:,} cells flagged ({100*n_pctl_fail/len(margins):.1f}%)")

    # Additional L6b-specific absolute margin filter
    if subclass_labels is not None and l6b_margin_threshold is not None:
        l6b_mask = subclass_labels == 'L6b'
        l6b_low = l6b_mask & (margins < l6b_margin_threshold) & corr_qc_pass
        n_l6b_extra = l6b_low.sum()
        corr_qc_pass[l6b_low] = False
        n_l6b_total = l6b_mask.sum()
        print(f"  L6b margin filter: {n_l6b_extra:,} additional L6b cells flagged "
              f"(margin < {l6b_margin_threshold}, {n_l6b_total:,} L6b total)")

    n_fail = (~corr_qc_pass).sum()
    print(f"  QC flagging total: {n_fail:,} cells flagged ({100*n_fail/len(margins):.1f}%)")

    return corr_qc_pass, thresholds


# ═══════════════════════════════════════════════════════════════════════
# 5. SPATIAL DOUBLET DETECTION
# ═══════════════════════════════════════════════════════════════════════

# ── Marker sets for doublet detection ──
# Panel-specific marker lists for spatial doublet detection.
# Each panel defines: GABAergic markers, Glutamatergic markers, and
# mutually exclusive interneuron subclass markers for within-GABA doublets.

# Xenium 300-gene panel (default):
# - 7 GABA markers, 4 Glut markers, 3 within-GABA triple
# - Validated FP rate: 0.098% in snRNAseq at gaba_threshold=4
DOUBLET_GABA_MARKERS_XENIUM = ['GAD1', 'GAD2', 'SLC32A1', 'SST', 'PVALB', 'VIP', 'LAMP5']
DOUBLET_GLUT_MARKERS_XENIUM = ['CUX2', 'RORB', 'GRIN2A', 'THEMIS']
DOUBLET_GABA_TRIPLE_XENIUM = ['SST', 'PVALB', 'LAMP5']

# MERFISH 180-gene panel:
# - GAD1 and SST are MISSING from MERFISH panel
# - 5 GABA markers (fewer → lower threshold needed)
# - SST missing → within-GABA uses PVALB+LAMP5 co-expression only (2-marker)
DOUBLET_GABA_MARKERS_MERFISH = ['GAD2', 'SLC32A1', 'PVALB', 'VIP', 'LAMP5']
DOUBLET_GLUT_MARKERS_MERFISH = ['CUX2', 'RORB', 'GRIN2A', 'THEMIS']
DOUBLET_GABA_TRIPLE_MERFISH = ['PVALB', 'LAMP5']  # SST missing → 2-marker

# Backward-compatible aliases (Xenium defaults)
DOUBLET_GABA_MARKERS = DOUBLET_GABA_MARKERS_XENIUM
DOUBLET_GLUT_MARKERS = DOUBLET_GLUT_MARKERS_XENIUM
DOUBLET_GABA_SUBCLASS_TRIPLE = DOUBLET_GABA_TRIPLE_XENIUM

# Default thresholds per panel
DOUBLET_THRESHOLD_XENIUM = 4   # requires 4/7 GABA markers (validated 0.098% FP)
DOUBLET_THRESHOLD_MERFISH = 3  # requires 3/5 GABA markers (calibrated below)


def flag_doublet_cells(adata, class_labels, subclass_to_class,
                       gaba_threshold=4, panel='xenium'):
    """Flag suspected spatial doublets using marker co-expression.

    Detects two types of doublets:
      1. Glutamatergic+GABAergic: cells classified as Glutamatergic that
         express >= gaba_threshold GABAergic markers AND more GABAergic
         markers than glutamatergic markers (diff_score > 0).
      2. Within-GABAergic: cells classified as GABAergic that co-express
         mutually exclusive subclass markers (SST+PVALB+LAMP5 for Xenium,
         PVALB+LAMP5 for MERFISH where SST is absent).

    Evidence (from 3-sample pilot + snRNAseq reference validation):
      - Type 1 cells have 1.7-1.9x UMI counts (two cells' worth)
      - Type 1 cells retain strong Glut markers (CUX2 at 69% vs 22% normal)
      - Triple co-expression (SST+PVALB+LAMP5) is <0.01% in snRNAseq
      - snRNAseq FP rate at threshold=4: 0.098% of Glut cells
      - Xenium rate: ~2.8% of Glut cells (29x amplification = true doublets)

    Parameters
    ----------
    adata : AnnData
        Expression data with raw counts in .X (NOT normalized).
    class_labels : np.ndarray of str
        Per-cell class labels ('Glutamatergic', 'GABAergic', 'Non-neuronal').
    subclass_to_class : dict
        Mapping from subclass name to class name.
    gaba_threshold : int
        Minimum number of GABAergic markers detected to flag a Glut cell.
        Default 4 (validated FP rate: 0.098% in snRNAseq at threshold=4
        with 7 Xenium markers). Use 3 for MERFISH (5 markers).
    panel : str
        Gene panel to use: 'xenium' (default) or 'merfish'.
        Selects appropriate marker gene lists and within-GABA detection
        strategy. MERFISH lacks GAD1 and SST, so uses adapted markers.

    Returns
    -------
    doublet_suspect : np.ndarray of bool
        True = suspected spatial doublet.
    doublet_type : np.ndarray of str
        '' for non-doublets, 'Glut+GABA' or 'GABA+GABA' for flagged cells.
    stats : dict
        Summary statistics for logging.
    """
    # Select panel-specific markers
    panel = panel.lower()
    if panel == 'merfish':
        gaba_markers = DOUBLET_GABA_MARKERS_MERFISH
        glut_markers = DOUBLET_GLUT_MARKERS_MERFISH
        gaba_triple = DOUBLET_GABA_TRIPLE_MERFISH
    else:
        gaba_markers = DOUBLET_GABA_MARKERS_XENIUM
        glut_markers = DOUBLET_GLUT_MARKERS_XENIUM
        gaba_triple = DOUBLET_GABA_TRIPLE_XENIUM

    gene_names = list(adata.var_names)

    X = adata.X
    if sparse.issparse(X):
        X = X.toarray()

    n_cells = X.shape[0]
    doublet_suspect = np.zeros(n_cells, dtype=bool)
    doublet_type = np.full(n_cells, '', dtype=object)

    # ── Helper: count detected markers ──
    def _marker_score(markers):
        score = np.zeros(n_cells, dtype=np.int8)
        for g in markers:
            if g in gene_names:
                score += (X[:, gene_names.index(g)] > 0).astype(np.int8)
        return score

    # ── Type 1: Glutamatergic + GABAergic doublets ──
    gaba_score = _marker_score(gaba_markers)
    glut_score = _marker_score(glut_markers)
    diff_score = gaba_score.astype(np.int16) - glut_score.astype(np.int16)

    is_glut = class_labels == 'Glutamatergic'
    glut_gaba_mask = is_glut & (gaba_score >= gaba_threshold) & (diff_score > 0)

    doublet_suspect[glut_gaba_mask] = True
    doublet_type[glut_gaba_mask] = 'Glut+GABA'

    # ── Type 2: Within-GABAergic doublets (co-expression of exclusive markers) ──
    is_gaba = class_labels == 'GABAergic'
    triple_mask = np.ones(n_cells, dtype=bool)
    for g in gaba_triple:
        if g in gene_names:
            triple_mask &= (X[:, gene_names.index(g)] > 0)
        else:
            triple_mask[:] = False
            break

    gaba_gaba_mask = is_gaba & triple_mask
    doublet_suspect[gaba_gaba_mask] = True
    doublet_type[gaba_gaba_mask] = 'GABA+GABA'

    # ── Summary stats ──
    n_glut_gaba = glut_gaba_mask.sum()
    n_gaba_gaba = gaba_gaba_mask.sum()
    n_total = doublet_suspect.sum()

    stats = {
        'n_glut_gaba': int(n_glut_gaba),
        'n_gaba_gaba': int(n_gaba_gaba),
        'n_total': int(n_total),
        'pct_of_glut': 100 * n_glut_gaba / is_glut.sum() if is_glut.sum() > 0 else 0,
        'pct_of_gaba': 100 * n_gaba_gaba / is_gaba.sum() if is_gaba.sum() > 0 else 0,
        'pct_of_all': 100 * n_total / n_cells,
    }

    print(f"  Doublet detection (panel={panel}):")
    print(f"    GABA markers: {gaba_markers}")
    print(f"    Glut markers: {glut_markers}")
    print(f"    Within-GABA:  {gaba_triple} (threshold={gaba_threshold})")
    print(f"    Glut+GABA: {n_glut_gaba:,} / {is_glut.sum():,} Glut "
          f"({stats['pct_of_glut']:.1f}%)")
    print(f"    GABA+GABA: {n_gaba_gaba:,} / {is_gaba.sum():,} GABA "
          f"({stats['pct_of_gaba']:.1f}%)")
    print(f"    Total:     {n_total:,} / {n_cells:,} ({stats['pct_of_all']:.1f}%)")

    return doublet_suspect, doublet_type, stats
