"""
Correlation-based label transfer from SEA-AD snRNAseq reference.

Reimplements the scrattch.mapping approach used by SEA-AD for MERFISH:
  1. Build reference centroids (mean log-normalized expression per cell type)
  2. Correlate each query cell against all centroids
  3. Assign best-matching type at each taxonomy level

This replaces the kNN-in-PCA approach with a more robust centroid-correlation
method that is specifically designed for gene-limited spatial platforms.
"""

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import sparse
from scipy.stats import pearsonr
import time


# SEA-AD taxonomy levels, ordered from coarse to fine
TAXONOMY_LEVELS = ["class_label", "subclass_label", "supertype_label"]


def build_hierarchy_table(ref):
    """
    Build a deterministic cluster → (subclass, class) lookup from reference.

    Each supertype_label maps to exactly one subclass_label and class_label.
    This enables hierarchical label transfer: assign at the finest level
    (cluster), then derive coarser labels from the hierarchy.

    Parameters
    ----------
    ref : anndata.AnnData
        Reference with 'supertype_label', 'subclass_label', 'class_label' in obs.

    Returns
    -------
    dict
        {supertype_label: {'subclass_label': str, 'class_label': str}}
    """
    hier = ref.obs[['supertype_label', 'subclass_label', 'class_label']].drop_duplicates()

    # Verify clean hierarchy (each cluster maps to exactly one subclass/class)
    cluster_counts = hier['supertype_label'].value_counts()
    ambiguous = cluster_counts[cluster_counts > 1]
    if len(ambiguous) > 0:
        print(f"  WARNING: {len(ambiguous)} clusters map to multiple parents!")
        for c in ambiguous.index:
            rows = hier[hier['supertype_label'] == c]
            print(f"    {c}: {rows[['subclass_label', 'class_label']].values.tolist()}")

    table = {}
    for _, row in hier.iterrows():
        table[row['supertype_label']] = {
            'subclass_label': row['subclass_label'],
            'class_label': row['class_label'],
        }

    print(f"  Hierarchy table: {len(table)} clusters → "
          f"{len(set(v['subclass_label'] for v in table.values()))} subclasses → "
          f"{len(set(v['class_label'] for v in table.values()))} classes")

    return table


def load_reference(reference_path, shared_genes=None):
    """
    Load the SEA-AD reference dataset.

    Parameters
    ----------
    reference_path : str
        Path to the reference .h5ad file.
    shared_genes : list of str, optional
        If provided, subset to these genes.

    Returns
    -------
    anndata.AnnData
    """
    print(f"Loading reference from {reference_path}...")
    ref = ad.read_h5ad(reference_path)
    print(f"  Reference: {ref.shape[0]:,} cells x {ref.shape[1]:,} genes")

    if shared_genes is not None:
        shared = sorted(set(shared_genes) & set(ref.var_names))
        ref = ref[:, shared].copy()
        print(f"  After subsetting to shared genes: {ref.shape[1]} genes")

    return ref


def normalize_log(adata, target_sum=1e4):
    """
    Normalize counts and log-transform in-place.
    """
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)


def build_centroids(ref, level, method="mean"):
    """
    Compute centroid expression profiles for each cell type at a given level.

    Parameters
    ----------
    ref : anndata.AnnData
        Reference AnnData, already log-normalized, subsetted to shared genes.
    level : str
        Taxonomy column name (e.g., 'subclass_label').
    method : str
        'mean' or 'median' for centroid computation.

    Returns
    -------
    centroids : pd.DataFrame
        Shape (n_types, n_genes). Index = cell type names, columns = gene names.
    counts : pd.Series
        Number of reference cells per type.
    """
    labels = ref.obs[level].astype(str)

    # Get dense matrix
    X = ref.X
    if sparse.issparse(X):
        X = X.toarray()

    df = pd.DataFrame(X, index=ref.obs_names, columns=ref.var_names)
    df["_label"] = labels.values

    if method == "mean":
        centroids = df.groupby("_label").mean()
    else:
        centroids = df.groupby("_label").median()

    # Clean NaN values in centroids
    centroids = centroids.fillna(0)

    counts = labels.value_counts()

    return centroids, counts


def correlate_to_centroids(query_expr, centroids):
    """
    Correlate each query cell against all reference centroids.

    Uses Pearson correlation, matching the scrattch.mapping approach.

    Parameters
    ----------
    query_expr : np.ndarray
        Dense matrix (n_cells, n_genes), log-normalized.
    centroids : pd.DataFrame
        (n_types, n_genes) centroid expression profiles.

    Returns
    -------
    corr_matrix : np.ndarray
        (n_cells, n_types) Pearson correlation values.
    type_names : list
        Cell type names corresponding to columns.
    """
    n_cells = query_expr.shape[0]
    type_names = list(centroids.index)
    n_types = len(type_names)
    centroid_arr = centroids.values.astype(np.float64)  # (n_types, n_genes)

    # Vectorized Pearson correlation:
    # Standardize both query and centroids, then dot product
    # Pearson r = mean( (x - mean_x)/std_x * (y - mean_y)/std_y )

    # Clean any NaN/inf in centroids
    centroid_arr = np.nan_to_num(centroid_arr, nan=0.0, posinf=0.0, neginf=0.0)

    # Standardize centroids
    c_mean = centroid_arr.mean(axis=1, keepdims=True)
    c_std = centroid_arr.std(axis=1, keepdims=True, ddof=0)
    c_std[c_std == 0] = 1.0
    c_norm = (centroid_arr - c_mean) / c_std  # (n_types, n_genes)

    # Standardize query in chunks to manage memory
    chunk_size = 5000
    corr_matrix = np.zeros((n_cells, n_types), dtype=np.float32)

    for start in range(0, n_cells, chunk_size):
        end = min(start + chunk_size, n_cells)
        chunk = query_expr[start:end].astype(np.float64)  # (chunk, n_genes)

        q_mean = chunk.mean(axis=1, keepdims=True)
        q_std = chunk.std(axis=1, keepdims=True, ddof=0)
        q_std[q_std == 0] = 1.0
        q_norm = (chunk - q_mean) / q_std  # (chunk, n_genes)

        # Clean NaN/inf
        q_norm = np.nan_to_num(q_norm, nan=0.0, posinf=0.0, neginf=0.0)

        # Dot product gives n_genes * pearson_r
        corr_matrix[start:end] = ((q_norm @ c_norm.T) / centroids.shape[1]).astype(np.float32)

    # Clean result
    corr_matrix = np.nan_to_num(corr_matrix, nan=0.0, posinf=0.0, neginf=0.0)

    return corr_matrix, type_names


def assign_labels_from_correlation(corr_matrix, type_names):
    """
    Assign best-matching label and confidence from correlation matrix.

    Parameters
    ----------
    corr_matrix : np.ndarray
        (n_cells, n_types)
    type_names : list
        Cell type names.

    Returns
    -------
    labels : np.ndarray
        Best-matching type for each cell.
    confidences : np.ndarray
        Max correlation value for each cell.
    """
    best_idx = np.argmax(corr_matrix, axis=1)
    labels = np.array([type_names[i] for i in best_idx])
    confidences = np.max(corr_matrix, axis=1)

    return labels, confidences


def correlation_label_transfer(query_adata, ref, levels=None,
                                centroid_method="mean", min_corr=0.0,
                                hierarchical=True):
    """
    Full correlation-based label transfer pipeline.

    Reimplements the scrattch.mapping approach used by SEA-AD:
    1. Subset both datasets to shared genes
    2. Normalize + log-transform both
    3. Build reference centroids per cell type
    4. Correlate each query cell to centroids
    5. Assign best match

    When hierarchical=True (default), assigns only at the finest level
    (supertype_label) and derives subclass_label and class_label from
    the reference taxonomy hierarchy. This guarantees internal consistency.

    When hierarchical=False (legacy), assigns independently at each level.

    Parameters
    ----------
    query_adata : anndata.AnnData
        Spatial data with raw counts in .X.
    ref : anndata.AnnData
        SEA-AD reference with raw counts and taxonomy labels.
    levels : list of str, optional
        Taxonomy levels to transfer. Default: TAXONOMY_LEVELS.
    centroid_method : str
        'mean' or 'median' for centroid computation.
    min_corr : float
        Minimum correlation to accept an assignment.
        Cells below this get 'Unassigned'.
    hierarchical : bool
        If True, assign only at supertype_label level and propagate
        subclass_label and class_label from hierarchy. Default: True.

    Returns
    -------
    anndata.AnnData
        Query AnnData with transferred labels, correlations, and
        PCA coordinates in .obsm['X_pca'].
    """
    if levels is None:
        levels = TAXONOMY_LEVELS

    t0 = time.time()

    # 1. Subset to shared genes
    shared_genes = sorted(set(query_adata.var_names) & set(ref.var_names))
    n_shared = len(shared_genes)
    print(f"  Shared genes: {n_shared}")

    query = query_adata[:, shared_genes].copy()
    ref_sub = ref[:, shared_genes].copy()

    # 2. Normalize
    already_normalized = query_adata.uns.get("already_normalized", False)
    if not already_normalized:
        query.layers["counts"] = query.X.copy()
        normalize_log(query)
    query.layers["lognorm"] = query.X.copy()
    normalize_log(ref_sub)

    # Get dense query expression
    X_query = query.X
    if sparse.issparse(X_query):
        X_query = X_query.toarray()
    X_query = np.asarray(X_query, dtype=np.float32)
    # Clean any NaN/inf values
    X_query = np.nan_to_num(X_query, nan=0.0, posinf=0.0, neginf=0.0)

    if hierarchical:
        # ── Hierarchical mode: assign at cluster level, propagate up ──
        print("  Mode: HIERARCHICAL (cluster → subclass → class)")
        hier_table = build_hierarchy_table(ref)

        finest_level = "supertype_label"
        if finest_level not in ref_sub.obs.columns:
            raise ValueError(f"'{finest_level}' not in reference for hierarchical mode")

        print(f"  Transferring {finest_level}...")
        centroids, counts = build_centroids(ref_sub, finest_level,
                                             method=centroid_method)
        print(f"    {len(centroids)} types, {n_shared} genes")

        corr_matrix, type_names = correlate_to_centroids(X_query, centroids)
        labels, confidences = assign_labels_from_correlation(corr_matrix, type_names)

        # Apply minimum correlation filter
        if min_corr > 0:
            labels[confidences < min_corr] = "Unassigned"

        # Assign cluster labels
        query.obs["supertype_label"] = labels
        query.obs["supertype_label_confidence"] = confidences

        # Store correlation matrix
        query.obsm["corr_matrix"] = corr_matrix

        # Derive subclass and class from hierarchy
        subclass_labels = np.array([
            hier_table.get(c, {}).get('subclass_label', 'Unknown')
            for c in labels
        ])
        class_labels = np.array([
            hier_table.get(c, {}).get('class_label', 'Unknown')
            for c in labels
        ])

        # Handle Unassigned cells
        subclass_labels[labels == "Unassigned"] = "Unassigned"
        class_labels[labels == "Unassigned"] = "Unassigned"

        query.obs["subclass_label"] = subclass_labels
        query.obs["subclass_label_confidence"] = confidences  # same as cluster
        query.obs["class_label"] = class_labels
        query.obs["class_label_confidence"] = confidences  # same as cluster

        # Print distributions
        for level_name in ["supertype_label", "subclass_label", "class_label"]:
            vc = query.obs[level_name].value_counts()
            print(f"    {level_name}: {len(vc)} types, "
                  f"top 5: {dict(vc.head(5))}")

        # Verify consistency
        n_unknown_sub = (subclass_labels == 'Unknown').sum()
        n_unknown_cls = (class_labels == 'Unknown').sum()
        if n_unknown_sub > 0 or n_unknown_cls > 0:
            print(f"    WARNING: {n_unknown_sub} cells got Unknown subclass, "
                  f"{n_unknown_cls} got Unknown class (cluster not in hierarchy)")

        print(f"    Median correlation: {np.median(confidences):.3f}, "
              f"range [{confidences.min():.3f}, {confidences.max():.3f}]")

    else:
        # ── Legacy mode: independent assignment at each level ──
        print("  Mode: INDEPENDENT (legacy, no hierarchy enforcement)")
        for level in levels:
            if level not in ref_sub.obs.columns:
                print(f"  Warning: '{level}' not in reference, skipping")
                continue

            print(f"  Transferring {level}...")
            centroids, counts = build_centroids(ref_sub, level,
                                                 method=centroid_method)
            print(f"    {len(centroids)} types, {n_shared} genes")

            corr_matrix, type_names = correlate_to_centroids(X_query, centroids)
            labels, confidences = assign_labels_from_correlation(
                corr_matrix, type_names
            )

            # Apply minimum correlation filter
            if min_corr > 0:
                labels[confidences < min_corr] = "Unassigned"

            query.obs[level] = labels
            query.obs[f"{level}_confidence"] = confidences

            # Store full correlation matrix for the finest level
            if level == levels[-1]:
                query.obsm["corr_matrix"] = corr_matrix

            # Print distribution
            vc = pd.Series(labels).value_counts()
            print(f"    Top 5: {dict(vc.head(5))}")
            print(f"    Median correlation: {np.median(confidences):.3f}, "
                  f"range [{confidences.min():.3f}, {confidences.max():.3f}]")

    # Also compute PCA for downstream use (e.g., UMAP)
    print("  Computing PCA...")
    try:
        query_pca = query.copy()
        sc.pp.scale(query_pca, max_value=10)
        # Replace NaN from scaling (constant genes)
        if sparse.issparse(query_pca.X):
            query_pca.X = query_pca.X.toarray()
        query_pca.X = np.nan_to_num(query_pca.X, nan=0.0)
        sc.tl.pca(query_pca, n_comps=min(50, n_shared - 1))
        query.obsm["X_pca"] = query_pca.obsm["X_pca"]
    except Exception as e:
        print(f"  Warning: PCA failed ({e}), skipping.")

    # Copy spatial coords and sample_id
    if "spatial" in query_adata.obsm:
        query.obsm["spatial"] = query_adata[:, shared_genes].obsm["spatial"]
    if "sample_id" in query_adata.obs.columns:
        query.obs["sample_id"] = query_adata.obs["sample_id"].values

    elapsed = time.time() - t0
    print(f"  Correlation label transfer complete in {elapsed:.1f}s")

    return query


def annotate_sample(query_adata, ref, method="correlation", **kwargs):
    """
    Unified interface: annotate a spatial sample using chosen method.

    Parameters
    ----------
    query_adata : anndata.AnnData
        Spatial data with raw counts.
    ref : anndata.AnnData
        SEA-AD reference.
    method : str
        'correlation' (scrattch-style) or 'knn' (PCA-based).
    **kwargs
        Passed to the chosen method.

    Returns
    -------
    anndata.AnnData
        Annotated query.
    """
    if method == "correlation":
        return correlation_label_transfer(query_adata, ref, **kwargs)
    elif method == "knn":
        return _knn_label_transfer(query_adata, ref, **kwargs)
    else:
        raise ValueError(f"Unknown method: {method}")


def _knn_label_transfer(query_adata, ref, n_comps=50, n_neighbors=15,
                        levels=None):
    """Legacy kNN-in-PCA label transfer (kept for comparison)."""
    from sklearn.neighbors import KNeighborsClassifier

    if levels is None:
        levels = TAXONOMY_LEVELS

    shared_genes = sorted(set(query_adata.var_names) & set(ref.var_names))
    query = query_adata[:, shared_genes].copy()
    ref_copy = ref[:, shared_genes].copy()

    query.layers["counts"] = query.X.copy()
    normalize_log(ref_copy)
    normalize_log(query)
    query.layers["lognorm"] = query.X.copy()

    # Build shared PCA
    ref_copy.obs["_source"] = "reference"
    query.obs["_source"] = "query"
    combined = ad.concat([ref_copy, query], join="inner")
    sc.pp.scale(combined, max_value=10)
    sc.tl.pca(combined, n_comps=n_comps)

    ref_mask = combined.obs["_source"] == "reference"
    query_mask = combined.obs["_source"] == "query"
    ref_pca = combined[ref_mask].obsm["X_pca"].copy()
    query_pca = combined[query_mask].obsm["X_pca"].copy()
    query.obsm["X_pca"] = query_pca

    # kNN transfer
    for level in levels:
        if level not in ref_copy.obs.columns:
            continue
        print(f"  Transferring {level} (kNN k={n_neighbors})...")
        knn = KNeighborsClassifier(
            n_neighbors=n_neighbors, weights="distance", n_jobs=-1
        )
        knn.fit(ref_pca, ref_copy.obs[level].values.astype(str))
        query.obs[level] = knn.predict(query_pca)

    if "spatial" in query_adata.obsm:
        query.obsm["spatial"] = query_adata[:, shared_genes].obsm["spatial"]
    if "sample_id" in query_adata.obs.columns:
        query.obs["sample_id"] = query_adata.obs["sample_id"].values

    return query


def get_seaad_colors(ref):
    """
    Extract SEA-AD color mappings from reference .obs.

    Returns
    -------
    dict
        Nested dict: {level: {label: hex_color}}.
    """
    colors = {}
    for level in TAXONOMY_LEVELS:
        color_col = level.replace("_label", "_color")
        if color_col in ref.obs.columns:
            pairs = ref.obs[[level, color_col]].drop_duplicates()
            colors[level] = dict(zip(
                pairs[level].astype(str),
                pairs[color_col].astype(str),
            ))
    return colors
