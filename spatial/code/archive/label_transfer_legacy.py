"""
Label transfer from SEA-AD snRNAseq reference to Xenium spatial data.

Uses kNN classification in a shared PCA space to transfer cell type
labels at multiple taxonomy levels (class, subclass, cluster).
"""

import numpy as np
import anndata as ad
import scanpy as sc
from sklearn.neighbors import KNeighborsClassifier


# SEA-AD taxonomy levels, ordered from coarse to fine
TAXONOMY_LEVELS = ["class_label", "subclass_label", "cluster_label"]


def load_reference(reference_path, shared_genes=None):
    """
    Load the SEA-AD reference dataset.

    Parameters
    ----------
    reference_path : str
        Path to the reference .h5ad file. Can be the full reference
        or the 10% subsample.
    shared_genes : list of str, optional
        If provided, subset to these genes. If None, uses all genes.

    Returns
    -------
    anndata.AnnData
        Reference AnnData, subsetted to shared genes if specified.
    """
    print(f"Loading reference from {reference_path}...")
    ref = ad.read_h5ad(reference_path)
    print(f"  Reference: {ref.shape[0]:,} cells x {ref.shape[1]:,} genes")

    if shared_genes is not None:
        shared = sorted(set(shared_genes) & set(ref.var_names))
        ref = ref[:, shared].copy()
        print(f"  After subsetting to shared genes: {ref.shape[1]} genes")

    return ref


def normalize_for_transfer(adata, target_sum=1e4):
    """
    Normalize an AnnData object for label transfer.

    Applies library-size normalization and log1p transformation.
    Operates in-place on .X.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData with raw counts in .X.
    target_sum : float
        Target sum for normalization.
    """
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)


def build_shared_pca(ref, query, n_comps=50, max_value=10):
    """
    Build a shared PCA space from reference and query data.

    Concatenates reference and query, scales, computes PCA, then
    splits back. Both inputs should already be log-normalized.

    Parameters
    ----------
    ref : anndata.AnnData
        Reference AnnData (log-normalized).
    query : anndata.AnnData
        Query AnnData (log-normalized).
    n_comps : int
        Number of PCA components.
    max_value : float
        Max value for scaling.

    Returns
    -------
    ref_pca : np.ndarray
        PCA coordinates for reference cells (n_ref x n_comps).
    query_pca : np.ndarray
        PCA coordinates for query cells (n_query x n_comps).
    """
    ref.obs["_source"] = "reference"
    query.obs["_source"] = "query"

    combined = ad.concat([ref, query], join="inner")
    sc.pp.scale(combined, max_value=max_value)
    sc.tl.pca(combined, n_comps=n_comps)

    ref_mask = combined.obs["_source"] == "reference"
    query_mask = combined.obs["_source"] == "query"

    ref_pca = combined[ref_mask].obsm["X_pca"].copy()
    query_pca = combined[query_mask].obsm["X_pca"].copy()

    # Store PCA in query for downstream use
    query.obsm["X_pca"] = query_pca

    # Clean up
    ref.obs.drop("_source", axis=1, inplace=True)
    query.obs.drop("_source", axis=1, inplace=True)

    return ref_pca, query_pca


def transfer_labels(ref_pca, query_pca, ref_obs, query_adata,
                    levels=None, n_neighbors=15):
    """
    Transfer cell type labels via kNN classification in PCA space.

    Parameters
    ----------
    ref_pca : np.ndarray
        Reference PCA coordinates.
    query_pca : np.ndarray
        Query PCA coordinates.
    ref_obs : pd.DataFrame
        Reference .obs with taxonomy columns.
    query_adata : anndata.AnnData
        Query AnnData to receive labels (modified in-place).
    levels : list of str, optional
        Taxonomy levels to transfer. Defaults to TAXONOMY_LEVELS.
    n_neighbors : int
        Number of neighbors for kNN.
    """
    if levels is None:
        levels = TAXONOMY_LEVELS

    for level in levels:
        if level not in ref_obs.columns:
            print(f"  Warning: '{level}' not in reference, skipping")
            continue
        print(f"  Transferring {level}...")
        knn = KNeighborsClassifier(
            n_neighbors=n_neighbors, weights="distance", n_jobs=-1
        )
        knn.fit(ref_pca, ref_obs[level].values.astype(str))
        query_adata.obs[level] = knn.predict(query_pca)

    print("  Label transfer complete.")


def annotate_sample(query_adata, ref, n_comps=50, n_neighbors=15):
    """
    Full label transfer pipeline for a single Xenium sample.

    Takes a raw-count query AnnData and a prepared reference,
    normalizes, builds shared PCA, and transfers labels.

    Parameters
    ----------
    query_adata : anndata.AnnData
        Xenium sample with raw counts in .X. Modified in-place.
    ref : anndata.AnnData
        SEA-AD reference, already subsetted to shared genes.
        Will be copied and normalized internally.
    n_comps : int
        Number of PCA components.
    n_neighbors : int
        Number of kNN neighbors.

    Returns
    -------
    anndata.AnnData
        The query AnnData with added taxonomy labels in .obs and
        PCA coordinates in .obsm['X_pca'].
    """
    # Subset query to shared genes
    shared_genes = sorted(set(query_adata.var_names) & set(ref.var_names))
    query = query_adata[:, shared_genes].copy()

    # Save raw counts
    query.layers["counts"] = query.X.copy()

    # Normalize both
    ref_copy = ref[:, shared_genes].copy()
    normalize_for_transfer(ref_copy)
    normalize_for_transfer(query)
    query.layers["lognorm"] = query.X.copy()

    # Build shared PCA and transfer
    print(f"  Building shared PCA ({n_comps} components)...")
    ref_pca, query_pca = build_shared_pca(ref_copy, query, n_comps=n_comps)

    print(f"  Running kNN (k={n_neighbors})...")
    transfer_labels(ref_pca, query_pca, ref_copy.obs, query)

    # Copy spatial coords and sample_id if present
    if "spatial" in query_adata.obsm:
        query.obsm["spatial"] = query_adata[:, shared_genes].obsm["spatial"]
    if "sample_id" in query_adata.obs.columns:
        query.obs["sample_id"] = query_adata.obs["sample_id"].values

    return query


def get_seaad_colors(ref):
    """
    Extract SEA-AD color mappings from reference .obs.

    Parameters
    ----------
    ref : anndata.AnnData
        SEA-AD reference with color columns in .obs.

    Returns
    -------
    dict
        Nested dict: {level: {label: hex_color}}.
        E.g., colors['subclass_label']['L2/3 IT'] = '#B1EC30'
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
