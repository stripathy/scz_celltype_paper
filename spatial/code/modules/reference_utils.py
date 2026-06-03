"""
snRNAseq reference loading, normalization, and subsampling utilities.

Used by: hierarchical_probe_selection.py, supertype_markers_panel_overlap.py,
         merscope_panel_assessment.py, nsforest_supertype_markers.py,
         compare_snrnaseq_merscope_expression.py, cross_platform_gene_corr.py,
         characterize_poor_genes.py
"""

import time
import numpy as np
import anndata as ad
import scanpy as sc
import scipy.sparse as sp


def load_and_normalize_reference(path, normalize=True, min_cells=10,
                                 target_sum=1e4):
    """Load SEA-AD snRNAseq reference with optional normalization.

    Automatically detects whether data is raw counts (integers) and
    normalizes accordingly.

    Parameters
    ----------
    path : str
        Path to the reference h5ad file.
    normalize : bool
        If True, apply normalize_total + log1p if data appears to be raw counts.
    min_cells : int
        Minimum cells for gene filtering. Set to 0 to skip.
    target_sum : float
        Target sum for normalize_total (default 1e4).

    Returns
    -------
    AnnData with (optionally) normalized expression in .X
    """
    print(f"Loading snRNAseq reference from {path}...")
    t0 = time.time()
    adata = ad.read_h5ad(path)
    print(f"  Shape: {adata.shape[0]} cells x {adata.shape[1]} genes "
          f"(loaded in {time.time()-t0:.0f}s)")

    if normalize:
        # Check if data looks like raw counts
        X_sample = adata.X[:100, :100]
        if sp.issparse(X_sample):
            X_sample = X_sample.toarray()
        if np.all(X_sample == X_sample.astype(int)):
            print("  Data appears to be raw counts. Normalizing...")
            sc.pp.normalize_total(adata, target_sum=target_sum)
            sc.pp.log1p(adata)
        else:
            print("  Data appears already normalized. Skipping normalization.")

    if min_cells > 0:
        n_before = adata.shape[1]
        sc.pp.filter_genes(adata, min_cells=min_cells)
        print(f"  Filtered genes: {n_before} -> {adata.shape[1]} "
              f"(min_cells={min_cells})")

    return adata


def subsample_by_group(adata, groupby, max_cells=500, min_cells=20, seed=42):
    """Subsample to max_cells per group, dropping groups below min_cells.

    Parameters
    ----------
    adata : AnnData
        Input data.
    groupby : str
        obs column to group by (e.g., 'Supertype', 'Subclass').
    max_cells : int
        Maximum cells per group.
    min_cells : int
        Groups with fewer cells are dropped.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    AnnData (subsampled copy)
    """
    type_counts = adata.obs[groupby].value_counts()
    valid_types = type_counts[type_counts >= min_cells].index.tolist()
    skipped = type_counts[type_counts < min_cells]
    if len(skipped) > 0:
        print(f"  Dropping {len(skipped)} groups with < {min_cells} cells")

    adata = adata[adata.obs[groupby].isin(valid_types)].copy()

    np.random.seed(seed)
    keep_idx = []
    for group in valid_types:
        idx = np.where(adata.obs[groupby] == group)[0]
        if len(idx) > max_cells:
            idx = np.random.choice(idx, max_cells, replace=False)
        keep_idx.extend(idx)

    adata = adata[sorted(keep_idx)].copy()

    # Clean up unused categories if categorical
    if hasattr(adata.obs[groupby], 'cat'):
        adata.obs[groupby] = adata.obs[groupby].cat.remove_unused_categories()

    print(f"  Subsampled to {adata.shape[0]} cells, "
          f"{adata.obs[groupby].nunique()} groups "
          f"(max {max_cells}/group)")

    return adata
