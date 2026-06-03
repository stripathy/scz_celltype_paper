"""
Pseudobulk computation utilities.

Used by: compare_snrnaseq_merscope_expression.py, cross_platform_gene_corr.py
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp


def compute_pseudobulk_mean(adata, groupby_col, gene_subset=None,
                            min_cells=10, verbose=True):
    """Compute pseudobulk mean expression per group.

    Expects adata.X to be log1p-normalized already.

    Parameters
    ----------
    adata : AnnData
        Annotated data with expression matrix in .X.
    groupby_col : str
        Column in adata.obs to group by (e.g., 'Subclass', 'corr_subclass').
    gene_subset : list of str, optional
        Subset of genes to compute pseudobulk for. If None, uses all genes.
    min_cells : int
        Minimum cells per group to include (default 10).
    verbose : bool
        Print group sizes.

    Returns
    -------
    DataFrame : rows = groups, columns = genes, values = mean expression
    """
    groups = adata.obs[groupby_col].unique()
    groups = [g for g in groups
              if g not in ("Unassigned", None, np.nan, "nan", "")]
    groups = sorted(groups)

    if gene_subset is not None:
        var_names = list(adata.var_names)
        gene_idx = [var_names.index(g) for g in gene_subset if g in var_names]
        gene_names = [g for g in gene_subset if g in var_names]
    else:
        gene_idx = list(range(adata.shape[1]))
        gene_names = list(adata.var_names)

    results = {}
    for g in groups:
        mask = adata.obs[groupby_col] == g
        n_cells = mask.sum()
        if n_cells < min_cells:
            if verbose:
                print(f"  Skipping {g}: only {n_cells} cells")
            continue
        X_sub = adata.X[mask.values][:, gene_idx]
        if sp.issparse(X_sub):
            X_sub = X_sub.toarray()
        results[g] = np.mean(X_sub, axis=0)
        if verbose:
            print(f"  {g}: {n_cells} cells")

    df = pd.DataFrame(results, index=gene_names).T
    return df
