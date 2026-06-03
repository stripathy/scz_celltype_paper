"""
specificity.py — Compute cell-type expression specificity from SEA-AD snRNA-seq.

Specificity measures what fraction of a gene's total expression (across all cell types)
falls within a given cell type. A gene with specificity=0.5 in Sst_2 has half of its
total (normalized) expression in that type.

Method:
  1. Load SEA-AD h5ad (137,303 cells × 36,601 genes)
  2. Normalize to target_sum=1e4, then log1p (standard scanpy preprocessing)
  3. Compute mean expression per supertype (~137 types)
  4. Column-normalize so each cell type sums to 1000 (equalizes library sizes)
  5. Row-normalize so each gene sums to 1.0 (converts to specificity fractions)

This follows the Bryois et al. / MAGMA gene property analysis convention.
"""
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path

from ..config import (
    NORMALIZE_TARGET_SUM,
    SPECIFICITY_COLUMN_SUM,
    SUPERTYPE_COLUMN,
)
from ..utils import Timer


def load_seaad_adata(h5ad_path):
    """
    Load the SEA-AD snRNA-seq reference h5ad file.

    Parameters
    ----------
    h5ad_path : str or Path
        Path to the h5ad file.

    Returns
    -------
    anndata.AnnData
        The loaded AnnData object.
    """
    with Timer("Loading SEA-AD h5ad"):
        adata = sc.read_h5ad(str(h5ad_path))
    print(f"    Shape: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
    n_types = adata.obs[SUPERTYPE_COLUMN].nunique()
    print(f"    Supertypes: {n_types}")
    return adata


def compute_mean_expression_per_supertype(adata, target_sum=NORMALIZE_TARGET_SUM):
    """
    Normalize expression and compute mean per supertype.

    Parameters
    ----------
    adata : anndata.AnnData
        Raw SEA-AD data.
    target_sum : float
        Target sum for normalize_total (default: 1e4).

    Returns
    -------
    pd.DataFrame
        Mean expression matrix, shape (n_genes, n_supertypes).
        Gene names as index, supertype names as columns.
    """
    with Timer("Normalizing (target_sum + log1p)"):
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)

    with Timer("Computing mean expression per supertype"):
        supertypes = adata.obs[SUPERTYPE_COLUMN].unique()
        means = {}
        for st in sorted(supertypes):
            mask = adata.obs[SUPERTYPE_COLUMN] == st
            # Use .X which may be sparse
            subset = adata[mask].X
            if hasattr(subset, "toarray"):
                means[st] = np.asarray(subset.mean(axis=0)).flatten()
            else:
                means[st] = subset.mean(axis=0).flatten()

        mean_expr = pd.DataFrame(means, index=adata.var_names)

    print(f"    Result: {mean_expr.shape[0]:,} genes × {mean_expr.shape[1]} supertypes")
    return mean_expr


def compute_specificity(mean_expr, column_sum=SPECIFICITY_COLUMN_SUM):
    """
    Convert mean expression matrix to specificity scores.

    Step 1: Column-normalize so each cell type sums to `column_sum`.
            This equalizes library sizes across cell types, preventing
            types with more detected genes from dominating.
    Step 2: Row-normalize so each gene's row sums to 1.0.
            This converts values to specificity (fraction of expression).

    Parameters
    ----------
    mean_expr : pd.DataFrame
        Mean expression matrix (genes × supertypes).
    column_sum : float
        Target column sum for normalization (default: 1000).

    Returns
    -------
    pd.DataFrame
        Specificity matrix (genes × supertypes). Each row sums to ~1.0.
        Genes with zero total expression are dropped.
    """
    with Timer("Computing specificity scores"):
        # Column normalize: each cell type sums to column_sum
        col_sums = mean_expr.sum(axis=0)
        col_normalized = mean_expr * (column_sum / col_sums)

        # Drop genes with zero total expression
        row_sums = col_normalized.sum(axis=1)
        nonzero = row_sums > 0
        n_dropped = (~nonzero).sum()
        if n_dropped > 0:
            print(f"    Dropping {n_dropped} genes with zero expression")

        # Row normalize: each gene sums to 1.0
        specificity = col_normalized.loc[nonzero].div(row_sums[nonzero], axis=0)

    print(f"    Result: {specificity.shape[0]:,} genes × {specificity.shape[1]} supertypes")

    # Validation: check rows sum to ~1.0
    row_sum_check = specificity.sum(axis=1)
    assert np.allclose(row_sum_check, 1.0, atol=1e-6), (
        f"Row sums not ~1.0: range [{row_sum_check.min():.6f}, {row_sum_check.max():.6f}]"
    )
    print("    Validated: all rows sum to 1.0")

    return specificity


def compute_and_save_specificity(h5ad_path, output_path, force=False):
    """
    Full pipeline: load h5ad → normalize → mean per supertype → specificity → save.

    Skips computation if output_path already exists (unless force=True).

    Parameters
    ----------
    h5ad_path : str or Path
        Path to SEA-AD h5ad.
    output_path : str or Path
        Path to save specificity CSV.
    force : bool
        If True, recompute even if output exists.

    Returns
    -------
    pd.DataFrame
        Specificity matrix.
    """
    output_path = Path(output_path)

    if output_path.exists() and not force:
        print(f"Specificity already cached at: {output_path}")
        print("  Use --force to recompute.")
        return pd.read_csv(output_path, index_col=0)

    print("=" * 60)
    print("Computing cell-type specificity from SEA-AD reference")
    print("=" * 60)

    adata = load_seaad_adata(h5ad_path)
    mean_expr = compute_mean_expression_per_supertype(adata)
    specificity = compute_specificity(mean_expr)

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    specificity.to_csv(output_path)
    print(f"\nSaved specificity matrix to: {output_path}")
    print(f"  Size: {output_path.stat().st_size / 1e6:.1f} MB")

    return specificity


def compute_and_save_mean_expression(h5ad_path, output_path, force=False):
    """
    Compute and save mean expression per supertype (not specificity-normalized).

    Unlike the specificity matrix (row-normalized fractions), this saves the
    log-normalized mean expression values directly. Useful for analyses where
    absolute expression level matters (e.g., correlating with depth across a
    subset of cell types).

    Parameters
    ----------
    h5ad_path : str or Path
        Path to SEA-AD h5ad.
    output_path : str or Path
        Path to save mean expression CSV.
    force : bool
        If True, recompute even if output exists.

    Returns
    -------
    pd.DataFrame
        Mean expression matrix (genes × supertypes).
    """
    output_path = Path(output_path)

    if output_path.exists() and not force:
        print(f"Mean expression already cached at: {output_path}")
        print("  Use --force to recompute.")
        return pd.read_csv(output_path, index_col=0)

    print("=" * 60)
    print("Computing mean expression per supertype from SEA-AD reference")
    print("=" * 60)

    adata = load_seaad_adata(h5ad_path)
    mean_expr = compute_mean_expression_per_supertype(adata)

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    mean_expr.to_csv(output_path)
    print(f"\nSaved mean expression matrix to: {output_path}")
    print(f"  Size: {output_path.stat().st_size / 1e6:.1f} MB")

    return mean_expr
