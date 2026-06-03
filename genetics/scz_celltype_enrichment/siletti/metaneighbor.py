"""
metaneighbor.py — Cross-dataset cell-type matching via MetaNeighbor.

Implements the canonical MetaNeighbor algorithm (Crow et al. 2018) adapted for
mean expression profiles:
  1. Select highly variable genes (HVGs) by coefficient of variation
  2. Compute Spearman correlation between all cell types across datasets
  3. Convert correlations to neighbor-voting AUROCs
  4. Identify reciprocal best hits (shared cell types)
"""
import numpy as np
import pandas as pd

from ..utils import Timer


def select_hvgs(mean_expr_a, mean_expr_b, n_hvgs=5000):
    """
    Select highly variable genes shared between two datasets.

    Computes coefficient of variation (CV = std/mean) across all cell types
    in both datasets jointly. Genes must have mean expression > 0.1.

    Parameters
    ----------
    mean_expr_a : pd.DataFrame
        Mean expression (genes × types) for dataset A.
    mean_expr_b : pd.DataFrame
        Mean expression (genes × types) for dataset B.
    n_hvgs : int
        Number of HVGs to select (default: 5000).

    Returns
    -------
    list of str
        Selected HVG gene names.
    """
    shared_genes = sorted(set(mean_expr_a.index) & set(mean_expr_b.index))
    combined = pd.concat([
        mean_expr_a.loc[shared_genes],
        mean_expr_b.loc[shared_genes],
    ], axis=1)

    gene_mean = combined.mean(axis=1)
    gene_std = combined.std(axis=1)

    expressed = gene_mean > 0.1
    cv = (gene_std[expressed] / gene_mean[expressed]).replace(
        [np.inf, -np.inf], np.nan
    ).dropna()

    hvgs = cv.nlargest(min(n_hvgs, len(cv))).index.tolist()
    print(f"  HVG selection: {len(shared_genes)} shared → {expressed.sum()} expressed "
          f"→ {len(hvgs)} HVGs (top by CV)")
    return hvgs


def compute_auroc_matrix(mean_expr_a, mean_expr_b, hvgs):
    """
    Compute cross-dataset AUROC matrix using neighbor-voting on HVGs.

    For each type in dataset A, computes its Spearman correlation with all types
    in dataset B. The AUROC for pair (i, j) measures: among all dataset-B types,
    what fraction have lower correlation with A-type i than B-type j does?

    Parameters
    ----------
    mean_expr_a : pd.DataFrame
        Mean expression for dataset A (genes × types).
    mean_expr_b : pd.DataFrame
        Mean expression for dataset B (genes × types).
    hvgs : list of str
        Highly variable genes to use.

    Returns
    -------
    auroc_df : pd.DataFrame
        AUROC matrix (A types × B types).
    corr_df : pd.DataFrame
        Spearman correlation matrix (A types × B types).
    """
    a_hvg = mean_expr_a.loc[hvgs]
    b_hvg = mean_expr_b.loc[hvgs]
    n_a = a_hvg.shape[1]
    n_b = b_hvg.shape[1]

    with Timer(f"Computing Spearman correlation ({n_a} × {n_b})"):
        a_ranks = a_hvg.rank(axis=0)
        b_ranks = b_hvg.rank(axis=0)
        corr_vals = np.corrcoef(a_ranks.values.T, b_ranks.values.T)[:n_a, n_a:]

    corr_df = pd.DataFrame(corr_vals, index=a_hvg.columns, columns=b_hvg.columns)
    print(f"    Correlation range: [{corr_vals.min():.3f}, {corr_vals.max():.3f}]")

    with Timer(f"Computing AUROC via neighbor voting ({n_a} × {n_b})"):
        auroc_vals = np.zeros((n_a, n_b))
        for i in range(n_a):
            corrs = corr_vals[i]
            for j in range(n_b):
                target = corrs[j]
                others = np.concatenate([corrs[:j], corrs[j + 1:]])
                auroc_vals[i, j] = (
                    np.sum(target > others) + 0.5 * np.sum(target == others)
                ) / len(others)

            if (i + 1) % 20 == 0:
                print(f"      {i + 1}/{n_a}...")

    auroc_df = pd.DataFrame(auroc_vals, index=a_hvg.columns, columns=b_hvg.columns)
    print(f"    AUROC range: [{auroc_vals.min():.4f}, {auroc_vals.max():.4f}]")

    return auroc_df, corr_df


def find_reciprocal_best_hits(auroc_matrix):
    """
    Find reciprocal best hits between two datasets.

    A reciprocal best hit (RBH) occurs when:
      - Dataset-A type i's best match in B is type j
      - Dataset-B type j's best match in A is type i

    Parameters
    ----------
    auroc_matrix : pd.DataFrame
        AUROC matrix (A types × B types).

    Returns
    -------
    rbh : pd.DataFrame
        Reciprocal best hits with columns: type_a, type_b,
        auroc_a_to_b, auroc_b_to_a, mean_auroc.
    all_matches : pd.DataFrame
        All A types with their top 3 B matches and reciprocal status.
    """
    a_best = auroc_matrix.idxmax(axis=1)
    a_best_auroc = auroc_matrix.max(axis=1)
    b_best = auroc_matrix.idxmax(axis=0)
    b_best_auroc = auroc_matrix.max(axis=0)

    rbh_rows = []
    all_rows = []

    for a_type in auroc_matrix.index:
        best_b = a_best[a_type]
        is_recip = b_best[best_b] == a_type
        top3 = auroc_matrix.loc[a_type].nlargest(3)

        all_rows.append({
            "seaad_type": a_type,
            "best_siletti": best_b,
            "best_auroc": a_best_auroc[a_type],
            "second_siletti": top3.index[1],
            "second_auroc": top3.iloc[1],
            "third_siletti": top3.index[2],
            "third_auroc": top3.iloc[2],
            "is_reciprocal": is_recip,
        })

        if is_recip:
            rbh_rows.append({
                "seaad_type": a_type,
                "siletti_cluster": best_b,
                "auroc_seaad_to_siletti": a_best_auroc[a_type],
                "auroc_siletti_to_seaad": b_best_auroc[best_b],
                "mean_auroc": (a_best_auroc[a_type] + b_best_auroc[best_b]) / 2,
            })

    rbh = pd.DataFrame(rbh_rows).sort_values("mean_auroc", ascending=False).reset_index(drop=True)
    all_matches = pd.DataFrame(all_rows)

    n_rbh = len(rbh)
    n_high = (rbh["mean_auroc"] > 0.9).sum() if n_rbh > 0 else 0
    print(f"\n  Reciprocal best hits: {n_rbh}/{len(auroc_matrix.index)}")
    print(f"  High-confidence (AUROC > 0.9): {n_high}")

    return rbh, all_matches


def build_rbh_combined_taxonomy(all_matches, seaad_cols, siletti_cols):
    """
    Determine which Siletti clusters are novel (not claimed by any SEA-AD type).

    Logic:
      - Every SEA-AD type is kept
      - Siletti clusters that are the best match for ANY SEA-AD type are excluded
        (already represented by that SEA-AD type)
      - Remaining Siletti clusters are included as novel types

    Parameters
    ----------
    all_matches : pd.DataFrame
        All SEA-AD types with best Siletti matches (from find_reciprocal_best_hits).
    seaad_cols : list of str
        SEA-AD type names.
    siletti_cols : list of str
        All Siletti cluster names (without "Siletti_" prefix).

    Returns
    -------
    novel_siletti : list of str
        Siletti clusters not claimed by any SEA-AD type.
    claimed_siletti : set of str
        Siletti clusters claimed as best match by ≥1 SEA-AD type.
    """
    claimed = set(all_matches["best_siletti"])
    novel = sorted(set(siletti_cols) - claimed)

    print(f"  SEA-AD types: {len(seaad_cols)}")
    print(f"  Siletti claimed by SEA-AD: {len(claimed)}")
    print(f"  Novel Siletti clusters: {len(novel)}")
    print(f"  Combined taxonomy: {len(seaad_cols)} + {len(novel)} = {len(seaad_cols) + len(novel)}")

    return novel, claimed
