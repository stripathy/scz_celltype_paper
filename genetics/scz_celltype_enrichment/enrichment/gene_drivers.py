"""
gene_drivers.py — Identify genes driving SCZ enrichment per cell type.

For each enriched cell type, the "contribution" of a gene is:
    contribution = specificity × GWAS_Z

This captures both how specific a gene is to the cell type AND how strongly
it is associated with schizophrenia. Genes with high contribution drive
the enrichment signal.

This module also decomposes the signal into:
  - Shared (pan-interneuron) genes: high contribution across many cell types
  - Type-specific genes: high contribution in one type but not others
  - Uniqueness scores: weight contribution by type-exclusivity
"""
import numpy as np
import pandas as pd
from itertools import combinations

from ..config import TOP_GENES_PER_TYPE
from ..utils import Timer


def compute_gene_contributions(specificity, gwas_df, supertype, scz_gene_set=None):
    """
    Compute per-gene contribution scores for a single cell type.

    contribution = specificity × GWAS_Z

    Parameters
    ----------
    specificity : pd.DataFrame
        Full specificity matrix (genes × supertypes).
    gwas_df : pd.DataFrame
        GWAS results with SYMBOL, ZSTAT, P columns.
    supertype : str
        Cell type name.
    scz_gene_set : set, optional
        If provided, add is_scz_gwas boolean column.

    Returns
    -------
    pd.DataFrame
        Sorted by contribution (descending), with columns:
        gene, specificity, gwas_zstat, gwas_p, contribution, is_scz_gwas
    """
    gwas_indexed = gwas_df.set_index("SYMBOL")
    genes = specificity.index

    spec_vals = specificity[supertype].values
    zstat_vals = gwas_indexed.loc[genes, "ZSTAT"].values
    pvals = gwas_indexed.loc[genes, "P"].values

    contributions = spec_vals * zstat_vals

    result = pd.DataFrame({
        "gene": genes,
        "specificity": spec_vals,
        "gwas_zstat": zstat_vals,
        "gwas_p": pvals,
        "contribution": contributions,
    })

    if scz_gene_set is not None:
        result["is_scz_gwas"] = result["gene"].isin(scz_gene_set)

    result = result.sort_values("contribution", ascending=False).reset_index(drop=True)
    return result


def compute_uniqueness_scores(specificity, gwas_df, independent_types, top_n=TOP_GENES_PER_TYPE):
    """
    For each independent type, compute uniqueness scores.

    uniqueness = contribution × type_fraction
    type_fraction = this_type_contribution / sum(all_independent_types_contribution)

    A gene with high uniqueness is both a strong contributor AND concentrated
    in this particular cell type rather than being shared across types.

    Parameters
    ----------
    specificity : pd.DataFrame
        Full specificity matrix.
    gwas_df : pd.DataFrame
        GWAS results with SYMBOL, ZSTAT columns.
    independent_types : list of str
        Independent cell type names.
    top_n : int
        Number of top genes to consider per type.

    Returns
    -------
    dict
        Mapping supertype → pd.DataFrame with uniqueness scores.
    """
    gwas_indexed = gwas_df.set_index("SYMBOL")
    genes = specificity.index
    zstat = gwas_indexed.loc[genes, "ZSTAT"].values

    # Compute contribution for all independent types
    all_contributions = {}
    for st in independent_types:
        all_contributions[st] = specificity[st].values * zstat

    # For each type, compute type_fraction
    results = {}
    for st in independent_types:
        this_contrib = all_contributions[st]

        # Sum of positive contributions across all independent types (per gene)
        total_contrib = np.zeros(len(genes))
        for other_st in independent_types:
            total_contrib += np.maximum(all_contributions[other_st], 0)

        # Avoid division by zero
        total_contrib = np.maximum(total_contrib, 1e-10)
        type_fraction = np.maximum(this_contrib, 0) / total_contrib

        uniqueness = this_contrib * type_fraction

        df = pd.DataFrame({
            "gene": genes,
            "specificity": specificity[st].values,
            "gwas_zstat": zstat,
            "contribution": this_contrib,
            "type_fraction": type_fraction,
            "uniqueness": uniqueness,
        })
        df = df.sort_values("uniqueness", ascending=False).reset_index(drop=True)
        results[st] = df

    return results


def get_top_gene_sets(specificity, gwas_df, types, top_n=TOP_GENES_PER_TYPE):
    """
    Get top-N contributing genes per cell type (by contribution score).

    Parameters
    ----------
    specificity : pd.DataFrame
        Specificity matrix.
    gwas_df : pd.DataFrame
        GWAS results with SYMBOL, ZSTAT columns.
    types : list of str
        Cell type names.
    top_n : int
        Number of top genes per type.

    Returns
    -------
    dict
        Mapping supertype → set of gene symbols (top N).
    """
    gwas_indexed = gwas_df.set_index("SYMBOL")
    genes = specificity.index
    zstat = gwas_indexed.loc[genes, "ZSTAT"].values

    gene_sets = {}
    for st in types:
        contrib = specificity[st].values * zstat
        top_idx = np.argsort(contrib)[::-1][:top_n]
        gene_sets[st] = set(genes[top_idx])

    return gene_sets


def compute_jaccard_matrix(gene_sets):
    """
    Compute pairwise Jaccard similarity of gene sets.

    Jaccard = |A ∩ B| / |A ∪ B|

    Parameters
    ----------
    gene_sets : dict
        Mapping cell type → set of gene symbols.

    Returns
    -------
    pd.DataFrame
        Square matrix of Jaccard indices.
    """
    types = sorted(gene_sets.keys())
    n = len(types)
    matrix = np.zeros((n, n))

    for i, t1 in enumerate(types):
        for j, t2 in enumerate(types):
            if i == j:
                matrix[i, j] = 1.0
            elif j > i:
                intersection = len(gene_sets[t1] & gene_sets[t2])
                union = len(gene_sets[t1] | gene_sets[t2])
                jacc = intersection / union if union > 0 else 0
                matrix[i, j] = jacc
                matrix[j, i] = jacc

    return pd.DataFrame(matrix, index=types, columns=types)


def identify_shared_genes(gene_sets, min_types=None):
    """
    Identify genes shared across multiple cell types.

    Parameters
    ----------
    gene_sets : dict
        Mapping cell type → set of gene symbols.
    min_types : int, optional
        Minimum number of types a gene must appear in.
        Default: half the number of types.

    Returns
    -------
    pd.DataFrame
        Columns: gene, n_types, types_list
        Sorted by n_types (descending).
    """
    if min_types is None:
        min_types = max(2, len(gene_sets) // 2)

    # Count how many types each gene appears in
    gene_counts = {}
    gene_types = {}
    for st, genes in gene_sets.items():
        for g in genes:
            gene_counts[g] = gene_counts.get(g, 0) + 1
            gene_types.setdefault(g, []).append(st)

    shared = []
    for gene, count in gene_counts.items():
        if count >= min_types:
            shared.append({
                "gene": gene,
                "n_types": count,
                "types_list": ", ".join(sorted(gene_types[gene])),
            })

    result = pd.DataFrame(shared)
    if len(result) > 0:
        result = result.sort_values("n_types", ascending=False).reset_index(drop=True)

    return result


def identify_type_specific_genes(gene_sets, target_types, other_types, top_n=TOP_GENES_PER_TYPE):
    """
    Find genes in top-N of target types but NOT in top-N of other types.

    Useful for identifying SST-specific vs pan-interneuron genes.

    Parameters
    ----------
    gene_sets : dict
        Mapping cell type → set of gene symbols.
    target_types : list of str
        Types of interest (e.g., SST supertypes).
    other_types : list of str
        Comparison types (e.g., non-SST interneurons).
    top_n : int
        Already applied in gene_sets.

    Returns
    -------
    pd.DataFrame
        Genes specific to target types: gene, n_target_types, target_types_list
    """
    # Genes in target types
    target_genes = {}
    for st in target_types:
        if st in gene_sets:
            for g in gene_sets[st]:
                target_genes.setdefault(g, []).append(st)

    # Genes in other types
    other_genes = set()
    for st in other_types:
        if st in gene_sets:
            other_genes |= gene_sets[st]

    # Specific = in ≥2 target types AND 0 other types
    specific = []
    for gene, types in target_genes.items():
        if len(types) >= 2 and gene not in other_genes:
            specific.append({
                "gene": gene,
                "n_target_types": len(types),
                "target_types_list": ", ".join(sorted(types)),
            })

    result = pd.DataFrame(specific)
    if len(result) > 0:
        result = result.sort_values("n_target_types", ascending=False).reset_index(drop=True)

    return result
