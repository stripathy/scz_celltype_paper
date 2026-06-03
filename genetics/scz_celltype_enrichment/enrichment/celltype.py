"""
celltype.py — MAGMA-style cell-type enrichment analysis.

For each cell type, tests whether genes with higher cell-type specificity also
have stronger GWAS associations, using OLS regression:

    GWAS_ZSTAT ~ specificity + log(gene_size) + log(nsnps) + log(N)

The specificity coefficient is tested one-sided (positive direction), matching
MAGMA's direction=greater option. This tests the hypothesis that genes
specifically expressed in a cell type are enriched for genetic associations.

Multiple testing correction uses both Bonferroni and FDR (Benjamini-Hochberg).
"""
import numpy as np
import pandas as pd

from ..utils import ols_regression, fdr_correct, bonferroni_correct, Timer


def run_enrichment_all_types(specificity, zstat, covariates):
    """
    Run gene property analysis for every cell type in the specificity matrix.

    Parameters
    ----------
    specificity : pd.DataFrame
        Specificity matrix (genes × supertypes). Index = gene symbols.
    zstat : np.ndarray, shape (n_genes,)
        GWAS Z-statistics for the same genes, in the same order.
    covariates : np.ndarray, shape (n_genes, k)
        Covariate matrix (intercept, log_gene_size, log_nsnps, log_N).

    Returns
    -------
    pd.DataFrame
        Enrichment results sorted by p_value (ascending), with columns:
        supertype, beta, se, t_stat, p_value, n_genes, p_bonferroni, p_fdr
    """
    supertypes = specificity.columns.tolist()
    n_types = len(supertypes)
    n_genes = len(zstat)

    print(f"Running gene property analysis for {n_types} cell types "
          f"({n_genes:,} genes)...")

    results = []
    for i, st in enumerate(supertypes):
        spec_values = specificity[st].values

        # Design matrix: [covariates, specificity]
        X = np.column_stack([covariates, spec_values])

        res = ols_regression(zstat, X, covariate_index=-1, one_sided=True)
        results.append({
            "supertype": st,
            "beta": res["beta"],
            "se": res["se"],
            "t_stat": res["t_stat"],
            "p_value": res["p_value"],
            "n_genes": res["n_obs"],
        })

        if (i + 1) % 50 == 0:
            print(f"  Processed {i + 1}/{n_types} cell types...")

    results_df = pd.DataFrame(results)

    # Multiple testing correction
    results_df["p_bonferroni"] = bonferroni_correct(
        results_df["p_value"].values, n_tests=n_types
    )
    results_df["p_fdr"] = fdr_correct(results_df["p_value"].values)

    # Sort by p-value
    results_df = results_df.sort_values("p_value").reset_index(drop=True)

    # Print summary
    n_bonf = (results_df["p_bonferroni"] < 0.05).sum()
    n_fdr = (results_df["p_fdr"] < 0.05).sum()
    print(f"\nResults summary:")
    print(f"  Bonferroni significant (p < 0.05/{n_types}): {n_bonf}")
    print(f"  FDR significant (q < 0.05): {n_fdr}")

    return results_df


def get_significant_types(results_df, fdr_threshold=0.05):
    """
    Extract cell types passing FDR threshold.

    Parameters
    ----------
    results_df : pd.DataFrame
        Enrichment results from run_enrichment_all_types.
    fdr_threshold : float
        FDR threshold (default: 0.05).

    Returns
    -------
    list of str
        Supertype names with FDR < threshold.
    """
    sig = results_df[results_df["p_fdr"] < fdr_threshold]
    return sig["supertype"].tolist()
