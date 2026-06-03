"""
conditional.py — Conditional analysis to identify independently enriched cell types.

Problem: Many cell types share gene expression programs (e.g., all SST interneurons
express similar genes). MAGMA-style marginal tests may identify multiple related
types that are really driven by the same underlying genetic signal.

Solution: Forward stepwise selection. Starting from all FDR-significant types,
iteratively add cell types that remain significant (p < 0.05) when conditioned
on ALL previously selected types. This yields a set of independently enriched
cell types whose signals are not redundant.

Implementation: Pure Python OLS (no MAGMA binary dependency). The conditional test
adds previously-selected cell types' specificity as additional covariates in the
regression, then tests whether the new type's specificity coefficient is still
significantly positive.

Note: This simplified forward selection differs from MAGMA's pairwise proportionality
score approach, but qualitatively identifies the same independent signals.
"""
import numpy as np
import pandas as pd

from ..utils import ols_regression, Timer, extract_subclass


def run_conditional_test(target, conditions, specificity, zstat, covariates):
    """
    Test whether 'target' cell type remains significant after conditioning
    on all cell types in 'conditions'.

    Design matrix: [covariates, cond_1_spec, ..., cond_k_spec, target_spec]
    Tests the last column (target specificity), one-sided positive.

    Parameters
    ----------
    target : str
        Supertype name to test.
    conditions : list of str
        Supertype names to condition on.
    specificity : pd.DataFrame
        Full specificity matrix.
    zstat : np.ndarray
        GWAS Z-statistics.
    covariates : np.ndarray
        Base covariate matrix (intercept + confounds).

    Returns
    -------
    dict
        Keys: beta, se, t_stat, p_value, n_obs, n_params
    """
    X_parts = [covariates]

    # Add conditioning cell types
    for ct in conditions:
        X_parts.append(specificity[ct].values.reshape(-1, 1))

    # Add target cell type (last column — this is what we test)
    X_parts.append(specificity[target].values.reshape(-1, 1))

    X = np.column_stack(X_parts)
    return ols_regression(zstat, X, covariate_index=-1, one_sided=True)


def forward_selection(
    significant_types,
    marginal_results,
    specificity,
    zstat,
    covariates,
    p_threshold=0.05,
):
    """
    Forward stepwise selection to identify independently enriched cell types.

    Algorithm:
      1. Sort significant types by marginal p-value (ascending)
      2. Accept the most significant type as the first independent type
      3. For each remaining type (in order of marginal significance):
         - Test it conditioned on all already-accepted independent types
         - If conditional p < p_threshold, accept it as independent
      4. Return the set of independent types with both marginal and conditional stats

    Parameters
    ----------
    significant_types : list of str
        Supertype names that passed FDR threshold in marginal analysis.
    marginal_results : pd.DataFrame
        Full marginal enrichment results (from enrichment.run_enrichment_all_types).
    specificity : pd.DataFrame
        Full specificity matrix.
    zstat : np.ndarray
        GWAS Z-statistics.
    covariates : np.ndarray
        Base covariate matrix.
    p_threshold : float
        Conditional p-value threshold (default: 0.05).

    Returns
    -------
    pd.DataFrame
        Independent supertypes with columns:
        supertype, subclass, marginal_p, marginal_t,
        conditional_p, conditional_t, conditional_beta
    """
    # Sort candidates by marginal p-value
    marginal_lookup = marginal_results.set_index("supertype")
    candidates = sorted(
        significant_types,
        key=lambda x: marginal_lookup.loc[x, "p_value"]
    )

    print(f"Forward selection: {len(candidates)} candidates, p_threshold={p_threshold}")
    print("-" * 70)

    independent = []
    independent_names = []

    for i, candidate in enumerate(candidates):
        marg_p = marginal_lookup.loc[candidate, "p_value"]
        marg_t = marginal_lookup.loc[candidate, "t_stat"]

        if len(independent_names) == 0:
            # First type: automatically accept
            cond_p = marg_p
            cond_t = marg_t
            cond_beta = marginal_lookup.loc[candidate, "beta"]
            accepted = True
        else:
            # Test conditioned on all previously accepted types
            res = run_conditional_test(
                candidate, independent_names, specificity, zstat, covariates
            )
            cond_p = res["p_value"]
            cond_t = res["t_stat"]
            cond_beta = res["beta"]
            accepted = cond_p < p_threshold

        status = "ACCEPTED" if accepted else "rejected"
        print(f"  [{i+1:2d}] {candidate:25s} | marginal p={marg_p:.2e} | "
              f"conditional p={cond_p:.2e} | {status}")

        if accepted:
            independent_names.append(candidate)
            independent.append({
                "supertype": candidate,
                "subclass": extract_subclass(candidate),
                "marginal_p": marg_p,
                "marginal_t": marg_t,
                "conditional_p": cond_p,
                "conditional_t": cond_t,
                "conditional_beta": cond_beta,
            })

    print("-" * 70)
    print(f"Result: {len(independent)} independent cell types "
          f"(from {len(candidates)} candidates)")

    return pd.DataFrame(independent)
