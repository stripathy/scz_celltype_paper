"""Mixed effects models for cell count analysis."""

import pandas as pd
import numpy as np
import statsmodels.formula.api as smf


def fit_mixed_model(df, formula, group_col="subject"):
    """Fit a linear mixed effects model.

    Parameters
    ----------
    df : DataFrame
        Must contain all variables referenced in the formula plus group_col.
    formula : str
        Wilkinson formula, e.g. "SST ~ C(diagnosis, Treatment('Control')) * layer + age + sex_binary"
    group_col : str
        Random intercept grouping variable.

    Returns
    -------
    statsmodels MixedLMResults
    """
    model = smf.mixedlm(formula, df, groups=df[group_col])
    return model.fit(reml=True)


def run_binary_layer_analysis(df_R, df_L, covariates=("age", "sex_binary"),
                              ref_diag="Control"):
    """Run diagnosis × layer mixed models for all 4 cell types.

    Parameters
    ----------
    df_R : DataFrame
        R-section data with SST, VIP columns.
    df_L : DataFrame
        L-section data with PV, PYR columns.
    covariates : tuple of str
        Covariate names to include in the model.
    ref_diag : str
        Reference diagnosis group.

    Returns
    -------
    DataFrame with columns: cell_type, term, coef, p
    """
    cov_str = " + ".join(covariates) if covariates else ""
    results = []

    cell_configs = [
        ("SST", df_R, "SST"),
        ("VIP", df_R, "VIP"),
        ("PV", df_L, "PV"),
        ("PYR", df_L, "PYR"),
    ]

    for cell_type, df, col in cell_configs:
        # Build formula
        diag_term = f"C(diagnosis, Treatment('{ref_diag}'))"
        formula = f"{col} ~ {diag_term} * C(layer, Treatment('L23'))"
        if cov_str:
            formula += f" + {cov_str}"

        # Fit model
        result = fit_mixed_model(df, formula)

        # Extract coefficients with cleaned term names
        for term, row in result.summary().tables[1].iterrows():
            # Clean up term names
            clean_term = (
                str(term)
                .replace(f"C(diagnosis, Treatment('{ref_diag}'))", "")
                .replace("C(layer, Treatment('L23'))", "")
                .replace("[T.", "")
                .replace("]", "")
                .replace(":", ":")
            )
            # Fix specific patterns
            if clean_term.startswith("Bipolar"):
                clean_term = clean_term.replace("Bipolar", "BP")
            if clean_term.startswith("SCHIZ"):
                pass  # already clean
            if "L56" in clean_term and ":" in clean_term:
                parts = clean_term.split(":")
                clean_term = f"{parts[0]}:{parts[1]}" if parts[0] else f"L56:{parts[1]}"

            results.append({
                "cell_type": cell_type,
                "term": clean_term,
                "coef": result.params[term],
                "p": result.pvalues[term],
            })

    stats_df = pd.DataFrame(results)
    return stats_df


def _clean_term(raw_term, ref_diag="Control"):
    """Clean up statsmodels term names to readable format."""
    s = str(raw_term)
    s = s.replace(f"C(diagnosis, Treatment('{ref_diag}'))", "")
    s = s.replace("C(layer, Treatment('L23'))", "")
    s = s.replace("[T.", "").replace("]", "")

    # Rename Bipolar -> BP for brevity
    s = s.replace("Bipolar", "BP")
    return s


def run_binary_layer_analysis_clean(df_R, df_L, covariates=("age", "sex_binary"),
                                     ref_diag="Control"):
    """Run diagnosis × layer mixed models, returning clean results.

    Returns DataFrame with columns: cell_type, term, coef, p
    """
    cov_str = " + ".join(covariates) if covariates else ""
    results = []

    cell_configs = [
        ("SST", df_R, "SST"),
        ("VIP", df_R, "VIP"),
        ("PV", df_L, "PV"),
        ("PYR", df_L, "PYR"),
    ]

    for cell_type, df, col in cell_configs:
        diag_term = f"C(diagnosis, Treatment('{ref_diag}'))"
        formula = f"{col} ~ {diag_term} * C(layer, Treatment('L23'))"
        if cov_str:
            formula += f" + {cov_str}"

        # Drop rows with NaN in outcome and reset index for statsmodels
        model_df = df.dropna(subset=[col]).reset_index(drop=True)
        result = fit_mixed_model(model_df, formula)

        for term in result.params.index:
            clean = _clean_term(term, ref_diag)
            results.append({
                "cell_type": cell_type,
                "term": clean,
                "coef": result.params[term],
                "p": result.pvalues[term],
            })

    return pd.DataFrame(results)
