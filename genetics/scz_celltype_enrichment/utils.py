"""
utils.py — Shared statistical and I/O utility functions.

Core OLS regression is implemented via np.linalg.lstsq with manual standard error
computation, matching the approach used in the original exploratory analysis and
staying close to MAGMA's internal gene property test implementation.
"""
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from pathlib import Path
import time


# ============================================================================
# Statistical functions
# ============================================================================

def ols_regression(y, X, covariate_index=-1, one_sided=True):
    """
    Run OLS regression and return inference for a specific covariate.

    Uses np.linalg.lstsq with manual SE computation (not statsmodels.OLS)
    to match the numerics from the exploratory analysis and MAGMA's approach.

    Parameters
    ----------
    y : np.ndarray, shape (n,)
        Response variable (e.g., GWAS Z-statistics).
    X : np.ndarray, shape (n, k)
        Design matrix. Should include intercept, covariates, and the
        variable of interest.
    covariate_index : int
        Column index of the covariate to test (default: -1 = last column).
    one_sided : bool
        If True, return one-sided p-value testing for POSITIVE effect.
        This matches MAGMA's direction=greater option.

    Returns
    -------
    dict with keys:
        beta : float — regression coefficient
        se : float — standard error
        t_stat : float — t-statistic
        p_value : float — p-value (one-sided or two-sided)
        n_obs : int — number of observations
        n_params : int — number of parameters
    """
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    residuals = y - X @ beta
    n_obs = len(y)
    n_params = X.shape[1]
    df = n_obs - n_params

    # Residual mean squared error
    mse = np.sum(residuals ** 2) / df

    # Standard error of the coefficient of interest
    XtX_inv = np.linalg.inv(X.T @ X)
    se = np.sqrt(mse * XtX_inv[covariate_index, covariate_index])

    t_stat = beta[covariate_index] / se

    if one_sided:
        # One-sided test: P(T > t_obs), testing positive direction
        p_value = 1.0 - stats.t.cdf(t_stat, df=df)
    else:
        p_value = 2.0 * (1.0 - stats.t.cdf(abs(t_stat), df=df))

    return {
        "beta": float(beta[covariate_index]),
        "se": float(se),
        "t_stat": float(t_stat),
        "p_value": float(p_value),
        "n_obs": n_obs,
        "n_params": n_params,
    }


def fdr_correct(pvalues, method="fdr_bh"):
    """
    Apply FDR (Benjamini-Hochberg) correction to an array of p-values.

    Parameters
    ----------
    pvalues : array-like
        Raw p-values.
    method : str
        Method for multipletests (default: 'fdr_bh').

    Returns
    -------
    np.ndarray
        FDR-corrected p-values.
    """
    pvals = np.asarray(pvalues, dtype=float)
    # Handle NaN p-values
    mask = ~np.isnan(pvals)
    corrected = np.full_like(pvals, np.nan)
    if mask.sum() > 0:
        _, corrected[mask], _, _ = multipletests(pvals[mask], method=method)
    return corrected


def bonferroni_correct(pvalues, n_tests=None):
    """
    Apply Bonferroni correction to an array of p-values.

    Parameters
    ----------
    pvalues : array-like
        Raw p-values.
    n_tests : int, optional
        Number of tests (default: length of pvalues).

    Returns
    -------
    np.ndarray
        Bonferroni-corrected p-values (capped at 1.0).
    """
    pvals = np.asarray(pvalues, dtype=float)
    if n_tests is None:
        n_tests = len(pvals)
    return np.minimum(pvals * n_tests, 1.0)


# ============================================================================
# I/O helpers
# ============================================================================

def ensure_dirs(*dirs):
    """Create directories if they don't exist."""
    for d in dirs:
        Path(d).mkdir(parents=True, exist_ok=True)


def load_or_compute(path, compute_fn, force=False):
    """
    Load cached result from path if it exists; otherwise run compute_fn,
    save the result, and return it.

    Parameters
    ----------
    path : str or Path
        Path to the cached CSV file.
    compute_fn : callable
        Function that returns a pd.DataFrame to save.
    force : bool
        If True, recompute even if cache exists.

    Returns
    -------
    pd.DataFrame
    """
    path = Path(path)
    if path.exists() and not force:
        print(f"  Loading cached: {path.name}")
        return pd.read_csv(path, index_col=0)
    else:
        print(f"  Computing (no cache at {path.name})...")
        t0 = time.time()
        result = compute_fn()
        result.to_csv(path)
        print(f"  Saved to {path} ({time.time() - t0:.1f}s)")
        return result


class Timer:
    """Simple context manager for timing code blocks."""

    def __init__(self, label=""):
        self.label = label

    def __enter__(self):
        self.t0 = time.time()
        if self.label:
            print(f"  {self.label}...", end=" ", flush=True)
        return self

    def __exit__(self, *args):
        elapsed = time.time() - self.t0
        print(f"done ({elapsed:.1f}s)")
        self.elapsed = elapsed


# ============================================================================
# Cell-type classification helpers
# ============================================================================

def classify_supertypes_by_class(supertypes):
    """
    Classify supertypes into broad cell classes based on naming conventions.

    Returns dict mapping supertype → class label.
    Class labels: 'GABAergic', 'Glutamatergic', 'Non-neuronal'
    """
    classes = {}
    gaba_prefixes = [
        "Lamp5", "Pvalb", "Sst", "Vip", "Sncg", "Pax6", "Chandelier"
    ]
    glut_prefixes = [
        "L2/3 IT", "L4 IT", "L5 ET", "L5 IT", "L5/6 NP",
        "L6 CT", "L6 IT", "L6b"
    ]

    for st in supertypes:
        if any(st.startswith(p) for p in gaba_prefixes):
            classes[st] = "GABAergic"
        elif any(st.startswith(p) for p in glut_prefixes):
            classes[st] = "Glutamatergic"
        else:
            classes[st] = "Non-neuronal"

    return classes


def extract_subclass(supertype):
    """
    Extract the subclass name from a supertype string.

    Examples: 'Sst_2' → 'Sst', 'L2/3 IT_1' → 'L2/3 IT', 'Lamp5_Lhx6_1' → 'Lamp5 Lhx6'
    """
    # Handle special cases
    if "Sst Chodl" in supertype:
        return "Sst Chodl"
    if "Lamp5_Lhx6" in supertype:
        return "Lamp5 Lhx6"

    # General: remove trailing _N (number)
    parts = supertype.rsplit("_", 1)
    if len(parts) == 2 and parts[1].replace("-SEAAD", "").isdigit():
        return parts[0]

    # Handle -SEAAD suffix
    if "-SEAAD" in supertype:
        return supertype.split("_")[0]

    return supertype
