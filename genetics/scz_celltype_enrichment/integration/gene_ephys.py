"""
gene_ephys.py — Gene-electrophysiology correlations from patch-seq data.

Patch-seq records simultaneous electrophysiology and transcriptomics from the
same cell. This module computes correlations between gene expression and
electrophysiological features (especially sag and tau) to identify genes
associated with specific membrane properties.

Key features analyzed:
  - SAG: voltage sag during hyperpolarization, driven by HCN channels.
    High sag = strong HCN activity. Varies across cortical layers.
  - TAU: membrane time constant. Reflects membrane capacitance/resistance.
    Higher tau = slower electrical integration.

The module supports:
  1. Cell-level correlations (gene expression × ephys within SST cells)
  2. Supertype-level correlations (SEA-AD mean expression × mean ephys)
  3. Residualized correlations controlling for cortical layer
  4. Discordant gene identification (sag+/tau− genes)
"""
import numpy as np
import pandas as pd
from scipy import stats

from ..utils import fdr_correct, Timer


def load_patchseq_data(metadata_path, ephys_path, scanvi_path):
    """
    Load and merge patch-seq metadata, electrophysiology, and scANVI labels.

    Parameters
    ----------
    metadata_path : str or Path
        Path to LeeDalley_manuscript_metadata_v2.csv.
    ephys_path : str or Path
        Path to LeeDalley_ephys_fx.csv.
    scanvi_path : str or Path
        Path to iterative_scANVI_results_patchseq_only.2022-11-22.csv.

    Returns
    -------
    pd.DataFrame
        Merged dataset with ephys features and cell type labels.
    """
    with Timer("Loading patch-seq data"):
        metadata = pd.read_csv(str(metadata_path))
        ephys = pd.read_csv(str(ephys_path))
        scanvi = pd.read_csv(str(scanvi_path), index_col=0)

    # Merge metadata and ephys on specimen_id
    # metadata has 'specimen_id_x', ephys has 'specimen_id'
    if "specimen_id_x" in metadata.columns:
        metadata = metadata.rename(columns={"specimen_id_x": "specimen_id"})

    merged = metadata.merge(ephys, on="specimen_id", how="inner")
    print(f"    Metadata + ephys: {len(merged)} cells")

    # Add scANVI labels (index is a different ID format, need to align)
    # scANVI index is like 'SM-GE4W8_S016_E1-50'
    # Match via patched_cell_container or other shared ID
    # For now, return merged data; scANVI labels will be matched in the script
    return merged, scanvi


def compute_gene_ephys_correlations(
    gene_values,
    ephys_values,
    gene_names=None,
    min_expressing=10,
):
    """
    Compute Spearman correlation between each gene's expression and an
    electrophysiology feature across cells.

    Parameters
    ----------
    gene_values : np.ndarray, shape (n_cells, n_genes)
        Gene expression values per cell.
    ephys_values : np.ndarray, shape (n_cells,)
        Electrophysiology feature values per cell.
    gene_names : list of str, optional
        Gene names corresponding to columns.
    min_expressing : int
        Minimum number of cells expressing a gene (>0) to compute correlation.

    Returns
    -------
    pd.DataFrame
        Columns: gene, spearman_rho, pval, n_expressing, pval_fdr, abs_rho
        Sorted by abs_rho (descending).
    """
    n_cells, n_genes = gene_values.shape

    # Mask: cells with valid ephys values
    valid = ~np.isnan(ephys_values)
    ephys_valid = ephys_values[valid]

    results = []
    for i in range(n_genes):
        expr = gene_values[valid, i]
        n_expressing = np.sum(expr > 0)

        if n_expressing < min_expressing:
            continue

        rho, pval = stats.spearmanr(expr, ephys_valid)

        name = gene_names[i] if gene_names is not None else f"gene_{i}"
        results.append({
            "gene": name,
            "spearman_rho": rho,
            "pval": pval,
            "n_expressing": int(n_expressing),
        })

    df = pd.DataFrame(results)

    if len(df) > 0:
        df["pval_fdr"] = fdr_correct(df["pval"].values)
        df["abs_rho"] = df["spearman_rho"].abs()
        df = df.sort_values("abs_rho", ascending=False).reset_index(drop=True)

    return df


def compute_supertype_correlations(
    specificity_or_mean_expr,
    ephys_by_supertype,
    feature_name="sag",
):
    """
    Compute Spearman correlation between gene expression (from SEA-AD reference)
    and mean electrophysiology per supertype.

    This is a supertype-level analysis: each data point is one supertype,
    and we correlate the gene's specificity with the type's mean ephys value.

    Parameters
    ----------
    specificity_or_mean_expr : pd.DataFrame
        Specificity or mean expression matrix (genes × supertypes).
    ephys_by_supertype : pd.Series
        Mean ephys feature per supertype. Index = supertype names.
    feature_name : str
        Name of the feature for labeling.

    Returns
    -------
    pd.DataFrame
        Columns: gene, spearman_rho, pval, pval_fdr, abs_rho
        Sorted by abs_rho (descending).
    """
    # Align supertypes
    shared_types = sorted(
        set(specificity_or_mean_expr.columns) & set(ephys_by_supertype.index)
    )
    print(f"  Computing {feature_name} correlations across {len(shared_types)} supertypes")

    expr = specificity_or_mean_expr[shared_types]
    ephys = ephys_by_supertype[shared_types].values

    results = []
    for gene in expr.index:
        gene_vals = expr.loc[gene].values
        if np.all(gene_vals == 0):
            continue
        rho, pval = stats.spearmanr(gene_vals, ephys)
        results.append({
            "gene": gene,
            "spearman_rho": rho,
            "pval": pval,
        })

    df = pd.DataFrame(results)
    if len(df) > 0:
        df["pval_fdr"] = fdr_correct(df["pval"].values)
        df["abs_rho"] = df["spearman_rho"].abs()
        df = df.sort_values("abs_rho", ascending=False).reset_index(drop=True)

    print(f"    Computed correlations for {len(df):,} genes")
    return df


def compute_residualized_correlations(
    gene_values,
    ephys_values,
    confound_values,
    gene_names=None,
    min_expressing=10,
):
    """
    Compute Spearman correlations after residualizing ephys on a confound.

    This controls for confounds like cortical layer: SAG varies with depth,
    so we want to find genes correlated with SAG *independent* of depth.

    Method: Regress ephys ~ confound, take residuals, then correlate with gene expr.

    Parameters
    ----------
    gene_values : np.ndarray, shape (n_cells, n_genes)
        Gene expression values.
    ephys_values : np.ndarray, shape (n_cells,)
        Electrophysiology feature values.
    confound_values : np.ndarray, shape (n_cells,)
        Confound variable (e.g., cortical layer).
    gene_names : list of str, optional
        Gene names.
    min_expressing : int
        Minimum cells expressing a gene.

    Returns
    -------
    pd.DataFrame
        Same format as compute_gene_ephys_correlations, but on residualized ephys.
    """
    # Remove NaN rows
    valid = ~(np.isnan(ephys_values) | np.isnan(confound_values))
    ephys_v = ephys_values[valid]
    confound_v = confound_values[valid]
    gene_v = gene_values[valid]

    # Residualize: regress ephys ~ confound, take residuals
    X = np.column_stack([np.ones(len(confound_v)), confound_v])
    beta = np.linalg.lstsq(X, ephys_v, rcond=None)[0]
    ephys_resid = ephys_v - X @ beta

    return compute_gene_ephys_correlations(
        gene_v, ephys_resid, gene_names=gene_names, min_expressing=min_expressing
    )


def find_discordant_genes(
    sag_corr_df,
    tau_corr_df,
    sag_direction="+",
    tau_direction="-",
    p_threshold=0.05,
):
    """
    Find genes with discordant sag/tau correlations (e.g., sag+/tau−).

    These genes are associated with cells that have strong HCN activity (high sag)
    but fast membrane dynamics (short tau). HCN1 is the prototypical example.

    Parameters
    ----------
    sag_corr_df : pd.DataFrame
        SAG correlation results (from compute_supertype_correlations).
    tau_corr_df : pd.DataFrame
        TAU correlation results.
    sag_direction : str
        '+' for positive sag correlation, '-' for negative.
    tau_direction : str
        '+' for positive tau correlation, '-' for negative.
    p_threshold : float
        Nominal p-value threshold for inclusion.

    Returns
    -------
    pd.DataFrame
        Discordant genes with columns:
        gene, spearman_rho_sag, pval_sag, pval_fdr_sag,
        spearman_rho_tau, pval_tau, pval_fdr_tau, combined_score
    """
    # Merge on gene
    sag = sag_corr_df.rename(columns={
        "spearman_rho": "spearman_rho_sag",
        "pval": "pval_sag",
        "pval_fdr": "pval_fdr_sag",
    })[["gene", "spearman_rho_sag", "pval_sag", "pval_fdr_sag"]]

    tau = tau_corr_df.rename(columns={
        "spearman_rho": "spearman_rho_tau",
        "pval": "pval_tau",
        "pval_fdr": "pval_fdr_tau",
    })[["gene", "spearman_rho_tau", "pval_tau", "pval_fdr_tau"]]

    merged = sag.merge(tau, on="gene", how="inner")

    # Filter by direction and significance
    if sag_direction == "+":
        sag_mask = merged["spearman_rho_sag"] > 0
    else:
        sag_mask = merged["spearman_rho_sag"] < 0

    if tau_direction == "-":
        tau_mask = merged["spearman_rho_tau"] < 0
    else:
        tau_mask = merged["spearman_rho_tau"] > 0

    p_mask = (merged["pval_sag"] < p_threshold) & (merged["pval_tau"] < p_threshold)

    result = merged[sag_mask & tau_mask & p_mask].copy()

    # Combined score: sum of absolute correlations
    result["combined_score"] = (
        result["spearman_rho_sag"].abs() + result["spearman_rho_tau"].abs()
    )

    result = result.sort_values("combined_score", ascending=False).reset_index(drop=True)

    print(f"  Found {len(result)} discordant genes "
          f"(sag{sag_direction}/tau{tau_direction}, p<{p_threshold})")

    return result


def annotate_with_gwas(gene_df, scz_gene_set):
    """
    Add is_scz_gwas boolean column to a gene DataFrame.

    Parameters
    ----------
    gene_df : pd.DataFrame
        Must have 'gene' column.
    scz_gene_set : set
        Set of SCZ GWAS gene symbols.

    Returns
    -------
    pd.DataFrame
        Input with added 'is_scz' column.
    """
    df = gene_df.copy()
    df["is_scz"] = df["gene"].isin(scz_gene_set)
    n_scz = df["is_scz"].sum()
    print(f"  {n_scz}/{len(df)} genes are in SCZ GWAS gene set "
          f"({n_scz/len(df)*100:.1f}%)")
    return df
