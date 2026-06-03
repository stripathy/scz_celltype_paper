"""
Cell QC filtering for Xenium spatial transcriptomics data.

Implements the QC approach from Kwon, Guo et al. (2026):
  1. Compute per-cell sums for 3 control feature types:
     - Negative Control Probes (20 features)
     - Negative Control Codewords (41 features)
     - Unassigned Codewords (180 features)
  2. Flag cells exceeding 99th percentile in any control type
  3. Flag cells with n_genes < 5 MADs below median
  4. Flag cells with total_counts outside 5 MADs from median
"""

import numpy as np
import pandas as pd
import h5py
from scipy.sparse import csc_matrix


# Feature type strings in the raw Xenium h5 files
GENE_EXPRESSION = "Gene Expression"
NEG_CONTROL_PROBE = "Negative Control Probe"
NEG_CONTROL_CODEWORD = "Negative Control Codeword"
UNASSIGNED_CODEWORD = "Unassigned Codeword"

CONTROL_TYPES = [NEG_CONTROL_PROBE, NEG_CONTROL_CODEWORD, UNASSIGNED_CODEWORD]


def compute_qc_metrics(h5_path):
    """
    Load a raw Xenium h5 file and compute per-cell QC metrics.

    Reads ALL features (genes + control probes + codewords), then computes
    per-cell summary statistics for each feature type.

    Parameters
    ----------
    h5_path : str
        Path to a GSM*-cell_feature_matrix.h5 file.

    Returns
    -------
    pd.DataFrame
        Index = cell barcodes, columns include:
        - total_counts: sum of gene expression features
        - n_genes: number of non-zero gene expression features
        - neg_probe_sum: sum of Negative Control Probe features
        - neg_codeword_sum: sum of Negative Control Codeword features
        - unassigned_sum: sum of Unassigned Codeword features
    """
    with h5py.File(h5_path, "r") as f:
        matrix = f["matrix"]
        barcodes = [b.decode() for b in matrix["barcodes"][:]]
        features = matrix["features"]
        feature_names = [g.decode() for g in features["name"][:]]
        feature_types = [t.decode() for t in features["feature_type"][:]]

        data = matrix["data"][:]
        indices = matrix["indices"][:]
        indptr = matrix["indptr"][:]
        shape = matrix["shape"][:]

        # CSC -> CSR (cells x features)
        X = csc_matrix((data, indices, indptr), shape=shape).T

    feature_types = np.array(feature_types)
    n_cells = len(barcodes)

    # Build masks for each feature type
    gene_mask = feature_types == GENE_EXPRESSION
    neg_probe_mask = feature_types == NEG_CONTROL_PROBE
    neg_codeword_mask = feature_types == NEG_CONTROL_CODEWORD
    unassigned_mask = feature_types == UNASSIGNED_CODEWORD

    # Compute per-cell metrics
    metrics = pd.DataFrame(index=barcodes)

    # Gene expression metrics
    X_genes = X[:, gene_mask]
    metrics["total_counts"] = np.asarray(X_genes.sum(axis=1)).ravel()
    metrics["n_genes"] = np.asarray((X_genes > 0).sum(axis=1)).ravel()

    # Control feature sums
    if neg_probe_mask.sum() > 0:
        metrics["neg_probe_sum"] = np.asarray(
            X[:, neg_probe_mask].sum(axis=1)
        ).ravel()
    else:
        metrics["neg_probe_sum"] = 0

    if neg_codeword_mask.sum() > 0:
        metrics["neg_codeword_sum"] = np.asarray(
            X[:, neg_codeword_mask].sum(axis=1)
        ).ravel()
    else:
        metrics["neg_codeword_sum"] = 0

    if unassigned_mask.sum() > 0:
        metrics["unassigned_sum"] = np.asarray(
            X[:, unassigned_mask].sum(axis=1)
        ).ravel()
    else:
        metrics["unassigned_sum"] = 0

    # Summary
    n_neg_probes = neg_probe_mask.sum()
    n_neg_codewords = neg_codeword_mask.sum()
    n_unassigned = unassigned_mask.sum()
    print(f"    Features: {gene_mask.sum()} genes, "
          f"{n_neg_probes} neg probes, "
          f"{n_neg_codewords} neg codewords, "
          f"{n_unassigned} unassigned codewords")

    return metrics


def _mad(x):
    """Compute Median Absolute Deviation."""
    med = np.median(x)
    return np.median(np.abs(x - med))


def flag_qc_failures(metrics, control_pctl=99, mad_threshold=5):
    """
    Flag cells failing QC criteria, matching Kwon et al. approach.

    Parameters
    ----------
    metrics : pd.DataFrame
        Output from compute_qc_metrics().
    control_pctl : int
        Percentile threshold for control features (default: 99).
    mad_threshold : float
        Number of MADs for outlier detection (default: 5).

    Returns
    -------
    pd.DataFrame
        Input DataFrame with added columns:
        - fail_neg_probe: exceeds control percentile for neg probes
        - fail_neg_codeword: exceeds control percentile for neg codewords
        - fail_unassigned: exceeds control percentile for unassigned codewords
        - fail_n_genes_low: n_genes < median - mad_threshold * MAD
        - fail_total_counts_low: total_counts < median - mad_threshold * MAD
        - fail_total_counts_high: total_counts > median + mad_threshold * MAD
        - qc_pass: True if no failures
    """
    df = metrics.copy()

    # --- Control feature thresholds (99th percentile) ---
    # Neg probes
    thresh_neg_probe = np.percentile(df["neg_probe_sum"], control_pctl)
    df["fail_neg_probe"] = df["neg_probe_sum"] > thresh_neg_probe

    # Neg codewords
    thresh_neg_codeword = np.percentile(df["neg_codeword_sum"], control_pctl)
    df["fail_neg_codeword"] = df["neg_codeword_sum"] > thresh_neg_codeword

    # Unassigned codewords
    thresh_unassigned = np.percentile(df["unassigned_sum"], control_pctl)
    df["fail_unassigned"] = df["unassigned_sum"] > thresh_unassigned

    # --- MAD-based outlier detection ---
    # n_genes: flag if below median - 5*MAD
    med_ngenes = np.median(df["n_genes"])
    mad_ngenes = _mad(df["n_genes"])
    ngenes_lower = med_ngenes - mad_threshold * mad_ngenes
    df["fail_n_genes_low"] = df["n_genes"] < ngenes_lower

    # total_counts: flag if outside median +/- 5*MAD
    med_counts = np.median(df["total_counts"])
    mad_counts = _mad(df["total_counts"])
    counts_lower = med_counts - mad_threshold * mad_counts
    counts_upper = med_counts + mad_threshold * mad_counts
    df["fail_total_counts_low"] = df["total_counts"] < counts_lower
    df["fail_total_counts_high"] = df["total_counts"] > counts_upper

    # --- Combined QC pass ---
    fail_cols = [c for c in df.columns if c.startswith("fail_")]
    df["qc_pass"] = ~df[fail_cols].any(axis=1)

    # Print summary
    n_total = len(df)
    n_pass = df["qc_pass"].sum()
    n_fail = n_total - n_pass

    print(f"    QC thresholds:")
    print(f"      Neg probe 99th pctl:     {thresh_neg_probe:.0f}")
    print(f"      Neg codeword 99th pctl:  {thresh_neg_codeword:.0f}")
    print(f"      Unassigned 99th pctl:    {thresh_unassigned:.0f}")
    print(f"      n_genes lower (med-5MAD): {ngenes_lower:.1f} "
          f"(median={med_ngenes:.0f}, MAD={mad_ngenes:.1f})")
    print(f"      total_counts range:       [{counts_lower:.1f}, {counts_upper:.1f}] "
          f"(median={med_counts:.0f}, MAD={mad_counts:.1f})")

    print(f"    Failures:")
    for col in fail_cols:
        n = df[col].sum()
        pct = 100 * n / n_total
        print(f"      {col:30s}: {n:>6,} ({pct:.2f}%)")
    print(f"    Overall: {n_pass:,}/{n_total:,} pass ({100*n_pass/n_total:.1f}%), "
          f"{n_fail:,} fail ({100*n_fail/n_total:.1f}%)")

    return df
