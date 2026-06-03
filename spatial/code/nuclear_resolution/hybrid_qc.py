"""
Hybrid QC utilities for nuclear doublet resolution.

Functions extracted from optional nuclear doublet resolution pipeline for
reusability and testability.

Functions:
  infer_class_from_markers()  — Assign Glut/GABA/NN class from marker genes
  compute_hybrid_qc_pass()    — Build hybrid_qc_pass from nuclear evidence
"""

import numpy as np
import pandas as pd


def infer_class_from_markers(adata, gaba_markers=None, glut_markers=None):
    """Infer broad cell class (Glut/GABA/NN) from marker gene expression.

    Used for cells that don't have a corr_class assignment (e.g., QC-fail cells
    that were skipped by step 02b). This is a simple thresholding approach —
    sufficient for doublet detection which only needs class-level identity.

    Parameters
    ----------
    adata : AnnData
        Expression data with raw counts in .X
    gaba_markers : list of str, optional
        GABAergic marker genes. Defaults to standard panel.
    glut_markers : list of str, optional
        Glutamatergic marker genes. Defaults to standard panel.

    Returns
    -------
    class_labels : np.ndarray of str
        Per-cell class labels ('Glutamatergic', 'GABAergic', 'Non-neuronal')
    """
    from scipy import sparse as sp

    if gaba_markers is None:
        gaba_markers = ['GAD1', 'GAD2', 'SLC32A1', 'SST', 'PVALB', 'VIP',
                        'LAMP5']
    if glut_markers is None:
        glut_markers = ['SLC17A7', 'CUX2', 'RORB', 'GRIN2A', 'THEMIS',
                        'FEZF2', 'FOXP2']

    gene_names = list(adata.var_names)
    X = adata.X
    if sp.issparse(X):
        X = X.toarray()

    def _count_detected(markers):
        score = np.zeros(adata.n_obs, dtype=np.int8)
        for g in markers:
            if g in gene_names:
                score += (X[:, gene_names.index(g)] > 0).astype(np.int8)
        return score

    gaba_score = _count_detected(gaba_markers)
    glut_score = _count_detected(glut_markers)

    # Assign: whichever has more markers detected; ties -> Non-neuronal
    labels = np.full(adata.n_obs, 'Non-neuronal', dtype=object)
    labels[(glut_score > gaba_score) & (glut_score >= 2)] = 'Glutamatergic'
    labels[(gaba_score > glut_score) & (gaba_score >= 2)] = 'GABAergic'

    n_glut = (labels == 'Glutamatergic').sum()
    n_gaba = (labels == 'GABAergic').sum()
    n_nn = (labels == 'Non-neuronal').sum()
    print(f"    Inferred classes: {n_glut:,} Glut, {n_gaba:,} GABA, {n_nn:,} NN",
          flush=True)
    return labels


def compute_hybrid_qc_pass(adata, status):
    """Build hybrid_qc_pass boolean array from nuclear doublet evidence.

    Starts from basic per-cell QC but WITHOUT the high-UMI filter, then applies
    nuclear-informed doublet filtering:
      - Rescued: cells that failed ONLY due to high UMI are re-evaluated
      - Persistent/nuclear_only doublets -> FAIL
      - Resolved doublets -> PASS (key benefit: cytoplasmic spillover rescued)
      - Insufficient evidence -> FAIL (conservative)

    Parameters
    ----------
    adata : AnnData
        Must have .obs columns: qc_pass, and optionally fail_total_counts_high,
        corr_qc_pass, corr_subclass_margin, doublet_suspect.
    status : np.ndarray of str
        Nuclear doublet status per cell: 'clean', 'resolved', 'persistent',
        'nuclear_only', 'insufficient'.

    Returns
    -------
    hybrid_qc : np.ndarray of bool
        Per-cell hybrid QC pass/fail array.
    basic_qc : np.ndarray of bool
        Intermediate QC mask (non-doublet QC with high-UMI rescue).
        Useful for computing summary statistics.
    n_high_umi_rescued : int
        Number of cells re-evaluated due to high-UMI-only failure.
    """
    # Start from basic QC, but undo the high-UMI exclusion
    basic_qc = adata.obs['qc_pass'].values.astype(bool).copy()

    # Rescue high-UMI cells that aren't true doublets
    n_high_umi_rescued = 0
    if 'fail_total_counts_high' in adata.obs.columns:
        high_umi_only = (
            adata.obs['fail_total_counts_high'].values.astype(bool) &
            ~adata.obs.get('fail_neg_probe',
                           pd.Series(False, index=adata.obs.index)).astype(bool) &
            ~adata.obs.get('fail_neg_codeword',
                           pd.Series(False, index=adata.obs.index)).astype(bool) &
            ~adata.obs.get('fail_unassigned',
                           pd.Series(False, index=adata.obs.index)).astype(bool) &
            ~adata.obs.get('fail_n_genes_low',
                           pd.Series(False, index=adata.obs.index)).astype(bool) &
            ~adata.obs.get('fail_total_counts_low',
                           pd.Series(False, index=adata.obs.index)).astype(bool)
        )
        n_high_umi_rescued = int(high_umi_only.sum())
        # Re-include these cells — they failed ONLY due to high UMI
        basic_qc[high_umi_only] = True

    # Apply correlation margin filter (keep this — catches low-quality assignments)
    if 'corr_qc_pass' in adata.obs.columns:
        margin_fail = basic_qc.copy()
        margin_fail[:] = False
        if 'corr_subclass_margin' in adata.obs.columns:
            # Cells that had corr_qc_pass=False and doublet_suspect=False
            # failed on margin -> keep them excluded
            old_doublet = adata.obs.get(
                'doublet_suspect',
                pd.Series(False, index=adata.obs.index)
            ).values.astype(bool)
            old_corr_qc = adata.obs['corr_qc_pass'].values.astype(bool)
            margin_fail = ~old_corr_qc & ~old_doublet & basic_qc
        hybrid_qc = basic_qc & ~margin_fail
    else:
        hybrid_qc = basic_qc.copy()

    # Apply nuclear-informed doublet filter
    hybrid_qc[status == 'persistent'] = False
    hybrid_qc[status == 'nuclear_only'] = False
    # Resolved doublets PASS (this is the key benefit)
    hybrid_qc[status == 'resolved'] = True
    # Insufficient evidence — keep excluded (conservative)
    hybrid_qc[status == 'insufficient'] = False

    return hybrid_qc, basic_qc, n_high_umi_rescued
