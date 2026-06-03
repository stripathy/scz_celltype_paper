"""
BANKSY-based spatial domain classification.

Uses BANKSY (Nature Genetics 2024) to cluster cells by gene expression +
spatial context, then classifies each cluster as:
  - Cortical: normal cortical tissue (including L1 border cells)
  - Vascular: >50% Endothelial + VLMC
  - WM: >40% Oligodendrocyte AND mean_depth > 0.80

L1 border detection: shallow non-neuronal-dominated BANKSY clusters are
recognized as L1 cortex (not meningeal/extra-cortical). Confirmed by
MERFISH comparison: L1 has ~81% non-neuronal composition. Flagged via
the `is_l1` output for downstream use.

This replaces the older K-NN composition → PCA → Leiden approach in
spatial_domains.py, which misclassified L1 cells as "Extra-cortical"
and lacked WM detection.

Requires: pybanksy (pip install pybanksy)

Typical usage:
    from banksy_domains import preprocess_for_banksy, run_banksy, classify_banksy_domains

    adata_b = preprocess_for_banksy(adata)
    labels = run_banksy(adata_b)
    domains, is_l1, info = classify_banksy_domains(adata, labels)
"""

import warnings
import numpy as np
import pandas as pd
import scanpy as sc

from spatial_domains import VASCULAR_TYPES, NON_NEURONAL_TYPES

# ── BANKSY clustering parameters ─────────────────────────────────────

BANKSY_LAMBDA = 0.8
BANKSY_RESOLUTION = 0.3
K_GEOM = 15
PCA_DIMS = [20]

# ── Domain classification thresholds ─────────────────────────────────

VASCULAR_THRESH = 0.50       # vasc_frac > 0.50 → Vascular
WM_OLIGO_THRESH = 0.40       # oligo_frac > 0.40 AND deep → WM
WM_DEPTH_THRESH = 0.80       # mean_depth > 0.80 for WM

# L1 border detection (formerly "Meningeal" — actually L1 cortex)
L1_NN_THRESH = 0.50           # non-neuronal fraction > 0.50
L1_DEPTH_THRESH = 0.20        # mean_depth < 0.20 (cortical surface)


def preprocess_for_banksy(adata):
    """Normalize, log-transform, and z-score for BANKSY input.

    Parameters
    ----------
    adata : anndata.AnnData
        Raw or count-based AnnData. Must have spatial coords in
        .obsm['spatial'].

    Returns
    -------
    anndata.AnnData
        Copy with normalized + scaled expression. Original counts
        stored in .layers['counts'].
    """
    adata_b = adata.copy()
    adata_b.layers["counts"] = adata_b.X.copy()
    sc.pp.normalize_total(adata_b, target_sum=1e4)
    sc.pp.log1p(adata_b)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.scale(adata_b)
    return adata_b


def run_banksy(adata_b, lambda_param=None, resolution=None, k_geom=None):
    """Run BANKSY clustering on preprocessed data.

    Parameters
    ----------
    adata_b : anndata.AnnData
        Preprocessed AnnData (from preprocess_for_banksy). Must have
        spatial coords in .obsm['spatial'].
    lambda_param : float, optional
        BANKSY spatial weighting parameter. Default: BANKSY_LAMBDA (0.8).
    resolution : float, optional
        Leiden clustering resolution. Default: BANKSY_RESOLUTION (0.3).
    k_geom : int, optional
        Number of spatial neighbors for BANKSY. Default: K_GEOM (15).

    Returns
    -------
    np.ndarray of int
        Cluster label per cell.
    """
    from banksy.initialize_banksy import initialize_banksy
    from banksy.embed_banksy import generate_banksy_matrix
    from banksy.cluster_methods import run_Leiden_partition
    from banksy_utils.umap_pca import pca_umap

    if lambda_param is None:
        lambda_param = BANKSY_LAMBDA
    if resolution is None:
        resolution = BANKSY_RESOLUTION
    if k_geom is None:
        k_geom = K_GEOM

    coord_keys = ("x", "y", "spatial")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        banksy_dict = initialize_banksy(
            adata_b, coord_keys, k_geom,
            nbr_weight_decay="scaled_gaussian", max_m=1,
            plt_edge_hist=False, plt_nbr_weights=False,
            plt_agf_angles=False,
        )
        banksy_dict, _ = generate_banksy_matrix(
            adata_b, banksy_dict, [lambda_param], max_m=1
        )
        pca_umap(banksy_dict, pca_dims=PCA_DIMS, add_umap=False)
        results_df, _ = run_Leiden_partition(
            banksy_dict, [resolution],
            num_nn=50, num_iterations=-1, partition_seed=42,
            match_labels=False,
        )

    labels = results_df.iloc[0]["labels"]
    if hasattr(labels, "dense"):
        labels = labels.dense
    return np.asarray(labels).astype(int)


def classify_banksy_domains(adata, banksy_labels,
                             subclass_col=None,
                             depth_col='predicted_norm_depth',
                             vascular_thresh=None,
                             wm_oligo_thresh=None,
                             wm_depth_thresh=None,
                             l1_nn_thresh=None,
                             l1_depth_thresh=None,
                             verbose=True):
    """Classify BANKSY clusters into spatial domains and flag L1 border.

    For each cluster, computes cell type composition and mean predicted
    depth, then applies rule-based classification:
      1. Vascular: vasc_frac > vascular_thresh
      2. White matter: oligo_frac > wm_oligo_thresh AND mean_depth > wm_depth_thresh
      3. L1 border: nn_frac > l1_nn_thresh AND mean_depth < l1_depth_thresh
         → classified as Cortical with is_l1=True
      4. Neuronal cortex: neuronal_frac > 0.20 AND 0 ≤ depth ≤ 0.90
      5. Deep WM: mean_depth > wm_depth_thresh
      6. Default: Cortical

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated sample with subclass labels and depth predictions.
    banksy_labels : np.ndarray of int
        Cluster label per cell from run_banksy().
    subclass_col : str, optional
        Column for subclass labels. Tries 'corr_subclass' then
        'subclass_label' if None.
    depth_col : str
        Column for predicted normalized depth.
    vascular_thresh, wm_oligo_thresh, wm_depth_thresh : float, optional
        Classification thresholds. Defaults to module constants.
    l1_nn_thresh, l1_depth_thresh : float, optional
        L1 detection thresholds. Defaults to module constants.
    verbose : bool
        Print per-cluster summary table.

    Returns
    -------
    domains : np.ndarray of str
        Per-cell domain label: 'Cortical', 'Vascular', or 'WM'.
    is_l1 : np.ndarray of bool
        Per-cell flag for L1 border cells.
    cluster_info : dict
        Per-cluster statistics (domain, is_l1, composition fractions, etc.).
    """
    # Defaults
    if vascular_thresh is None:
        vascular_thresh = VASCULAR_THRESH
    if wm_oligo_thresh is None:
        wm_oligo_thresh = WM_OLIGO_THRESH
    if wm_depth_thresh is None:
        wm_depth_thresh = WM_DEPTH_THRESH
    if l1_nn_thresh is None:
        l1_nn_thresh = L1_NN_THRESH
    if l1_depth_thresh is None:
        l1_depth_thresh = L1_DEPTH_THRESH

    # Resolve subclass column
    if subclass_col is None:
        if 'corr_subclass' in adata.obs.columns:
            subclass_col = 'corr_subclass'
        else:
            subclass_col = 'subclass_label'

    subclass = adata.obs[subclass_col].values.astype(str)
    pred_depth = pd.to_numeric(
        adata.obs[depth_col], errors="coerce"
    ).values

    unique_cl = np.unique(banksy_labels)
    cluster_info = {}

    for cl in unique_cl:
        mask = banksy_labels == cl
        n_cl = mask.sum()
        cl_sub = subclass[mask]
        cl_depth = pred_depth[mask]

        # Cell type composition
        sub_counts = {}
        for s in cl_sub:
            sub_counts[s] = sub_counts.get(s, 0) + 1

        vasc_frac = sum(sub_counts.get(v, 0) for v in VASCULAR_TYPES) / n_cl
        nn_frac = sum(sub_counts.get(v, 0) for v in NON_NEURONAL_TYPES) / n_cl
        oligo_frac = sub_counts.get("Oligodendrocyte", 0) / n_cl
        neuronal_frac = 1.0 - nn_frac

        valid_depth = cl_depth[~np.isnan(cl_depth)]
        mean_depth = (float(np.mean(valid_depth))
                      if len(valid_depth) > 0 else np.nan)

        # Top 3 types for diagnostics
        top3 = sorted(sub_counts.items(), key=lambda x: -x[1])[:3]
        top3_str = ", ".join(f"{s}({100*c/n_cl:.0f}%)" for s, c in top3)

        # ── Classification ────────────────────────────────────
        is_l1_cluster = False

        if vasc_frac > vascular_thresh:
            domain = "Vascular"
        elif (oligo_frac > wm_oligo_thresh
              and mean_depth > wm_depth_thresh):
            domain = "WM"
        elif (nn_frac > l1_nn_thresh
              and (np.isnan(mean_depth) or mean_depth < l1_depth_thresh)):
            # Shallow non-neuronal cluster = L1 cortex (not meningeal)
            domain = "Cortical"
            is_l1_cluster = True
        elif (neuronal_frac > 0.20
              and not np.isnan(mean_depth)
              and 0.0 <= mean_depth <= 0.90):
            domain = "Cortical"
        elif mean_depth > wm_depth_thresh:
            domain = "WM"
        else:
            domain = "Cortical"

        cluster_info[cl] = {
            "domain": domain,
            "is_l1": is_l1_cluster,
            "n_cells": n_cl,
            "vasc_frac": vasc_frac,
            "nn_frac": nn_frac,
            "oligo_frac": oligo_frac,
            "mean_depth": mean_depth,
            "top3": top3_str,
        }

    # Print summary table
    if verbose:
        print(f"  {'Cl':>4} | {'N':>6} | {'Domain':<9} | {'L1?':>3} "
              f"| {'Vasc%':>6} | {'NN%':>6} | {'Oligo%':>6} "
              f"| {'Depth':>6} | Top types")
        print(f"  {'-'*100}")
        for cl in sorted(cluster_info.keys()):
            info = cluster_info[cl]
            l1_str = "YES" if info["is_l1"] else ""
            print(f"  {cl:>4} | {info['n_cells']:>6,} "
                  f"| {info['domain']:<9} | {l1_str:>3} "
                  f"| {info['vasc_frac']*100:>5.1f}% "
                  f"| {info['nn_frac']*100:>5.1f}% "
                  f"| {info['oligo_frac']*100:>5.1f}% "
                  f"| {info['mean_depth']:>6.3f} "
                  f"| {info['top3']}")

    # Map to per-cell arrays
    domains = np.array([cluster_info[cl]["domain"] for cl in banksy_labels])
    is_l1 = np.array([cluster_info[cl]["is_l1"] for cl in banksy_labels])

    return domains, is_l1, cluster_info
