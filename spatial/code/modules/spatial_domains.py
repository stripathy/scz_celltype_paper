"""
Spatial domain classification — cell type constants and legacy classifier.

SUPERSEDED: The classify_spatial_domains() function in this module has been
replaced by banksy_domains.py, which uses BANKSY clustering for more
accurate domain classification:
  - Correctly identifies L1 border cells as Cortical (not "Extra-cortical")
  - Detects white matter via oligo + depth thresholds
  - Lowered vascular threshold (0.50 vs 0.80) enabled by spatially coherent clusters

This module is retained for:
  - VASCULAR_TYPES, NON_NEURONAL_TYPES constants (used throughout codebase)
  - classify_spatial_domains_legacy() for reference/comparison

For new code, use banksy_domains.py instead.
"""

import numpy as np
import anndata as ad
import scanpy as sc

VASCULAR_TYPES = {'Endothelial', 'VLMC'}
NON_NEURONAL_TYPES = {'Endothelial', 'VLMC', 'Astrocyte', 'Microglia-PVM',
                      'Oligodendrocyte', 'OPC'}


def classify_spatial_domains_legacy(adata, subclass_names, K=50, resolution=0.8,
                              subclass_col='subclass_label',
                              depth_col='predicted_norm_depth',
                              pia_nn_thresh=0.60, pia_depth_thresh=0.15,
                              vascular_thresh=0.80):
    """
    Classify cells into spatial domains (Extra-cortical, Vascular, Cortical).

    Steps:
      1. Build K-NN subclass composition features from spatial coordinates
      2. PCA on composition features
      3. Build KNN graph in PCA space
      4. Leiden clustering
      5. Classify each cluster by cell type composition thresholds

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated Xenium sample with subclass labels, depth predictions,
        and spatial coordinates in .obsm['spatial'].
    subclass_names : list of str
        Ordered list of all subclass names (from model bundle or reference).
    K : int
        Number of spatial nearest neighbors for composition features.
    resolution : float
        Leiden clustering resolution.
    subclass_col : str
        Column name for subclass labels.
    depth_col : str
        Column name for predicted depth.
    pia_nn_thresh : float
        Non-neuronal fraction threshold for pia classification.
    pia_depth_thresh : float
        Mean depth threshold for pia classification.
    vascular_thresh : float
        Vascular fraction threshold for vascular classification.

    Returns
    -------
    np.ndarray of str
        Domain label per cell: 'Extra-cortical', 'Vascular', or 'Cortical'.
    dict
        Per-cluster statistics for diagnostics.
    """
    from depth_model import build_neighborhood_features

    n_cells = adata.shape[0]
    coords = adata.obsm['spatial']
    subclass = adata.obs[subclass_col].values.astype(str)
    pred_depth = adata.obs[depth_col].values if depth_col in adata.obs.columns else np.full(n_cells, np.nan)

    # Step 1: Build neighborhood composition features
    features = build_neighborhood_features(coords, subclass, subclass_names,
                                            K=K, sections=None)
    n_sub = len(subclass_names)
    neigh_fracs = features[:, :n_sub]

    # Step 2-4: PCA → KNN graph → Leiden
    adata_temp = ad.AnnData(X=neigh_fracs)
    adata_temp.var_names = [f'neigh_{s}' for s in subclass_names]
    sc.pp.pca(adata_temp, n_comps=min(20, n_sub - 1))
    sc.pp.neighbors(adata_temp, n_neighbors=30, use_rep='X_pca')
    sc.tl.leiden(adata_temp, resolution=resolution, flavor='igraph',
                 n_iterations=2)
    clusters = adata_temp.obs['leiden'].values.astype(str)

    # Step 5: Classify each cluster
    unique_clusters = sorted(set(clusters), key=lambda c: int(c))
    cluster_stats = {}

    for cl in unique_clusters:
        cl_mask = clusters == cl
        n_cl = cl_mask.sum()
        cl_sub = subclass[cl_mask]
        cl_depth = pred_depth[cl_mask]

        # Count cell types
        sc_counts = {}
        for s in cl_sub:
            sc_counts[s] = sc_counts.get(s, 0) + 1

        vasc_frac = sum(sc_counts.get(v, 0) for v in VASCULAR_TYPES) / n_cl
        nn_frac = sum(sc_counts.get(v, 0) for v in NON_NEURONAL_TYPES) / n_cl

        valid_depth = cl_depth[~np.isnan(cl_depth)]
        mean_depth = float(np.mean(valid_depth)) if len(valid_depth) > 0 else np.nan

        # Top 3 cell types
        top3 = sorted(sc_counts.items(), key=lambda x: -x[1])[:3]
        top3_str = ', '.join(f'{s}({100*c/n_cl:.0f}%)' for s, c in top3)

        # Classification
        if vasc_frac > vascular_thresh:
            domain = 'Vascular'
        elif nn_frac > pia_nn_thresh and (np.isnan(mean_depth) or mean_depth < pia_depth_thresh):
            domain = 'Extra-cortical'
        else:
            domain = 'Cortical'

        cluster_stats[cl] = {
            'n_cells': n_cl, 'domain': domain,
            'vasc_frac': vasc_frac, 'nn_frac': nn_frac,
            'mean_depth': mean_depth, 'top3': top3_str,
        }

    # Map clusters to domains
    domain_labels = np.array([cluster_stats[cl]['domain'] for cl in clusters])

    return domain_labels, cluster_stats
