#!/usr/bin/env python3
"""
Test whether spatial domain clustering can directly assign cortical layers
from Xenium data, rather than transferring depth from MERFISH.

Approach:
  1. Cluster cells by K=50 neighborhood composition (same as before)
  2. Order clusters by their mean spatial position along the pia→WM axis
  3. Characterize each cluster's cell type composition
  4. Compare to MERFISH-transferred layer assignments

Usage:
    python3 -u spatial_domain_layers.py [sample_id] [resolution]
"""

import os
import sys
import time
import numpy as np
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sklearn.decomposition import PCA

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from depth_model import (
    load_model, build_neighborhood_features, assign_discrete_layers,
    LAYER_BINS, LAYER_COLORS
)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5AD_DIR = os.path.join(BASE_DIR, "output", "h5ad")
OUTPUT_DIR = os.path.join(BASE_DIR, "output")

VASCULAR_TYPES = {'Endothelial', 'VLMC'}
NON_NEURONAL_TYPES = {'Endothelial', 'VLMC', 'Astrocyte', 'Microglia-PVM',
                      'Oligodendrocyte', 'OPC'}

# Expected layer markers (from SEA-AD taxonomy)
LAYER_MARKER_TYPES = {
    'L1': ['Lamp5', 'Vip', 'Lamp5 Lhx6'],
    'L2/3': ['L2/3 IT'],
    'L4': ['L4 IT'],
    'L5': ['L5 IT', 'L5 ET', 'L5/6 NP'],
    'L6': ['L6 IT', 'L6 CT', 'L6b', 'L6 IT Car3'],
    'WM': ['Oligodendrocyte'],
}


def estimate_pia_wm_axis(coords, subclass, pred_depth):
    """
    Estimate the pia→WM axis direction using predicted depth as a guide.
    Returns axis as a unit vector and the projection of each cell onto it.
    """
    # Use PCA on spatial coords weighted by depth to find the laminar axis
    # Simple approach: the direction that maximizes correlation with depth
    valid = ~np.isnan(pred_depth)
    c = coords[valid] - coords[valid].mean(axis=0)
    d = pred_depth[valid]

    # Find direction that best correlates with depth
    # This is just the regression coefficient direction
    # proj = c @ w, maximize corr(proj, d)
    # Solution: w = cov(c, d) / ||cov(c, d)||
    cov_cd = c.T @ d / len(d)
    axis = cov_cd / np.linalg.norm(cov_cd)

    # Project all cells
    projections = (coords - coords.mean(axis=0)) @ axis
    return axis, projections


def main():
    t0 = time.time()
    sample_id = sys.argv[1] if len(sys.argv) > 1 else "Br8667"
    resolution = float(sys.argv[2]) if len(sys.argv) > 2 else 0.8

    model_bundle = load_model(os.path.join(OUTPUT_DIR, "depth_model_normalized.pkl"))
    subclass_names = model_bundle['subclass_names']

    h5ad_path = os.path.join(H5AD_DIR, f"{sample_id}_annotated.h5ad")
    adata = ad.read_h5ad(h5ad_path)
    qc_mask = adata.obs['qc_pass'].values.astype(bool)
    adata_pass = adata[qc_mask].copy()
    n_total = adata_pass.shape[0]
    print(f"Loaded {sample_id}: {n_total:,} QC-pass cells")

    coords = adata_pass.obsm['spatial']
    subclass = adata_pass.obs['subclass_label'].values.astype(str)
    pred_depth = adata_pass.obs['predicted_norm_depth'].values
    x, y = coords[:, 0], coords[:, 1]

    # Cluster by neighborhood composition
    print("Building neighborhood features...")
    features = build_neighborhood_features(coords, subclass, subclass_names,
                                            K=50, sections=None)
    n_sub = len(subclass_names)
    neigh_fracs = features[:, :n_sub]

    adata_temp = ad.AnnData(X=neigh_fracs, obs=adata_pass.obs.copy())
    adata_temp.var_names = [f'neigh_{s}' for s in subclass_names]
    adata_temp.obsm['spatial'] = coords

    sc.pp.pca(adata_temp, n_comps=min(20, n_sub - 1))
    sc.pp.neighbors(adata_temp, n_neighbors=30, use_rep='X_pca')
    sc.tl.leiden(adata_temp, resolution=resolution, flavor='igraph',
                 n_iterations=2)
    clusters = adata_temp.obs['leiden'].values.astype(str)
    unique_cl = sorted(set(clusters), key=lambda c: int(c))
    n_cl = len(unique_cl)
    print(f"  Found {n_cl} clusters at resolution {resolution}")

    # Estimate pia→WM axis
    print("Estimating pia→WM axis...")
    axis, projections = estimate_pia_wm_axis(coords, subclass, pred_depth)
    print(f"  Axis direction: [{axis[0]:.3f}, {axis[1]:.3f}]")

    # ── Characterize each cluster ──────────────────────────────────
    cluster_info = []
    for cl in unique_cl:
        cl_mask = clusters == cl
        n_cl_cells = cl_mask.sum()
        cl_sub = subclass[cl_mask]
        cl_depth = pred_depth[cl_mask]

        # Composition
        sc_counts = {}
        for s in cl_sub:
            sc_counts[s] = sc_counts.get(s, 0) + 1

        vasc_frac = sum(sc_counts.get(v, 0) for v in VASCULAR_TYPES) / n_cl_cells
        nn_frac = sum(sc_counts.get(v, 0) for v in NON_NEURONAL_TYPES) / n_cl_cells

        # Layer marker enrichment
        layer_scores = {}
        for layer, markers in LAYER_MARKER_TYPES.items():
            marker_count = sum(sc_counts.get(m, 0) for m in markers)
            layer_scores[layer] = marker_count / n_cl_cells

        # Mean depth and projection
        valid_depth = cl_depth[~np.isnan(cl_depth)]
        mean_depth = np.mean(valid_depth) if len(valid_depth) > 0 else np.nan
        mean_proj = np.mean(projections[cl_mask])

        # Top 3 cell types
        top3 = sorted(sc_counts.items(), key=lambda x: -x[1])[:3]
        top3_str = ', '.join(f'{s}({100*c/n_cl_cells:.0f}%)' for s, c in top3)

        # Best matching layer
        best_layer = max(layer_scores, key=layer_scores.get)
        best_score = layer_scores[best_layer]

        cluster_info.append({
            'cluster': cl, 'n_cells': n_cl_cells,
            'mean_depth': mean_depth, 'mean_proj': mean_proj,
            'vasc_frac': vasc_frac, 'nn_frac': nn_frac,
            'layer_scores': layer_scores, 'best_layer': best_layer,
            'best_score': best_score, 'top3': top3_str,
            'sc_counts': sc_counts,
        })

    # Sort by mean projection (pia→WM)
    cluster_info.sort(key=lambda x: x['mean_proj'])

    print(f"\n{'='*100}")
    print(f"Clusters ordered by pia→WM axis (resolution={resolution}):")
    print(f"{'='*100}")
    print(f"{'Cl':>3} {'N':>6} {'Depth':>6} {'Vasc%':>5} {'NN%':>4} "
          f"{'L1':>5} {'L2/3':>5} {'L4':>5} {'L5':>5} {'L6':>5} {'WM':>5} "
          f"{'Best':>5} {'Top types'}")
    print('-' * 100)
    for info in cluster_info:
        ls = info['layer_scores']
        domain = 'PIA' if info['vasc_frac'] > 0.8 or (info['nn_frac'] > 0.6 and info['mean_depth'] < 0.15) else ''
        print(f"{info['cluster']:>3} {info['n_cells']:>6,} {info['mean_depth']:>6.3f} "
              f"{100*info['vasc_frac']:>5.0f} {100*info['nn_frac']:>4.0f} "
              f"{100*ls['L1']:>5.1f} {100*ls['L2/3']:>5.1f} {100*ls['L4']:>5.1f} "
              f"{100*ls['L5']:>5.1f} {100*ls['L6']:>5.1f} {100*ls['WM']:>5.1f} "
              f"{info['best_layer']:>5} {info['top3']}  {domain}")

    # ── Assign layers based on cluster identity ────────────────────
    # Strategy: for each cluster, assign the layer whose marker types
    # are most enriched, with special handling for pia/vascular
    def assign_cluster_layer(info):
        if info['vasc_frac'] > 0.80:
            return 'Vascular'
        if info['nn_frac'] > 0.60 and info['mean_depth'] < 0.15:
            return 'Extra-cortical'
        # Use marker enrichment
        return info['best_layer']

    cluster_to_layer = {}
    for info in cluster_info:
        cluster_to_layer[info['cluster']] = assign_cluster_layer(info)

    # Apply
    cluster_layers = np.array([cluster_to_layer[cl] for cl in clusters])

    print(f"\n\nCluster-based layer assignment:")
    for lname in ['Extra-cortical', 'L1', 'L2/3', 'L4', 'L5', 'L6', 'WM', 'Vascular']:
        n = (cluster_layers == lname).sum()
        if n > 0:
            print(f"  {lname:20s}: {n:>6,} ({100*n/n_total:5.1f}%)")

    # MERFISH-transferred layers for comparison
    merfish_layers = assign_discrete_layers(pred_depth)

    print(f"\nMERFISH-transferred layer assignment:")
    for lname in list(LAYER_BINS.keys()):
        n = (merfish_layers == lname).sum()
        print(f"  {lname:20s}: {n:>6,} ({100*n/n_total:5.1f}%)")

    # Concordance
    # Map Vascular→closest layer for fair comparison
    comparable_cluster_layers = cluster_layers.copy()
    # Don't count vascular or extra-cortical in concordance
    comparable_mask = ~np.isin(cluster_layers, ['Vascular', 'Extra-cortical'])
    agree = (comparable_cluster_layers[comparable_mask] ==
             merfish_layers[comparable_mask]).sum()
    print(f"\nConcordance (excl. Vasc/Extra): {agree:,}/{comparable_mask.sum():,} "
          f"({100*agree/comparable_mask.sum():.1f}%)")

    # ── Figure ─────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 3, figsize=(36, 24))

    # Cluster-based layers
    ax = axes[0, 0]
    ax.set_facecolor('black')
    layer_color_map = {**LAYER_COLORS, 'Vascular': (0.2, 0.9, 0.2)}
    all_labels = ['Extra-cortical', 'L1', 'L2/3', 'L4', 'L5', 'L6', 'WM', 'Vascular']
    for lname in all_labels:
        mask = cluster_layers == lname
        if mask.sum() > 0:
            c = layer_color_map.get(lname, (0.5, 0.5, 0.5))
            ax.scatter(x[mask], y[mask], c=[c], s=0.1, alpha=0.4, rasterized=True)
    patches = [mpatches.Patch(color=layer_color_map.get(l, (0.5, 0.5, 0.5)),
                              label=f'{l} ({(cluster_layers==l).sum():,})')
               for l in all_labels if (cluster_layers == l).sum() > 0]
    ax.legend(handles=patches, fontsize=11, loc='upper right')
    ax.set_title(f'Cluster-based Layers', fontsize=22)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    # MERFISH-transferred layers
    ax = axes[0, 1]
    ax.set_facecolor('black')
    for lname in LAYER_BINS.keys():
        mask = merfish_layers == lname
        if mask.sum() > 0:
            c = LAYER_COLORS.get(lname, (0.5, 0.5, 0.5))
            ax.scatter(x[mask], y[mask], c=[c], s=0.1, alpha=0.4, rasterized=True)
    patches = [mpatches.Patch(color=LAYER_COLORS.get(l, (0.5, 0.5, 0.5)),
                              label=f'{l} ({(merfish_layers==l).sum():,})')
               for l in LAYER_BINS.keys()]
    ax.legend(handles=patches, fontsize=11, loc='upper right')
    ax.set_title('MERFISH-transferred Layers', fontsize=22)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    # Agreement/disagreement
    ax = axes[0, 2]
    ax.set_facecolor('black')
    # Compare only non-vascular/extra-cortical cells
    agree_mask = (cluster_layers == merfish_layers) & comparable_mask
    disagree_mask = (cluster_layers != merfish_layers) & comparable_mask
    special_mask = ~comparable_mask

    ax.scatter(x[agree_mask], y[agree_mask], c='#33cc33', s=0.05, alpha=0.3,
               rasterized=True)
    ax.scatter(x[disagree_mask], y[disagree_mask], c='#ff3333', s=0.3, alpha=0.5,
               rasterized=True)
    ax.scatter(x[special_mask], y[special_mask], c='#ffcc00', s=0.3, alpha=0.5,
               rasterized=True)
    patches = [
        mpatches.Patch(color='#33cc33', label=f'Agree ({agree_mask.sum():,})'),
        mpatches.Patch(color='#ff3333', label=f'Disagree ({disagree_mask.sum():,})'),
        mpatches.Patch(color='#ffcc00', label=f'Pia/Vasc ({special_mask.sum():,})'),
    ]
    ax.legend(handles=patches, fontsize=14, loc='upper right')
    ax.set_title('Layer Agreement', fontsize=22)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    # Confusion matrix (heatmap)
    ax = axes[1, 0]
    layer_order_cl = ['Extra-cortical', 'L1', 'L2/3', 'L4', 'L5', 'L6', 'WM', 'Vascular']
    layer_order_mf = list(LAYER_BINS.keys())
    conf = np.zeros((len(layer_order_cl), len(layer_order_mf)))
    for i, cl_l in enumerate(layer_order_cl):
        for j, mf_l in enumerate(layer_order_mf):
            conf[i, j] = ((cluster_layers == cl_l) & (merfish_layers == mf_l)).sum()
    # Normalize by row (cluster layer)
    row_sums = conf.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    conf_norm = conf / row_sums

    im = ax.imshow(conf_norm, cmap='Blues', aspect='auto', vmin=0, vmax=1)
    ax.set_xticks(range(len(layer_order_mf)))
    ax.set_xticklabels(layer_order_mf, fontsize=12, rotation=45)
    ax.set_yticks(range(len(layer_order_cl)))
    ax.set_yticklabels(layer_order_cl, fontsize=12)
    ax.set_xlabel('MERFISH-transferred Layer', fontsize=14)
    ax.set_ylabel('Cluster-based Layer', fontsize=14)
    ax.set_title('Confusion Matrix (row-normalized)', fontsize=18)
    # Add text annotations
    for i in range(len(layer_order_cl)):
        for j in range(len(layer_order_mf)):
            val = conf_norm[i, j]
            count = int(conf[i, j])
            if count > 0:
                color = 'white' if val > 0.5 else 'black'
                ax.text(j, i, f'{val:.2f}\n({count:,})', ha='center',
                        va='center', fontsize=8, color=color)
    plt.colorbar(im, ax=ax, shrink=0.7)

    # Cluster composition heatmap
    ax = axes[1, 1]
    # Show subclass fractions per cluster (ordered by depth)
    key_types = ['L2/3 IT', 'L4 IT', 'L5 IT', 'L5 ET', 'L6 IT', 'L6 CT',
                 'L6b', 'Oligodendrocyte', 'Astrocyte', 'Endothelial',
                 'VLMC', 'Microglia-PVM', 'Pvalb', 'Sst', 'Vip', 'Lamp5']
    comp_matrix = np.zeros((len(cluster_info), len(key_types)))
    for i, info in enumerate(cluster_info):
        for j, kt in enumerate(key_types):
            comp_matrix[i, j] = info['sc_counts'].get(kt, 0) / info['n_cells']

    im2 = ax.imshow(comp_matrix, cmap='YlOrRd', aspect='auto', vmin=0, vmax=0.4)
    ax.set_xticks(range(len(key_types)))
    ax.set_xticklabels(key_types, fontsize=9, rotation=60, ha='right')
    cl_labels = [f"Cl{info['cluster']} ({cluster_to_layer[info['cluster']]})"
                 for info in cluster_info]
    ax.set_yticks(range(len(cluster_info)))
    ax.set_yticklabels(cl_labels, fontsize=9)
    ax.set_title('Subclass Composition per Cluster\n(ordered pia→WM)', fontsize=16)
    plt.colorbar(im2, ax=ax, shrink=0.7, label='Fraction')

    # Depth distribution per cluster-assigned layer
    ax = axes[1, 2]
    layer_order = ['Extra-cortical', 'L1', 'L2/3', 'L4', 'L5', 'L6', 'WM']
    data_for_box = []
    labels_for_box = []
    for lname in layer_order:
        mask = cluster_layers == lname
        if mask.sum() > 0:
            d = pred_depth[mask]
            d = d[~np.isnan(d)]
            data_for_box.append(d)
            labels_for_box.append(f'{lname}\n(n={mask.sum():,})')
    bp = ax.boxplot(data_for_box, vert=True, patch_artist=True,
                    showfliers=False, tick_labels=labels_for_box)
    colors_box = [LAYER_COLORS.get(l.split('\n')[0], (0.5, 0.5, 0.5))
                  for l in labels_for_box]
    for patch, color in zip(bp['boxes'], colors_box):
        patch.set_facecolor(color)
    ax.set_ylabel('MERFISH-predicted Depth', fontsize=14)
    ax.set_title('Depth Distribution per Cluster-Layer', fontsize=18)
    ax.tick_params(axis='x', labelsize=10)

    plt.suptitle(f'{sample_id}: Cluster-based vs MERFISH-transferred Layers (res={resolution})',
                 fontsize=24)
    plt.tight_layout()

    fig_path = os.path.join(OUTPUT_DIR, f"layer_comparison_{sample_id}.png")
    plt.savefig(fig_path, dpi=150, facecolor='white', bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {fig_path}")
    print(f"Total: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
