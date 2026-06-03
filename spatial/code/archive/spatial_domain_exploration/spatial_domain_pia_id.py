#!/usr/bin/env python3
"""
Refine spatial domain clustering to identify pia/meningeal zone.

Focus on separating:
  - Pia/meningeal: superficial, non-neuronal enriched, contiguous surface
  - Vascular: scattered throughout, nearly pure Endothelial+VLMC
  - Cortical layers: normal laminar structure

Usage:
    python3 -u spatial_domain_pia_id.py [sample_id]
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


def classify_clusters(clusters, subclass, pred_depth, coords):
    """
    Classify each cluster as 'pia', 'vascular', or 'cortical'.

    Criteria:
      - vascular: >80% Endothelial+VLMC
      - pia: >60% non-neuronal, mean_depth < 0.15, NOT >80% vascular
      - cortical: everything else
    """
    unique_cls = sorted(set(clusters), key=lambda c: int(c))
    classifications = {}

    for cl in unique_cls:
        cl_mask = clusters == cl
        cl_sub = subclass[cl_mask]
        n_cl = cl_mask.sum()

        vasc_frac = sum(1 for s in cl_sub if s in VASCULAR_TYPES) / n_cl
        nn_frac = sum(1 for s in cl_sub if s in NON_NEURONAL_TYPES) / n_cl

        cl_depths = pred_depth[cl_mask]
        valid = cl_depths[~np.isnan(cl_depths)]
        mean_depth = np.mean(valid) if len(valid) > 0 else 0.5

        if vasc_frac > 0.80:
            classifications[cl] = 'vascular'
        elif nn_frac > 0.60 and mean_depth < 0.15:
            classifications[cl] = 'pia'
        else:
            classifications[cl] = 'cortical'

    return classifications


def main():
    t0 = time.time()
    sample_id = sys.argv[1] if len(sys.argv) > 1 else "Br8667"

    model_bundle = load_model(os.path.join(OUTPUT_DIR, "depth_model_normalized.pkl"))
    subclass_names = model_bundle['subclass_names']

    h5ad_path = os.path.join(H5AD_DIR, f"{sample_id}_annotated.h5ad")
    adata = ad.read_h5ad(h5ad_path)
    qc_mask = adata.obs['qc_pass'].values.astype(bool)
    adata_pass = adata[qc_mask].copy()
    n_total = adata_pass.shape[0]
    print(f"Loaded {sample_id}: {n_total:,} QC-pass cells")

    # Build neighborhood features
    coords = adata_pass.obsm['spatial']
    subclass = adata_pass.obs['subclass_label'].values.astype(str)
    pred_depth = adata_pass.obs['predicted_norm_depth'].values
    x, y = coords[:, 0], coords[:, 1]

    print("Building neighborhood features...")
    features = build_neighborhood_features(coords, subclass, subclass_names,
                                            K=50, sections=None)
    n_sub = len(subclass_names)
    neigh_fracs = features[:, :n_sub]

    # Cluster at multiple resolutions to find stable pia domains
    resolutions = [0.3, 0.5, 0.8, 1.0]

    # Create temp adata for clustering
    adata_temp = ad.AnnData(X=neigh_fracs, obs=adata_pass.obs.copy())
    adata_temp.var_names = [f'neigh_{s}' for s in subclass_names]
    sc.pp.pca(adata_temp, n_comps=min(20, n_sub - 1))
    sc.pp.neighbors(adata_temp, n_neighbors=30, use_rep='X_pca')

    fig, axes = plt.subplots(3, len(resolutions), figsize=(10*len(resolutions), 30))

    for ri, res in enumerate(resolutions):
        print(f"\nLeiden resolution={res}...")
        sc.tl.leiden(adata_temp, resolution=res, flavor='igraph',
                     n_iterations=2, key_added=f'leiden_{res}')
        clusters = adata_temp.obs[f'leiden_{res}'].values.astype(str)
        n_cl = len(set(clusters))

        # Classify clusters
        classifications = classify_clusters(clusters, subclass, pred_depth, coords)

        # Map to domain labels
        domain_labels = np.array([classifications[cl] for cl in clusters])
        n_pia = (domain_labels == 'pia').sum()
        n_vasc = (domain_labels == 'vascular').sum()
        n_cort = (domain_labels == 'cortical').sum()

        print(f"  {n_cl} clusters -> pia:{n_pia:,} vasc:{n_vasc:,} cort:{n_cort:,}")

        # Top row: All clusters colored
        ax = axes[0, ri]
        ax.set_facecolor('black')
        cmap = plt.colormaps.get_cmap('tab20')
        unique_cls = sorted(set(clusters), key=lambda c: int(c))
        for i, cl in enumerate(unique_cls):
            mask = clusters == cl
            c = cmap(i % 20)
            ax.scatter(x[mask], y[mask], c=[c], s=0.1, alpha=0.4, rasterized=True)
        ax.set_title(f'res={res}: {n_cl} clusters', fontsize=18)
        ax.set_aspect('equal')
        ax.invert_yaxis()

        # Middle row: Domain classification
        ax = axes[1, ri]
        ax.set_facecolor('black')
        domain_colors = {'cortical': '#666666', 'pia': '#ff3333', 'vascular': '#33ff33'}
        for domain in ['cortical', 'vascular', 'pia']:  # pia on top
            mask = domain_labels == domain
            if mask.sum() > 0:
                ax.scatter(x[mask], y[mask], c=domain_colors[domain],
                          s=0.1 if domain == 'cortical' else 2,
                          alpha=0.3 if domain == 'cortical' else 0.8,
                          rasterized=True)
        patches = [
            mpatches.Patch(color='#ff3333', label=f'Pia ({n_pia:,}, {100*n_pia/n_total:.1f}%)'),
            mpatches.Patch(color='#33ff33', label=f'Vascular ({n_vasc:,}, {100*n_vasc/n_total:.1f}%)'),
            mpatches.Patch(color='#666666', label=f'Cortical ({n_cort:,})'),
        ]
        ax.legend(handles=patches, fontsize=14, loc='upper right')
        ax.set_title(f'Domain Classification (res={res})', fontsize=18)
        ax.set_aspect('equal')
        ax.invert_yaxis()

        # Bottom row: Layers with pia/vascular marked
        ax = axes[2, ri]
        ax.set_facecolor('black')
        layers = assign_discrete_layers(pred_depth)
        layers[domain_labels == 'pia'] = 'Extra-cortical'
        all_layers = list(LAYER_BINS.keys()) + ['Extra-cortical']
        for lname in all_layers:
            mask = layers == lname
            if mask.sum() > 0:
                c = LAYER_COLORS.get(lname, (0.5, 0.5, 0.5))
                ax.scatter(x[mask], y[mask], c=[c], s=0.1, alpha=0.4,
                           rasterized=True)
        patches = [mpatches.Patch(color=LAYER_COLORS.get(l, (0.5, 0.5, 0.5)),
                                  label=f'{l} ({(layers==l).sum():,})')
                   for l in all_layers if (layers == l).sum() > 0]
        ax.legend(handles=patches, fontsize=10, loc='upper right')
        ax.set_title(f'Layers + Pia (res={res})', fontsize=18)
        ax.set_aspect('equal')
        ax.invert_yaxis()

    plt.suptitle(f'{sample_id}: Pia/Vascular Domain Identification', fontsize=26)
    plt.tight_layout()

    fig_path = os.path.join(OUTPUT_DIR, f"spatial_pia_id_{sample_id}.png")
    plt.savefig(fig_path, dpi=150, facecolor='white', bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {fig_path}")
    print(f"Total: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
