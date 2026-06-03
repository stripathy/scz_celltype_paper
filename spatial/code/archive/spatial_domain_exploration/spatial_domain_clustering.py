#!/usr/bin/env python3
"""
Spatial domain clustering using neighborhood cell type composition.

BANKSY-inspired approach: cluster cells based on their local neighborhood
composition (K=50 NN subclass fractions). Cells in the pia/meninges should
form a distinct cluster because their neighborhoods are dominated by
vascular/non-neuronal types and lack the characteristic cortical layering.

Usage:
    python3 -u spatial_domain_clustering.py [sample_id] [resolution]
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


def cluster_by_neighborhood(adata, subclass_names, K=50, resolution=0.5,
                             n_neighbors_graph=30):
    """
    Cluster cells by spatial neighborhood composition.

    1. Build K-NN subclass composition features
    2. PCA on composition features
    3. Build KNN graph in PCA space
    4. Leiden clustering

    Returns cluster labels and the composition features.
    """
    coords = adata.obsm['spatial']
    subclass = adata.obs['subclass_label'].values.astype(str)

    print(f"  Building K={K} neighborhood features...")
    t0 = time.time()
    features = build_neighborhood_features(coords, subclass, subclass_names,
                                            K=K, sections=None)
    # Use only neighbor fractions (first n_sub columns)
    n_sub = len(subclass_names)
    neigh_fracs = features[:, :n_sub]
    print(f"  Features built: {time.time()-t0:.0f}s")

    # Create a temporary AnnData for scanpy operations
    adata_temp = ad.AnnData(X=neigh_fracs,
                            obs=adata.obs.copy())
    adata_temp.var_names = [f'neigh_{s}' for s in subclass_names]
    adata_temp.obsm['spatial'] = coords

    # PCA
    print(f"  Running PCA...")
    sc.pp.pca(adata_temp, n_comps=min(20, n_sub - 1))

    # KNN graph + Leiden
    print(f"  Building KNN graph (k={n_neighbors_graph})...")
    sc.pp.neighbors(adata_temp, n_neighbors=n_neighbors_graph, use_rep='X_pca')

    print(f"  Leiden clustering (resolution={resolution})...")
    sc.tl.leiden(adata_temp, resolution=resolution, flavor='igraph',
                 n_iterations=2)

    clusters = adata_temp.obs['leiden'].values.astype(str)
    n_clusters = len(set(clusters))
    print(f"  Found {n_clusters} clusters")

    return clusters, neigh_fracs, adata_temp


def main():
    t0 = time.time()
    sample_id = sys.argv[1] if len(sys.argv) > 1 else "Br8667"
    resolution = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5

    # Load model for subclass names and depth predictions
    model_bundle = load_model(os.path.join(OUTPUT_DIR, "depth_model_normalized.pkl"))
    subclass_names = model_bundle['subclass_names']

    # Load sample
    h5ad_path = os.path.join(H5AD_DIR, f"{sample_id}_annotated.h5ad")
    adata = ad.read_h5ad(h5ad_path)
    qc_mask = adata.obs['qc_pass'].values.astype(bool)
    adata_pass = adata[qc_mask].copy()
    n_total = adata_pass.shape[0]
    print(f"Loaded {sample_id}: {n_total:,} QC-pass cells")

    # Cluster
    clusters, neigh_fracs, adata_temp = cluster_by_neighborhood(
        adata_pass, subclass_names, K=50, resolution=resolution
    )

    coords = adata_pass.obsm['spatial']
    x, y = coords[:, 0], coords[:, 1]
    subclass = adata_pass.obs['subclass_label'].values
    pred_depth = adata_pass.obs['predicted_norm_depth'].values

    # ── Characterize each cluster ──────────────────────────────────
    unique_clusters = sorted(set(clusters), key=lambda c: int(c))
    print(f"\n{'='*80}")
    print(f"Cluster characterization for {sample_id} (resolution={resolution}):")
    print(f"{'='*80}")

    cluster_stats = []
    for cl in unique_clusters:
        cl_mask = clusters == cl
        n_cl = cl_mask.sum()

        # Mean depth
        cl_depths = pred_depth[cl_mask]
        valid_depths = cl_depths[~np.isnan(cl_depths)]
        mean_depth = np.mean(valid_depths) if len(valid_depths) > 0 else np.nan

        # Spatial extent
        cl_x, cl_y = x[cl_mask], y[cl_mask]
        y_range = cl_y.max() - cl_y.min()
        x_range = cl_x.max() - cl_x.min()

        # Top subclasses
        cl_subclass = subclass[cl_mask]
        sc_counts = {}
        for s in cl_subclass:
            sc_counts[s] = sc_counts.get(s, 0) + 1
        top3 = sorted(sc_counts.items(), key=lambda x: -x[1])[:3]
        top3_str = ', '.join(f'{s}({100*c/n_cl:.0f}%)' for s, c in top3)

        # Vascular fraction
        vasc_types = {'Endothelial', 'VLMC'}
        n_vasc = sum(1 for s in cl_subclass if s in vasc_types)
        vasc_pct = 100 * n_vasc / n_cl

        # Non-neuronal fraction
        nn_types = {'Endothelial', 'VLMC', 'Astrocyte', 'Microglia-PVM',
                    'Oligodendrocyte', 'OPC'}
        n_nn = sum(1 for s in cl_subclass if s in nn_types)
        nn_pct = 100 * n_nn / n_cl

        cluster_stats.append({
            'cluster': cl, 'n_cells': n_cl, 'mean_depth': mean_depth,
            'vasc_pct': vasc_pct, 'nn_pct': nn_pct, 'top3': top3_str
        })

        flag = " <<< PIA?" if vasc_pct > 30 and mean_depth < 0.2 else ""
        flag = " <<< PIA?" if vasc_pct > 40 else flag
        print(f"\n  Cluster {cl}: {n_cl:>6,} cells, "
              f"mean_depth={mean_depth:.3f}, vasc={vasc_pct:.0f}%, "
              f"non-neur={nn_pct:.0f}%{flag}")
        print(f"    Top: {top3_str}")

    # ── Figure ─────────────────────────────────────────────────────
    n_clusters = len(unique_clusters)
    # Use a qualitative colormap
    cmap = plt.cm.get_cmap('tab20', n_clusters)
    cluster_colors = {cl: cmap(i) for i, cl in enumerate(unique_clusters)}

    fig, axes = plt.subplots(2, 3, figsize=(36, 24))

    # Panel 1: Spatial clusters
    ax = axes[0, 0]
    ax.set_facecolor('black')
    for cl in unique_clusters:
        mask = clusters == cl
        c = cluster_colors[cl]
        ax.scatter(x[mask], y[mask], c=[c], s=0.1, alpha=0.4, rasterized=True)
    ax.set_title(f'{sample_id}: Spatial Domain Clusters (n={n_clusters})',
                 fontsize=20)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    # Panel 2: Current layer assignments
    ax = axes[0, 1]
    ax.set_facecolor('black')
    layers = assign_discrete_layers(pred_depth)
    all_layers = list(LAYER_BINS.keys())
    for lname in all_layers:
        mask = layers == lname
        if mask.sum() > 0:
            c = LAYER_COLORS.get(lname, (0.5, 0.5, 0.5))
            ax.scatter(x[mask], y[mask], c=[c], s=0.1, alpha=0.4,
                       rasterized=True)
    patches = [mpatches.Patch(color=LAYER_COLORS.get(l, (0.5, 0.5, 0.5)),
                              label=l) for l in all_layers]
    ax.legend(handles=patches, fontsize=12, loc='upper right')
    ax.set_title('Current Layer Assignment', fontsize=20)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    # Panel 3: Vascular fraction per cluster (bar)
    ax = axes[0, 2]
    cl_labels = [s['cluster'] for s in cluster_stats]
    vasc_pcts = [s['vasc_pct'] for s in cluster_stats]
    n_cells = [s['n_cells'] for s in cluster_stats]
    colors_bar = [cluster_colors[cl] for cl in cl_labels]
    bars = ax.bar(range(len(cl_labels)), vasc_pcts, color=colors_bar)
    ax.set_xticks(range(len(cl_labels)))
    ax.set_xticklabels(cl_labels, fontsize=12)
    ax.set_ylabel('% Vascular (Endo + VLMC)', fontsize=14)
    ax.set_xlabel('Cluster', fontsize=14)
    ax.set_title('Vascular Fraction by Cluster', fontsize=20)
    # Add cell count labels
    for i, (n, v) in enumerate(zip(n_cells, vasc_pcts)):
        ax.text(i, v + 1, f'{n:,}', ha='center', fontsize=9, rotation=45)

    # Panel 4: Mean depth per cluster
    ax = axes[1, 0]
    mean_depths = [s['mean_depth'] for s in cluster_stats]
    bars = ax.bar(range(len(cl_labels)), mean_depths, color=colors_bar)
    ax.set_xticks(range(len(cl_labels)))
    ax.set_xticklabels(cl_labels, fontsize=12)
    ax.set_ylabel('Mean Predicted Depth', fontsize=14)
    ax.set_xlabel('Cluster', fontsize=14)
    ax.set_title('Mean Depth by Cluster', fontsize=20)
    ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
    ax.axhline(1, color='gray', linestyle='--', alpha=0.5)

    # Panel 5: Highlight candidate pia clusters
    # Flag clusters with >30% vascular OR mean_depth < 0 as candidates
    ax = axes[1, 1]
    ax.set_facecolor('black')
    pia_candidates = set()
    for s in cluster_stats:
        if s['vasc_pct'] > 35 or (s['vasc_pct'] > 20 and s['mean_depth'] < 0.05):
            pia_candidates.add(s['cluster'])

    is_pia = np.isin(clusters, list(pia_candidates))
    ax.scatter(x[~is_pia], y[~is_pia], c='#333333', s=0.02, alpha=0.3,
               rasterized=True)
    for cl in pia_candidates:
        mask = clusters == cl
        c = cluster_colors[cl]
        ax.scatter(x[mask], y[mask], c=[c], s=1.5, alpha=0.8, rasterized=True,
                   label=f'Cluster {cl} ({mask.sum():,})')
    ax.legend(fontsize=14, loc='upper right', markerscale=5)
    n_pia = is_pia.sum()
    ax.set_title(f'Candidate Pia Clusters: {n_pia:,} ({100*n_pia/n_total:.1f}%)',
                 fontsize=20)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    # Panel 6: UMAP of composition space, colored by cluster
    ax = axes[1, 2]
    print("\n  Computing UMAP of composition space...")
    sc.tl.umap(adata_temp)
    umap = adata_temp.obsm['X_umap']
    for cl in unique_clusters:
        mask = clusters == cl
        c = cluster_colors[cl]
        ax.scatter(umap[mask, 0], umap[mask, 1], c=[c], s=0.3, alpha=0.3,
                   rasterized=True, label=cl)
    ax.set_title('UMAP of Neighborhood Composition', fontsize=20)
    ax.set_xlabel('UMAP1', fontsize=14)
    ax.set_ylabel('UMAP2', fontsize=14)

    plt.suptitle(f'{sample_id}: Spatial Domain Clustering (res={resolution})',
                 fontsize=26)
    plt.tight_layout()

    fig_path = os.path.join(OUTPUT_DIR, f"spatial_domains_{sample_id}.png")
    plt.savefig(fig_path, dpi=150, facecolor='white', bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {fig_path}")
    print(f"Total: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
