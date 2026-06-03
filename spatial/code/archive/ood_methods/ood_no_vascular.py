#!/usr/bin/env python3
"""
Test OOD scoring with vascular types excluded from composition features.

The idea: Endothelial and VLMC have unusual neighborhoods everywhere,
not just at the pia. By excluding them from the composition vector,
we focus on whether the *non-vascular* cellular neighborhood is cortex-like.

Usage:
    python3 -u ood_no_vascular.py [sample_id]
"""

import os
import sys
import time
import numpy as np
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from depth_model import (
    load_model, build_neighborhood_features, assign_discrete_layers,
    LAYER_BINS, LAYER_COLORS
)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5AD_DIR = os.path.join(BASE_DIR, "output", "h5ad")
OUTPUT_DIR = os.path.join(BASE_DIR, "output")
MERFISH_PATH = os.path.join(BASE_DIR, "sea-ad_reference",
                            "SEAAD_MTG_MERFISH.2024-12-11.h5ad")

# Types to exclude from OOD composition features
VASCULAR_TYPES = {'Endothelial', 'VLMC'}


def build_ood_features_no_vascular(coords, subclass_labels, subclass_names, K=50):
    """
    Build neighborhood composition features EXCLUDING vascular types.
    Returns only the neighbor-fraction columns (no own-type one-hot).
    """
    # Get indices of non-vascular subclasses
    non_vasc_names = [s for s in subclass_names if s not in VASCULAR_TYPES]
    non_vasc_idx = {s: i for i, s in enumerate(non_vasc_names)}
    n_nv = len(non_vasc_names)

    n_cells = len(subclass_labels)

    # Build KNN
    nn = NearestNeighbors(n_neighbors=K+1, algorithm='ball_tree')
    nn.fit(coords)
    _, nn_idx = nn.kneighbors(coords)
    nn_idx = nn_idx[:, 1:]  # exclude self

    # Map subclass labels to non-vascular indices (-1 for vascular)
    sub_idx = np.array([non_vasc_idx.get(s, -1) for s in subclass_labels])

    # Compute fractions using only non-vascular neighbors
    neighbor_types = sub_idx[nn_idx]  # (n_cells, K)
    valid = neighbor_types >= 0  # True for non-vascular neighbors
    neighbor_types_safe = np.where(valid, neighbor_types, 0)

    fractions = np.zeros((n_cells, n_nv), dtype=np.float32)
    for k in range(K):
        col = neighbor_types_safe[:, k]
        mask = valid[:, k]
        np.add.at(fractions, (np.arange(n_cells)[mask], col[mask]), 1)

    n_valid = valid.sum(axis=1, keepdims=True).astype(np.float32)
    n_valid = np.maximum(n_valid, 1)
    fractions /= n_valid

    return fractions, non_vasc_names


def main():
    t0 = time.time()
    sample_id = sys.argv[1] if len(sys.argv) > 1 else "Br8667"

    # Load main model for depth predictions and standard OOD
    model_bundle = load_model(os.path.join(OUTPUT_DIR, "depth_model_normalized.pkl"))
    thresh_99_orig = model_bundle['ood_threshold_99']

    # Load MERFISH for building vascular-excluded reference
    print("Loading MERFISH reference...")
    merfish = ad.read_h5ad(MERFISH_PATH)
    depth_col = 'Normalized depth from pia'
    has_depth = ~merfish.obs[depth_col].isna()
    merfish_depth = merfish[has_depth].copy()
    print(f"  {merfish_depth.shape[0]:,} depth-annotated cells")

    # Build vascular-excluded features for MERFISH training data
    print("Building vascular-excluded features for MERFISH...")
    t1 = time.time()
    merfish_coords = merfish_depth.obsm['X_spatial_raw']
    merfish_subclass = merfish_depth.obs['Subclass'].values.astype(str)
    merfish_sections = merfish_depth.obs['Section'].values.astype(str)
    all_subclass_names = sorted(set(merfish_subclass))

    # For MERFISH, build per-section (but skip vascular in composition)
    non_vasc_names = [s for s in all_subclass_names if s not in VASCULAR_TYPES]
    non_vasc_idx = {s: i for i, s in enumerate(non_vasc_names)}
    n_nv = len(non_vasc_names)
    sub_idx_merfish = np.array([non_vasc_idx.get(s, -1) for s in merfish_subclass])

    n_merfish = len(merfish_subclass)
    merfish_feats = np.zeros((n_merfish, n_nv), dtype=np.float32)
    K = 50

    unique_sections = np.unique(merfish_sections)
    for sec_i, sec in enumerate(unique_sections):
        if sec_i % 20 == 0:
            print(f"  Section {sec_i+1}/{len(unique_sections)}")
        sec_mask = merfish_sections == sec
        sec_indices = np.where(sec_mask)[0]
        sec_coords = merfish_coords[sec_indices]

        k_actual = min(K + 1, len(sec_indices))
        nn = NearestNeighbors(n_neighbors=k_actual, algorithm='ball_tree')
        nn.fit(sec_coords)
        _, nn_idx_local = nn.kneighbors(sec_coords)
        nn_neighbors_global = sec_indices[nn_idx_local[:, 1:K+1]]

        # Compute non-vascular fractions
        neighbor_types = sub_idx_merfish[nn_neighbors_global]
        valid = neighbor_types >= 0
        neighbor_types_safe = np.where(valid, neighbor_types, 0)

        sec_fracs = np.zeros((len(sec_indices), n_nv), dtype=np.float32)
        for k in range(min(K, k_actual - 1)):
            col = neighbor_types_safe[:, k]
            mask = valid[:, k]
            np.add.at(sec_fracs, (np.arange(len(sec_indices))[mask], col[mask]), 1)

        n_valid = valid.sum(axis=1, keepdims=True).astype(np.float32)
        n_valid = np.maximum(n_valid, 1)
        sec_fracs /= n_valid
        merfish_feats[sec_indices] = sec_fracs

    print(f"  MERFISH features built: {time.time()-t1:.0f}s")

    # Split train/test same as depth model
    donors = merfish_depth.obs.get('Donor ID',
             merfish_depth.obs.get('Specimen Barcode')).values.astype(str)
    train_donors = model_bundle['train_donors']
    test_donors = model_bundle['test_donors']
    train_mask = np.isin(donors, train_donors)
    test_mask = np.isin(donors, test_donors)

    X_train_nv = merfish_feats[train_mask]
    X_test_nv = merfish_feats[test_mask]
    print(f"  Train: {X_train_nv.shape[0]:,}, Test: {X_test_nv.shape[0]:,}")

    # Fit 1-NN on non-vascular features
    print("Fitting 1-NN on non-vascular features...")
    t2 = time.time()
    ood_nn_nv = NearestNeighbors(n_neighbors=1, algorithm='ball_tree', metric='euclidean')
    ood_nn_nv.fit(X_train_nv)

    # Calibrate on test set
    test_dists, _ = ood_nn_nv.kneighbors(X_test_nv)
    test_dists = test_dists.ravel()
    thresh_99_nv = float(np.percentile(test_dists, 99))
    thresh_95_nv = float(np.percentile(test_dists, 95))
    print(f"  1-NN fit + calibration: {time.time()-t2:.0f}s")
    print(f"  Non-vascular thresholds: 95th={thresh_95_nv:.4f}, 99th={thresh_99_nv:.4f}")
    print(f"  (Original thresholds:   95th={model_bundle['ood_threshold_95']:.4f}, "
          f"99th={thresh_99_orig:.4f})")

    # ── Load Xenium sample ─────────────────────────────────────────
    print(f"\nLoading {sample_id}...")
    h5ad_path = os.path.join(H5AD_DIR, f"{sample_id}_annotated.h5ad")
    adata = ad.read_h5ad(h5ad_path)
    qc_mask = adata.obs['qc_pass'].values.astype(bool)
    adata_pass = adata[qc_mask].copy()
    n_total = adata_pass.shape[0]

    # Standard depth + OOD
    from depth_model import predict_depth
    pred_depth, ood_scores_orig = predict_depth(adata_pass, model_bundle,
                                                 subclass_col='subclass_label',
                                                 compute_ood=True)

    # Non-vascular OOD features
    print("Building non-vascular features for Xenium...")
    coords = adata_pass.obsm['spatial']
    subclass = adata_pass.obs['subclass_label'].values.astype(str)
    xenium_feats_nv, _ = build_ood_features_no_vascular(
        coords, subclass, all_subclass_names, K=K
    )

    # Query the non-vascular 1-NN
    ood_dists_nv, _ = ood_nn_nv.kneighbors(xenium_feats_nv)
    ood_scores_nv = ood_dists_nv.ravel()

    # ── Compare approaches ─────────────────────────────────────────
    x, y = coords[:, 0], coords[:, 1]

    ood_orig = ood_scores_orig > thresh_99_orig
    ood_nv = ood_scores_nv > thresh_99_nv

    is_vascular = np.isin(subclass, list(VASCULAR_TYPES))

    print(f"\n{'='*60}")
    print(f"Comparison for {sample_id} ({n_total:,} cells):")
    print(f"  Original OOD (99th): {ood_orig.sum():,} ({100*ood_orig.sum()/n_total:.1f}%)")
    print(f"  No-vasc OOD (99th):  {ood_nv.sum():,} ({100*ood_nv.sum()/n_total:.1f}%)")
    print(f"  Overlap:             {(ood_orig & ood_nv).sum():,}")
    print(f"  Original-only:       {(ood_orig & ~ood_nv).sum():,}")
    print(f"  No-vasc-only:        {(~ood_orig & ood_nv).sum():,}")

    print(f"\n  Vascular cells (Endo+VLMC): {is_vascular.sum():,}")
    print(f"  Original OOD & vascular: {(ood_orig & is_vascular).sum():,}")
    print(f"  No-vasc OOD & vascular:  {(ood_nv & is_vascular).sum():,}")

    # Layer impact
    print(f"\n  No-vasc OOD by original layer:")
    for lname, (lo, hi) in LAYER_BINS.items():
        lmask = (pred_depth >= lo) & (pred_depth < hi)
        n_layer = lmask.sum()
        n_ood = (lmask & ood_nv).sum()
        print(f"    {lname:10s}: {n_ood:>5,}/{n_layer:>6,} "
              f"({100*n_ood/n_layer:.1f}%)" if n_layer > 0 else "")

    # Subclass breakdown
    print(f"\n  No-vasc OOD rate by subclass:")
    for sc in sorted(set(subclass)):
        sc_m = subclass == sc
        n_sc = sc_m.sum()
        n_ood = (sc_m & ood_nv).sum()
        pct = 100 * n_ood / n_sc if n_sc > 0 else 0
        flag = " <<<" if pct > 5 else ""
        print(f"    {sc:25s}: {n_ood:>5,}/{n_sc:>6,} ({pct:5.1f}%){flag}")

    # ── Figure ─────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 3, figsize=(36, 24))

    # Row 1: Spatial comparison
    # Panel 1: Original OOD
    ax = axes[0, 0]
    ax.set_facecolor('black')
    ax.scatter(x[~ood_orig], y[~ood_orig], c='#222222', s=0.02, alpha=0.3,
               rasterized=True)
    ax.scatter(x[ood_orig], y[ood_orig], c='#ff3333', s=1, alpha=0.7,
               rasterized=True)
    ax.set_title(f'Original OOD: {ood_orig.sum():,} ({100*ood_orig.sum()/n_total:.1f}%)',
                 fontsize=20)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    # Panel 2: No-vascular OOD
    ax = axes[0, 1]
    ax.set_facecolor('black')
    ax.scatter(x[~ood_nv], y[~ood_nv], c='#222222', s=0.02, alpha=0.3,
               rasterized=True)
    ax.scatter(x[ood_nv], y[ood_nv], c='#ffcc00', s=1, alpha=0.7,
               rasterized=True)
    ax.set_title(f'No-vascular OOD: {ood_nv.sum():,} ({100*ood_nv.sum()/n_total:.1f}%)',
                 fontsize=20)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    # Panel 3: Layers with no-vasc Extra-cortical
    ax = axes[0, 2]
    ax.set_facecolor('black')
    layers = assign_discrete_layers(pred_depth)
    layers[ood_nv] = 'Extra-cortical'
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
    ax.legend(handles=patches, fontsize=11, loc='upper right')
    ax.set_title('Layers + Extra-cortical (no-vasc)', fontsize=20)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    # Row 2: Diagnostics
    # Panel 4: OOD score distributions
    ax = axes[1, 0]
    ax.hist(ood_scores_nv, bins=200, color='#ffcc00', alpha=0.7,
            label='No-vascular', edgecolor='none')
    ax.axvline(thresh_99_nv, color='red', linestyle='--', linewidth=2,
               label=f'99th = {thresh_99_nv:.4f}')
    ax.set_xlabel('OOD Distance', fontsize=14)
    ax.set_ylabel('# Cells', fontsize=14)
    ax.set_title('No-vascular OOD Score Distribution', fontsize=18)
    ax.legend(fontsize=14)
    ax.set_xlim(0, min(np.percentile(ood_scores_nv, 99.9) * 2, ood_scores_nv.max()))

    # Panel 5: Venn-like comparison
    ax = axes[1, 1]
    ax.set_facecolor('black')
    both = ood_orig & ood_nv
    orig_only = ood_orig & ~ood_nv
    nv_only = ~ood_orig & ood_nv
    ax.scatter(x[~(ood_orig | ood_nv)], y[~(ood_orig | ood_nv)], c='#222222',
               s=0.02, alpha=0.3, rasterized=True)
    ax.scatter(x[orig_only], y[orig_only], c='#ff3333', s=1, alpha=0.7,
               rasterized=True)
    ax.scatter(x[nv_only], y[nv_only], c='#33ff33', s=1, alpha=0.7,
               rasterized=True)
    ax.scatter(x[both], y[both], c='#ffcc00', s=1, alpha=0.7,
               rasterized=True)
    patches = [mpatches.Patch(color='#ff3333',
                              label=f'Original-only ({orig_only.sum():,})'),
               mpatches.Patch(color='#33ff33',
                              label=f'No-vasc-only ({nv_only.sum():,})'),
               mpatches.Patch(color='#ffcc00',
                              label=f'Both ({both.sum():,})')]
    ax.legend(handles=patches, fontsize=14, loc='upper right')
    ax.set_title('OOD comparison', fontsize=18)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    # Panel 6: No-vasc OOD vs depth
    ax = axes[1, 2]
    ax.set_facecolor('#111111')
    n_plot = min(30000, n_total)
    idx = np.random.choice(n_total, n_plot, replace=False)
    vasc_plot = is_vascular[idx]
    ax.scatter(pred_depth[idx][~vasc_plot], ood_scores_nv[idx][~vasc_plot],
               s=0.3, alpha=0.3, c='#4ecdc4', rasterized=True, label='Non-vasc')
    ax.scatter(pred_depth[idx][vasc_plot], ood_scores_nv[idx][vasc_plot],
               s=0.5, alpha=0.5, c='#ff6633', rasterized=True, label='Vascular')
    ax.axhline(thresh_99_nv, color='red', linestyle='--', linewidth=2)
    ax.set_xlabel('Predicted Depth', fontsize=14)
    ax.set_ylabel('No-vasc OOD Distance', fontsize=14)
    ax.set_title('No-vasc OOD vs Depth', fontsize=18)
    ax.legend(fontsize=14)

    plt.suptitle(f'{sample_id}: Vascular-excluded OOD Comparison', fontsize=26)
    plt.tight_layout()

    fig_path = os.path.join(OUTPUT_DIR, f"ood_novasc_{sample_id}.png")
    plt.savefig(fig_path, dpi=150, facecolor='white', bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {fig_path}")
    print(f"Total: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
