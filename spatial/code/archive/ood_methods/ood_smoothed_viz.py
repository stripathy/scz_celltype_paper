#!/usr/bin/env python3
"""
Smoothed OOD approach: compute average OOD score in each cell's spatial
neighborhood. Cells in regions where MANY neighbors are OOD (contiguous
pia/meninges) get flagged, but isolated vascular cells don't.

Usage:
    python3 -u ood_smoothed_viz.py [sample_id]
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
    load_model, predict_depth, assign_discrete_layers,
    LAYER_BINS, LAYER_COLORS
)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5AD_DIR = os.path.join(BASE_DIR, "output", "h5ad")
OUTPUT_DIR = os.path.join(BASE_DIR, "output")


def smooth_ood_scores(coords, ood_scores, ood_threshold, K=50):
    """
    Compute the fraction of K spatial neighbors that are OOD.

    Returns a value between 0 and 1 for each cell. High values
    indicate the cell is in a contiguous OOD region (e.g., pia).
    Low values indicate isolated OOD cells (e.g., scattered vascular).
    """
    nn = NearestNeighbors(n_neighbors=K+1, algorithm='ball_tree')
    nn.fit(coords)
    _, nn_idx = nn.kneighbors(coords)
    nn_idx = nn_idx[:, 1:]  # exclude self

    is_ood = (ood_scores > ood_threshold).astype(float)
    # For each cell, fraction of neighbors that are OOD
    neighbor_ood_frac = is_ood[nn_idx].mean(axis=1)
    return neighbor_ood_frac


def main():
    t0 = time.time()
    sample_id = sys.argv[1] if len(sys.argv) > 1 else "Br8667"

    model_path = os.path.join(OUTPUT_DIR, "depth_model_normalized.pkl")
    model_bundle = load_model(model_path)
    thresh_99 = model_bundle['ood_threshold_99']

    h5ad_path = os.path.join(H5AD_DIR, f"{sample_id}_annotated.h5ad")
    adata = ad.read_h5ad(h5ad_path)
    qc_mask = adata.obs['qc_pass'].values.astype(bool)
    adata_pass = adata[qc_mask].copy()
    n_total = adata_pass.shape[0]
    print(f"Loaded {sample_id}: {n_total:,} QC-pass cells")

    pred_depth, ood_scores = predict_depth(adata_pass, model_bundle,
                                           subclass_col='subclass_label',
                                           compute_ood=True)

    coords = adata_pass.obsm['spatial']
    x, y = coords[:, 0], coords[:, 1]
    subclass = adata_pass.obs['subclass_label'].values

    # Compute spatially smoothed OOD
    print("Computing spatially smoothed OOD (K=50)...")
    t1 = time.time()
    neighbor_ood_frac = smooth_ood_scores(coords, ood_scores, thresh_99, K=50)
    print(f"  Smoothing took {time.time()-t1:.0f}s")

    # Try different thresholds for the smoothed score
    frac_thresholds = [0.3, 0.4, 0.5, 0.6, 0.7]

    print("\nSmoothed Extra-cortical (OOD + spatial context):")
    for ft in frac_thresholds:
        extra = (neighbor_ood_frac >= ft).sum()
        print(f"  neighbor_ood_frac >= {ft:.1f}: {extra:>5,} ({100*extra/n_total:.1f}%)")

    # ── Figure ─────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 3, figsize=(36, 24))

    # Panel 1: Raw OOD binary
    ax = axes[0, 0]
    ax.set_facecolor('black')
    ood_mask = ood_scores > thresh_99
    ax.scatter(x[~ood_mask], y[~ood_mask], c='#222222', s=0.02, alpha=0.3,
               rasterized=True)
    ax.scatter(x[ood_mask], y[ood_mask], c='#ff3333', s=1, alpha=0.7,
               rasterized=True)
    ax.set_title(f'Raw OOD: {ood_mask.sum():,} ({100*ood_mask.sum()/n_total:.1f}%)',
                 fontsize=20)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    # Panel 2: Smoothed OOD score (continuous)
    ax = axes[0, 1]
    ax.set_facecolor('black')
    sc = ax.scatter(x, y, c=neighbor_ood_frac, cmap='hot', s=0.1, alpha=0.5,
                    vmin=0, vmax=0.8, rasterized=True)
    plt.colorbar(sc, ax=ax, label='Fraction of K=50 neighbors that are OOD',
                 shrink=0.7)
    ax.set_title('Smoothed OOD (neighbor fraction)', fontsize=20)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    # Panel 3: Smoothed OOD binary (threshold = 0.5)
    best_thresh = 0.5
    extra_smoothed = neighbor_ood_frac >= best_thresh
    ax = axes[0, 2]
    ax.set_facecolor('black')
    ax.scatter(x[~extra_smoothed], y[~extra_smoothed], c='#222222', s=0.02,
               alpha=0.3, rasterized=True)
    ax.scatter(x[extra_smoothed], y[extra_smoothed], c='#ffcc00', s=1.5,
               alpha=0.8, rasterized=True)
    ax.set_title(f'Smoothed OOD >= {best_thresh}: {extra_smoothed.sum():,} '
                 f'({100*extra_smoothed.sum()/n_total:.1f}%)', fontsize=20)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    # Panel 4: Layer assignment with smoothed extra-cortical
    ax = axes[1, 0]
    ax.set_facecolor('black')
    layers = assign_discrete_layers(pred_depth)
    layers[extra_smoothed] = 'Extra-cortical'
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
    ax.set_title('Layers + Extra-cortical (smoothed)', fontsize=20)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    # Panel 5: Subclass composition of extra-cortical
    ax = axes[1, 1]
    if extra_smoothed.sum() > 0:
        sc_counts = {}
        for sc in set(subclass[extra_smoothed]):
            sc_counts[sc] = (subclass[extra_smoothed] == sc).sum()
        sorted_sc = sorted(sc_counts.items(), key=lambda x: -x[1])
        names = [s[0] for s in sorted_sc[:12]]
        counts = [s[1] for s in sorted_sc[:12]]
        pcts = [100 * c / extra_smoothed.sum() for c in counts]
        bars = ax.barh(range(len(names)), pcts, color='#e94560')
        for i, (c, p) in enumerate(zip(counts, pcts)):
            ax.text(p + 0.5, i, f'{c:,}', va='center', fontsize=11)
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels(names, fontsize=12)
        ax.set_xlabel('% of Extra-cortical cells', fontsize=14)
        ax.invert_yaxis()
    ax.set_title('Subclass composition of Extra-cortical', fontsize=18)

    # Panel 6: Compare smoothed vs raw on multiple thresholds
    ax = axes[1, 2]
    for ft in frac_thresholds:
        extra = neighbor_ood_frac >= ft
        layers_temp = assign_discrete_layers(pred_depth)
        layers_temp[extra] = 'Extra-cortical'
        # Count by original layer assignment
        orig_layers = assign_discrete_layers(pred_depth)
        layer_names = list(LAYER_BINS.keys())
        frac_stolen = []
        for ln in layer_names:
            in_layer = orig_layers == ln
            n_layer = in_layer.sum()
            n_stolen = (in_layer & extra).sum()
            frac_stolen.append(100 * n_stolen / n_layer if n_layer > 0 else 0)
        ax.plot(layer_names, frac_stolen, 'o-', label=f'frac >= {ft:.1f}')
    ax.set_xlabel('Original Layer', fontsize=14)
    ax.set_ylabel('% reassigned to Extra-cortical', fontsize=14)
    ax.set_title('Layer impact by threshold', fontsize=18)
    ax.legend(fontsize=12)
    ax.set_ylim(0, None)

    plt.suptitle(f'{sample_id}: Smoothed OOD for Extra-cortical Detection',
                 fontsize=26)
    plt.tight_layout()

    fig_path = os.path.join(OUTPUT_DIR, f"ood_smoothed_{sample_id}.png")
    plt.savefig(fig_path, dpi=150, facecolor='white', bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {fig_path}")
    print(f"Total: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
