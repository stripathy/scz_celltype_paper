#!/usr/bin/env python3
"""
Quick OOD visualization — loads pre-trained model, generates spatial plots.
Much faster than test_ood_one_sample.py (no retraining).

Usage:
    python3 -u quick_ood_viz.py [sample_id]
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from depth_model import (
    load_model, predict_depth, assign_layers_with_ood,
    LAYER_BINS, LAYER_COLORS
)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5AD_DIR = os.path.join(BASE_DIR, "output", "h5ad")
OUTPUT_DIR = os.path.join(BASE_DIR, "output")


def main():
    t0 = time.time()
    sample_id = sys.argv[1] if len(sys.argv) > 1 else "Br8667"

    # Load model
    model_path = os.path.join(OUTPUT_DIR, "depth_model_normalized.pkl")
    print(f"Loading model from {model_path}...")
    model_bundle = load_model(model_path)
    thresh_99 = model_bundle['ood_threshold_99']
    thresh_95 = model_bundle['ood_threshold_95']
    print(f"  OOD thresholds: 95th={thresh_95:.4f}, 99th={thresh_99:.4f}")

    # Load sample
    h5ad_path = os.path.join(H5AD_DIR, f"{sample_id}_annotated.h5ad")
    print(f"Loading {sample_id}...")
    adata = ad.read_h5ad(h5ad_path)
    qc_mask = adata.obs['qc_pass'].values.astype(bool)
    adata_pass = adata[qc_mask].copy()
    n_total = adata_pass.shape[0]
    print(f"  {n_total:,} QC-pass cells")

    # Predict
    print("Predicting depth + OOD scores...")
    t1 = time.time()
    pred_depth, ood_scores = predict_depth(adata_pass, model_bundle,
                                           subclass_col='subclass_label',
                                           compute_ood=True)
    print(f"  Prediction took {time.time()-t1:.0f}s")

    # Layer assignments
    layers = assign_layers_with_ood(pred_depth, ood_scores,
                                     model_bundle=model_bundle)
    ood_mask = ood_scores > thresh_99

    # Print summary
    print(f"\nOOD summary for {sample_id}:")
    print(f"  OOD cells: {ood_mask.sum():,} / {n_total:,} ({100*ood_mask.sum()/n_total:.1f}%)")
    print(f"  Layer distribution:")
    all_layers = list(LAYER_BINS.keys()) + ['Extra-cortical']
    for lname in all_layers:
        n = (layers == lname).sum()
        print(f"    {lname:20s}: {n:>6,} ({100*n/n_total:5.1f}%)")

    # ── Generate figures ───────────────────────────────────────────
    coords = adata_pass.obsm['spatial']
    x, y = coords[:, 0], coords[:, 1]

    fig, axes = plt.subplots(1, 3, figsize=(36, 12))

    # Panel 1: OOD score continuous
    ax = axes[0]
    ax.set_facecolor('black')
    sc = ax.scatter(x, y, c=np.log10(ood_scores + 1e-6), cmap='hot',
                    s=0.1, alpha=0.5, rasterized=True)
    plt.colorbar(sc, ax=ax, label='log10(OOD distance)', shrink=0.7)
    ax.set_title(f'{sample_id}: OOD Score', fontsize=22)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    # Panel 2: OOD binary
    ax = axes[1]
    ax.set_facecolor('black')
    ax.scatter(x[~ood_mask], y[~ood_mask], c='#33cccc', s=0.05, alpha=0.3,
               rasterized=True)
    ax.scatter(x[ood_mask], y[ood_mask], c='#ff3333', s=1.0, alpha=0.8,
               rasterized=True)
    # Legend
    patches = [mpatches.Patch(color='#33cccc',
                              label=f'In-dist ({(~ood_mask).sum():,})'),
               mpatches.Patch(color='#ff3333',
                              label=f'OOD ({ood_mask.sum():,})')]
    ax.legend(handles=patches, fontsize=16, loc='upper right')
    ax.set_title(f'{sample_id}: OOD Cells (99th pctl)', fontsize=22)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    # Panel 3: Layers with Extra-cortical
    ax = axes[2]
    ax.set_facecolor('black')
    for lname in all_layers:
        mask = layers == lname
        if mask.sum() > 0:
            c = LAYER_COLORS.get(lname, (0.5, 0.5, 0.5))
            ax.scatter(x[mask], y[mask], c=[c], s=0.1, alpha=0.4,
                       rasterized=True)
    # Legend
    patches = [mpatches.Patch(color=LAYER_COLORS.get(l, (0.5, 0.5, 0.5)),
                              label=f'{l} ({(layers==l).sum():,})')
               for l in all_layers if (layers == l).sum() > 0]
    ax.legend(handles=patches, fontsize=12, loc='upper right')
    ax.set_title(f'{sample_id}: Layers + Extra-cortical', fontsize=22)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    plt.tight_layout()
    fig_path = os.path.join(OUTPUT_DIR, f"ood_test_{sample_id}.png")
    plt.savefig(fig_path, dpi=150, facecolor='white', bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {fig_path}")
    print(f"Total time: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
