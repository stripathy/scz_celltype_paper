#!/usr/bin/env python3
"""
Refined OOD visualization with spatial context.

Strategy: Extra-cortical = OOD AND predicted_depth < threshold
This avoids flagging vascular cells deep in the cortex as extra-cortical.

Usage:
    python3 -u ood_refined_viz.py [sample_id]
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
    load_model, predict_depth, assign_discrete_layers,
    LAYER_BINS, LAYER_COLORS
)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H5AD_DIR = os.path.join(BASE_DIR, "output", "h5ad")
OUTPUT_DIR = os.path.join(BASE_DIR, "output")


def assign_layers_refined(pred_depth, ood_scores, ood_threshold,
                           depth_cutoff=0.05):
    """
    Assign layers with Extra-cortical for cells that are BOTH:
      1. OOD (ood_score > threshold)
      2. Near pial surface (predicted_depth < depth_cutoff)

    This avoids flagging vascular cells deep in cortex as extra-cortical.
    """
    layers = assign_discrete_layers(pred_depth)
    # Extra-cortical: OOD AND superficial
    extra_mask = (ood_scores > ood_threshold) & (pred_depth < depth_cutoff)
    layers[extra_mask] = 'Extra-cortical'
    return layers


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

    ood_mask = ood_scores > thresh_99
    coords = adata_pass.obsm['spatial']
    x, y = coords[:, 0], coords[:, 1]
    subclass = adata_pass.obs['subclass_label'].values

    # Try multiple depth cutoffs to find the sweet spot
    cutoffs = [-0.05, 0.0, 0.05, 0.10, 0.15]

    fig, axes = plt.subplots(2, len(cutoffs), figsize=(8*len(cutoffs), 16))

    for ci, depth_cut in enumerate(cutoffs):
        layers = assign_layers_refined(pred_depth, ood_scores, thresh_99,
                                        depth_cutoff=depth_cut)
        extra_mask = layers == 'Extra-cortical'
        n_extra = extra_mask.sum()

        # Top row: spatial
        ax = axes[0, ci]
        ax.set_facecolor('black')
        all_layers = list(LAYER_BINS.keys()) + ['Extra-cortical']
        for lname in all_layers:
            mask = layers == lname
            if mask.sum() > 0:
                c = LAYER_COLORS.get(lname, (0.5, 0.5, 0.5))
                ax.scatter(x[mask], y[mask], c=[c], s=0.1, alpha=0.4,
                           rasterized=True)
        ax.set_title(f'depth < {depth_cut}\nExtra-cortical: {n_extra:,} '
                     f'({100*n_extra/n_total:.1f}%)', fontsize=16)
        ax.set_aspect('equal')
        ax.invert_yaxis()
        if ci == 0:
            patches = [mpatches.Patch(color=LAYER_COLORS.get(l, (0.5, 0.5, 0.5)),
                                      label=l) for l in all_layers]
            ax.legend(handles=patches, fontsize=9, loc='upper right')

        # Bottom row: subclass breakdown of extra-cortical
        ax = axes[1, ci]
        if n_extra > 0:
            sc_counts = {}
            for sc in set(subclass[extra_mask]):
                sc_counts[sc] = (subclass[extra_mask] == sc).sum()
            sorted_sc = sorted(sc_counts.items(), key=lambda x: -x[1])
            names = [s[0] for s in sorted_sc[:10]]
            counts = [s[1] for s in sorted_sc[:10]]
            ax.barh(range(len(names)), counts, color='#e94560')
            ax.set_yticks(range(len(names)))
            ax.set_yticklabels(names, fontsize=11)
            ax.set_xlabel('# Extra-cortical cells', fontsize=12)
            ax.invert_yaxis()
        ax.set_title(f'Top subclasses (depth < {depth_cut})', fontsize=14)

    plt.suptitle(f'{sample_id}: Refined Extra-cortical (OOD + depth cutoff)',
                 fontsize=24)
    plt.tight_layout()

    fig_path = os.path.join(OUTPUT_DIR, f"ood_refined_{sample_id}.png")
    plt.savefig(fig_path, dpi=150, facecolor='white', bbox_inches='tight')
    plt.close()

    # Also print a clean summary for the best cutoff
    print("\nRefined Extra-cortical summary (varying depth cutoff):")
    for depth_cut in cutoffs:
        layers = assign_layers_refined(pred_depth, ood_scores, thresh_99,
                                        depth_cutoff=depth_cut)
        n_extra = (layers == 'Extra-cortical').sum()
        print(f"  depth < {depth_cut:+.2f}: {n_extra:>5,} cells "
              f"({100*n_extra/n_total:.1f}%)")

    print(f"\nSaved: {fig_path}")
    print(f"Total: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
