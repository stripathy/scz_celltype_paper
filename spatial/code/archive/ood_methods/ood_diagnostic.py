#!/usr/bin/env python3
"""
Diagnostic: examine spatial distribution of OOD cells by cell type.
Are they concentrated at pia or scattered (vascular)?

Usage:
    python3 -u ood_diagnostic.py [sample_id]
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
    subclass = adata_pass.obs['subclass_label'].values
    class_label = adata_pass.obs['class_label'].values
    coords = adata_pass.obsm['spatial']
    x, y = coords[:, 0], coords[:, 1]

    # Split OOD into vascular vs non-vascular
    vascular_types = {'Endothelial', 'VLMC'}
    is_vascular = np.isin(subclass, list(vascular_types))
    is_non_neuronal = np.isin(class_label, ['Non-neuronal'])

    ood_vascular = ood_mask & is_vascular
    ood_non_vasc = ood_mask & ~is_vascular
    ood_neuronal = ood_mask & ~is_non_neuronal
    ood_nonneuronal = ood_mask & is_non_neuronal

    print(f"\nOOD breakdown:")
    print(f"  Total OOD: {ood_mask.sum():,} ({100*ood_mask.sum()/n_total:.1f}%)")
    print(f"  Vascular OOD (Endothelial+VLMC): {ood_vascular.sum():,}")
    print(f"  Non-vascular OOD: {ood_non_vasc.sum():,}")
    print(f"  Neuronal OOD: {ood_neuronal.sum():,}")
    print(f"  Non-neuronal OOD: {ood_nonneuronal.sum():,}")

    # ── Figure: 4 panels ──────────────────────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(24, 24))

    # Panel 1: OOD vascular vs non-vascular
    ax = axes[0, 0]
    ax.set_facecolor('black')
    ax.scatter(x[~ood_mask], y[~ood_mask], c='#222222', s=0.02, alpha=0.3,
               rasterized=True)
    ax.scatter(x[ood_non_vasc], y[ood_non_vasc], c='#ff3333', s=1.5, alpha=0.8,
               rasterized=True)
    ax.scatter(x[ood_vascular], y[ood_vascular], c='#33ff33', s=1.5, alpha=0.8,
               rasterized=True)
    patches = [mpatches.Patch(color='#ff3333',
                              label=f'OOD non-vasc ({ood_non_vasc.sum():,})'),
               mpatches.Patch(color='#33ff33',
                              label=f'OOD vascular ({ood_vascular.sum():,})')]
    ax.legend(handles=patches, fontsize=16, loc='upper right')
    ax.set_title(f'{sample_id}: OOD by vascular status', fontsize=22)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    # Panel 2: OOD neuronal vs non-neuronal
    ax = axes[0, 1]
    ax.set_facecolor('black')
    ax.scatter(x[~ood_mask], y[~ood_mask], c='#222222', s=0.02, alpha=0.3,
               rasterized=True)
    ax.scatter(x[ood_neuronal], y[ood_neuronal], c='#ff6600', s=1.5, alpha=0.8,
               rasterized=True)
    ax.scatter(x[ood_nonneuronal], y[ood_nonneuronal], c='#6666ff', s=1.5, alpha=0.8,
               rasterized=True)
    patches = [mpatches.Patch(color='#ff6600',
                              label=f'OOD neuronal ({ood_neuronal.sum():,})'),
               mpatches.Patch(color='#6666ff',
                              label=f'OOD non-neuronal ({ood_nonneuronal.sum():,})')]
    ax.legend(handles=patches, fontsize=16, loc='upper right')
    ax.set_title(f'{sample_id}: OOD by class', fontsize=22)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    # Panel 3: OOD score by subclass (box plot of top 10)
    ax = axes[1, 0]
    sc_names = sorted(set(subclass))
    medians = [(sc, np.median(ood_scores[subclass == sc])) for sc in sc_names]
    medians.sort(key=lambda x: -x[1])
    top10 = [m[0] for m in medians[:12]]
    data = [ood_scores[subclass == sc] for sc in top10]
    bp = ax.boxplot(data, labels=top10, vert=True, patch_artist=True,
                    showfliers=False)
    for patch in bp['boxes']:
        patch.set_facecolor('#4ecdc4')
    ax.axhline(thresh_99, color='red', linestyle='--', linewidth=2,
               label=f'99th pctl = {thresh_99:.4f}')
    ax.set_ylabel('OOD distance', fontsize=16)
    ax.set_title('OOD Score by Subclass (top 12 median)', fontsize=20)
    ax.tick_params(axis='x', rotation=45, labelsize=11)
    ax.legend(fontsize=14)

    # Panel 4: OOD score vs predicted depth, colored by vascular
    ax = axes[1, 1]
    ax.set_facecolor('#111111')
    n_plot = min(30000, n_total)
    idx = np.random.choice(n_total, n_plot, replace=False)
    vasc_plot = is_vascular[idx]
    ax.scatter(pred_depth[idx][~vasc_plot], ood_scores[idx][~vasc_plot],
               s=0.3, alpha=0.3, c='#4ecdc4', rasterized=True, label='Non-vascular')
    ax.scatter(pred_depth[idx][vasc_plot], ood_scores[idx][vasc_plot],
               s=0.5, alpha=0.5, c='#ff6633', rasterized=True, label='Vascular')
    ax.axhline(thresh_99, color='red', linestyle='--', linewidth=2)
    ax.set_xlabel('Predicted Depth', fontsize=16)
    ax.set_ylabel('OOD Distance', fontsize=16)
    ax.set_title('OOD vs Depth (vascular highlighted)', fontsize=20)
    ax.legend(fontsize=14)

    plt.suptitle(f'OOD Diagnostic: {sample_id}', fontsize=26)
    plt.tight_layout()

    fig_path = os.path.join(OUTPUT_DIR, f"ood_diagnostic_{sample_id}.png")
    plt.savefig(fig_path, dpi=150, facecolor='white', bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {fig_path}")
    print(f"Total: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
