#!/usr/bin/env python3
"""
Test adaptive OOD thresholding: instead of using MERFISH test-set percentiles,
use the Xenium sample's OWN OOD distribution to find outliers.

Rationale: There's a platform shift between MERFISH and Xenium, so all Xenium
cells are somewhat more distant from MERFISH training data. We want to find
the Xenium cells that are unusually far even by Xenium standards.

Approach: Fit a mixture model or just use high percentiles of the Xenium
OOD distribution itself.

Usage:
    python3 -u ood_adaptive_thresh.py [sample_id]
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


def main():
    t0 = time.time()
    sample_id = sys.argv[1] if len(sys.argv) > 1 else "Br8667"

    model_bundle = load_model(os.path.join(OUTPUT_DIR, "depth_model_normalized.pkl"))
    thresh_99_merfish = model_bundle['ood_threshold_99']

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

    # Xenium-adaptive thresholds
    pctls = [95, 97, 98, 99, 99.5]
    print(f"\nOOD score distribution for {sample_id}:")
    print(f"  MERFISH 99th threshold: {thresh_99_merfish:.4f}")
    print(f"  Xenium stats: median={np.median(ood_scores):.4f}, "
          f"mean={np.mean(ood_scores):.4f}")
    for p in pctls:
        t = np.percentile(ood_scores, p)
        n_ood = (ood_scores > t).sum()
        print(f"  Xenium {p}th pctl = {t:.4f}: {n_ood:,} OOD ({100*n_ood/n_total:.1f}%)")

    # MAD-based threshold
    med = np.median(ood_scores)
    mad = np.median(np.abs(ood_scores - med))
    for n_mad in [3, 4, 5, 7, 10]:
        t_mad = med + n_mad * 1.4826 * mad  # 1.4826 scales MAD to std
        n_ood = (ood_scores > t_mad).sum()
        print(f"  {n_mad}-MAD = {t_mad:.4f}: {n_ood:,} OOD ({100*n_ood/n_total:.1f}%)")

    # ── Figure: sweep Xenium-adaptive thresholds ───────────────────
    test_thresholds = {
        f'MERFISH 99th ({thresh_99_merfish:.3f})': thresh_99_merfish,
        f'Xenium 99th ({np.percentile(ood_scores, 99):.3f})': np.percentile(ood_scores, 99),
        f'Xenium 99.5th ({np.percentile(ood_scores, 99.5):.3f})': np.percentile(ood_scores, 99.5),
        f'5-MAD ({med + 5 * 1.4826 * mad:.3f})': med + 5 * 1.4826 * mad,
        f'7-MAD ({med + 7 * 1.4826 * mad:.3f})': med + 7 * 1.4826 * mad,
    }

    n_thresh = len(test_thresholds)
    fig, axes = plt.subplots(2, n_thresh, figsize=(8*n_thresh, 16))

    for ti, (tname, tval) in enumerate(test_thresholds.items()):
        ood_mask = ood_scores > tval
        layers = assign_discrete_layers(pred_depth)
        layers[ood_mask] = 'Extra-cortical'
        n_extra = ood_mask.sum()

        # Top row: spatial
        ax = axes[0, ti]
        ax.set_facecolor('black')
        all_layers = list(LAYER_BINS.keys()) + ['Extra-cortical']
        for lname in all_layers:
            mask = layers == lname
            if mask.sum() > 0:
                c = LAYER_COLORS.get(lname, (0.5, 0.5, 0.5))
                ax.scatter(x[mask], y[mask], c=[c], s=0.1, alpha=0.4,
                           rasterized=True)
        ax.set_title(f'{tname}\nExtra: {n_extra:,} ({100*n_extra/n_total:.1f}%)',
                     fontsize=14)
        ax.set_aspect('equal')
        ax.invert_yaxis()
        if ti == 0:
            patches = [mpatches.Patch(color=LAYER_COLORS.get(l, (0.5, 0.5, 0.5)),
                                      label=l) for l in all_layers]
            ax.legend(handles=patches, fontsize=9, loc='upper right')

        # Bottom row: subclass composition
        ax = axes[1, ti]
        if n_extra > 0:
            is_vasc = np.isin(subclass, ['Endothelial', 'VLMC'])
            n_vasc_ood = (ood_mask & is_vasc).sum()
            n_nonvasc_ood = (ood_mask & ~is_vasc).sum()

            sc_counts = {}
            for sc in sorted(set(subclass[ood_mask])):
                sc_counts[sc] = (subclass[ood_mask] == sc).sum()
            sorted_sc = sorted(sc_counts.items(), key=lambda x: -x[1])[:10]
            names = [s[0] for s in sorted_sc]
            counts = [s[1] for s in sorted_sc]
            colors = ['#ff6633' if n in ['Endothelial', 'VLMC'] else '#e94560'
                      for n in names]
            ax.barh(range(len(names)), counts, color=colors)
            ax.set_yticks(range(len(names)))
            ax.set_yticklabels(names, fontsize=10)
            ax.set_xlabel('# cells', fontsize=12)
            ax.invert_yaxis()
            ax.set_title(f'Vasc: {n_vasc_ood} / Non-vasc: {n_nonvasc_ood}',
                         fontsize=12)

    plt.suptitle(f'{sample_id}: Adaptive OOD Threshold Comparison', fontsize=24)
    plt.tight_layout()

    fig_path = os.path.join(OUTPUT_DIR, f"ood_adaptive_{sample_id}.png")
    plt.savefig(fig_path, dpi=150, facecolor='white', bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {fig_path}")
    print(f"Total: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
