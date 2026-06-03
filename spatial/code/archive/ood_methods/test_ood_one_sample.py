#!/usr/bin/env python3
"""
Test OOD (out-of-distribution) scoring on one example Xenium section.

Steps:
  1. Load MERFISH reference and retrain depth model (with OOD 1-NN)
  2. Load one Xenium sample
  3. Predict depth + compute OOD scores
  4. Generate diagnostic figures:
     - Spatial plot colored by OOD score
     - Spatial plot colored by predicted layer (with Extra-cortical)
     - OOD score distribution by cell type
     - OOD score vs predicted depth scatter
  5. Print summary stats

Usage:
    python test_ood_one_sample.py [sample_id]
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patches as mpatches

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from depth_model import (
    train_depth_model, save_model, predict_depth,
    compute_ood_scores, assign_layers_with_ood,
    build_neighborhood_features, LAYER_BINS, LAYER_COLORS
)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MERFISH_PATH = os.path.join(BASE_DIR, "sea-ad_reference",
                            "SEAAD_MTG_MERFISH.2024-12-11.h5ad")
H5AD_DIR = os.path.join(BASE_DIR, "output", "h5ad")
OUTPUT_DIR = os.path.join(BASE_DIR, "output")


def main():
    t_start = time.time()

    # Pick sample
    sample_id = sys.argv[1] if len(sys.argv) > 1 else "Br8667"
    print(f"Testing OOD scoring on sample: {sample_id}")
    print(f"{'='*60}")

    # ── Step 1: Retrain depth model with OOD 1-NN ──────────────────
    print("\n[Step 1] Retraining depth model with OOD scoring...")
    t0 = time.time()
    merfish = ad.read_h5ad(MERFISH_PATH)
    print(f"  MERFISH loaded: {merfish.shape}")
    print(f"  Loading took {time.time()-t0:.0f}s")

    model_bundle = train_depth_model(merfish, K=50)

    # Save the new model
    model_path = os.path.join(OUTPUT_DIR, "depth_model_normalized.pkl")
    save_model(model_bundle, model_path)
    print(f"  Model saved: {model_path}")
    print(f"  OOD threshold (99th pctl): {model_bundle['ood_threshold_99']:.4f}")
    print(f"  OOD threshold (95th pctl): {model_bundle['ood_threshold_95']:.4f}")

    # ── Step 2: Load one Xenium sample ─────────────────────────────
    print(f"\n[Step 2] Loading {sample_id}...")
    h5ad_path = os.path.join(H5AD_DIR, f"{sample_id}_annotated.h5ad")
    adata = ad.read_h5ad(h5ad_path)
    print(f"  Shape: {adata.shape}")
    print(f"  QC pass: {adata.obs['qc_pass'].sum():,} / {len(adata):,}")

    # Subset to QC-pass cells
    qc_mask = adata.obs['qc_pass'].values.astype(bool)
    adata_pass = adata[qc_mask].copy()
    print(f"  Working with {adata_pass.shape[0]:,} QC-pass cells")

    # ── Step 3: Predict depth + OOD scores ─────────────────────────
    print(f"\n[Step 3] Predicting depth and OOD scores...")
    t1 = time.time()
    pred_depth, ood_scores = predict_depth(adata_pass, model_bundle,
                                           subclass_col='subclass_label',
                                           compute_ood=True)
    print(f"  Prediction took {time.time()-t1:.0f}s")

    # Add to adata
    adata_pass.obs['predicted_norm_depth'] = pred_depth
    adata_pass.obs['ood_score'] = ood_scores

    # ── Step 4: Summary statistics ─────────────────────────────────
    print(f"\n[Step 4] OOD Score Summary")
    print(f"  {'='*50}")
    thresh_99 = model_bundle['ood_threshold_99']
    thresh_95 = model_bundle['ood_threshold_95']

    n_ood_99 = (ood_scores > thresh_99).sum()
    n_ood_95 = (ood_scores > thresh_95).sum()
    n_total = len(ood_scores)

    print(f"  Total QC-pass cells: {n_total:,}")
    print(f"  OOD (>99th pctl, {thresh_99:.4f}): {n_ood_99:,} ({100*n_ood_99/n_total:.1f}%)")
    print(f"  OOD (>95th pctl, {thresh_95:.4f}): {n_ood_95:,} ({100*n_ood_95/n_total:.1f}%)")
    print(f"  OOD score: median={np.median(ood_scores):.4f}, "
          f"mean={np.mean(ood_scores):.4f}, max={np.max(ood_scores):.4f}")

    # OOD by cell type
    print(f"\n  OOD rate by subclass (>99th pctl threshold):")
    ood_mask_99 = ood_scores > thresh_99
    subclass = adata_pass.obs['subclass_label'].values
    for sc in sorted(set(subclass)):
        sc_mask = subclass == sc
        n_sc = sc_mask.sum()
        n_sc_ood = (sc_mask & ood_mask_99).sum()
        if n_sc > 0:
            pct = 100 * n_sc_ood / n_sc
            flag = " <<<" if pct > 5 else ""
            print(f"    {sc:25s}: {n_sc_ood:>5,}/{n_sc:>6,} ({pct:5.1f}%){flag}")

    # OOD by current layer assignment
    print(f"\n  OOD rate by current layer:")
    old_depth = adata_pass.obs['predicted_norm_depth'].values
    for lname, (lo, hi) in LAYER_BINS.items():
        layer_mask = (old_depth >= lo) & (old_depth < hi)
        n_layer = layer_mask.sum()
        n_layer_ood = (layer_mask & ood_mask_99).sum()
        if n_layer > 0:
            pct = 100 * n_layer_ood / n_layer
            print(f"    {lname:15s}: {n_layer_ood:>5,}/{n_layer:>6,} ({pct:5.1f}%)")

    # Assign layers with OOD
    layers_new = assign_layers_with_ood(
        pred_depth, ood_scores, model_bundle=model_bundle
    )
    adata_pass.obs['layer_with_ood'] = layers_new

    print(f"\n  Layer distribution (with Extra-cortical):")
    for lname in list(LAYER_BINS.keys()) + ['Extra-cortical']:
        n = (layers_new == lname).sum()
        pct = 100 * n / n_total
        print(f"    {lname:20s}: {n:>6,} ({pct:5.1f}%)")

    # ── Step 5: Diagnostic figures ─────────────────────────────────
    print(f"\n[Step 5] Generating diagnostic figures...")

    coords = adata_pass.obsm['spatial']
    x, y = coords[:, 0], coords[:, 1]

    fig, axes = plt.subplots(2, 3, figsize=(30, 20))

    # Panel 1: Spatial OOD score (continuous)
    ax = axes[0, 0]
    ax.set_facecolor('black')
    # Use log scale for better contrast
    ood_plot = np.log10(ood_scores + 1e-6)
    sc = ax.scatter(x, y, c=ood_plot, cmap='hot', s=0.1, alpha=0.5, rasterized=True)
    plt.colorbar(sc, ax=ax, label='log10(OOD distance)', shrink=0.7)
    ax.axhline(y=0, color='white', alpha=0.1)
    ax.set_title(f'{sample_id}: OOD Score (log10)', fontsize=20)
    ax.set_xlabel('X', fontsize=14)
    ax.set_ylabel('Y', fontsize=14)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    # Panel 2: Spatial OOD binary (99th threshold)
    ax = axes[0, 1]
    ax.set_facecolor('black')
    ax.scatter(x[~ood_mask_99], y[~ood_mask_99], c='#33cccc', s=0.05, alpha=0.3,
               rasterized=True, label=f'In-distribution ({(~ood_mask_99).sum():,})')
    ax.scatter(x[ood_mask_99], y[ood_mask_99], c='#ff3333', s=1.0, alpha=0.8,
               rasterized=True, label=f'OOD ({ood_mask_99.sum():,})')
    ax.legend(fontsize=14, loc='upper right', markerscale=10)
    ax.set_title(f'{sample_id}: OOD Cells (99th pctl)', fontsize=20)
    ax.set_xlabel('X', fontsize=14)
    ax.set_ylabel('Y', fontsize=14)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    # Panel 3: Spatial layers WITH Extra-cortical
    ax = axes[0, 2]
    ax.set_facecolor('black')
    all_layers = list(LAYER_BINS.keys()) + ['Extra-cortical']
    for lname in all_layers:
        mask = layers_new == lname
        if mask.sum() > 0:
            c = LAYER_COLORS.get(lname, (0.5, 0.5, 0.5))
            ax.scatter(x[mask], y[mask], c=[c], s=0.1, alpha=0.4,
                       rasterized=True, label=f'{lname} ({mask.sum():,})')
    ax.legend(fontsize=11, loc='upper right', markerscale=10)
    ax.set_title(f'{sample_id}: Layers with Extra-cortical', fontsize=20)
    ax.set_xlabel('X', fontsize=14)
    ax.set_ylabel('Y', fontsize=14)
    ax.set_aspect('equal')
    ax.invert_yaxis()

    # Panel 4: OOD score distribution (histogram)
    ax = axes[1, 0]
    ax.hist(ood_scores, bins=200, color='#4ecdc4', alpha=0.8, edgecolor='none')
    ax.axvline(thresh_99, color='red', linestyle='--', linewidth=2,
               label=f'99th pctl = {thresh_99:.4f}')
    ax.axvline(thresh_95, color='orange', linestyle='--', linewidth=2,
               label=f'95th pctl = {thresh_95:.4f}')
    ax.set_xlabel('OOD Distance (1-NN to MERFISH training)', fontsize=14)
    ax.set_ylabel('# Cells', fontsize=14)
    ax.set_title('OOD Score Distribution', fontsize=20)
    ax.legend(fontsize=14)
    ax.set_xlim(0, min(np.percentile(ood_scores, 99.9) * 2, ood_scores.max()))

    # Panel 5: OOD score vs predicted depth
    ax = axes[1, 1]
    ax.set_facecolor('#111111')
    # Subsample for plotting
    n_plot = min(50000, n_total)
    idx = np.random.choice(n_total, n_plot, replace=False)
    ax.scatter(pred_depth[idx], ood_scores[idx], s=0.2, alpha=0.3,
               c='#4ecdc4', rasterized=True)
    ax.axhline(thresh_99, color='red', linestyle='--', linewidth=2,
               label=f'99th pctl threshold')
    ax.set_xlabel('Predicted Normalized Depth', fontsize=14)
    ax.set_ylabel('OOD Distance', fontsize=14)
    ax.set_title('OOD Score vs Predicted Depth', fontsize=20)
    ax.legend(fontsize=14)

    # Panel 6: OOD rate by subclass (bar chart)
    ax = axes[1, 2]
    sc_names = sorted(set(subclass))
    ood_rates = []
    sc_counts = []
    for sc_name in sc_names:
        sc_m = subclass == sc_name
        n_sc = sc_m.sum()
        n_ood = (sc_m & ood_mask_99).sum()
        ood_rates.append(100 * n_ood / n_sc if n_sc > 0 else 0)
        sc_counts.append(n_sc)

    # Sort by OOD rate
    sort_idx = np.argsort(ood_rates)[::-1]
    sc_sorted = [sc_names[i] for i in sort_idx]
    rates_sorted = [ood_rates[i] for i in sort_idx]

    bars = ax.barh(range(len(sc_sorted)), rates_sorted, color='#e94560', alpha=0.8)
    ax.set_yticks(range(len(sc_sorted)))
    ax.set_yticklabels(sc_sorted, fontsize=10)
    ax.set_xlabel('% Cells OOD (>99th pctl)', fontsize=14)
    ax.set_title('OOD Rate by Subclass', fontsize=20)
    ax.invert_yaxis()

    plt.suptitle(f'OOD Analysis: {sample_id}', fontsize=24, y=1.01)
    plt.tight_layout()

    fig_path = os.path.join(OUTPUT_DIR, f"ood_test_{sample_id}.png")
    plt.savefig(fig_path, dpi=150, facecolor='white', bbox_inches='tight')
    plt.close()
    print(f"  Saved: {fig_path}")

    # ── Also save a zoomed-in view of a region with high OOD ──────
    if ood_mask_99.sum() > 10:
        # Find the densest OOD cluster
        ood_x = x[ood_mask_99]
        ood_y = y[ood_mask_99]
        # Use median of OOD cells as center
        cx, cy = np.median(ood_x), np.median(ood_y)
        # Zoom window: ±500 microns around the OOD cluster center
        margin = 1000
        zoom_mask = (np.abs(x - cx) < margin) & (np.abs(y - cy) < margin)

        if zoom_mask.sum() > 100:
            fig2, ax2 = plt.subplots(1, 1, figsize=(14, 14))
            ax2.set_facecolor('black')

            zoom_layers = layers_new[zoom_mask]
            zoom_x = x[zoom_mask]
            zoom_y = y[zoom_mask]

            for lname in all_layers:
                mask = zoom_layers == lname
                if mask.sum() > 0:
                    c = LAYER_COLORS.get(lname, (0.5, 0.5, 0.5))
                    ax2.scatter(zoom_x[mask], zoom_y[mask], c=[c], s=3,
                                alpha=0.7, rasterized=True,
                                label=f'{lname} ({mask.sum():,})')

            ax2.legend(fontsize=14, loc='upper right', markerscale=5)
            ax2.set_title(f'{sample_id}: Zoomed Region (layers + OOD)', fontsize=20)
            ax2.set_aspect('equal')
            ax2.invert_yaxis()

            zoom_path = os.path.join(OUTPUT_DIR, f"ood_test_{sample_id}_zoom.png")
            plt.savefig(zoom_path, dpi=150, facecolor='white', bbox_inches='tight')
            plt.close()
            print(f"  Saved: {zoom_path}")

    total_time = time.time() - t_start
    print(f"\nTotal time: {total_time:.0f}s")


if __name__ == "__main__":
    main()
