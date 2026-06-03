#!/usr/bin/env python3
"""
Apply depth prediction model to ALL SEA-AD MERFISH cells.

The original MERFISH h5ad has manual 'Normalized depth from pia' annotations
for only ~19.5% of cells. This script uses the trained depth model (based on
KNN subclass composition features) to predict depth for ALL cells, then
assigns layer labels using the same hybrid approach as Xenium (depth bins +
Vascular isn't relevant here since MERFISH doesn't have vascular tissue).

Adds to the MERFISH h5ad:
  - predicted_norm_depth: model-predicted depth for all cells
  - predicted_layer: discrete layer from depth bins (L1, L2/3, L4, L5, L6, WM)
  - ood_score: OOD distance score per cell

Usage:
    python3 -u annotate_merfish_depth.py
"""

import os
import sys
import time
import numpy as np
import anndata as ad

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Navigate from code/archive/one_time_utils/ up to project root
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(SCRIPT_DIR)))
MODULES_DIR = os.path.join(BASE_DIR, "code", "modules")
sys.path.insert(0, MODULES_DIR)

from depth_model import (
    load_model, build_neighborhood_features, assign_discrete_layers
)

MERFISH_PATH = os.path.join(BASE_DIR, "data", "reference",
                            "SEAAD_MTG_MERFISH.2024-12-11.h5ad")
MODEL_PATH = os.path.join(BASE_DIR, "output", "depth_model_normalized.pkl")


def main():
    t0 = time.time()

    # Load model
    print("Loading depth model...")
    model_bundle = load_model(MODEL_PATH)
    subclass_names = model_bundle['subclass_names']
    K = model_bundle['K']
    model = model_bundle['model']
    print(f"  K={K}, {len(subclass_names)} subclasses")

    # Load MERFISH data (fully into memory since we'll modify it)
    print(f"\nLoading MERFISH data: {MERFISH_PATH}")
    merfish = ad.read_h5ad(MERFISH_PATH)
    n_cells = merfish.n_obs
    print(f"  {n_cells:,} cells, {merfish.obs['Section'].nunique()} sections")

    # Extract data
    coords = merfish.obsm['X_spatial_raw']
    subclass = merfish.obs['Subclass'].values.astype(str)
    sections = merfish.obs['Section'].values.astype(str)

    # Check existing depth annotations
    manual_depth = merfish.obs['Normalized depth from pia'].values.astype(float)
    n_manual = (~np.isnan(manual_depth)).sum()
    print(f"  Manual depth annotations: {n_manual:,} ({n_manual/n_cells*100:.1f}%)")

    # Build KNN features per section (this is the slow step)
    print(f"\nBuilding K={K} neighborhood features per section...")
    t_feat = time.time()
    features = build_neighborhood_features(
        coords, subclass, subclass_names, K=K, sections=sections
    )
    print(f"  Feature building: {time.time() - t_feat:.0f}s")

    # Predict depth for all cells
    print("\nPredicting depth for all cells...")
    t_pred = time.time()
    pred_depth = model.predict(features).astype(np.float32)
    print(f"  Prediction: {time.time() - t_pred:.0f}s")

    # Assign discrete layers (depth bins only — no Vascular/OOD for MERFISH)
    pred_layers = assign_discrete_layers(pred_depth)

    # Stats
    print(f"\n  Predicted depth range: [{pred_depth.min():.3f}, {pred_depth.max():.3f}]")

    # Compare predicted vs manual where available
    has_manual = ~np.isnan(manual_depth)
    if has_manual.sum() > 0:
        from scipy.stats import pearsonr
        from sklearn.metrics import r2_score, mean_absolute_error
        r, _ = pearsonr(manual_depth[has_manual], pred_depth[has_manual])
        r2 = r2_score(manual_depth[has_manual], pred_depth[has_manual])
        mae = mean_absolute_error(manual_depth[has_manual], pred_depth[has_manual])
        print(f"\n  Predicted vs manual depth (where available):")
        print(f"    r = {r:.4f}, R² = {r2:.4f}, MAE = {mae:.4f}")

    # Layer distribution
    from collections import Counter
    layer_counts = Counter(pred_layers)
    print(f"\n  Predicted layer distribution:")
    for l in ['L1', 'L2/3', 'L4', 'L5', 'L6', 'WM']:
        n = layer_counts.get(l, 0)
        print(f"    {l:6s}: {n:>9,} ({n/n_cells*100:5.1f}%)")

    # Add columns to MERFISH h5ad
    merfish.obs['predicted_norm_depth'] = pred_depth
    merfish.obs['predicted_layer'] = pred_layers

    # Save
    print(f"\nSaving updated MERFISH h5ad...")
    t_save = time.time()
    merfish.write_h5ad(MERFISH_PATH)
    print(f"  Save: {time.time() - t_save:.0f}s")

    total = time.time() - t0
    print(f"\nDone in {total:.0f}s")


if __name__ == "__main__":
    main()
