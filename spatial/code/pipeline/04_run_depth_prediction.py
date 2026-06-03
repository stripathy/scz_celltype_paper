#!/usr/bin/env python3
"""
Step 4: Retrain depth model from MERFISH and predict depth for all Xenium samples.

This retrains the GradientBoostingRegressor depth model from the SEA-AD MERFISH
reference (using K=50 neighborhood composition features), then applies it to all
24 Xenium samples. Uses corr_qc_pass (from step 02b) when available to exclude
low-margin and doublet-suspect cells from neighborhood features.

Saves:
  - output/depth_model_normalized.pkl (model bundle)
  - Updates predicted_norm_depth column in each h5ad file

Note: Domain classification (Vascular, WM, L1 border) is handled by step 05
using BANKSY spatial clustering (see banksy_domains.py).

Usage:
    python3 -u 04_run_depth_prediction.py
"""

import os
import sys
import time
import numpy as np
import anndata as ad
from multiprocessing import Pool, current_process

# Pipeline config
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipeline_config import (
    MERFISH_PATH, H5AD_DIR, DEPTH_MODEL_PATH, MODULES_DIR, N_WORKERS,
)

# Modules
sys.path.insert(0, MODULES_DIR)
from depth_model import (
    train_depth_model, save_model, load_model, predict_depth
)

MODEL_PATH = DEPTH_MODEL_PATH

_model_bundle = None


def _init_worker():
    global _model_bundle
    pid = current_process().pid
    print(f"  [Worker {pid}] Loading depth model...")
    _model_bundle = load_model(MODEL_PATH)
    print(f"  [Worker {pid}] Ready.")


def _process_sample(h5ad_path):
    global _model_bundle
    pid = current_process().pid
    t0 = time.time()
    sample_id = os.path.basename(h5ad_path).replace("_annotated.h5ad", "")

    try:
        adata = ad.read_h5ad(h5ad_path)
        # Use corr_qc_pass (spatial QC + margin filter + doublet exclusion);
        # fall back to qc_pass if step 02b not run
        if 'corr_qc_pass' in adata.obs.columns:
            qc_mask = adata.obs['corr_qc_pass'].values.astype(bool)
        else:
            qc_mask = adata.obs.get("qc_pass", np.ones(adata.shape[0], dtype=bool))
            if hasattr(qc_mask, 'values'):
                qc_mask = qc_mask.values.astype(bool)
        adata_pass = adata[qc_mask].copy()

        # Predict depth
        # Use correlation-derived subclass labels if available (from step 02b)
        subclass_col = ('corr_subclass' if 'corr_subclass' in adata_pass.obs.columns
                        else 'subclass_label')
        pred_depth = predict_depth(adata_pass, _model_bundle,
                                    subclass_col=subclass_col)

        # Store back
        adata.obs['predicted_norm_depth'] = np.nan
        adata.obs.iloc[np.where(qc_mask)[0],
                       adata.obs.columns.get_loc('predicted_norm_depth')] = pred_depth

        adata.write_h5ad(h5ad_path)

        elapsed = time.time() - t0
        print(f"  [{pid}] {sample_id}: {qc_mask.sum():,} cells, "
              f"depth [{pred_depth.min():.3f}, {pred_depth.max():.3f}], {elapsed:.0f}s")
        return {"sample_id": sample_id, "status": "success",
                "n_cells": int(qc_mask.sum()), "time": elapsed}
    except Exception as e:
        import traceback
        traceback.print_exc()
        return {"sample_id": sample_id, "status": "error", "error": str(e)}


def main():
    t_start = time.time()

    # Step 1: Retrain depth model
    print("="*60)
    print("Step 1: Retraining depth model from MERFISH reference")
    print("="*60)

    if not os.path.exists(MERFISH_PATH):
        print(f"ERROR: MERFISH reference not found at {MERFISH_PATH}")
        sys.exit(1)

    t0 = time.time()
    merfish = ad.read_h5ad(MERFISH_PATH)
    print(f"  MERFISH loaded: {merfish.shape} ({time.time()-t0:.0f}s)")

    model_bundle = train_depth_model(merfish, K=50)
    save_model(model_bundle, MODEL_PATH)
    print(f"  Model saved: {MODEL_PATH}")

    # Step 2: Predict depth for all samples
    print(f"\n{'='*60}")
    print(f"Step 2: Predicting depth for all samples ({N_WORKERS} workers)")
    print(f"{'='*60}")

    h5ad_files = sorted(
        os.path.join(H5AD_DIR, f)
        for f in os.listdir(H5AD_DIR)
        if f.endswith("_annotated.h5ad")
    )
    print(f"  Found {len(h5ad_files)} samples")

    with Pool(processes=N_WORKERS, initializer=_init_worker) as pool:
        results = pool.map(_process_sample, h5ad_files)

    n_ok = sum(1 for r in results if r["status"] == "success")
    total_cells = sum(r.get("n_cells", 0) for r in results if r["status"] == "success")

    print(f"\n{'='*60}")
    print(f"SUMMARY: {n_ok}/{len(results)} samples, {total_cells:,} cells")
    print(f"Total time: {time.time()-t_start:.0f}s")


if __name__ == "__main__":
    main()
