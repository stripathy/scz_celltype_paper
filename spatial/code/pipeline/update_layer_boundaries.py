#!/usr/bin/env python3
"""
Re-assign discrete layers using updated LAYER_BINS boundaries.

This is a lightweight alternative to re-running 05_run_spatial_domains.py:
it reuses the existing BANKSY domain classifications (banksy_domain,
banksy_is_l1) and only re-runs assign_discrete_layers() +
smooth_layers_spatial() with the new boundary values.

Avoids re-running BANKSY clustering (~15-20 min/sample), reducing
total runtime from ~6 hours to ~30-45 minutes.

Usage:
    python3 -u update_layer_boundaries.py
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import anndata as ad
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipeline_config import H5AD_DIR, OUTPUT_DIR, MODULES_DIR

sys.path.insert(0, MODULES_DIR)
from depth_model import assign_discrete_layers, smooth_layers_spatial, LAYER_BINS


def update_one_sample(h5ad_path):
    """Re-assign layers for one sample using existing BANKSY domains."""
    t0 = time.time()
    sample_id = os.path.basename(h5ad_path).replace("_annotated.h5ad", "")

    try:
        adata = ad.read_h5ad(h5ad_path)
        n_total = adata.shape[0]

        # Check that BANKSY columns exist
        required = ["banksy_domain", "banksy_is_l1", "predicted_norm_depth"]
        missing = [c for c in required if c not in adata.obs.columns]
        if missing:
            raise ValueError(f"Missing columns: {missing}. Run 05_run_spatial_domains.py first.")

        # Get QC mask (same logic as step 05)
        if "corr_qc_pass" in adata.obs.columns:
            qc_mask = adata.obs["corr_qc_pass"].values.astype(bool)
        elif "qc_pass" in adata.obs.columns:
            qc_mask = adata.obs["qc_pass"].values.astype(bool)
        else:
            qc_mask = np.ones(n_total, dtype=bool)

        n_pass = qc_mask.sum()

        # Extract QC-pass data
        pred_depth = adata.obs["predicted_norm_depth"].values[qc_mask]
        domains = adata.obs["banksy_domain"].values[qc_mask].astype(str)
        is_l1 = adata.obs["banksy_is_l1"].values[qc_mask].astype(bool)
        coords = adata.obsm["spatial"][qc_mask]

        # Re-assign layers with new boundaries
        depth_layers = assign_discrete_layers(pred_depth)

        # Hybrid: depth bins for all, override Vascular
        combined_layers = depth_layers.copy()
        combined_layers[domains == "Vascular"] = "Vascular"

        # Spatial smoothing (same 3-step pipeline as step 05)
        smoothed_layers = smooth_layers_spatial(
            coords=coords,
            layers=combined_layers,
            domains=domains,
            is_l1_banksy=is_l1,
            depths=pred_depth,
            verbose=False,
        )

        # spatial_domain: update WM based on smoothed layers
        spatial_domain = domains.copy()
        spatial_domain[smoothed_layers == "WM"] = "WM"

        # Write back to full adata
        full_spatial_domain = np.full(n_total, "Unassigned", dtype=object)
        full_spatial_domain[qc_mask] = spatial_domain

        full_layer = np.full(n_total, "Unassigned", dtype=object)
        full_layer[qc_mask] = smoothed_layers

        # Impute for QC-fail cells via KNN (same as step 05)
        n_fail = n_total - n_pass
        if n_fail > 0:
            coords_all = adata.obsm["spatial"][:, :2]
            pass_idx = np.where(qc_mask)[0]
            fail_idx = np.where(~qc_mask)[0]
            tree = cKDTree(coords_all[pass_idx])
            _, nn_idx = tree.query(coords_all[fail_idx], k=min(15, len(pass_idx)))
            for col_full, col_pass in [
                (full_spatial_domain, spatial_domain),
                (full_layer, smoothed_layers),
            ]:
                for j, fi in enumerate(fail_idx):
                    neighbors = nn_idx[j] if nn_idx.ndim > 1 else [nn_idx[j]]
                    votes = {}
                    for ni in neighbors:
                        v = col_pass[ni]
                        votes[v] = votes.get(v, 0) + 1
                    col_full[fi] = max(votes, key=votes.get)

        adata.obs["spatial_domain"] = full_spatial_domain
        adata.obs["layer"] = full_layer

        # Save
        adata.write_h5ad(h5ad_path)

        elapsed = time.time() - t0

        # Layer distribution summary
        layer_dist = {}
        for lname in list(LAYER_BINS.keys()) + ["Vascular"]:
            layer_dist[lname] = int((smoothed_layers == lname).sum())

        print(f"  [{sample_id}] {n_pass:,} cells, "
              + ", ".join(f"{k}={v:,}" for k, v in layer_dist.items())
              + f" ({elapsed:.0f}s)")

        return {"sample_id": sample_id, "status": "success",
                "n_pass": int(n_pass), "layer_dist": layer_dist, "time": elapsed}

    except Exception as e:
        import traceback
        elapsed = time.time() - t0
        print(f"  [{sample_id}] FAILED - {e}")
        traceback.print_exc()
        return {"sample_id": sample_id, "status": "error",
                "error": str(e), "time": elapsed}


def main():
    t_start = time.time()

    print(f"Updating layer boundaries using:")
    for lname, (lo, hi) in LAYER_BINS.items():
        print(f"  {lname}: ({lo}, {hi})")
    print()

    h5ad_files = sorted(
        os.path.join(H5AD_DIR, f)
        for f in os.listdir(H5AD_DIR)
        if f.endswith("_annotated.h5ad")
    )
    print(f"Found {len(h5ad_files)} samples\n")

    results = []
    for i, h5ad_path in enumerate(h5ad_files):
        sid = os.path.basename(h5ad_path).replace("_annotated.h5ad", "")
        print(f"[{i+1}/{len(h5ad_files)}] {sid}")
        result = update_one_sample(h5ad_path)
        results.append(result)

    # Summary
    total_time = time.time() - t_start
    success = [r for r in results if r["status"] == "success"]
    print(f"\n{'='*60}")
    print(f"SUMMARY - {total_time:.0f}s total")
    print(f"  Success: {len(success)}/{len(results)}")

    if success:
        # Aggregate layer distribution
        agg = {}
        total_pass = 0
        for r in success:
            total_pass += r["n_pass"]
            for lname, n in r["layer_dist"].items():
                agg[lname] = agg.get(lname, 0) + n

        print(f"\n  Combined layer distribution ({total_pass:,} QC-pass cells):")
        for lname in ["L1", "L2/3", "L4", "L5", "L6", "WM", "Vascular"]:
            n = agg.get(lname, 0)
            print(f"    {lname:20s}: {n:>8,} ({100*n/total_pass:5.1f}%)")

    # Merge into combined h5ad
    print("\nMerging all samples into combined h5ad...")
    t_merge = time.time()
    adatas = []
    for r in success:
        path = os.path.join(H5AD_DIR, f"{r['sample_id']}_annotated.h5ad")
        adatas.append(ad.read_h5ad(path))
    combined = ad.concat(adatas, join="outer")
    combined_path = os.path.join(OUTPUT_DIR, "all_samples_annotated.h5ad")
    combined.write_h5ad(combined_path)
    file_size = os.path.getsize(combined_path) / 1e9
    print(f"  Saved: {combined_path} ({file_size:.2f} GB, {combined.shape})")
    print(f"  Merge took {time.time()-t_merge:.0f}s")

    print(f"\nTotal time: {time.time()-t_start:.0f}s")


if __name__ == "__main__":
    main()
