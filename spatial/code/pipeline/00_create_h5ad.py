#!/usr/bin/env python3
"""
Step 0: Create initial h5ad files from raw Xenium data.

For each sample, loads the raw .h5 expression matrix and cell boundary
CSV, computes cell centroids, and saves an AnnData (.h5ad) file with:
  - .X: sparse count matrix (cells × genes, Gene Expression only)
  - .obsm['spatial']: (x, y) centroid coordinates
  - .obs['sample_id']: sample identifier (e.g. 'Br8667')

This is the first step in the pipeline. Subsequent steps add QC columns
(step 01), cell type labels (step 02), transcript export (step 03),
depth predictions (step 04), and spatial domain/layer assignments (step 05).

Usage:
    python3 -u 00_create_h5ad.py
"""

import os
import sys
import time

# Pipeline config
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipeline_config import RAW_DIR, H5AD_DIR, MODULES_DIR

# Modules
sys.path.insert(0, MODULES_DIR)
from loading import discover_samples, load_xenium_sample


def main():
    t_start = time.time()

    os.makedirs(H5AD_DIR, exist_ok=True)

    # Discover all raw samples
    samples = discover_samples(RAW_DIR)
    print(f"Found {len(samples)} samples in {RAW_DIR}")
    print(f"Output directory: {H5AD_DIR}")
    print(f"{'='*60}")

    results = []
    for i, info in enumerate(samples):
        sid = info["sample_id"]
        t0 = time.time()
        print(f"\n[{i+1}/{len(samples)}] {sid}")

        try:
            # Load raw expression + spatial coordinates
            adata = load_xenium_sample(info["h5_path"], info["boundaries_path"])
            n_cells = adata.shape[0]
            n_genes = adata.shape[1]

            # Save initial h5ad
            out_path = os.path.join(H5AD_DIR, f"{sid}_annotated.h5ad")
            adata.write_h5ad(out_path)

            elapsed = time.time() - t0
            file_mb = os.path.getsize(out_path) / 1e6
            print(f"  {n_cells:,} cells × {n_genes} genes -> "
                  f"{out_path} ({file_mb:.0f} MB, {elapsed:.1f}s)")
            results.append({"sample_id": sid, "n_cells": n_cells,
                            "n_genes": n_genes, "status": "success"})

        except Exception as e:
            print(f"  ERROR: {e}")
            results.append({"sample_id": sid, "status": "error", "error": str(e)})

    # Summary
    total_time = time.time() - t_start
    n_ok = sum(1 for r in results if r["status"] == "success")
    total_cells = sum(r.get("n_cells", 0) for r in results)

    print(f"\n{'='*60}")
    print(f"SUMMARY")
    print(f"{'='*60}")
    print(f"  Success: {n_ok}/{len(results)} samples")
    print(f"  Total cells: {total_cells:,}")
    print(f"  Output: {H5AD_DIR}")
    print(f"  Time: {total_time:.0f}s")


if __name__ == "__main__":
    main()
