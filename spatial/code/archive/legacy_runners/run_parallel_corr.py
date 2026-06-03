"""
Parallel pipeline: correlation-based label transfer, depth prediction, and plotting.

Full pipeline for all Xenium samples using multiprocessing:
  1. Load raw Xenium data
  2. Correlation-based label transfer from SEA-AD snRNAseq reference
  3. Predict normalized cortical depth (MERFISH-trained model)
  4. Generate 2x3 summary figure per sample
  5. Save annotated h5ad per sample
  6. Merge all into combined h5ad

Usage:
    python run_parallel_corr.py \
        --data_dir /path/to/xenium_data \
        --reference /path/to/Reference_MTG_RNAseq.h5ad \
        --depth_model /path/to/depth_model_normalized.pkl \
        --output_dir /path/to/output \
        --n_workers 4 --save_combined
"""

import os
import sys
import time
import argparse
import traceback
import numpy as np
import anndata as ad
from multiprocessing import Pool, current_process

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from loading import load_xenium_sample, discover_samples
from label_transfer_corr import (
    load_reference, correlation_label_transfer,
    get_seaad_colors, TAXONOMY_LEVELS
)
from depth_model import load_model, predict_depth
from plotting import plot_summary


# Globals (loaded once per worker via initializer)
_ref = None
_colors = None
_depth_model = None
_output_dir = None


def _init_worker(reference_path, depth_model_path, output_dir):
    """Initialize worker: load reference and depth model."""
    global _ref, _colors, _depth_model, _output_dir
    _output_dir = output_dir
    pid = current_process().pid

    print(f"  [Worker {pid}] Loading reference...")
    _ref = load_reference(reference_path)
    _colors = get_seaad_colors(_ref)

    if depth_model_path and os.path.exists(depth_model_path):
        print(f"  [Worker {pid}] Loading depth model...")
        _depth_model = load_model(depth_model_path)
    else:
        _depth_model = None
        print(f"  [Worker {pid}] No depth model, skipping depth prediction.")

    print(f"  [Worker {pid}] Ready.")


def _process_one_sample(sample_info):
    """Process a single sample: load -> annotate -> depth -> plot -> save."""
    global _ref, _colors, _depth_model, _output_dir
    sid = sample_info["sample_id"]
    pid = current_process().pid
    t0 = time.time()

    try:
        # Load
        adata = load_xenium_sample(
            sample_info["h5_path"], sample_info["boundaries_path"]
        )
        n_cells = adata.shape[0]
        print(f"  [{pid}] {sid}: {n_cells:,} cells, annotating (correlation)...")

        # Correlation-based label transfer
        annotated = correlation_label_transfer(adata, _ref)

        # Depth prediction
        if _depth_model is not None:
            print(f"  [{pid}] {sid}: Predicting depth...")
            pred_depth = predict_depth(annotated, _depth_model)
            annotated.obs['predicted_norm_depth'] = pred_depth
            depth_range = f"[{pred_depth.min():.3f}, {pred_depth.max():.3f}]"
        else:
            depth_range = "N/A"

        # Save h5ad
        h5ad_dir = os.path.join(_output_dir, "h5ad")
        os.makedirs(h5ad_dir, exist_ok=True)
        h5ad_path = os.path.join(h5ad_dir, f"{sid}_annotated.h5ad")
        annotated.write_h5ad(h5ad_path)

        # Summary figure
        plot_dir = os.path.join(_output_dir, "plots")
        os.makedirs(plot_dir, exist_ok=True)
        plot_path = os.path.join(plot_dir, f"{sid}_summary.png")
        plot_summary(annotated, plot_path, model_bundle=_depth_model)

        elapsed = time.time() - t0
        print(f"  [{pid}] {sid}: Done in {elapsed:.0f}s "
              f"(depth {depth_range})")
        return {
            "sample_id": sid, "status": "success",
            "n_cells": n_cells, "time": elapsed,
            "depth_range": depth_range,
        }

    except Exception as e:
        elapsed = time.time() - t0
        print(f"  [{pid}] {sid}: FAILED after {elapsed:.0f}s - {e}")
        traceback.print_exc()
        return {
            "sample_id": sid, "status": "error",
            "error": str(e), "time": elapsed,
        }


def main():
    parser = argparse.ArgumentParser(
        description="Parallel Xenium pipeline: correlation label transfer + depth"
    )
    parser.add_argument("--data_dir", required=True,
                        help="Directory with Xenium .h5 and boundary files")
    parser.add_argument("--reference", required=True,
                        help="SEA-AD full snRNAseq reference h5ad path")
    parser.add_argument("--depth_model", default=None,
                        help="Path to trained depth model .pkl")
    parser.add_argument("--output_dir", default=None,
                        help="Output directory (default: data_dir/output)")
    parser.add_argument("--n_workers", type=int, default=4,
                        help="Number of parallel workers")
    parser.add_argument("--samples", nargs="*", default=None,
                        help="Specific sample IDs (default: all)")
    parser.add_argument("--save_combined", action="store_true",
                        help="Merge all samples into combined h5ad")
    args = parser.parse_args()

    if args.output_dir is None:
        args.output_dir = os.path.join(args.data_dir, "output")
    os.makedirs(args.output_dir, exist_ok=True)

    # Discover samples
    all_samples = discover_samples(args.data_dir)
    if args.samples:
        all_samples = [s for s in all_samples
                       if s["sample_id"] in args.samples]
    print(f"Found {len(all_samples)} samples, using {args.n_workers} workers")
    print(f"Using CORRELATION-based label transfer (scrattch.mapping style)")

    t_start = time.time()

    # Process in parallel
    with Pool(
        processes=args.n_workers,
        initializer=_init_worker,
        initargs=(args.reference, args.depth_model, args.output_dir),
    ) as pool:
        results = pool.map(_process_one_sample, all_samples)

    # Summary
    total_time = time.time() - t_start
    print(f"\n{'='*60}")
    print(f"SUMMARY - {total_time:.0f}s total")
    print(f"{'='*60}")
    n_ok = sum(1 for r in results if r["status"] == "success")
    total_cells = sum(r.get("n_cells", 0) for r in results)
    print(f"  Success: {n_ok}/{len(results)}")
    print(f"  Total cells: {total_cells:,}")
    for r in results:
        status = "✓" if r["status"] == "success" else "✗"
        time_s = r.get("time", 0)
        n = r.get("n_cells", 0)
        depth = r.get("depth_range", "N/A")
        print(f"  {status} {r['sample_id']:10s}: {n:>6,} cells, {time_s:>4.0f}s, depth {depth}")
    for r in results:
        if r["status"] == "error":
            print(f"  FAILED: {r['sample_id']}: {r['error']}")

    # Optionally merge
    if args.save_combined:
        print("\nMerging all samples into combined h5ad...")
        t_merge = time.time()
        h5ad_dir = os.path.join(args.output_dir, "h5ad")
        adatas = []
        for r in results:
            if r["status"] == "success":
                path = os.path.join(h5ad_dir, f"{r['sample_id']}_annotated.h5ad")
                adatas.append(ad.read_h5ad(path))
        combined = ad.concat(adatas, join='outer')
        combined_path = os.path.join(args.output_dir, "all_samples_annotated.h5ad")
        combined.write_h5ad(combined_path)
        file_size = os.path.getsize(combined_path) / 1e9
        print(f"  Saved: {combined_path} ({file_size:.2f} GB)")
        print(f"  Shape: {combined.shape}")
        print(f"  Merge took {time.time()-t_merge:.0f}s")

    print(f"\nOutput: {args.output_dir}")


if __name__ == "__main__":
    main()
