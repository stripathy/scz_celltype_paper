#!/usr/bin/env python3
"""
Build input data for depth-stratified proportion and density analyses.

Loads all 24 Xenium samples and produces two output CSVs:

1. Cell-level table (for GAM continuous-depth models):
   One row per cell with: sample_id, subclass_label, predicted_norm_depth,
   x, y, diagnosis, sex, age, pmi

2. Layer-level count table (for per-layer mixed models):
   One row per (sample × layer × subclass) with: count, total, proportion,
   area_um2, density (cells/mm²), plus metadata covariates.

Both tables use cortical-only QC-pass cells (consistent with crumblr).

Output:
  output/depth_proportions/cell_level_data.csv
  output/depth_proportions/layer_counts.csv

Usage:
    python3 -u build_depth_proportion_input.py [--qc-mode corr|hybrid]
"""

import os
import sys
import glob
import time
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import BASE_DIR, H5AD_DIR, EXCLUDE_SAMPLES, load_cells

sys.path.insert(0, os.path.join(BASE_DIR, "code"))
from modules.metadata import get_subject_info
from modules.depth_model import LAYER_BINS

METADATA_PATH = os.path.join(BASE_DIR, "data", "sample_metadata.xlsx")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "depth_proportions")

# Cortical layers only (exclude WM)
CORTICAL_LAYERS = ["L1", "L2/3", "L4", "L5", "L6"]


def estimate_layer_area(obs_layer, min_cells=20):
    """Estimate spatial area of a layer using convex hull of cell positions.

    Returns area in µm² (Xenium coords are in µm).
    Returns NaN if too few cells for reliable hull estimation.
    """
    if len(obs_layer) < min_cells:
        return np.nan
    coords = obs_layer[["x", "y"]].values
    try:
        hull = ConvexHull(coords)
        return hull.volume  # in 2D, .volume gives area
    except Exception:
        return np.nan


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Build depth proportion input data")
    parser.add_argument("--qc-mode", default="corr", choices=["corr", "hybrid"])
    args = parser.parse_args()

    t0 = time.time()
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load metadata
    print("Loading donor metadata...")
    meta = get_subject_info(METADATA_PATH).set_index("sample_id")
    print(f"  {len(meta)} subjects")

    # Discover samples
    h5ad_files = sorted(glob.glob(os.path.join(H5AD_DIR, "*_annotated.h5ad")))
    print(f"Found {len(h5ad_files)} h5ad files\n")

    all_cells = []
    all_layer_counts = []

    for fpath in h5ad_files:
        sid = os.path.basename(fpath).replace("_annotated.h5ad", "")
        if sid in EXCLUDE_SAMPLES:
            print(f"  Skipping {sid} (excluded)")
            continue

        # Load with depth and spatial coords
        obs = load_cells(sid, cortical_only=True, qc_mode=args.qc_mode,
                         extra_obs_columns=["predicted_norm_depth"])

        n_total = len(obs)
        n_with_depth = obs["predicted_norm_depth"].notna().sum()
        print(f"  {sid}: {n_total:,} cortical cells ({n_with_depth:,} with depth)")

        # ── Cell-level data ──
        cell_df = obs[["sample_id", "subclass_label", "supertype_label",
                        "predicted_norm_depth", "layer", "x", "y"]].copy()
        all_cells.append(cell_df)

        # ── Layer-level counts ──
        for layer in CORTICAL_LAYERS:
            layer_mask = obs["layer"] == layer
            obs_layer = obs[layer_mask]
            n_layer = len(obs_layer)

            if n_layer == 0:
                continue

            # Area estimate for density
            area_um2 = estimate_layer_area(obs_layer)
            area_mm2 = area_um2 / 1e6 if not np.isnan(area_um2) else np.nan

            # Subclass counts within this layer
            sc_counts = obs_layer.groupby("subclass_label").size().reset_index(name="count")
            sc_counts["total"] = n_layer
            sc_counts["proportion"] = sc_counts["count"] / sc_counts["total"]
            sc_counts["area_um2"] = area_um2
            sc_counts["density_per_mm2"] = sc_counts["count"] / area_mm2 if not np.isnan(area_mm2) else np.nan
            sc_counts["sample_id"] = sid
            sc_counts["layer"] = layer

            # Depth bin midpoint for reference
            lo, hi = LAYER_BINS[layer]
            sc_counts["depth_midpoint"] = (max(lo, 0) + min(hi, 1)) / 2

            all_layer_counts.append(sc_counts)

    # ── Assemble cell-level table ──
    print("\nAssembling cell-level table...")
    cells_df = pd.concat(all_cells, ignore_index=True)

    # Attach metadata
    cells_df = cells_df.merge(
        meta[["diagnosis", "sex", "age", "pmi"]].reset_index(),
        on="sample_id", how="left"
    )
    print(f"  {len(cells_df):,} cells from {cells_df['sample_id'].nunique()} samples")
    print(f"  Subclasses: {cells_df['subclass_label'].nunique()}")
    print(f"  Depth range: {cells_df['predicted_norm_depth'].min():.3f} - "
          f"{cells_df['predicted_norm_depth'].max():.3f}")

    cell_path = os.path.join(OUTPUT_DIR, "cell_level_data.csv")
    cells_df.to_csv(cell_path, index=False)
    print(f"  Saved: {cell_path}")

    # ── Assemble layer-level count table ──
    print("\nAssembling layer-level count table...")
    layer_df = pd.concat(all_layer_counts, ignore_index=True)

    # Attach metadata
    layer_df = layer_df.merge(
        meta[["diagnosis", "sex", "age", "pmi"]].reset_index(),
        on="sample_id", how="left"
    )

    # Sort for readability
    layer_df = layer_df.sort_values(
        ["sample_id", "layer", "subclass_label"]
    ).reset_index(drop=True)

    n_combos = len(layer_df)
    n_types = layer_df["subclass_label"].nunique()
    print(f"  {n_combos:,} rows ({layer_df['sample_id'].nunique()} samples × "
          f"{layer_df['layer'].nunique()} layers × {n_types} subclasses)")

    layer_path = os.path.join(OUTPUT_DIR, "layer_counts.csv")
    layer_df.to_csv(layer_path, index=False)
    print(f"  Saved: {layer_path}")

    # ── Summary statistics ──
    print("\n── Layer cell counts (median across samples) ──")
    layer_totals = layer_df.groupby(["sample_id", "layer"])["total"].first()
    for layer in CORTICAL_LAYERS:
        vals = layer_totals.xs(layer, level="layer") if layer in layer_totals.index.get_level_values("layer") else pd.Series()
        if len(vals) > 0:
            print(f"  {layer:5s}: median {vals.median():,.0f} cells  "
                  f"(range {vals.min():,.0f} - {vals.max():,.0f})")

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.0f}s")


if __name__ == "__main__":
    main()
