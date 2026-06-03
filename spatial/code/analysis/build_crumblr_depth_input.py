#!/usr/bin/env python3
"""
Build crumblr input for depth x diagnosis interaction analysis.

Bins cells by density-adaptive (quantile-based) depth bins — each bin
contains approximately equal numbers of cells across the full cohort.
This gives more resolution in cell-dense regions (L2/3, L5) and less
in sparse regions (L1, deep L6).

Produces a long-format count table where each row is one
(sample x depth_bin x celltype) observation. This lets crumblr model:

    CLR(proportion) ~ depth * diagnosis + sex + scale(age) + scale(pmi) + (1|donor)

The depth term is the bin midpoint (continuous, 0=pia, 1=WM), so the
interaction tests whether the compositional depth profile differs by diagnosis.

Output:
  output/crumblr/crumblr_depth_input_subclass.csv

Each row: donor, depth_bin, depth_midpoint, celltype, count, total,
          diagnosis, sex, age, pmi
"""

import os
import sys
import glob
import time
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import BASE_DIR, H5AD_DIR, CRUMBLR_DIR, EXCLUDE_SAMPLES, load_cells

sys.path.insert(0, os.path.join(BASE_DIR, "code"))
from modules.metadata import get_subject_info

METADATA_PATH = os.path.join(BASE_DIR, "data", "sample_metadata.xlsx")

# 15 density-adaptive depth bins: edges set so each bin contains
# approximately the same number of cells across the full cohort.
# This gives more resolution in cell-dense regions (L2/3, L5) and
# less in sparse regions (L1, deep L6).
N_BINS = 15
# Edges will be computed from data at runtime (see main()).
# Initialized as None; set in main() before use.
DEPTH_EDGES = None
DEPTH_MIDPOINTS = None


def compute_quantile_edges(all_depths, n_bins=15):
    """Compute depth bin edges so each bin has equal cell count.

    Parameters
    ----------
    all_depths : array-like
        All predicted_norm_depth values (pooled across samples, clipped to [0,1]).
    n_bins : int
        Number of bins.

    Returns
    -------
    edges : np.ndarray of shape (n_bins + 1,)
    midpoints : np.ndarray of shape (n_bins,)
    """
    quantiles = np.linspace(0, 1, n_bins + 1)
    edges = np.quantile(all_depths, quantiles)
    # Force exact boundaries at 0 and 1
    edges[0] = 0.0
    edges[-1] = 1.0
    midpoints = (edges[:-1] + edges[1:]) / 2
    return edges, midpoints


def estimate_bin_area(obs_bin, min_cells=20):
    """Estimate spatial area of a depth bin using convex hull of cell positions.

    Returns area in µm². Returns NaN if too few cells.
    """
    if len(obs_bin) < min_cells:
        return np.nan
    coords = obs_bin[["x", "y"]].values
    try:
        hull = ConvexHull(coords)
        return hull.volume  # in 2D, .volume gives area
    except Exception:
        return np.nan


def main():
    global DEPTH_EDGES, DEPTH_MIDPOINTS

    t0 = time.time()
    os.makedirs(CRUMBLR_DIR, exist_ok=True)

    # Load metadata
    print("Loading donor metadata...")
    meta = get_subject_info(METADATA_PATH).set_index("sample_id")
    print(f"  {len(meta)} subjects")

    # Discover samples
    h5ad_files = sorted(glob.glob(os.path.join(H5AD_DIR, "*_annotated.h5ad")))
    sample_ids = []
    for fpath in h5ad_files:
        sid = os.path.basename(fpath).replace("_annotated.h5ad", "")
        if sid not in EXCLUDE_SAMPLES:
            sample_ids.append(sid)
    print(f"Found {len(h5ad_files)} h5ad files ({len(sample_ids)} after exclusions)\n")

    # ── Pass 1: collect all depths to compute quantile bin edges ──
    print("Pass 1: collecting depths for quantile bin edges...")
    all_depths = []
    for sid in sample_ids:
        obs = load_cells(sid, cortical_only=True, qc_mode="corr",
                         extra_obs_columns=["predicted_norm_depth"])
        depths = obs["predicted_norm_depth"].dropna().clip(0.0, 1.0).values
        all_depths.append(depths)
        print(f"  {sid}: {len(depths):,} cells with depth")
    all_depths = np.concatenate(all_depths)
    print(f"  Total: {len(all_depths):,} cells pooled")

    DEPTH_EDGES, DEPTH_MIDPOINTS = compute_quantile_edges(all_depths, N_BINS)
    print(f"  Quantile edges: {np.round(DEPTH_EDGES, 4)}")
    print(f"  Bin midpoints:  {np.round(DEPTH_MIDPOINTS, 4)}")
    print(f"  ~{len(all_depths) // N_BINS:,} cells per bin\n")

    # ── Pass 2: bin cells and build count table ──
    print("Pass 2: binning cells and building count table...")
    all_rows = []

    for sid in sample_ids:
        obs = load_cells(sid, cortical_only=True, qc_mode="corr",
                         extra_obs_columns=["predicted_norm_depth"])

        n_total = len(obs)
        n_depth = obs["predicted_norm_depth"].notna().sum()

        # Drop cells without depth prediction
        obs = obs[obs["predicted_norm_depth"].notna()].copy()

        # Assign depth bins (clip to 0-1 range for binning)
        depth_clipped = obs["predicted_norm_depth"].clip(0.0, 1.0 - 1e-9)
        obs["depth_bin"] = np.digitize(depth_clipped, DEPTH_EDGES) - 1
        obs["depth_bin"] = obs["depth_bin"].clip(0, N_BINS - 1)
        obs["depth_midpoint"] = DEPTH_MIDPOINTS[obs["depth_bin"]]

        print(f"  {sid}: {n_total:,} cortical ({n_depth:,} with depth) -> "
              f"{len(obs):,} binned")

        # Count cells per (sample x depth_bin x subclass)
        counts = (obs.groupby(["sample_id", "depth_bin", "subclass_label"],
                               observed=True)
                  .size()
                  .reset_index(name="count"))

        # Add depth_midpoint (deterministic from depth_bin, not a groupby key)
        counts["depth_midpoint"] = DEPTH_MIDPOINTS[counts["depth_bin"].values]

        # Total cells per (sample x depth_bin) for composition denominator
        totals = (obs.groupby(["sample_id", "depth_bin"], observed=True)
                  .size()
                  .reset_index(name="total"))
        counts = counts.merge(totals, on=["sample_id", "depth_bin"])

        # Estimate area per depth bin (convex hull of cell positions)
        bin_areas = {}
        for b in range(N_BINS):
            obs_bin = obs[obs["depth_bin"] == b]
            bin_areas[b] = estimate_bin_area(obs_bin)
        counts["area_um2"] = counts["depth_bin"].map(bin_areas)
        counts["area_mm2"] = counts["area_um2"] / 1e6
        counts["density_per_mm2"] = counts["count"] / counts["area_mm2"]

        counts = counts.rename(columns={
            "sample_id": "donor",
            "subclass_label": "celltype",
        })
        all_rows.append(counts)

    # Assemble
    df = pd.concat(all_rows, ignore_index=True)

    # Attach metadata
    df = df.merge(
        meta[["diagnosis", "sex", "age", "pmi"]].reset_index(),
        left_on="donor", right_on="sample_id", how="left"
    ).drop(columns=["sample_id"])

    # Sort
    df = df.sort_values(["donor", "depth_bin", "celltype"]).reset_index(drop=True)

    n_donors = df["donor"].nunique()
    n_types = df["celltype"].nunique()
    n_bins = df["depth_bin"].nunique()
    print(f"\nTotal: {len(df):,} rows "
          f"({n_donors} donors x {n_bins} depth bins x {n_types} subclasses)")
    print(f"Depth bins: {n_bins} (edges: {DEPTH_EDGES})")

    # Save
    outpath = os.path.join(CRUMBLR_DIR, "crumblr_depth_input_subclass.csv")
    df.to_csv(outpath, index=False)
    print(f"\nSaved: {outpath}")

    # Quick summary: cells per bin
    print("\n── Cells per depth bin (median across samples) ──")
    bin_totals = df.groupby(["donor", "depth_bin"])["total"].first()
    for b in range(N_BINS):
        vals = bin_totals.xs(b, level="depth_bin")
        print(f"  Bin {b} (depth {DEPTH_MIDPOINTS[b]:.2f}): "
              f"median {vals.median():,.0f} cells "
              f"(range {vals.min():,.0f} - {vals.max():,.0f})")

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.0f}s")


if __name__ == "__main__":
    main()
