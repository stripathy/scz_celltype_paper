#!/usr/bin/env python3
"""
Build crumblr input CSVs with the L6b margin filter applied.

Strategy: Drop L6b cells with corr_subclass_margin < 0.02.
All other cell types unaffected.

Outputs to output/crumblr/ with suffix _l6b_margin for R analysis.
"""

import os
import sys
import time
import glob
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    BASE_DIR, H5AD_DIR, CRUMBLR_DIR, EXCLUDE_SAMPLES,
    load_sample_adata,
)

sys.path.insert(0, os.path.join(BASE_DIR, "code"))
from modules.metadata import get_subject_info

METADATA_PATH = os.path.join(BASE_DIR, "data", "sample_metadata.xlsx")
MARGIN_THRESH = 0.02
SUFFIX = "_l6b_margin"


def main():
    t0 = time.time()
    os.makedirs(CRUMBLR_DIR, exist_ok=True)

    # Load donor metadata
    print("Loading donor metadata...")
    meta = get_subject_info(METADATA_PATH)
    meta = meta.set_index("sample_id")
    print(f"  {len(meta)} Xenium subjects")

    # Discover samples
    h5ad_files = sorted(glob.glob(os.path.join(H5AD_DIR, "*_annotated.h5ad")))
    print(f"Found {len(h5ad_files)} h5ad files")

    all_obs = []
    for fpath in h5ad_files:
        sid = os.path.basename(fpath).replace("_annotated.h5ad", "")
        if sid in EXCLUDE_SAMPLES:
            print(f"  Skipping {sid} (excluded)")
            continue

        adata = load_sample_adata(sid)
        obs = adata.obs.copy()

        # Apply base QC (same as current pipeline)
        mask = obs['corr_qc_pass'].values.astype(bool)
        obs = obs[mask].copy()

        # Filter to cortical
        obs = obs[obs['spatial_domain'] == 'Cortical'].copy()

        # Count L6b before/after
        n_l6b_before = (obs['corr_subclass'] == 'L6b').sum()

        # Apply L6b margin filter: drop L6b cells with margin < threshold
        l6b_low_margin = (
            (obs['corr_subclass'] == 'L6b') &
            (obs['corr_subclass_margin'] < MARGIN_THRESH)
        )
        obs = obs[~l6b_low_margin].copy()

        n_l6b_after = (obs['corr_subclass'] == 'L6b').sum()
        n_dropped = n_l6b_before - n_l6b_after

        print(f"  {sid}: {len(obs):,} cortical cells "
              f"(L6b: {n_l6b_before} -> {n_l6b_after}, dropped {n_dropped})")

        # Use corr_subclass for subclass level, corr_supertype for supertype
        obs['sample_id'] = sid
        all_obs.append(obs)

    df = pd.concat(all_obs, ignore_index=True)
    n_samples = df['sample_id'].nunique()
    print(f"\nTotal: {len(df):,} cortical cells from {n_samples} samples")

    # Build count tables at both levels
    for level_col, level_name in [("corr_subclass", "subclass"),
                                   ("corr_supertype", "supertype")]:
        # Build counts
        counts = df.groupby(["sample_id", level_col]).size().reset_index(name="count")
        totals = df.groupby("sample_id").size().reset_index(name="total")
        counts = counts.merge(totals, on="sample_id")
        counts = counts.rename(columns={"sample_id": "donor", level_col: "celltype"})

        # Attach metadata
        counts = counts.merge(
            meta[["diagnosis", "sex", "age"]].reset_index(),
            left_on="donor", right_on="sample_id", how="left"
        ).drop(columns=["sample_id"])

        counts = counts.sort_values(["donor", "celltype"]).reset_index(drop=True)

        n_types = counts["celltype"].nunique()
        outpath = os.path.join(CRUMBLR_DIR, f"crumblr_input_{level_name}{SUFFIX}.csv")
        counts.to_csv(outpath, index=False)
        print(f"  {level_name}: {n_samples} donors x {n_types} types -> {outpath}")

    # Also build baseline (no L6b filter) for fair comparison using same labels
    print("\nBuilding baseline (no L6b filter) with same label columns...")
    all_obs_base = []
    for fpath in h5ad_files:
        sid = os.path.basename(fpath).replace("_annotated.h5ad", "")
        if sid in EXCLUDE_SAMPLES:
            continue
        adata = load_sample_adata(sid)
        obs = adata.obs.copy()
        mask = obs['corr_qc_pass'].values.astype(bool)
        obs = obs[mask].copy()
        obs = obs[obs['spatial_domain'] == 'Cortical'].copy()
        obs['sample_id'] = sid
        all_obs_base.append(obs)

    df_base = pd.concat(all_obs_base, ignore_index=True)

    for level_col, level_name in [("corr_subclass", "subclass"),
                                   ("corr_supertype", "supertype")]:
        counts = df_base.groupby(["sample_id", level_col]).size().reset_index(name="count")
        totals = df_base.groupby("sample_id").size().reset_index(name="total")
        counts = counts.merge(totals, on="sample_id")
        counts = counts.rename(columns={"sample_id": "donor", level_col: "celltype"})
        counts = counts.merge(
            meta[["diagnosis", "sex", "age"]].reset_index(),
            left_on="donor", right_on="sample_id", how="left"
        ).drop(columns=["sample_id"])
        counts = counts.sort_values(["donor", "celltype"]).reset_index(drop=True)

        outpath = os.path.join(CRUMBLR_DIR, f"crumblr_input_{level_name}_baseline.csv")
        counts.to_csv(outpath, index=False)
        print(f"  {level_name} baseline: {outpath}")

    print(f"\nDone in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
