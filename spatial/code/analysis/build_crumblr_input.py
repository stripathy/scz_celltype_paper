#!/usr/bin/env python3
"""
Build crumblr input CSVs for SCZ vs Control compositional analysis.

Uses cortical-only cells (spatial_domain=='Cortical', via load_cells())
from all 24 Xenium samples. Generates:
  1. Whole-composition inputs (neurons + non-neurons together)
  2. Stratified inputs (neuronal-only and non-neuronal-only)

For stratified inputs, "total" is the total within each stratum
(e.g., neuronal total = all neuronal cells in that sample), so
proportions are within-class (e.g., Sst proportion out of all neurons).

Output:
  output/crumblr/crumblr_input_subclass.csv
  output/crumblr/crumblr_input_supertype.csv
  output/crumblr/crumblr_input_subclass_neuronal.csv
  output/crumblr/crumblr_input_supertype_neuronal.csv
  output/crumblr/crumblr_input_subclass_nonneuronal.csv
  output/crumblr/crumblr_input_supertype_nonneuronal.csv

Each CSV is long-format: donor, celltype, count, total, diagnosis, sex, age, pmi
"""

import os
import sys
import glob
import time
import argparse
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    BASE_DIR, H5AD_DIR, CRUMBLR_DIR, EXCLUDE_SAMPLES, load_cells,
    GABA_PREFIXES, GLUT_PREFIXES, NN_PREFIXES,
)

sys.path.insert(0, os.path.join(BASE_DIR, "code"))
from modules.metadata import get_subject_info

METADATA_PATH = os.path.join(BASE_DIR, "data", "sample_metadata.xlsx")

NEURONAL_PREFIXES = GABA_PREFIXES + GLUT_PREFIXES


def is_neuronal(celltype):
    """Check if a cell type is neuronal based on subclass/supertype name."""
    for p in NEURONAL_PREFIXES:
        if celltype.startswith(p):
            return True
    return False


def build_counts(obs_df, level_col):
    """Build long-format count table: one row per (sample, celltype)."""
    counts = obs_df.groupby(["sample_id", level_col]).size().reset_index(name="count")
    totals = obs_df.groupby("sample_id").size().reset_index(name="total")
    counts = counts.merge(totals, on="sample_id")
    counts = counts.rename(columns={"sample_id": "donor", level_col: "celltype"})
    return counts


def main():
    parser = argparse.ArgumentParser(description="Build crumblr input CSVs")
    parser.add_argument("--qc-mode", default="corr", choices=["corr", "hybrid"],
                        help="QC mode: 'corr' (default, spatial QC + margin filter) or 'hybrid' (nuclear doublet-resolved)")
    args = parser.parse_args()

    qc_mode = args.qc_mode
    # Output suffix: default (corr) writes to standard filenames; hybrid gets suffix
    suffix = "_hybrid" if qc_mode == "hybrid" else ""

    t0 = time.time()
    os.makedirs(CRUMBLR_DIR, exist_ok=True)

    print(f"QC mode: {qc_mode}")

    # Load donor metadata
    print("Loading donor metadata...")
    meta = get_subject_info(METADATA_PATH)
    meta = meta.set_index("sample_id")
    print(f"  {len(meta)} Xenium subjects in metadata")

    # Discover h5ad files
    h5ad_files = sorted(glob.glob(os.path.join(H5AD_DIR, "*_annotated.h5ad")))
    print(f"\nFound {len(h5ad_files)} h5ad files")

    # Load and concatenate
    all_obs = []
    for fpath in h5ad_files:
        sid = os.path.basename(fpath).replace("_annotated.h5ad", "")
        if sid in EXCLUDE_SAMPLES:
            print(f"  Skipping {sid} (excluded)")
            continue

        obs = load_cells(sid, cortical_only=True, qc_mode=qc_mode)
        print(f"  {sid}: {len(obs):,} cortical cells")
        all_obs.append(obs)

    df = pd.concat(all_obs, ignore_index=True)
    n_samples = df["sample_id"].nunique()
    print(f"\nTotal: {len(df):,} cortical cells from {n_samples} samples")

    # Classify cells as neuronal vs non-neuronal using subclass label
    df["is_neuronal"] = df["subclass_label"].apply(is_neuronal)
    n_neuronal = df["is_neuronal"].sum()
    n_nonneuronal = (~df["is_neuronal"]).sum()
    print(f"  Neuronal: {n_neuronal:,} ({100*n_neuronal/len(df):.1f}%)")
    print(f"  Non-neuronal: {n_nonneuronal:,} ({100*n_nonneuronal/len(df):.1f}%)")

    # Define strata: (stratum_name, filter_func)
    strata = [
        ("", None),                                    # whole composition
        ("_neuronal", lambda x: x["is_neuronal"]),      # neurons only
        ("_nonneuronal", lambda x: ~x["is_neuronal"]),  # non-neurons only
    ]

    # Build count tables at both levels, for each stratum
    for level_col, level_name in [("subclass_label", "subclass"),
                                   ("supertype_label", "supertype")]:
        for stratum_suffix, filter_fn in strata:
            if filter_fn is not None:
                sub_df = df[filter_fn(df)].copy()
                stratum_label = stratum_suffix.strip("_")
            else:
                sub_df = df
                stratum_label = "whole"

            counts = build_counts(sub_df, level_col)

            # Attach metadata
            counts = counts.merge(
                meta[["diagnosis", "sex", "age", "pmi"]].reset_index(),
                left_on="donor", right_on="sample_id", how="left"
            ).drop(columns=["sample_id"])

            # Sort for readability
            counts = counts.sort_values(["donor", "celltype"]).reset_index(drop=True)

            n_types = counts["celltype"].nunique()
            outpath = os.path.join(CRUMBLR_DIR,
                                   f"crumblr_input_{level_name}{stratum_suffix}{suffix}.csv")
            counts.to_csv(outpath, index=False)
            print(f"\n  {level_name} ({stratum_label}): {n_samples} donors x "
                  f"{n_types} types -> {os.path.basename(outpath)}")

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.0f}s")


if __name__ == "__main__":
    main()
