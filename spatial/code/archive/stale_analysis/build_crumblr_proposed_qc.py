#!/usr/bin/env python3
"""
Build crumblr input CSVs with the PROPOSED final QC filter stack:

  1. qc_pass (spatial QC — neg probes, unassigned, low counts)
  2. 5th percentile margin filter (per-sample bottom 5% of corr_subclass_margin)
  3. doublet_suspect exclusion
  4. Cortical cells only (spatial_domain == 'Cortical')

Also builds a "current default" baseline with the same label columns
(corr_subclass / corr_supertype) for fair comparison.

Outputs to output/crumblr/ with suffix _proposed.
"""

import os
import sys
import time
import glob
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    BASE_DIR, H5AD_DIR, CRUMBLR_DIR, EXCLUDE_SAMPLES,
    load_sample_adata,
)

sys.path.insert(0, os.path.join(BASE_DIR, "code"))
from modules.metadata import get_subject_info

METADATA_PATH = os.path.join(BASE_DIR, "data", "sample_metadata.xlsx")
MARGIN_PCTL = 5  # 5th percentile


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

    all_proposed = []
    all_current = []

    for fpath in h5ad_files:
        sid = os.path.basename(fpath).replace("_annotated.h5ad", "")
        if sid in EXCLUDE_SAMPLES:
            print(f"  Skipping {sid} (excluded)")
            continue

        adata = load_sample_adata(sid)
        obs = adata.obs.copy()

        # ── Step 1: Spatial QC ──
        mask_qc = obs['qc_pass'].values.astype(bool)
        obs_qc = obs[mask_qc].copy()

        # ── Cortical only ──
        obs_cortical = obs_qc[obs_qc['spatial_domain'] == 'Cortical'].copy()

        # ── CURRENT DEFAULT: just qc_pass + cortical ──
        n_current = len(obs_cortical)
        obs_cortical_copy = obs_cortical.copy()
        obs_cortical_copy['sample_id'] = sid
        all_current.append(obs_cortical_copy)

        # ── Step 2: 5th percentile margin filter (per-sample) ──
        margins = obs_cortical['corr_subclass_margin'].values
        thresh = np.nanpercentile(margins, MARGIN_PCTL)
        low_margin = obs_cortical['corr_subclass_margin'] < thresh
        n_low_margin = low_margin.sum()

        # ── Step 3: Doublet suspect exclusion ──
        is_doublet = obs_cortical['doublet_suspect'].values.astype(bool)
        n_doublet = is_doublet.sum()

        # Combined filter: remove low margin OR doublet suspect
        remove = low_margin | is_doublet
        # Some cells may be both
        n_both = (low_margin & is_doublet).sum()
        obs_proposed = obs_cortical[~remove].copy()

        n_proposed = len(obs_proposed)
        n_dropped = n_current - n_proposed
        pct_dropped = 100 * n_dropped / n_current if n_current > 0 else 0

        print(f"  {sid}: {n_current:,} -> {n_proposed:,} cortical cells "
              f"(dropped {n_dropped}: {n_low_margin} low-margin, "
              f"{n_doublet} doublet, {n_both} both) "
              f"[margin thresh={thresh:.4f}, {pct_dropped:.1f}% dropped]")

        obs_proposed['sample_id'] = sid
        all_proposed.append(obs_proposed)

    # ── Concatenate ──
    df_proposed = pd.concat(all_proposed, ignore_index=True)
    df_current = pd.concat(all_current, ignore_index=True)

    n_samples = df_proposed['sample_id'].nunique()
    print(f"\nProposed: {len(df_proposed):,} cortical cells from {n_samples} samples")
    print(f"Current:  {len(df_current):,} cortical cells from {n_samples} samples")
    print(f"Dropped:  {len(df_current) - len(df_proposed):,} cells "
          f"({100*(len(df_current)-len(df_proposed))/len(df_current):.1f}%)")

    # ── Build count tables ──
    for tag, df in [("_proposed", df_proposed), ("_current_baseline", df_current)]:
        for level_col, level_name in [("corr_subclass", "subclass"),
                                       ("corr_supertype", "supertype")]:
            counts = df.groupby(["sample_id", level_col]).size().reset_index(name="count")
            totals = df.groupby("sample_id").size().reset_index(name="total")
            counts = counts.merge(totals, on="sample_id")
            counts = counts.rename(columns={"sample_id": "donor", level_col: "celltype"})

            counts = counts.merge(
                meta[["diagnosis", "sex", "age"]].reset_index(),
                left_on="donor", right_on="sample_id", how="left"
            ).drop(columns=["sample_id"])

            counts = counts.sort_values(["donor", "celltype"]).reset_index(drop=True)
            n_types = counts["celltype"].nunique()
            outpath = os.path.join(CRUMBLR_DIR, f"crumblr_input_{level_name}{tag}.csv")
            counts.to_csv(outpath, index=False)
            print(f"  {level_name}{tag}: {n_samples} donors x {n_types} types -> {outpath}")

    # ── Summary: what gets dropped by type ──
    print("\n── Per-subclass impact of proposed filters ──")
    for sc_name in sorted(df_current['corr_subclass'].unique()):
        n_cur = (df_current['corr_subclass'] == sc_name).sum()
        n_pro = (df_proposed['corr_subclass'] == sc_name).sum()
        n_drop = n_cur - n_pro
        pct = 100 * n_drop / n_cur if n_cur > 0 else 0
        flag = " ***" if pct > 10 else ""
        print(f"  {sc_name:25s}: {n_cur:7,} -> {n_pro:7,} (dropped {n_drop:5,}, {pct:5.1f}%){flag}")

    print(f"\nDone in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
