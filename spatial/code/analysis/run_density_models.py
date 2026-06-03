#!/usr/bin/env python3
"""
Supertype-level cell density analysis (cells per mm²).

Computes density for each supertype in each sample, then tests SCZ vs Control
using the same covariates as crumblr:

    log(density) ~ diagnosis + sex + scale(age) + scale(pmi)

This is NOT compositional — each supertype is modeled independently, so
changes in one type do not mechanically affect others (unlike crumblr/CLR).

Area is estimated per sample as the convex hull of cells in the region.

Runs two versions:
  1. "cortical" — all cells with spatial_domain == Cortical (matches crumblr)
  2. "L1-L6" — restricted to cells assigned to layers L1 through L6

Output:
  output/density_analysis/density_results_supertype_cortical.csv
  output/density_analysis/density_results_supertype_L1-L6.csv
  output/density_analysis/density_per_sample_supertype.csv  (raw data, both)

Usage:
    python3 -u run_density_models.py
"""

import os
import sys
import glob
import time
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import BASE_DIR, H5AD_DIR, EXCLUDE_SAMPLES, load_cells, infer_class

sys.path.insert(0, os.path.join(BASE_DIR, "code"))
from modules.metadata import get_subject_info

METADATA_PATH = os.path.join(BASE_DIR, "data", "sample_metadata.xlsx")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "density_analysis")
CORTICAL_LAYERS = {"L1", "L2/3", "L4", "L5", "L6"}


def estimate_area_mm2(obs):
    """Estimate area from convex hull of cell positions."""
    coords = obs[["x", "y"]].values
    if len(coords) < 10:
        return np.nan
    try:
        return ConvexHull(coords).volume / 1e6  # µm² -> mm²
    except Exception:
        return np.nan


def build_density_table(all_obs, meta, region_label, layer_filter=None):
    """Build per-sample, per-supertype density table.

    Parameters
    ----------
    all_obs : dict
        {sample_id: obs DataFrame}
    meta : DataFrame
        Donor metadata indexed by sample_id.
    region_label : str
        Label for this region version (e.g., "cortical", "L1-L6").
    layer_filter : set or None
        If provided, restrict to cells in these layers.
    """
    records = []

    for sid, obs in all_obs.items():
        if layer_filter is not None:
            obs = obs[obs["layer"].isin(layer_filter)]

        if len(obs) == 0:
            continue

        area_mm2 = estimate_area_mm2(obs)
        st_counts = obs["supertype_label"].value_counts()
        total_cells = len(obs)

        for st, count in st_counts.items():
            records.append({
                "sample_id": sid,
                "supertype": st,
                "region": region_label,
                "count": count,
                "total_cells": total_cells,
                "area_mm2": area_mm2,
                "density_per_mm2": count / area_mm2 if area_mm2 > 0 else np.nan,
                "proportion": count / total_cells,
            })

    df = pd.DataFrame(records)
    df = df.merge(
        meta[["diagnosis", "sex", "age", "pmi"]].reset_index(),
        on="sample_id", how="left"
    )
    return df


def fit_density_models(df, region_label):
    """Fit OLS density models per supertype, return results DataFrame."""
    n_samples = df["sample_id"].nunique()
    n_types = df["supertype"].nunique()

    # Filter to supertypes present in >= 50% of samples (matching crumblr)
    presence = df.groupby("supertype")["sample_id"].nunique()
    min_samples = n_samples * 0.5
    keep_types = presence[presence >= min_samples].index.tolist()
    print(f"  Testing {len(keep_types)} / {n_types} supertypes "
          f"(present in >= {min_samples:.0f} samples)")

    results = []

    for st in sorted(keep_types):
        st_data = df[df["supertype"] == st].copy()

        n_ctrl = st_data[st_data["diagnosis"] == "Control"]["sample_id"].nunique()
        n_scz = st_data[st_data["diagnosis"] == "SCZ"]["sample_id"].nunique()
        if n_ctrl < 3 or n_scz < 3:
            continue

        st_data["y"] = np.log1p(st_data["density_per_mm2"])
        st_data["age_z"] = (st_data["age"] - st_data["age"].mean()) / st_data["age"].std()
        st_data["pmi_z"] = (st_data["pmi"] - st_data["pmi"].mean()) / st_data["pmi"].std()

        try:
            formula = "y ~ C(diagnosis, Treatment('Control')) + C(sex) + age_z + pmi_z"
            fit = smf.ols(formula, data=st_data).fit()

            scz_key = [k for k in fit.params.index if "SCZ" in k]
            if scz_key:
                coef = fit.params[scz_key[0]]
                se = fit.bse[scz_key[0]]
                tstat = fit.tvalues[scz_key[0]]
                pval = fit.pvalues[scz_key[0]]
            else:
                coef = se = tstat = pval = np.nan

            ctrl_density = st_data[st_data["diagnosis"] == "Control"]["density_per_mm2"].mean()
            scz_density = st_data[st_data["diagnosis"] == "SCZ"]["density_per_mm2"].mean()
            logFC = np.log(scz_density / ctrl_density) if ctrl_density > 0 else np.nan

        except Exception as e:
            print(f"    {st}: failed ({e})")
            coef = se = tstat = pval = logFC = np.nan
            ctrl_density = scz_density = np.nan

        results.append({
            "supertype": st,
            "region": region_label,
            "n_samples": len(st_data),
            "n_ctrl": n_ctrl,
            "n_scz": n_scz,
            "ctrl_mean_density": ctrl_density,
            "scz_mean_density": scz_density,
            "logFC": logFC,
            "coef": coef,
            "se": se,
            "tstat": tstat,
            "pval": pval,
            "class": infer_class(st),
        })

    results_df = pd.DataFrame(results)

    # FDR correction
    valid = results_df["pval"].notna()
    if valid.any():
        _, fdr_vals, _, _ = multipletests(results_df.loc[valid, "pval"],
                                           method="fdr_bh")
        results_df.loc[valid, "fdr"] = fdr_vals

    results_df = results_df.sort_values("pval").reset_index(drop=True)
    return results_df


def print_top_hits(results_df, n=20):
    """Print top hits from results."""
    for _, row in results_df.head(n).iterrows():
        sig = "***" if row.get("fdr", 1) < 0.05 else (
              "**" if row.get("fdr", 1) < 0.10 else (
              "*" if row["pval"] < 0.05 else ""))
        print(f"    {row['supertype']:30s} logFC={row['logFC']:+.3f}  "
              f"coef={row['coef']:+.3f}  p={row['pval']:.4g}  "
              f"FDR={row.get('fdr', np.nan):.4g}  {sig}")

    n_fdr05 = (results_df.get("fdr", pd.Series()) < 0.05).sum()
    n_fdr10 = (results_df.get("fdr", pd.Series()) < 0.10).sum()
    n_nom = (results_df["pval"] < 0.05).sum()
    print(f"    => {n_fdr05} FDR<0.05, {n_fdr10} FDR<0.10, {n_nom} nom p<0.05 "
          f"(of {len(results_df)} tested)")


def main():
    t0 = time.time()
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load metadata
    print("Loading donor metadata...")
    meta = get_subject_info(METADATA_PATH).set_index("sample_id")
    print(f"  {len(meta)} subjects")

    # Load all samples once
    h5ad_files = sorted(glob.glob(os.path.join(H5AD_DIR, "*_annotated.h5ad")))
    print(f"Found {len(h5ad_files)} h5ad files\n")

    all_obs = {}
    for fpath in h5ad_files:
        sid = os.path.basename(fpath).replace("_annotated.h5ad", "")
        if sid in EXCLUDE_SAMPLES:
            print(f"  Skipping {sid} (excluded)")
            continue

        obs = load_cells(sid, cortical_only=True, qc_mode="corr")
        all_obs[sid] = obs

        n_l16 = obs["layer"].isin(CORTICAL_LAYERS).sum()
        print(f"  {sid}: {len(obs):,} cortical cells ({n_l16:,} in L1-L6)")

    # ── Build density tables for both versions ──
    print("\n" + "=" * 70)
    print("Building density tables...")
    print("=" * 70)

    df_cortical = build_density_table(all_obs, meta, "cortical")
    df_l16 = build_density_table(all_obs, meta, "L1-L6",
                                  layer_filter=CORTICAL_LAYERS)

    print(f"\n  Cortical: {len(df_cortical):,} rows, "
          f"{df_cortical['sample_id'].nunique()} samples, "
          f"{df_cortical['supertype'].nunique()} supertypes")
    print(f"  L1-L6:    {len(df_l16):,} rows, "
          f"{df_l16['sample_id'].nunique()} samples, "
          f"{df_l16['supertype'].nunique()} supertypes")

    # Save combined raw data
    raw_df = pd.concat([df_cortical, df_l16], ignore_index=True)
    raw_path = os.path.join(OUTPUT_DIR, "density_per_sample_supertype.csv")
    raw_df.to_csv(raw_path, index=False)
    print(f"\nSaved raw data: {raw_path}")

    # ── Fit models ──
    for label, df in [("cortical", df_cortical), ("L1-L6", df_l16)]:
        print(f"\n{'=' * 70}")
        print(f"Fitting density models: {label}")
        print(f"  log(density) ~ diagnosis + sex + scale(age) + scale(pmi)")
        print(f"{'=' * 70}")

        results = fit_density_models(df, label)

        out_path = os.path.join(OUTPUT_DIR,
                                 f"density_results_supertype_{label}.csv")
        results.to_csv(out_path, index=False)
        print(f"\n  Saved: {out_path}")

        print(f"\n  Top hits ({label}):")
        print_top_hits(results)

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.0f}s")


if __name__ == "__main__":
    main()
