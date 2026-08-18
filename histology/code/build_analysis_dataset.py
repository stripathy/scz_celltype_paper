#!/usr/bin/env python3
"""Rebuild the consolidated RNAscope analysis dataset from the raw count sheets.

This is the provenance for data/sst_analysis_data.csv, the single input to the
manuscript analysis (code/fit_sst_mixed_models.py). Pipeline:

    1. Per-site SST / VIP counts for the R sections (data/Cell_counts_NU.csv),
       subject-ID corrections and exclusions from code/config.py, diagnosis
       from data/full cell counts(Excel).xlsx, demographics (age, sex, PMI)
       from data/pTable with correct med info.csv  -> code/data_loading.py
    2. Layer assignment: sites 1-10 = L2/3, sites 11-20 = L5/6 (binary layers).
    3. VIP laminar quality filter: VIP interneurons are concentrated in the
       superficial layers, so for each subject we test VIP(L2/3) > VIP(L5/6)
       (one-sided two-sample t-test, sites with missing counts dropped).
       vip_pass = P < 0.05. Subjects failing the test have layer annotations
       that do not reflect laminar biology and are excluded from the primary
       analysis (54 of 68 subjects retained).
    4. Sanity check against the original study's per-subject summary counts
       (Dwight Newton's sheet in the same workbook).

Usage (from histology/):
    python3 code/build_analysis_dataset.py            # verify against the tracked CSV
    python3 code/build_analysis_dataset.py --write    # (re)write data/sst_analysis_data.csv
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from code.data_loading import build_analysis_df, load_demographics, load_dwight_summary  # noqa: E402

OUT = ROOT / "data" / "sst_analysis_data.csv"
LAYER_LABEL = {"L23": "L2/3", "L56": "L5/6"}
COLS = ["subject", "section", "site", "layer", "SST", "VIP", "diagnosis",
        "age", "sex", "sex_binary", "PMI", "vip_ttest_p", "vip_pass"]


def vip_filter(df_R: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for subj, g in df_R.groupby("subject"):
        a = g.loc[g["layer"] == "L23", "VIP"].dropna()
        b = g.loc[g["layer"] == "L56", "VIP"].dropna()
        p = stats.ttest_ind(a, b, alternative="greater").pvalue if len(a) > 1 and len(b) > 1 else np.nan
        rows.append({"subject": subj, "vip_ttest_p": p, "vip_pass": bool(p < 0.05)})
    return pd.DataFrame(rows)


def build() -> pd.DataFrame:
    df_R, _ = build_analysis_df()                       # merges age and sex
    pmi = load_demographics()[["subject", "PMI"]]
    df = (df_R.merge(pmi, on="subject", how="left")
              .merge(vip_filter(df_R), on="subject", how="left"))
    df["layer"] = df["layer"].map(LAYER_LABEL)
    df["subject"] = df["subject"].astype(int)
    df = df[COLS].sort_values(["subject", "section", "site"]).reset_index(drop=True)
    return df


def check_against_original_summary(df: pd.DataFrame) -> None:
    """Per-subject mean SST/VIP per site should match the original study's summary."""
    ours = df.groupby("subject")[["SST", "VIP"]].mean()
    orig = load_dwight_summary()
    orig["subject"] = orig["subject"].astype(int)
    orig = orig.set_index("subject")[["SST", "VIP"]]
    both = ours.join(orig, lsuffix="_ours", rsuffix="_orig", how="inner")
    for ct in ["SST", "VIP"]:
        r = np.corrcoef(both[f"{ct}_ours"], both[f"{ct}_orig"])[0, 1]
        mad = (both[f"{ct}_ours"] - both[f"{ct}_orig"]).abs().max()
        print(f"  {ct}: {len(both)} subjects vs original summary, r = {r:.4f}, max |diff| = {mad:.3f}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true", help="write data/sst_analysis_data.csv")
    args = ap.parse_args()

    df = build()
    n_by_dx = df.groupby("diagnosis")["subject"].nunique().to_dict()
    n_pass = df[df.vip_pass].groupby("diagnosis")["subject"].nunique().to_dict()
    print(f"subjects: {df.subject.nunique()} ({n_by_dx}); VIP-QC pass: "
          f"{df[df.vip_pass].subject.nunique()} ({n_pass}); site-rows: {len(df)}")
    print("validation against the original study's per-subject summary:")
    check_against_original_summary(df)

    if OUT.exists():
        ref = pd.read_csv(OUT)
        m = ref.merge(df, on=["subject", "section", "site"], suffixes=("_ref", ""))
        same_counts = np.allclose(m.SST_ref.fillna(-1), m.SST.fillna(-1)) and np.allclose(m.VIP_ref.fillna(-1), m.VIP.fillna(-1))
        same_p = np.allclose(m.vip_ttest_p_ref, m.vip_ttest_p, atol=1e-12)
        same_pass = (m.vip_pass_ref == m.vip_pass).all()
        print(f"vs tracked CSV: rows matched {len(m)}/{len(ref)}; counts identical: {same_counts}; "
              f"VIP P identical: {same_p}; vip_pass identical: {same_pass}")
    if args.write:
        df.to_csv(OUT, index=False)
        print(f"wrote {OUT.relative_to(ROOT)} ({len(df)} rows, {df.shape[1]} columns)")


if __name__ == "__main__":
    main()
