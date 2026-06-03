"""Subject filtering funnel — tracks where subjects get included or dropped.

Builds a per-subject DataFrame documenting status at each pipeline step,
making it easy to audit how many subjects survive each filter and why
specific subjects were excluded.
"""

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from . import config


def build_subject_funnel(counts, diag_map, coords_df, vip_val=None):
    """Build a per-subject filtering status DataFrame.

    Parameters
    ----------
    counts : DataFrame
        Output of load_cell_counts(). Columns: subject, section, site, ...
    diag_map : dict
        {subject_id_str: diagnosis_str} from load_diagnosis_mapping()
    coords_df : DataFrame
        Output of compute_all_depths(). Columns: subject, section, site, depth, depth_raw
    vip_val : DataFrame, optional
        Output of validate_depth_with_vip(). Columns: subject, vip_depth_r, quality.
        If None, VIP validation columns are omitted.

    Returns
    -------
    DataFrame with one row per subject and boolean/status columns for each step.
    """
    # Start from all subjects in the diagnosis mapping
    all_subjects = sorted(diag_map.keys(), key=lambda x: int(x) if x.isdigit() else x)
    funnel = pd.DataFrame({"subject": all_subjects})
    funnel["diagnosis"] = funnel["subject"].map(diag_map)

    # Step 1: Present in cell counts?
    r_subjects = set(counts[counts["section"] == "R"]["subject"].unique())
    l_subjects = set(counts[counts["section"] == "L"]["subject"].unique())
    funnel["in_counts_R"] = funnel["subject"].isin(r_subjects)
    funnel["in_counts_L"] = funnel["subject"].isin(l_subjects)

    # Step 2: Not excluded?
    funnel["not_excluded_R"] = ~funnel["subject"].isin(config.EXCLUDE_R)
    funnel["not_excluded_L"] = ~funnel["subject"].isin(config.EXCLUDE_L)

    # Step 3: Has coordinates + valid depth?
    r_depth_subjects = set(
        coords_df[(coords_df["section"] == "R") & coords_df["depth"].notna()]["subject"].unique()
    )
    l_depth_subjects = set(
        coords_df[(coords_df["section"] == "L") & coords_df["depth"].notna()]["subject"].unique()
    )
    funnel["has_depth_R"] = funnel["subject"].isin(r_depth_subjects)
    funnel["has_depth_L"] = funnel["subject"].isin(l_depth_subjects)

    # Combined: in the "all subjects" analysis?
    funnel["in_final_all_R"] = (
        funnel["in_counts_R"] & funnel["not_excluded_R"] &
        funnel["has_depth_R"] & funnel["diagnosis"].notna()
    )
    funnel["in_final_all_L"] = (
        funnel["in_counts_L"] & funnel["not_excluded_L"] &
        funnel["has_depth_L"] & funnel["diagnosis"].notna()
    )

    # Step 4: VIP validation (R-section only)
    if vip_val is not None:
        vip_map_r = dict(zip(vip_val["subject"].astype(str), vip_val["vip_depth_r"]))
        vip_map_q = dict(zip(vip_val["subject"].astype(str), vip_val["quality"]))
        funnel["vip_r"] = funnel["subject"].map(vip_map_r)
        funnel["vip_quality"] = funnel["subject"].map(vip_map_q)
        funnel["in_final_strong_R"] = funnel["in_final_all_R"] & (funnel["vip_quality"] == "strong")
        # L-section constrained to R strong-VIP subjects
        strong_R_subjects = set(funnel[funnel["in_final_strong_R"]]["subject"])
        funnel["in_final_strong_L"] = funnel["in_final_all_L"] & funnel["subject"].isin(strong_R_subjects)
    else:
        funnel["vip_r"] = np.nan
        funnel["vip_quality"] = None
        funnel["in_final_strong_R"] = False
        funnel["in_final_strong_L"] = False

    # Assign drop reason (first reason encountered)
    reasons = []
    for _, row in funnel.iterrows():
        r_reasons = []
        if not row["in_counts_R"] and not row["in_counts_L"]:
            r_reasons.append("no cell counts")
        if row["in_counts_R"] and not row["not_excluded_R"]:
            r_reasons.append(f"excluded_R ({row['subject']} in EXCLUDE_R)")
        if row["in_counts_L"] and not row["not_excluded_L"]:
            r_reasons.append(f"excluded_L ({row['subject']} in EXCLUDE_L)")
        if row["in_counts_R"] and row["not_excluded_R"] and not row["has_depth_R"]:
            r_reasons.append("no coordinates/depth (R)")
        if row["in_counts_L"] and row["not_excluded_L"] and not row["has_depth_L"]:
            r_reasons.append("no coordinates/depth (L)")
        if row["in_final_all_R"] and not row["in_final_strong_R"]:
            q = row.get("vip_quality", "?")
            r_val = row.get("vip_r", np.nan)
            r_reasons.append(f"VIP validation {q} (r={r_val:.2f})" if pd.notna(r_val) else f"VIP validation {q}")

        if row["in_final_strong_R"] or row["in_final_strong_L"]:
            reasons.append("included (strong VIP)")
        elif row["in_final_all_R"] or row["in_final_all_L"]:
            reasons.append("included (all); " + "; ".join(r_reasons) if r_reasons else "included (all)")
        else:
            reasons.append("; ".join(r_reasons) if r_reasons else "no diagnosis")

    funnel["status"] = reasons
    return funnel


def print_funnel(funnel_df):
    """Pretty-print the subject filtering funnel to stdout."""
    n_total = len(funnel_df)

    # R-section funnel
    print("\nSUBJECT FILTERING FUNNEL")
    print("=" * 60)

    print(f"\n  R-section (SST / VIP)")
    print(f"  {'Diagnosis file:':<30s} {n_total:>4d} subjects")
    n_counts = funnel_df["in_counts_R"].sum()
    print(f"  {'In cell counts:':<30s} {n_counts:>4d}")
    n_not_excl = (funnel_df["in_counts_R"] & funnel_df["not_excluded_R"]).sum()
    print(f"  {'After exclusions:':<30s} {n_not_excl:>4d}")
    n_depth = funnel_df["in_final_all_R"].sum()
    print(f"  {'With valid depth:':<30s} {n_depth:>4d}")
    n_strong = funnel_df["in_final_strong_R"].sum()
    print(f"  {'Strong VIP (r < -0.3):':<30s} {n_strong:>4d}")

    # Diagnosis breakdown for all-subjects
    print(f"\n  Diagnosis breakdown (all R, N={n_depth}):")
    for d in config.DIAG_ORDER:
        nd = (funnel_df["in_final_all_R"] & (funnel_df["diagnosis"] == d)).sum()
        print(f"    {d:<12s} {nd:>3d}")

    # Diagnosis breakdown for strong-VIP
    print(f"\n  Diagnosis breakdown (strong VIP R, N={n_strong}):")
    for d in config.DIAG_ORDER:
        nd = (funnel_df["in_final_strong_R"] & (funnel_df["diagnosis"] == d)).sum()
        print(f"    {d:<12s} {nd:>3d}")

    # L-section funnel
    print(f"\n  L-section (PV / PYR)")
    print(f"  {'Diagnosis file:':<30s} {n_total:>4d} subjects")
    n_counts_l = funnel_df["in_counts_L"].sum()
    print(f"  {'In cell counts:':<30s} {n_counts_l:>4d}")
    n_not_excl_l = (funnel_df["in_counts_L"] & funnel_df["not_excluded_L"]).sum()
    print(f"  {'After exclusions:':<30s} {n_not_excl_l:>4d}")
    n_depth_l = funnel_df["in_final_all_L"].sum()
    print(f"  {'With valid depth:':<30s} {n_depth_l:>4d}")
    n_strong_l = funnel_df["in_final_strong_L"].sum()
    print(f"  {'Constrained to strong VIP R:':<30s} {n_strong_l:>4d}")

    print(f"\n  Diagnosis breakdown (all L, N={n_depth_l}):")
    for d in config.DIAG_ORDER:
        nd = (funnel_df["in_final_all_L"] & (funnel_df["diagnosis"] == d)).sum()
        print(f"    {d:<12s} {nd:>3d}")

    # Dropped subjects
    dropped_R = funnel_df[~funnel_df["in_final_all_R"] & funnel_df["in_counts_R"]]
    if len(dropped_R) > 0:
        print(f"\n  Dropped from R ({len(dropped_R)}):")
        for _, row in dropped_R.iterrows():
            print(f"    {row['subject']:>6s}  {row['diagnosis'] or 'no diagnosis':<12s}  {row['status']}")

    dropped_L = funnel_df[~funnel_df["in_final_all_L"] & funnel_df["in_counts_L"]]
    if len(dropped_L) > 0:
        print(f"\n  Dropped from L ({len(dropped_L)}):")
        for _, row in dropped_L.iterrows():
            print(f"    {row['subject']:>6s}  {row['diagnosis'] or 'no diagnosis':<12s}  {row['status']}")

    # Not-strong VIP subjects (in all but not in strong)
    not_strong = funnel_df[funnel_df["in_final_all_R"] & ~funnel_df["in_final_strong_R"]]
    if len(not_strong) > 0:
        print(f"\n  In all-R but not strong-VIP ({len(not_strong)}):")
        for _, row in not_strong.iterrows():
            r_val = f"r={row['vip_r']:.3f}" if pd.notna(row["vip_r"]) else "r=NaN"
            print(f"    {row['subject']:>6s}  {row['diagnosis']:<12s}  {row['vip_quality']:<12s}  {r_val}")

    print("=" * 60)
