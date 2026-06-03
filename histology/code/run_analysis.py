#!/usr/bin/env python3
"""Main analysis script: reproduces the full binary-layer cell density analysis.

Usage:
    python -m code.run_analysis

Outputs:
    results/fig_final_matched_exclusions.png
    results/mixed_model_stats_final.csv
"""

import sys
from pathlib import Path

# Ensure project root is on path
PROJECT_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_DIR))

from code import config
from code.data_loading import build_analysis_df, load_dwight_summary, load_cell_counts
from code.models import run_binary_layer_analysis_clean
from code.plotting import plot_diagnosis_barplots, plot_validation_scatter

import numpy as np
import pandas as pd


def compute_subject_means(counts_df):
    """Compute per-subject mean cell count per site (matching Dwight's raw count format).

    NOTE: This returns raw counts (not density) for validation against Dwight's summary.
    """
    # SST/VIP from R sections
    r_data = counts_df[counts_df["section"] == "R"]
    sst = r_data.groupby("subject")["ch488"].mean().rename("SST")
    vip = r_data.groupby("subject")["ch568"].mean().rename("VIP")

    # PV/PYR from L sections, split by layer
    l_data = counts_df[counts_df["section"] == "L"]
    pv = l_data.groupby("subject")["ch568"].mean().rename("PV")
    pyr_23 = l_data[l_data["layer"] == "L23"].groupby("subject")["ch488"].mean().rename("PYR_23")
    pyr_56 = l_data[l_data["layer"] == "L56"].groupby("subject")["ch488"].mean().rename("PYR_56")

    totals = pd.concat([sst, vip, pv, pyr_23, pyr_56], axis=1)
    totals.index.name = "subject"
    return totals.reset_index()


def main():
    print("=" * 60)
    print("Dwight Densiometry Analysis — Binary Layer Mixed Models")
    print("=" * 60)

    # ── 1. Load data ──────────────────────────────────────────────
    print("\n--- Loading data ---")
    config.RESULTS_DIR.mkdir(exist_ok=True)
    df_R, df_L = build_analysis_df()

    # Print diagnosis breakdown
    for label, df in [("R-section (SST/VIP)", df_R), ("L-section (PV/PYR)", df_L)]:
        diag_counts = df.groupby("diagnosis")["subject"].nunique()
        print(f"\n{label} diagnosis breakdown:")
        for d in config.DIAG_ORDER:
            n = diag_counts.get(d, 0)
            print(f"  {d}: {n}")

    # ── 2. Run mixed effects models ───────────────────────────────
    print("\n--- Running mixed effects models ---")
    print("Formula: cell_count ~ diagnosis * layer + age + sex + (1|subject)")
    stats_df = run_binary_layer_analysis_clean(df_R, df_L)

    # Save stats
    stats_path = config.RESULTS_DIR / "mixed_model_stats_final.csv"
    stats_df.to_csv(stats_path, index=False)
    print(f"\nSaved statistics: {stats_path}")

    # Print key findings
    print("\n--- Key findings ---")
    for cell_type in ["SST", "VIP", "PV", "PYR"]:
        ct_stats = stats_df[stats_df["cell_type"] == cell_type]
        sig = ct_stats[ct_stats["p"] < 0.05]
        if len(sig) > 0:
            print(f"\n{cell_type}:")
            for _, row in sig.iterrows():
                if row["term"] == "Intercept":
                    continue
                print(f"  {row['term']}: coef={row['coef']:.3f}, p={row['p']:.4f}")

    # ── 3. Generate main figure ───────────────────────────────────
    print("\n--- Generating figures ---")
    fig_path = config.RESULTS_DIR / "fig_final_matched_exclusions.png"
    plot_diagnosis_barplots(df_R, df_L, stats_df, fig_path)

    # ── 4. Validation against Dwight's summary ────────────────────
    print("\n--- Validating against Dwight's summary ---")
    counts_raw = load_cell_counts()
    our_means = compute_subject_means(counts_raw)
    dwight_totals = load_dwight_summary()

    merged = our_means.merge(dwight_totals, on="subject", suffixes=("_ours", "_dwight"))
    for ct in ["SST", "VIP", "PV", "PYR_23", "PYR_56"]:
        x, y = merged[f"{ct}_ours"], merged[f"{ct}_dwight"]
        mask = x.notna() & y.notna()
        r = np.corrcoef(x[mask], y[mask])[0, 1]
        print(f"  {ct}: r = {r:.4f} (n={mask.sum()})")

    val_path = config.RESULTS_DIR / "fig_validation_vs_dwight.png"
    plot_validation_scatter(our_means, dwight_totals, val_path)

    print("\n" + "=" * 60)
    print("Analysis complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
