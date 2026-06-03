#!/usr/bin/env python3
"""
Run cell QC on all 24 Xenium samples.

For each sample:
  1. Load the raw h5 file to compute QC metrics (including control features)
  2. Flag QC failures using Kwon et al. approach
  3. Add QC columns to the existing annotated h5ad
  4. Save a per-sample summary

Outputs:
  - Updated h5ad files with QC columns in obs
  - output/qc_summary.csv — per-sample QC statistics
  - output/qc_validation.png — validation figures
"""

import os
import sys
import time
import re
import glob
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Pipeline config
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipeline_config import RAW_DIR, H5AD_DIR, OUTPUT_DIR, MODULES_DIR

# Modules
sys.path.insert(0, MODULES_DIR)
from cell_qc import compute_qc_metrics, flag_qc_failures


def find_raw_h5_for_sample(sample_id, base_dir):
    """Find the raw h5 file matching a sample ID (e.g., Br8667)."""
    pattern = os.path.join(base_dir, f"GSM*{sample_id}-cell_feature_matrix.h5")
    matches = glob.glob(pattern)
    if len(matches) == 1:
        return matches[0]
    elif len(matches) > 1:
        print(f"  WARNING: Multiple h5 files for {sample_id}, using first")
        return matches[0]
    return None


def main():
    t_start = time.time()

    # Find all annotated h5ad files
    h5ad_files = sorted(glob.glob(os.path.join(H5AD_DIR, "*_annotated.h5ad")))
    print(f"Found {len(h5ad_files)} annotated h5ad files")

    summary_rows = []

    for h5ad_path in h5ad_files:
        # Extract sample ID
        basename = os.path.basename(h5ad_path)
        sample_id = basename.replace("_annotated.h5ad", "")
        print(f"\n{'='*60}")
        print(f"Processing {sample_id}")
        print(f"{'='*60}")

        t0 = time.time()

        # Find raw h5
        raw_h5_path = find_raw_h5_for_sample(sample_id, RAW_DIR)
        if raw_h5_path is None:
            print(f"  ERROR: No raw h5 found for {sample_id}, skipping")
            continue

        print(f"  Raw h5: {os.path.basename(raw_h5_path)}")

        # Compute QC metrics from raw h5
        print(f"  Computing QC metrics...")
        metrics = compute_qc_metrics(raw_h5_path)
        print(f"    {len(metrics):,} cells in raw h5")

        # Flag QC failures
        print(f"  Flagging QC failures...")
        qc_df = flag_qc_failures(metrics)

        # Load existing annotated h5ad
        print(f"  Loading annotated h5ad...")
        adata = ad.read_h5ad(h5ad_path)
        print(f"    {adata.shape[0]:,} cells in h5ad")

        # Align QC metrics to h5ad cell barcodes
        # The h5ad barcodes should be a subset of (or same as) the raw h5 barcodes
        shared_barcodes = list(set(adata.obs_names) & set(qc_df.index))
        print(f"    {len(shared_barcodes):,} shared barcodes")

        if len(shared_barcodes) < adata.shape[0]:
            print(f"    WARNING: {adata.shape[0] - len(shared_barcodes)} h5ad cells "
                  f"not found in raw h5!")

        # Reindex QC metrics to match h5ad
        qc_aligned = qc_df.reindex(adata.obs_names)

        # Add QC columns to adata.obs
        qc_columns = [
            "total_counts", "n_genes",
            "neg_probe_sum", "neg_codeword_sum", "unassigned_sum",
            "fail_neg_probe", "fail_neg_codeword", "fail_unassigned",
            "fail_n_genes_low", "fail_total_counts_low", "fail_total_counts_high",
            "qc_pass"
        ]

        for col in qc_columns:
            if col in qc_aligned.columns:
                adata.obs[col] = qc_aligned[col].values

        # Handle any NaN from mismatched barcodes (mark as QC fail)
        if adata.obs["qc_pass"].isna().any():
            n_na = adata.obs["qc_pass"].isna().sum()
            print(f"    Setting {n_na} cells with missing QC to fail")
            adata.obs["qc_pass"] = adata.obs["qc_pass"].fillna(False)
            for col in qc_columns:
                if col.startswith("fail_"):
                    adata.obs[col] = adata.obs[col].fillna(True)

        # Save updated h5ad
        print(f"  Saving updated h5ad...")
        adata.write_h5ad(h5ad_path)

        # Collect summary stats
        n_total = len(adata)
        n_pass = adata.obs["qc_pass"].sum()
        n_fail = n_total - n_pass
        elapsed = time.time() - t0

        row = {
            "sample_id": sample_id,
            "n_cells": n_total,
            "n_pass": int(n_pass),
            "n_fail": int(n_fail),
            "pct_pass": round(100 * n_pass / n_total, 2),
            "fail_neg_probe": int(adata.obs["fail_neg_probe"].sum()),
            "fail_neg_codeword": int(adata.obs["fail_neg_codeword"].sum()),
            "fail_unassigned": int(adata.obs["fail_unassigned"].sum()),
            "fail_n_genes_low": int(adata.obs["fail_n_genes_low"].sum()),
            "fail_total_counts_low": int(adata.obs["fail_total_counts_low"].sum()),
            "fail_total_counts_high": int(adata.obs["fail_total_counts_high"].sum()),
            "median_total_counts": float(adata.obs["total_counts"].median()),
            "median_n_genes": float(adata.obs["n_genes"].median()),
            "time_s": round(elapsed, 1),
        }
        summary_rows.append(row)
        print(f"  Done in {elapsed:.1f}s: {n_pass:,}/{n_total:,} pass "
              f"({100*n_pass/n_total:.1f}%)")

    # Save summary CSV
    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(OUTPUT_DIR, "qc_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    print(f"\nSummary saved: {summary_path}")

    # Print overall summary
    total_cells = summary_df["n_cells"].sum()
    total_pass = summary_df["n_pass"].sum()
    total_fail = summary_df["n_fail"].sum()
    print(f"\n{'='*60}")
    print(f"OVERALL QC SUMMARY")
    print(f"{'='*60}")
    print(f"  Total cells:  {total_cells:,}")
    print(f"  QC pass:      {total_pass:,} ({100*total_pass/total_cells:.1f}%)")
    print(f"  QC fail:      {total_fail:,} ({100*total_fail/total_cells:.1f}%)")
    print(f"  Per-sample pass rates: "
          f"{summary_df['pct_pass'].min():.1f}% - {summary_df['pct_pass'].max():.1f}%")

    # Generate validation figures
    print("\nGenerating validation figures...")
    generate_qc_figures(summary_df, H5AD_DIR, OUTPUT_DIR)

    total_time = time.time() - t_start
    print(f"\nTotal time: {total_time:.0f}s")


def generate_qc_figures(summary_df, h5ad_dir, output_dir):
    """Generate QC validation figures."""

    fig, axes = plt.subplots(2, 3, figsize=(20, 14))

    # --- Panel 1: Per-sample QC pass rates ---
    ax = axes[0, 0]
    colors = ['#4ecdc4' if 'Br' in sid else '#e94560'
              for sid in summary_df['sample_id']]
    bars = ax.bar(range(len(summary_df)), summary_df['pct_pass'], color='#4ecdc4')
    ax.set_xticks(range(len(summary_df)))
    ax.set_xticklabels(summary_df['sample_id'], rotation=90, fontsize=10)
    ax.set_ylabel('% QC Pass', fontsize=14)
    ax.set_title('Per-Sample QC Pass Rate', fontsize=16)
    ax.axhline(y=summary_df['pct_pass'].mean(), color='red', linestyle='--',
               label=f"Mean: {summary_df['pct_pass'].mean():.1f}%")
    ax.legend(fontsize=12)
    ax.set_ylim(80, 100)

    # --- Panel 2: Failure reasons (stacked) ---
    ax = axes[0, 1]
    fail_cols = ['fail_neg_probe', 'fail_neg_codeword', 'fail_unassigned',
                 'fail_n_genes_low', 'fail_total_counts_low', 'fail_total_counts_high']
    fail_labels = ['Neg probe', 'Neg codeword', 'Unassigned', 'Low n_genes',
                   'Low counts', 'High counts']
    fail_colors = ['#e94560', '#ff6b6b', '#ffa07a', '#4ecdc4', '#45b7d1', '#96ceb4']

    x = range(len(summary_df))
    bottom = np.zeros(len(summary_df))
    for col, label, color in zip(fail_cols, fail_labels, fail_colors):
        pcts = 100 * summary_df[col] / summary_df['n_cells']
        ax.bar(x, pcts, bottom=bottom, label=label, color=color, alpha=0.8)
        bottom += pcts.values
    ax.set_xticks(range(len(summary_df)))
    ax.set_xticklabels(summary_df['sample_id'], rotation=90, fontsize=10)
    ax.set_ylabel('% Cells Failing', fontsize=14)
    ax.set_title('QC Failure Breakdown', fontsize=16)
    ax.legend(fontsize=10, loc='upper right')

    # --- Panels 3-6: Load one sample to show distributions ---
    # Pick a representative sample (median pass rate)
    median_idx = (summary_df['pct_pass'] - summary_df['pct_pass'].median()).abs().idxmin()
    rep_sid = summary_df.loc[median_idx, 'sample_id']

    try:
        rep_adata = ad.read_h5ad(os.path.join(h5ad_dir, f"{rep_sid}_annotated.h5ad"))

        # Panel 3: total_counts distribution
        ax = axes[0, 2]
        counts = rep_adata.obs["total_counts"]
        qc = rep_adata.obs["qc_pass"]
        ax.hist(counts[qc], bins=100, alpha=0.7, color='#4ecdc4', label='Pass')
        ax.hist(counts[~qc], bins=100, alpha=0.7, color='#e94560', label='Fail')
        ax.set_xlabel('Total Counts', fontsize=14)
        ax.set_ylabel('# Cells', fontsize=14)
        ax.set_title(f'Total Counts ({rep_sid})', fontsize=16)
        ax.legend(fontsize=12)

        # Panel 4: n_genes distribution
        ax = axes[1, 0]
        ngenes = rep_adata.obs["n_genes"]
        ax.hist(ngenes[qc], bins=100, alpha=0.7, color='#4ecdc4', label='Pass')
        ax.hist(ngenes[~qc], bins=100, alpha=0.7, color='#e94560', label='Fail')
        ax.set_xlabel('# Detected Genes', fontsize=14)
        ax.set_ylabel('# Cells', fontsize=14)
        ax.set_title(f'Detected Genes ({rep_sid})', fontsize=16)
        ax.legend(fontsize=12)

        # Panel 5: Control feature sums
        ax = axes[1, 1]
        for col, label, color in [
            ("neg_probe_sum", "Neg probe", "#e94560"),
            ("neg_codeword_sum", "Neg codeword", "#ff6b6b"),
            ("unassigned_sum", "Unassigned", "#ffa07a"),
        ]:
            vals = rep_adata.obs[col]
            ax.hist(vals, bins=50, alpha=0.6, color=color, label=label)
        ax.set_xlabel('Sum of Control Features', fontsize=14)
        ax.set_ylabel('# Cells', fontsize=14)
        ax.set_title(f'Control Feature Sums ({rep_sid})', fontsize=16)
        ax.legend(fontsize=12)

        # Panel 6: Cell type composition of QC fail vs pass
        ax = axes[1, 2]
        if "subclass_label" in rep_adata.obs.columns:
            # Get proportions
            pass_types = rep_adata.obs.loc[qc, "subclass_label"].value_counts(normalize=True)
            fail_types = rep_adata.obs.loc[~qc, "subclass_label"].value_counts(normalize=True)
            # Combine
            comp_df = pd.DataFrame({
                "Pass": pass_types,
                "Fail": fail_types,
            }).fillna(0)
            # Sort by absolute difference
            comp_df["diff"] = comp_df["Fail"] - comp_df["Pass"]
            comp_df = comp_df.sort_values("diff", ascending=True)
            # Plot
            y = range(len(comp_df))
            ax.barh(y, comp_df["diff"] * 100, color=[
                '#e94560' if d > 0 else '#4ecdc4' for d in comp_df["diff"]
            ])
            ax.set_yticks(y)
            ax.set_yticklabels(comp_df.index, fontsize=9)
            ax.set_xlabel('Fail% - Pass% (pp)', fontsize=14)
            ax.set_title(f'Cell Type Enrichment in QC Failures ({rep_sid})', fontsize=16)
            ax.axvline(0, color='white', linewidth=0.5)
        else:
            ax.text(0.5, 0.5, "No subclass labels", transform=ax.transAxes,
                    ha='center', va='center', fontsize=14)

    except Exception as e:
        print(f"  Warning: Could not load {rep_sid} for figure: {e}")
        for ax in [axes[0, 2], axes[1, 0], axes[1, 1], axes[1, 2]]:
            ax.text(0.5, 0.5, "Could not load data", transform=ax.transAxes,
                    ha='center', va='center', fontsize=14)

    plt.tight_layout()
    fig_path = os.path.join(output_dir, "qc_validation.png")
    plt.savefig(fig_path, dpi=150, facecolor='white')
    plt.close()
    print(f"  Saved: {fig_path}")


if __name__ == "__main__":
    main()
