#!/usr/bin/env python3
"""
Cell type-specific differential expression analysis using edgePython.

For each cell type (subclass), aggregates Xenium counts into pseudobulk
(summing across cells within each sample × cell type), then runs edgeR-style
GLM quasi-likelihood testing for SCZ vs Control, adjusting for sex and age.

This is the DE analogue of the crumblr differential abundance analysis.

Usage:
    /Users/shreejoy/venv/bin/python3 -u run_de_edgepython.py [--level subclass|supertype]

Output:
    output/de/de_results_{level}.csv         — all DE results
    output/de/de_summary_{level}.csv         — summary: n_sig per cell type
    output/de/de_volcano_{celltype}.png      — volcano plot per cell type
    output/de/de_summary_heatmap_{level}.png — heatmap of top DE genes
"""

import os
import sys
import time
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import edgepython as ep

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    BASE_DIR, H5AD_DIR, EXCLUDE_SAMPLES, SUBCLASS_TO_CLASS,
    load_sample_adata,
)

sys.path.insert(0, os.path.join(BASE_DIR, "code"))
from modules.metadata import get_subject_info

METADATA_PATH = os.path.join(BASE_DIR, "data", "sample_metadata.xlsx")
OUT_DIR = os.path.join(BASE_DIR, "output", "de")

# Minimum cells per sample × cell type to include in pseudobulk
MIN_CELLS_PER_PB = 10
# Minimum samples per group (SCZ and Control) for a cell type to be tested
MIN_SAMPLES_PER_GROUP = 3
# FDR threshold for significance
FDR_THRESH = 0.20


# ──────────────────────────────────────────────────────────────────────
# Data loading
# ──────────────────────────────────────────────────────────────────────

def load_sample_counts(sample_id, qc_mode='corr'):
    """
    Load one sample, filter to QC-pass cortical cells, return raw counts
    and metadata.

    Uses load_sample_adata() from config for standardized loading.

    Returns (counts_matrix, gene_names, obs_df) where counts_matrix is
    cells × genes (sparse).
    """
    adata = load_sample_adata(sample_id, cortical_only=True, qc_mode=qc_mode)
    obs = adata.obs[["sample_id", "subclass_label", "supertype_label"]].copy()
    return adata.X, list(adata.var_names), obs


def build_pseudobulk(all_counts, all_genes, all_obs, min_cells=MIN_CELLS_PER_PB,
                     group_col="subclass_label"):
    """
    Aggregate single-cell counts into pseudobulk per (sample, group).

    Parameters
    ----------
    group_col : str
        Column in obs to group by (e.g. "subclass_label" or "supertype_label").

    Returns a dict: {group_name: (counts_matrix, gene_names, sample_info_df)}
    where counts_matrix is genes × samples (numpy array).
    """
    import scipy.sparse as sp

    # Combine all obs
    obs_combined = pd.concat(all_obs, ignore_index=True)

    # Verify all samples share the same gene set
    assert all(g == all_genes[0] for g in all_genes), \
        "Gene sets differ across samples!"
    gene_names = all_genes[0]
    n_genes = len(gene_names)

    # Stack all count matrices
    if sp.issparse(all_counts[0]):
        X = sp.vstack(all_counts)
    else:
        X = np.vstack(all_counts)

    print(f"  Combined: {X.shape[0]:,} cells × {X.shape[1]:,} genes")

    # Group by (sample, group_col) and sum
    subclasses = sorted(obs_combined[group_col].unique())
    samples = sorted(obs_combined["sample_id"].unique())

    result = {}
    for sc in subclasses:
        sc_mask = obs_combined[group_col] == sc
        sc_samples = []
        sc_counts = []

        for sid in samples:
            sample_mask = (obs_combined["sample_id"] == sid) & sc_mask
            n_cells = sample_mask.sum()

            if n_cells < min_cells:
                continue

            # Sum counts across cells for this sample × cell type
            if sp.issparse(X):
                pb = np.asarray(X[sample_mask.values].sum(axis=0)).ravel()
            else:
                pb = X[sample_mask.values].sum(axis=0).ravel()

            sc_samples.append({"sample_id": sid, "n_cells": n_cells})
            sc_counts.append(pb)

        if len(sc_samples) < 2:
            continue

        # genes × samples matrix
        counts_matrix = np.column_stack(sc_counts).astype(np.float64)
        sample_info = pd.DataFrame(sc_samples)

        result[sc] = (counts_matrix, gene_names, sample_info)

    return result


# ──────────────────────────────────────────────────────────────────────
# DE analysis
# ──────────────────────────────────────────────────────────────────────

def run_de_for_celltype(counts, gene_names, sample_info, meta_df):
    """
    Run edgePython GLM-QL DE analysis for one cell type.

    Parameters
    ----------
    counts : np.ndarray (genes × samples)
    gene_names : list of str
    sample_info : DataFrame with sample_id, n_cells
    meta_df : DataFrame with sample_id, diagnosis, sex, age

    Returns
    -------
    DataFrame with DE results (one row per gene), or None if insufficient data.
    """
    # Merge metadata
    info = sample_info.merge(meta_df, on="sample_id", how="left")

    # Check we have both groups
    n_scz = (info["diagnosis"] == "SCZ").sum()
    n_ctrl = (info["diagnosis"] == "Control").sum()
    if n_scz < MIN_SAMPLES_PER_GROUP or n_ctrl < MIN_SAMPLES_PER_GROUP:
        return None

    n_samples = len(info)

    # Build design matrix: intercept + SCZ + sex + age
    design = pd.DataFrame({
        "Intercept": np.ones(n_samples),
        "SCZ": (info["diagnosis"].values == "SCZ").astype(float),
        "Male": (info["sex"].values == "M").astype(float),
        "Age": info["age"].values.astype(float),
    })

    # Center age
    design["Age"] = design["Age"] - design["Age"].mean()

    # Create DGEList
    genes_df = pd.DataFrame({"gene": gene_names})
    samples_df = info[["sample_id", "diagnosis", "n_cells"]].copy()
    samples_df.index = range(n_samples)

    try:
        # Create temporary DGEList just for filtering
        dge_tmp = ep.make_dgelist(counts=counts, genes=genes_df,
                                   samples=samples_df)

        # Filter low-expressed genes
        keep = ep.filter_by_expr(dge_tmp, design=design)
        if isinstance(keep, np.ndarray):
            n_keep = int(keep.sum())
        else:
            n_keep = sum(keep)

        if n_keep < 5:
            return None

        # Subset counts and genes BEFORE creating final DGEList
        # (DGEList indexing has issues, so we filter manually)
        counts_filt = counts[keep]
        genes_filt = genes_df[keep].reset_index(drop=True)

        dge = ep.make_dgelist(counts=counts_filt, genes=genes_filt,
                               samples=samples_df)

        # TMM normalization
        dge = ep.calc_norm_factors(dge, method="TMM")

        # Estimate dispersion
        dge = ep.estimate_disp(dge, design=design, robust=True)

        # GLM quasi-likelihood fit
        fit = ep.glm_ql_fit(dge, design=design, robust=True)

        # Test SCZ coefficient (column index 1)
        results = ep.glm_ql_ftest(fit, coef=1)  # SCZ column

        # Extract results — use n=n_keep (not None, which is buggy)
        tt = ep.top_tags(results, n=n_keep, sort_by="PValue")

        # Convert to DataFrame
        if isinstance(tt, dict) and "table" in tt:
            df = tt["table"].copy()
        elif hasattr(tt, 'table'):
            df = tt.table.copy()
        elif isinstance(tt, pd.DataFrame):
            df = tt.copy()
        else:
            df = pd.DataFrame(tt)

        df["n_scz"] = n_scz
        df["n_ctrl"] = n_ctrl
        df["n_samples"] = n_samples
        df["n_genes_tested"] = n_keep

        return df

    except Exception as e:
        print(f"    ERROR: {e}")
        import traceback
        traceback.print_exc()
        return None


# ──────────────────────────────────────────────────────────────────────
# Plotting
# ──────────────────────────────────────────────────────────────────────

def plot_volcano(df, celltype, out_path, fdr_thresh=FDR_THRESH):
    """Simple volcano plot for one cell type."""
    fig, ax = plt.subplots(figsize=(8, 6))

    # Compute -log10 p
    pvals = df["PValue"].values
    pvals = np.maximum(pvals, 1e-300)  # avoid log(0)
    neg_log_p = -np.log10(pvals)
    logfc = df["logFC"].values

    # Color by significance
    if "FDR" in df.columns:
        sig = df["FDR"].values < fdr_thresh
    else:
        sig = np.zeros(len(df), dtype=bool)

    ax.scatter(logfc[~sig], neg_log_p[~sig], s=8, c="#999999", alpha=0.5,
               linewidths=0, rasterized=True)
    ax.scatter(logfc[sig], neg_log_p[sig], s=20, c="#e74c3c", alpha=0.8,
               linewidths=0, rasterized=True)

    # Label top genes
    if sig.sum() > 0:
        top_idx = np.argsort(pvals)[:10]
        gene_col = "gene" if "gene" in df.columns else df.columns[0]
        for idx in top_idx:
            if sig[idx]:
                ax.annotate(df[gene_col].iloc[idx],
                            (logfc[idx], neg_log_p[idx]),
                            fontsize=8, xytext=(5, 5),
                            textcoords="offset points")

    n_sig = sig.sum()
    n_up = ((logfc > 0) & sig).sum()
    n_down = ((logfc < 0) & sig).sum()

    ax.axhline(-np.log10(0.05), color="#bbbbbb", linestyle="--", linewidth=0.8)
    ax.set_xlabel("log₂ Fold Change (SCZ vs Control)", fontsize=14)
    ax.set_ylabel("-log₁₀ p-value", fontsize=14)
    ax.set_title(f"{celltype}\n{n_sig} DE genes (FDR<{fdr_thresh}): "
                  f"{n_up} up, {n_down} down",
                  fontsize=16, fontweight="bold")

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()


def plot_summary(all_results, out_path):
    """Summary heatmap of DE results across cell types."""
    # Compute number of significant genes per cell type
    summary = []
    for ct, df in sorted(all_results.items()):
        if df is None or len(df) == 0:
            continue
        n_tested = len(df)
        n_sig = (df["FDR"] < FDR_THRESH).sum() if "FDR" in df.columns else 0
        n_up = ((df["logFC"] > 0) & (df["FDR"] < FDR_THRESH)).sum() \
            if "FDR" in df.columns else 0
        n_down = ((df["logFC"] < 0) & (df["FDR"] < FDR_THRESH)).sum() \
            if "FDR" in df.columns else 0
        summary.append({
            "celltype": ct,
            "n_tested": n_tested,
            "n_sig": n_sig,
            "n_up": n_up,
            "n_down": n_down,
            "class": SUBCLASS_TO_CLASS.get(ct, "Other"),
        })

    if not summary:
        print("  No results to summarize")
        return

    sum_df = pd.DataFrame(summary).sort_values("n_sig", ascending=True)

    fig, ax = plt.subplots(figsize=(10, max(6, len(sum_df) * 0.4)))

    y_pos = range(len(sum_df))
    ax.barh(y_pos, sum_df["n_up"], color="#e74c3c", alpha=0.7, label="Up in SCZ")
    ax.barh(y_pos, -sum_df["n_down"], color="#3498db", alpha=0.7,
            label="Down in SCZ")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(sum_df["celltype"], fontsize=12)
    ax.set_xlabel(f"Number of DE genes (FDR < {FDR_THRESH})", fontsize=14)
    ax.set_title("Differential Expression Summary: SCZ vs Control",
                  fontsize=18, fontweight="bold")
    ax.legend(fontsize=12)
    ax.axvline(0, color="black", linewidth=0.5)

    # Add text labels
    for i, (_, row) in enumerate(sum_df.iterrows()):
        if row["n_up"] > 0:
            ax.text(row["n_up"] + 0.3, i, str(int(row["n_up"])),
                    va="center", fontsize=9, color="#e74c3c")
        if row["n_down"] > 0:
            ax.text(-row["n_down"] - 0.3, i, str(int(row["n_down"])),
                    va="center", ha="right", fontsize=9, color="#3498db")

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--level", default="subclass",
                        choices=["subclass", "supertype"],
                        help="Taxonomy level for pseudobulk grouping")
    parser.add_argument("--qc-mode", default="corr", choices=["corr", "hybrid"],
                        help="QC mode: 'corr' (default) or 'hybrid' (nuclear doublet-resolved)")
    args = parser.parse_args()

    level = args.level
    qc_mode = args.qc_mode
    suffix = "_hybrid" if qc_mode == "hybrid" else ""
    group_col = f"{level}_label"
    print(f"Running DE at {level} level (grouping by {group_col})")
    print(f"QC mode: {qc_mode}")
    print(f"FDR threshold: {FDR_THRESH}")

    t_start = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)

    # Load metadata
    print("Loading donor metadata...")
    meta = get_subject_info(METADATA_PATH)
    print(f"  {len(meta)} subjects")
    print(f"  SCZ: {(meta['diagnosis']=='SCZ').sum()}, "
          f"Control: {(meta['diagnosis']=='Control').sum()}")

    # Load all samples
    h5ad_files = sorted(glob.glob(os.path.join(H5AD_DIR, "*_annotated.h5ad")))
    print(f"\nLoading {len(h5ad_files)} h5ad files...")

    all_counts = []
    all_genes = []
    all_obs = []

    for fpath in h5ad_files:
        sid = os.path.basename(fpath).replace("_annotated.h5ad", "")
        if sid in EXCLUDE_SAMPLES:
            print(f"  Skipping {sid} (excluded)")
            continue

        t0 = time.time()
        counts, genes, obs = load_sample_counts(sid, qc_mode=qc_mode)
        all_counts.append(counts)
        all_genes.append(genes)
        all_obs.append(obs)
        print(f"  {sid}: {counts.shape[0]:,} cells, {time.time()-t0:.1f}s")

    # Build pseudobulk
    print(f"\nBuilding pseudobulk (min {MIN_CELLS_PER_PB} cells per group)...")
    pb_data = build_pseudobulk(all_counts, all_genes, all_obs,
                               group_col=group_col)
    print(f"  Cell types with sufficient data: {len(pb_data)}")

    for ct, (counts, genes, sinfo) in sorted(pb_data.items()):
        print(f"  {ct}: {counts.shape[0]} genes × {counts.shape[1]} samples "
              f"(cells per sample: {sinfo['n_cells'].median():.0f} median)")

    # Run DE for each cell type
    print(f"\n{'='*60}")
    print("Running differential expression analysis...")
    print(f"{'='*60}")

    all_results = {}
    for ct, (counts, genes, sinfo) in sorted(pb_data.items()):
        print(f"\n  {ct}...")
        t_ct = time.time()

        result = run_de_for_celltype(counts, genes, sinfo, meta)

        if result is not None:
            result["celltype"] = ct
            result["cell_class"] = SUBCLASS_TO_CLASS.get(ct, "Other")
            all_results[ct] = result

            n_sig = (result["FDR"] < FDR_THRESH).sum() if "FDR" in result.columns else 0
            print(f"    {len(result)} genes tested, {n_sig} significant "
                  f"(FDR<{FDR_THRESH}), {time.time()-t_ct:.1f}s")
        else:
            all_results[ct] = None
            print(f"    Skipped (insufficient data)")

    # Combine results
    valid_results = [df for df in all_results.values() if df is not None]
    if valid_results:
        combined = pd.concat(valid_results, ignore_index=True)
        csv_path = os.path.join(OUT_DIR, f"de_results_{level}{suffix}.csv")
        combined.to_csv(csv_path, index=False)
        print(f"\nSaved: {csv_path} ({len(combined):,} rows)")

    # Generate plots
    print(f"\n{'='*60}")
    print("Generating plots...")
    print(f"{'='*60}")

    for ct, df in all_results.items():
        if df is not None and len(df) > 0:
            safe_name = ct.replace("/", "_").replace(" ", "_")
            plot_volcano(df, ct,
                         os.path.join(OUT_DIR, f"de_volcano_{level}{suffix}_{safe_name}.png"))

    # Summary
    plot_summary(all_results,
                 os.path.join(OUT_DIR, f"de_summary_heatmap_{level}{suffix}.png"))

    # Summary table
    summary_rows = []
    for ct, df in sorted(all_results.items()):
        if df is None:
            summary_rows.append({
                "celltype": ct, "n_tested": 0, "n_sig_fdr05": 0,
                "n_up": 0, "n_down": 0, "status": "skipped"
            })
        else:
            n_sig = (df["FDR"] < FDR_THRESH).sum() if "FDR" in df.columns else 0
            n_up = ((df["logFC"] > 0) & (df["FDR"] < FDR_THRESH)).sum() \
                if "FDR" in df.columns else 0
            n_down = ((df["logFC"] < 0) & (df["FDR"] < FDR_THRESH)).sum() \
                if "FDR" in df.columns else 0
            summary_rows.append({
                "celltype": ct,
                "n_tested": len(df),
                "n_sig_fdr05": n_sig,
                "n_up": n_up,
                "n_down": n_down,
                "status": "ok",
            })

    sum_df = pd.DataFrame(summary_rows)
    sum_path = os.path.join(OUT_DIR, f"de_summary_{level}{suffix}.csv")
    sum_df.to_csv(sum_path, index=False)
    print(f"Saved: {sum_path}")

    # Print summary
    total_time = time.time() - t_start
    print(f"\n{'='*60}")
    print(f"SUMMARY — {total_time:.0f}s total")
    print(f"{'='*60}")
    total_sig = sum(r["n_sig_fdr05"] for r in summary_rows)
    n_tested_ct = sum(1 for r in summary_rows if r["status"] == "ok")
    print(f"  Cell types tested: {n_tested_ct}")
    print(f"  Total DE genes (FDR<{FDR_THRESH}): {total_sig}")
    for r in sorted(summary_rows, key=lambda x: -x["n_sig_fdr05"]):
        if r["status"] == "ok":
            print(f"    {r['celltype']:20s}: {r['n_sig_fdr05']:3d} sig "
                  f"({r['n_up']} up, {r['n_down']} down) "
                  f"of {r['n_tested']} tested")


if __name__ == "__main__":
    main()
