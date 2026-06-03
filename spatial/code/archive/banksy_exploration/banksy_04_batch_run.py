#!/usr/bin/env python3
"""
BANKSY Phase 5: Batch run on all 24 samples.

For each sample:
  1. Run BANKSY (λ=0.8, res=0.3)
  2. Classify domains (Cortical / Vascular / WM)
  3. Flag L1 border cells (shallow non-neuronal BANKSY clusters = L1, not meningeal)
  4. Identify cortical strips (complete + partial tiers)
  5. Add new columns to h5ad and save

New h5ad columns (coexist with existing columns):
  - banksy_cluster:       int   — raw BANKSY cluster ID
  - banksy_domain:        str   — Cortical / Vascular / WM
  - banksy_is_l1:         bool  — True if cell in L1-enriched BANKSY cluster (pia border)
  - cortical_strip_id:    int   — strip ID (-1 if not in strip)
  - cortical_strip_tier:  str   — "complete" / "partial" / ""
  - in_cortical_strip:    bool  — True if cell is in any selected strip

Note: Previous "Meningeal" domain was found to be L1 cortex (confirmed by
MERFISH comparison: 81% non-neuronal composition matches L1). These cells
are now classified as Cortical with banksy_is_l1=True. The L1 border will
serve as input to the curved cortex strip algorithm.

Usage:
    python3 -u code/analysis/banksy_04_batch_run.py           # all 24 samples
    python3 -u code/analysis/banksy_04_batch_run.py --dry-run # list samples only
    python3 -u code/analysis/banksy_04_batch_run.py --samples Br6437 Br6389
"""

import os
import sys
import time
import argparse
import warnings
import traceback
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")

# BANKSY imports
from banksy.initialize_banksy import initialize_banksy
from banksy.embed_banksy import generate_banksy_matrix
from banksy.cluster_methods import run_Leiden_partition
from banksy_utils.umap_pca import pca_umap

# Project imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import H5AD_DIR, SAMPLE_TO_DX

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "modules"))
from depth_model import LAYER_BINS
from spatial_domains import VASCULAR_TYPES, NON_NEURONAL_TYPES

# ── Constants (matching banksy_02_domains_and_strips.py) ─────────────

# BANKSY parameters
BANKSY_LAMBDA = 0.8
BANKSY_RESOLUTION = 0.3
K_GEOM = 15
PCA_DIMS = [20]

# Domain classification thresholds
VASCULAR_THRESH = 0.50
WM_OLIGO_THRESH = 0.40
WM_DEPTH_THRESH = 0.80

# L1 identification thresholds (formerly "Meningeal" — actually L1 cortex)
L1_NN_THRESH = 0.50       # non-neuronal fraction
L1_DEPTH_THRESH = 0.20    # mean depth < 0.20 (cortical surface)

# Cortical strip parameters
STRIP_WIDTH_UM = 750
MIN_CELLS_PER_LAYER = 15
DEPTH_GRID_BINS = 50
DEPTH_SMOOTH_SIGMA = 2.0

REQUIRED_LAYERS = ["L1", "L2/3", "L4", "L5", "L6"]
LAYER_DEPTH_ORDER = {"L1": 0, "L2/3": 1, "L4": 2, "L5": 3, "L6": 4}

# Output
OUT_DIR = os.path.expanduser("~/Github/SCZ_Xenium/output/banksy")
os.makedirs(OUT_DIR, exist_ok=True)


# ── Processing Functions ─────────────────────────────────────────────

def preprocess_for_banksy(adata):
    """Normalize → log1p → z-score for BANKSY."""
    adata_b = adata.copy()
    adata_b.layers["counts"] = adata_b.X.copy()
    sc.pp.normalize_total(adata_b, target_sum=1e4)
    sc.pp.log1p(adata_b)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.scale(adata_b)
    return adata_b


def run_banksy(adata_b):
    """Run BANKSY, return cluster labels."""
    coord_keys = ("x", "y", "spatial")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        banksy_dict = initialize_banksy(
            adata_b, coord_keys, K_GEOM,
            nbr_weight_decay="scaled_gaussian", max_m=1,
            plt_edge_hist=False, plt_nbr_weights=False, plt_agf_angles=False,
        )
        banksy_dict, _ = generate_banksy_matrix(
            adata_b, banksy_dict, [BANKSY_LAMBDA], max_m=1
        )
        pca_umap(banksy_dict, pca_dims=PCA_DIMS, add_umap=False)
        results_df, _ = run_Leiden_partition(
            banksy_dict, [BANKSY_RESOLUTION],
            num_nn=50, num_iterations=-1, partition_seed=42,
            match_labels=False,
        )

    labels = results_df.iloc[0]["labels"]
    if hasattr(labels, "dense"):
        labels = labels.dense
    return np.asarray(labels).astype(int)


def classify_banksy_domains(adata, banksy_labels):
    """Classify BANKSY clusters into spatial domains + flag L1 border clusters.

    Domains: Cortical / Vascular / WM (no more "Meningeal" — those are L1).

    Returns
    -------
    domains : np.ndarray of str
        Per-cell domain labels.
    is_l1 : np.ndarray of bool
        Per-cell flag for L1 border cells (shallow, non-neuronal-dominated clusters).
    cluster_info : dict
        Per-cluster statistics including domain and l1 status.
    """
    subclass = adata.obs.get("corr_subclass",
                              adata.obs.get("subclass_label")).values.astype(str)
    pred_depth = pd.to_numeric(adata.obs["predicted_norm_depth"],
                                errors="coerce").values

    unique_cl = np.unique(banksy_labels)
    cluster_info = {}

    for cl in unique_cl:
        mask = banksy_labels == cl
        n_cl = mask.sum()
        cl_sub = subclass[mask]
        cl_depth = pred_depth[mask]

        sub_counts = {}
        for s in cl_sub:
            sub_counts[s] = sub_counts.get(s, 0) + 1

        vasc_frac = sum(sub_counts.get(v, 0) for v in VASCULAR_TYPES) / n_cl
        nn_frac = sum(sub_counts.get(v, 0) for v in NON_NEURONAL_TYPES) / n_cl
        oligo_frac = sub_counts.get("Oligodendrocyte", 0) / n_cl
        neuronal_frac = 1.0 - nn_frac

        valid_depth = cl_depth[~np.isnan(cl_depth)]
        mean_depth = float(np.mean(valid_depth)) if len(valid_depth) > 0 else np.nan

        # Top 3 types for diagnostics
        top3 = sorted(sub_counts.items(), key=lambda x: -x[1])[:3]
        top3_str = ", ".join(f"{s}({100*c/n_cl:.0f}%)" for s, c in top3)

        # Domain classification (Cortical / Vascular / WM)
        # L1 detection: shallow non-neuronal clusters are L1, not meningeal
        is_l1_cluster = False
        if vasc_frac > VASCULAR_THRESH:
            domain = "Vascular"
        elif oligo_frac > WM_OLIGO_THRESH and mean_depth > WM_DEPTH_THRESH:
            domain = "WM"
        elif nn_frac > L1_NN_THRESH and (np.isnan(mean_depth) or mean_depth < L1_DEPTH_THRESH):
            # Previously "Meningeal" — actually L1 cortex
            domain = "Cortical"
            is_l1_cluster = True
        elif neuronal_frac > 0.20 and not np.isnan(mean_depth) and 0.0 <= mean_depth <= 0.90:
            domain = "Cortical"
        elif mean_depth > WM_DEPTH_THRESH:
            domain = "WM"
        else:
            domain = "Cortical"

        cluster_info[cl] = {
            "domain": domain,
            "is_l1": is_l1_cluster,
            "n_cells": n_cl,
            "vasc_frac": vasc_frac,
            "nn_frac": nn_frac,
            "oligo_frac": oligo_frac,
            "mean_depth": mean_depth,
            "top3": top3_str,
        }

    # Print per-cluster summary
    print(f"  {'Cl':>4} | {'N':>6} | {'Domain':<9} | {'L1?':>3} | {'Vasc%':>6} | {'NN%':>6} "
          f"| {'Oligo%':>6} | {'Depth':>6} | Top types")
    print(f"  {'-'*100}")
    for cl in sorted(cluster_info.keys()):
        info = cluster_info[cl]
        l1_str = "YES" if info["is_l1"] else ""
        print(f"  {cl:>4} | {info['n_cells']:>6,} | {info['domain']:<9} | {l1_str:>3} "
              f"| {info['vasc_frac']*100:>5.1f}% | {info['nn_frac']*100:>5.1f}% "
              f"| {info['oligo_frac']*100:>5.1f}% | {info['mean_depth']:>6.3f} "
              f"| {info['top3']}")

    # Map to per-cell arrays
    domains = np.array([cluster_info[cl]["domain"] for cl in banksy_labels])
    is_l1 = np.array([cluster_info[cl]["is_l1"] for cl in banksy_labels])
    return domains, is_l1, cluster_info


def identify_cortical_strips(adata, domains, pred_depth):
    """Identify cortical strips using depth gradient orientation."""
    from scipy.ndimage import gaussian_filter, generic_filter
    from scipy.stats import kendalltau

    coords = adata.obsm["spatial"]
    x, y = coords[:, 0], coords[:, 1]
    cortical_mask = domains == "Cortical"

    # Estimate depth gradient orientation
    valid = ~np.isnan(pred_depth)
    depth_sum, x_edges, y_edges = np.histogram2d(
        x[valid], y[valid], bins=DEPTH_GRID_BINS, weights=pred_depth[valid])
    count, _, _ = np.histogram2d(x[valid], y[valid], bins=DEPTH_GRID_BINS)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mean_depth_grid = np.where(count > 5, depth_sum / count, np.nan)

    def fill_nan(values):
        center = values[len(values) // 2]
        if np.isnan(center):
            valid_v = values[~np.isnan(values)]
            return np.mean(valid_v) if len(valid_v) > 0 else np.nan
        return center

    filled = generic_filter(mean_depth_grid, fill_nan, size=5,
                             mode='constant', cval=np.nan)
    smoothed = gaussian_filter(np.nan_to_num(filled, nan=0.5), sigma=DEPTH_SMOOTH_SIGMA)
    grad_y, grad_x = np.gradient(smoothed)

    valid_mask = count > 10
    if valid_mask.sum() > 0:
        angles = np.arctan2(grad_y[valid_mask], grad_x[valid_mask])
        angle = np.median(angles)
    else:
        angle = np.pi / 2

    # Tile into strips perpendicular to depth gradient
    perp_angle = angle + np.pi / 2
    proj = x * np.cos(perp_angle) + y * np.sin(perp_angle)
    proj_min, proj_max = proj.min(), proj.max()
    n_strips = max(1, int((proj_max - proj_min) / STRIP_WIDTH_UM))
    strip_edges = np.linspace(proj_min, proj_max, n_strips + 1)
    strip_assignment = np.clip(np.digitize(proj, strip_edges) - 1, 0, n_strips - 1)

    # Depth-bin layers for scoring
    depth_layers = np.full(len(pred_depth), "Unknown", dtype=object)
    for layer_name, (lo, hi) in LAYER_BINS.items():
        mask = (pred_depth >= lo) & (pred_depth < hi)
        depth_layers[mask] = layer_name

    # Score strips
    complete_ids = set()
    partial_ids = set()

    for s in range(n_strips):
        in_strip = (strip_assignment == s) & cortical_mask
        n_total = in_strip.sum()
        if n_total < MIN_CELLS_PER_LAYER * 3:
            continue

        # Layer completeness
        strip_layers = depth_layers[in_strip]
        strip_depths = pred_depth[in_strip]
        layer_counts = {l: (strip_layers == l).sum() for l in REQUIRED_LAYERS}
        layers_present = [l for l in REQUIRED_LAYERS
                         if layer_counts.get(l, 0) >= MIN_CELLS_PER_LAYER]
        completeness = len(layers_present) / len(REQUIRED_LAYERS)

        # Laminar order (Kendall tau)
        order_score = 0.0
        if len(layers_present) >= 3:
            expected_order = sorted(layers_present, key=lambda l: LAYER_DEPTH_ORDER[l])
            actual_medians = {l: np.nanmedian(strip_depths[strip_layers == l])
                             for l in layers_present
                             if (strip_layers == l).sum() >= MIN_CELLS_PER_LAYER}
            if len(actual_medians) >= 3:
                actual_order = sorted(actual_medians.keys(),
                                     key=lambda l: actual_medians[l])
                expected_ranks = [expected_order.index(l) for l in actual_order]
                tau, _ = kendalltau(range(len(expected_ranks)), expected_ranks)
                order_score = max(0, tau)

        # Purity
        all_in_strip = strip_assignment == s
        purity = cortical_mask[all_in_strip].sum() / max(1, all_in_strip.sum())

        # Tier classification
        if completeness == 1.0 and order_score > 0.8 and purity > 0.75:
            complete_ids.add(s)
        elif completeness >= 0.80 and order_score > 0.6 and purity > 0.60:
            partial_ids.add(s)

    # Build per-cell outputs
    selected = complete_ids | partial_ids
    strip_ids = np.full(len(pred_depth), -1, dtype=int)
    strip_tiers = np.full(len(pred_depth), "", dtype=object)
    in_strip = np.zeros(len(pred_depth), dtype=bool)

    for i, s in enumerate(strip_assignment):
        if s in selected and cortical_mask[i]:
            strip_ids[i] = s
            strip_tiers[i] = "complete" if s in complete_ids else "partial"
            in_strip[i] = True

    return strip_ids, strip_tiers, in_strip, len(complete_ids), len(partial_ids)


def process_one_sample(h5ad_path):
    """Process one sample: run BANKSY, classify domains, find strips, save to h5ad."""
    t0 = time.time()
    sample_id = os.path.basename(h5ad_path).replace("_annotated.h5ad", "")
    dx = SAMPLE_TO_DX.get(sample_id, "?")

    print(f"\n  [{sample_id}] ({dx}) Loading...")
    adata = ad.read_h5ad(h5ad_path)
    n_total = adata.n_obs

    # QC mask
    if "hybrid_qc_pass" in adata.obs.columns:
        qc_mask = adata.obs["hybrid_qc_pass"].values.astype(bool)
    elif "qc_pass" in adata.obs.columns:
        qc_mask = adata.obs["qc_pass"].values.astype(bool)
    else:
        qc_mask = np.ones(n_total, dtype=bool)

    n_pass = qc_mask.sum()
    adata_pass = adata[qc_mask].copy()
    print(f"  [{sample_id}] {n_pass:,} / {n_total:,} QC-pass cells")

    # Preprocess and run BANKSY
    t_banksy = time.time()
    adata_b = preprocess_for_banksy(adata_pass)
    banksy_labels = run_banksy(adata_b)
    del adata_b  # Free memory
    n_clusters = len(np.unique(banksy_labels))
    print(f"  [{sample_id}] BANKSY: {n_clusters} clusters ({time.time()-t_banksy:.0f}s)")

    # Domain classification
    domains, is_l1, cluster_info = classify_banksy_domains(adata_pass, banksy_labels)
    domain_counts = {d: (domains == d).sum() for d in ["Cortical", "Vascular", "WM"]}
    n_l1 = is_l1.sum()
    print(f"  [{sample_id}] Domains: " + ", ".join(f"{d}={n:,}" for d, n in domain_counts.items())
          + f", L1_border={n_l1:,}")

    # Cortical strips
    pred_depth = pd.to_numeric(adata_pass.obs["predicted_norm_depth"],
                                errors="coerce").values
    strip_ids, strip_tiers, in_strip, n_complete, n_partial = \
        identify_cortical_strips(adata_pass, domains, pred_depth)
    n_in_strip = in_strip.sum()
    n_cortical = domain_counts.get("Cortical", 0)
    coverage = n_in_strip / max(1, n_cortical) * 100
    print(f"  [{sample_id}] Strips: {n_complete} complete + {n_partial} partial, "
          f"{n_in_strip:,} cells ({coverage:.1f}% of cortical)")

    # Write back to full adata (NaN/-1/""/False for non-QC cells)
    full_banksy_cluster = np.full(n_total, -1, dtype=int)
    full_banksy_cluster[qc_mask] = banksy_labels
    adata.obs["banksy_cluster"] = full_banksy_cluster

    full_banksy_domain = np.full(n_total, "", dtype=object)
    full_banksy_domain[qc_mask] = domains
    adata.obs["banksy_domain"] = full_banksy_domain

    full_banksy_is_l1 = np.zeros(n_total, dtype=bool)
    full_banksy_is_l1[qc_mask] = is_l1
    adata.obs["banksy_is_l1"] = full_banksy_is_l1

    full_strip_ids = np.full(n_total, -1, dtype=int)
    full_strip_ids[qc_mask] = strip_ids
    adata.obs["cortical_strip_id"] = full_strip_ids

    full_strip_tiers = np.full(n_total, "", dtype=object)
    full_strip_tiers[qc_mask] = strip_tiers
    adata.obs["cortical_strip_tier"] = full_strip_tiers

    full_in_strip = np.zeros(n_total, dtype=bool)
    full_in_strip[qc_mask] = in_strip
    adata.obs["in_cortical_strip"] = full_in_strip

    # Save
    adata.write_h5ad(h5ad_path)
    elapsed = time.time() - t0
    print(f"  [{sample_id}] Saved h5ad ({elapsed:.0f}s)")

    return {
        "sample_id": sample_id,
        "diagnosis": dx,
        "n_cells": n_total,
        "n_qc_pass": n_pass,
        "n_clusters": n_clusters,
        "n_cortical": domain_counts.get("Cortical", 0),
        "n_vascular": domain_counts.get("Vascular", 0),
        "n_wm": domain_counts.get("WM", 0),
        "n_l1_border": n_l1,
        "pct_cortical": domain_counts.get("Cortical", 0) / n_pass * 100,
        "pct_vascular": domain_counts.get("Vascular", 0) / n_pass * 100,
        "pct_l1_border": n_l1 / n_pass * 100,
        "n_complete_strips": n_complete,
        "n_partial_strips": n_partial,
        "n_cells_in_strips": n_in_strip,
        "strip_coverage_pct": coverage,
        "time_sec": elapsed,
        "status": "success",
    }


# ── Main ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", nargs="+", default=None,
                        help="Specific sample IDs (default: all 24)")
    parser.add_argument("--dry-run", action="store_true",
                        help="List samples without processing")
    args = parser.parse_args()

    # Find all h5ad files
    all_h5ad = sorted([
        os.path.join(H5AD_DIR, f)
        for f in os.listdir(H5AD_DIR)
        if f.endswith("_annotated.h5ad")
    ])

    if args.samples:
        all_h5ad = [p for p in all_h5ad
                    if os.path.basename(p).replace("_annotated.h5ad", "") in args.samples]

    sample_ids = [os.path.basename(p).replace("_annotated.h5ad", "") for p in all_h5ad]

    print(f"BANKSY Batch Run: {len(all_h5ad)} samples")
    print(f"  Parameters: λ={BANKSY_LAMBDA}, res={BANKSY_RESOLUTION}")
    print(f"  Strips: width={STRIP_WIDTH_UM}μm, min_cells/layer={MIN_CELLS_PER_LAYER}")
    print(f"  Samples: {sample_ids}")

    if args.dry_run:
        for p in all_h5ad:
            sid = os.path.basename(p).replace("_annotated.h5ad", "")
            size_mb = os.path.getsize(p) / 1e6
            print(f"    {sid} ({SAMPLE_TO_DX.get(sid, '?')}): {size_mb:.0f} MB")
        return

    t_start = time.time()
    results = []

    for i, h5ad_path in enumerate(all_h5ad):
        sid = sample_ids[i]
        print(f"\n{'='*60}")
        print(f"[{i+1}/{len(all_h5ad)}] {sid}")
        print(f"{'='*60}")

        try:
            result = process_one_sample(h5ad_path)
            results.append(result)
        except Exception as e:
            print(f"\n  ERROR processing {sid}: {e}")
            traceback.print_exc()
            results.append({
                "sample_id": sid,
                "status": "error",
                "error": str(e),
            })

    # Summary
    print(f"\n\n{'='*80}")
    print("BATCH SUMMARY")
    print(f"{'='*80}")

    success = [r for r in results if r.get("status") == "success"]
    failed = [r for r in results if r.get("status") != "success"]

    if success:
        df = pd.DataFrame(success)
        cols = ["sample_id", "diagnosis", "n_qc_pass", "n_clusters",
                "n_cortical", "n_vascular", "n_wm", "n_l1_border",
                "n_complete_strips", "n_partial_strips",
                "n_cells_in_strips", "strip_coverage_pct", "time_sec"]
        print(df[cols].to_string(index=False))
        print(f"\nSucceeded: {len(success)}/{len(results)}")

        # Aggregates
        total_cortical = df["n_cortical"].sum()
        total_vascular = df["n_vascular"].sum()
        total_wm = df["n_wm"].sum()
        total_l1 = df["n_l1_border"].sum()
        total_cells = df["n_qc_pass"].sum()
        mean_coverage = df["strip_coverage_pct"].mean()
        total_complete = df["n_complete_strips"].sum()
        total_partial = df["n_partial_strips"].sum()

        print(f"\nAggregates across {len(success)} samples:")
        print(f"  Total QC-pass cells: {total_cells:,}")
        print(f"  Cortical: {total_cortical:,} ({total_cortical/total_cells*100:.1f}%)")
        print(f"    of which L1 border: {total_l1:,} ({total_l1/total_cells*100:.1f}%)")
        print(f"  Vascular: {total_vascular:,} ({total_vascular/total_cells*100:.1f}%)")
        print(f"  WM: {total_wm:,} ({total_wm/total_cells*100:.1f}%)")
        print(f"  Complete strips: {total_complete}, Partial strips: {total_partial}")
        print(f"  Mean strip coverage: {mean_coverage:.1f}%")

        csv_path = os.path.join(OUT_DIR, "batch_summary.csv")
        df.to_csv(csv_path, index=False)
        print(f"\nSaved: {csv_path}")

    if failed:
        print(f"\nFailed ({len(failed)}):")
        for r in failed:
            print(f"  {r['sample_id']}: {r.get('error', 'unknown')}")

    total_time = time.time() - t_start
    print(f"\nTotal time: {total_time:.0f}s ({total_time/60:.1f} min)")


if __name__ == "__main__":
    main()
