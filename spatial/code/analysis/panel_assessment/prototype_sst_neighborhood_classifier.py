#!/usr/bin/env python3
"""
Prototype: Neighborhood-enhanced Sst supertype classification.

The ~300 gene Xenium v1 panel contains zero within-Sst supertype-discriminating
markers. SCZ shifts 54-112 panel genes per Sst subtype, causing diagnosis-dependent
misclassification (e.g., Sst_3 → Sst_20 confusion). This script tests whether
spatial neighborhood information can stabilize classification.

Two methods:
  A) Expression smoothing: average each Sst cell's expression with its K nearest
     Sst neighbors, then re-correlate against centroids
  B) Label majority voting: smooth supertype labels using spatial neighbors

Validation:
  - Reclassification should be asymmetric by diagnosis (more SCZ changes)
  - Margins should increase after smoothing
  - Sst_20 proportion should shift toward snRNAseq expectation (decreased in SCZ)

Usage:
    python3 -u prototype_sst_neighborhood_classifier.py
"""

import os
import sys
import glob
import time
import numpy as np
import pandas as pd
from scipy.sparse import issparse
from sklearn.neighbors import NearestNeighbors
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    BASE_DIR, H5AD_DIR, EXCLUDE_SAMPLES, SST_TYPES, SST_COLORS,
    REPRESENTATIVE_SAMPLES, load_sample_adata,
)

sys.path.insert(0, os.path.join(BASE_DIR, "code"))
from modules.metadata import get_subject_info
from modules.correlation_classifier import correlate, assign_labels

METADATA_PATH = os.path.join(BASE_DIR, "data", "sample_metadata.xlsx")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "density_analysis")


def normalize_expression(X):
    """Counts-per-10K + log1p normalization (matches correlation classifier)."""
    if issparse(X):
        X = X.toarray()
    X = X.astype(np.float64)
    total = X.sum(axis=1, keepdims=True)
    total[total == 0] = 1
    return np.log1p(X * 10000 / total)


def find_sst_neighbors(all_coords, all_is_sst, sst_indices, K=30):
    """Find K nearest Sst neighbors for each Sst cell within a sample.

    Builds KNN on ALL cells, then for each Sst cell filters to Sst-only
    neighbors. Uses a larger initial K to ensure enough Sst neighbors.

    Returns
    -------
    neighbor_idx : np.ndarray (n_sst, K)
        Indices into sst_indices array (not global indices).
        -1 for missing neighbors.
    n_found : np.ndarray (n_sst,)
        Number of Sst neighbors found per cell.
    """
    n_all = len(all_coords)
    n_sst = len(sst_indices)

    # Map global index → sst local index
    global_to_sst = np.full(n_all, -1, dtype=np.int32)
    global_to_sst[sst_indices] = np.arange(n_sst)

    # Use larger K to compensate for non-Sst neighbors being filtered out
    # Sst cells are ~5-10% of cortical cells, so we need ~10x K to reliably
    # find K Sst neighbors
    k_search = min(K * 15, n_all)

    nn = NearestNeighbors(n_neighbors=k_search, algorithm='ball_tree')
    nn.fit(all_coords)
    _, nn_global = nn.kneighbors(all_coords[sst_indices])

    # For each Sst cell, filter to Sst neighbors (exclude self)
    neighbor_idx = np.full((n_sst, K), -1, dtype=np.int32)
    n_found = np.zeros(n_sst, dtype=np.int32)

    for i in range(n_sst):
        global_self = sst_indices[i]
        sst_neighs = []
        for g_idx in nn_global[i]:
            if g_idx == global_self:
                continue
            sst_local = global_to_sst[g_idx]
            if sst_local >= 0:
                sst_neighs.append(sst_local)
                if len(sst_neighs) >= K:
                    break
        n_found[i] = len(sst_neighs)
        neighbor_idx[i, :len(sst_neighs)] = sst_neighs

    return neighbor_idx, n_found


def smooth_expression(X_norm, neighbor_idx, n_found, alpha=0.5, min_neighbors=5):
    """Smooth expression: alpha * own + (1-alpha) * mean(Sst neighbors).

    Cells with fewer than min_neighbors Sst neighbors keep their own expression.
    """
    n_cells, n_genes = X_norm.shape
    X_smooth = X_norm.copy()

    for i in range(n_cells):
        if n_found[i] < min_neighbors:
            continue  # keep original
        neighs = neighbor_idx[i, :n_found[i]]
        neigh_mean = X_norm[neighs].mean(axis=0)
        X_smooth[i] = alpha * X_norm[i] + (1 - alpha) * neigh_mean

    return X_smooth


def label_voting(original_labels, original_margins, neighbor_idx, n_found,
                 min_neighbors=5):
    """Majority vote among Sst neighbors' supertype labels.

    Only reclassifies cells with below-median margin for their type.
    """
    n_cells = len(original_labels)
    new_labels = original_labels.copy()

    # Compute per-type median margin
    type_median_margin = {}
    for st in np.unique(original_labels):
        mask = original_labels == st
        type_median_margin[st] = np.median(original_margins[mask])

    n_reclassified = 0
    for i in range(n_cells):
        if n_found[i] < min_neighbors:
            continue
        # Only consider low-confidence cells
        own_type = original_labels[i]
        if original_margins[i] >= type_median_margin.get(own_type, 0):
            continue

        neighs = neighbor_idx[i, :n_found[i]]
        neigh_labels = original_labels[neighs]
        vote_counts = Counter(neigh_labels)
        winner = vote_counts.most_common(1)[0][0]
        if winner != own_type:
            new_labels[i] = winner
            n_reclassified += 1

    return new_labels, n_reclassified


def build_control_centroids(X_norm, sup_labels, diagnoses, gene_names):
    """Build Sst supertype centroids from Control cells only (avoid circularity)."""
    ctrl_mask = diagnoses == "Control"
    unique_sups = sorted(set(sup_labels))
    centroids = {}
    for sup in unique_sups:
        mask = (sup_labels == sup) & ctrl_mask
        if mask.sum() >= 10:
            centroids[sup] = X_norm[mask].mean(axis=0)
        else:
            # Fall back to all cells if too few controls
            mask_all = sup_labels == sup
            centroids[sup] = X_norm[mask_all].mean(axis=0)
    return pd.DataFrame(centroids, index=gene_names).T


def main():
    t0 = time.time()
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    meta = get_subject_info(METADATA_PATH).set_index("sample_id")

    # ── Load all samples, keeping Sst cells with spatial coords ──
    print("Loading samples...")
    h5ad_files = sorted(glob.glob(os.path.join(H5AD_DIR, "*_annotated.h5ad")))

    # Store per-sample data for KNN (need spatial coords within each sample)
    sample_data = {}  # sid -> {adata_sst, sst_indices_in_sample, all_coords}

    for fpath in h5ad_files:
        sid = os.path.basename(fpath).replace("_annotated.h5ad", "")
        if sid in EXCLUDE_SAMPLES:
            continue

        adata = load_sample_adata(sid, cortical_only=True, qc_mode="corr")
        sst_mask = adata.obs["subclass_label"] == "Sst"
        n_sst = sst_mask.sum()
        if n_sst < 10:
            continue

        dx = meta.loc[sid, "diagnosis"]
        all_coords = adata.obsm["spatial"][:, :2].copy()
        sst_indices = np.where(sst_mask)[0]

        adata_sst = adata[sst_mask].copy()
        adata_sst.obs["diagnosis"] = dx

        sample_data[sid] = {
            "adata_sst": adata_sst,
            "sst_indices": sst_indices,
            "all_coords": all_coords,
        }
        print(f"  {sid}: {n_sst:,} Sst cells ({dx})")

    # ── Concatenate Sst cells, normalize ──
    import anndata as ad
    all_sst = [sd["adata_sst"] for sd in sample_data.values()]
    adata_all = ad.concat(all_sst, merge="same")
    print(f"\nTotal: {adata_all.n_obs:,} Sst cells from {len(sample_data)} samples")
    print(f"Supertypes: {adata_all.obs['supertype_label'].value_counts().to_dict()}")

    X_norm = normalize_expression(adata_all.X)
    gene_names = list(adata_all.var_names)
    sup_labels = adata_all.obs["supertype_label"].values.astype(str)
    diagnoses = adata_all.obs["diagnosis"].values.astype(str)
    sample_ids = adata_all.obs["sample_id"].values.astype(str)

    # ── Build centroids from Control cells ──
    print("\nBuilding Sst supertype centroids (Control-only)...")
    centroids = build_control_centroids(X_norm, sup_labels, diagnoses, gene_names)
    print(f"  {len(centroids)} supertype centroids")

    # ── Original classification (baseline) ──
    print("\nBaseline: re-correlating against Control centroids...")
    corr_matrix_orig, type_names = correlate(X_norm, centroids)
    labels_orig, corr_orig, margin_orig = assign_labels(corr_matrix_orig, type_names)

    # ── Build per-sample KNN and smooth ──
    # We need to map cells back to their sample for spatial KNN
    K_DEFAULT = 30
    ALPHA_DEFAULT = 0.5

    # Build sample offsets
    sample_order = list(sample_data.keys())
    sample_sizes = [sample_data[sid]["adata_sst"].n_obs for sid in sample_order]
    sample_offsets = np.cumsum([0] + sample_sizes)

    print(f"\nBuilding per-sample KNN graphs (K={K_DEFAULT})...")
    # Pre-allocate neighbor arrays for all cells
    n_total = adata_all.n_obs
    all_neighbor_idx = np.full((n_total, K_DEFAULT), -1, dtype=np.int32)
    all_n_found = np.zeros(n_total, dtype=np.int32)

    for i, sid in enumerate(sample_order):
        sd = sample_data[sid]
        offset = sample_offsets[i]
        n_sst_sample = sample_sizes[i]

        # Find Sst neighbors within this sample
        neigh_idx, n_found = find_sst_neighbors(
            sd["all_coords"], None, sd["sst_indices"], K=K_DEFAULT
        )

        # Shift local Sst indices to global (offset within concatenated array)
        valid_mask = neigh_idx >= 0
        neigh_idx_global = neigh_idx.copy()
        neigh_idx_global[valid_mask] += offset
        neigh_idx_global[~valid_mask] = -1

        all_neighbor_idx[offset:offset + n_sst_sample] = neigh_idx_global
        all_n_found[offset:offset + n_sst_sample] = n_found

        median_found = np.median(n_found)
        print(f"  {sid}: median {median_found:.0f} Sst neighbors found")

    # ── Method A: Expression smoothing ──
    print(f"\nMethod A: Smoothing expression (alpha={ALPHA_DEFAULT})...")
    X_smooth = smooth_expression(X_norm, all_neighbor_idx, all_n_found,
                                  alpha=ALPHA_DEFAULT)

    corr_matrix_smooth, _ = correlate(X_smooth, centroids)
    labels_smooth, corr_smooth, margin_smooth = assign_labels(
        corr_matrix_smooth, type_names
    )

    n_changed_a = (labels_smooth != labels_orig).sum()
    pct_changed_a = 100 * n_changed_a / len(labels_orig)
    print(f"  Reclassified: {n_changed_a:,} / {len(labels_orig):,} "
          f"({pct_changed_a:.1f}%)")

    # ── Method B: Label voting ──
    print("\nMethod B: Label majority voting...")
    labels_vote, n_reclass_b = label_voting(
        labels_orig, margin_orig, all_neighbor_idx, all_n_found
    )
    pct_changed_b = 100 * (labels_vote != labels_orig).sum() / len(labels_orig)
    print(f"  Reclassified: {(labels_vote != labels_orig).sum():,} / "
          f"{len(labels_orig):,} ({pct_changed_b:.1f}%)")

    # ══════════════════════════════════════════════════════════════════
    # VALIDATION OUTPUTS
    # ══════════════════════════════════════════════════════════════════

    # ── 1. Reclassification matrix by diagnosis ──
    print("\n" + "=" * 70)
    print("Reclassification analysis (Method A: expression smoothing)")
    print("=" * 70)

    for dx in ["Control", "SCZ"]:
        dx_mask = diagnoses == dx
        changed = labels_smooth[dx_mask] != labels_orig[dx_mask]
        n_dx = dx_mask.sum()
        n_changed_dx = changed.sum()
        print(f"\n  {dx}: {n_changed_dx:,} / {n_dx:,} reclassified "
              f"({100*n_changed_dx/n_dx:.1f}%)")

        # Show flow for key Sst types
        for st in SST_TYPES:
            st_mask = dx_mask & (labels_orig == st)
            if st_mask.sum() == 0:
                continue
            new_labels_st = labels_smooth[st_mask]
            changed_st = new_labels_st != st
            if changed_st.sum() > 0:
                flows = Counter(new_labels_st[changed_st])
                flow_str = ", ".join(f"→{k}: {v}" for k, v in
                                     flows.most_common(5))
                print(f"    {st} (n={st_mask.sum()}): "
                      f"{changed_st.sum()} left ({flow_str})")

    # ── 2. Margin comparison ──
    print("\n" + "=" * 70)
    print("Margin comparison: original vs smoothed")
    print("=" * 70)

    for dx in ["Control", "SCZ"]:
        dx_mask = diagnoses == dx
        print(f"\n  {dx}:")
        print(f"    Original margin: {margin_orig[dx_mask].mean():.4f} "
              f"± {margin_orig[dx_mask].std():.4f}")
        print(f"    Smoothed margin: {margin_smooth[dx_mask].mean():.4f} "
              f"± {margin_smooth[dx_mask].std():.4f}")
        print(f"    Improvement: {(margin_smooth[dx_mask] - margin_orig[dx_mask]).mean():+.4f}")

    # Per Sst type
    print("\n  Per supertype:")
    for st in SST_TYPES:
        for dx in ["Control", "SCZ"]:
            mask = (sup_labels == st) & (diagnoses == dx)
            if mask.sum() < 10:
                continue
            orig_m = margin_orig[mask].mean()
            smooth_m = margin_smooth[mask].mean()
            print(f"    {st} ({dx}): {orig_m:.4f} → {smooth_m:.4f} "
                  f"(Δ={smooth_m - orig_m:+.4f})")

    # ── 3. Proportion changes ──
    print("\n" + "=" * 70)
    print("Proportion changes: original vs smoothed")
    print("=" * 70)

    prop_records = []
    for sid in sample_order:
        i = sample_order.index(sid)
        offset = sample_offsets[i]
        end = sample_offsets[i + 1]
        dx = meta.loc[sid, "diagnosis"]
        n_sst = end - offset

        for method_name, method_labels in [("original", labels_orig),
                                             ("smoothed_A", labels_smooth),
                                             ("voted_B", labels_vote)]:
            sample_labels = method_labels[offset:end]
            counts = Counter(sample_labels)
            for st in SST_TYPES:
                prop_records.append({
                    "sample_id": sid,
                    "diagnosis": dx,
                    "method": method_name,
                    "supertype": st,
                    "count": counts.get(st, 0),
                    "n_sst": n_sst,
                    "proportion": counts.get(st, 0) / n_sst,
                })

    prop_df = pd.DataFrame(prop_records)
    prop_path = os.path.join(OUTPUT_DIR, "sst_neighborhood_proportions.csv")
    prop_df.to_csv(prop_path, index=False)
    print(f"\nSaved: {prop_path}")

    # Summary: mean proportion by diagnosis and method
    print("\n  Mean Sst supertype proportions (within Sst cells):")
    for st in SST_TYPES:
        print(f"\n  {st}:")
        for method in ["original", "smoothed_A", "voted_B"]:
            for dx in ["Control", "SCZ"]:
                vals = prop_df[(prop_df["supertype"] == st) &
                               (prop_df["method"] == method) &
                               (prop_df["diagnosis"] == dx)]["proportion"]
                print(f"    {method:15s} {dx:8s}: {vals.mean():.4f} "
                      f"± {vals.std():.4f}")

    # ── 4. Direction of SCZ effect ──
    print("\n" + "=" * 70)
    print("SCZ effect direction: does smoothing change it?")
    print("=" * 70)

    for st in SST_TYPES:
        print(f"\n  {st}:")
        for method in ["original", "smoothed_A", "voted_B"]:
            ctrl = prop_df[(prop_df["supertype"] == st) &
                           (prop_df["method"] == method) &
                           (prop_df["diagnosis"] == "Control")]["proportion"]
            scz = prop_df[(prop_df["supertype"] == st) &
                          (prop_df["method"] == method) &
                          (prop_df["diagnosis"] == "SCZ")]["proportion"]
            diff = scz.mean() - ctrl.mean()
            direction = "↑SCZ" if diff > 0 else "↓SCZ"
            print(f"    {method:15s}: ctrl={ctrl.mean():.4f} "
                  f"scz={scz.mean():.4f} diff={diff:+.4f} {direction}")

    # ── 5. Parameter sweep ──
    print("\n" + "=" * 70)
    print("Parameter sweep")
    print("=" * 70)

    sweep_results = []
    for K in [15, 30, 50]:
        # Rebuild KNN if K differs from default
        if K != K_DEFAULT:
            print(f"\n  Rebuilding KNN with K={K}...")
            sweep_neigh = np.full((n_total, K), -1, dtype=np.int32)
            sweep_nfound = np.zeros(n_total, dtype=np.int32)
            for i, sid in enumerate(sample_order):
                sd = sample_data[sid]
                offset = sample_offsets[i]
                n_sst_s = sample_sizes[i]
                ni, nf = find_sst_neighbors(
                    sd["all_coords"], None, sd["sst_indices"], K=K
                )
                valid = ni >= 0
                ni_g = ni.copy()
                ni_g[valid] += offset
                ni_g[~valid] = -1
                sweep_neigh[offset:offset + n_sst_s] = ni_g
                sweep_nfound[offset:offset + n_sst_s] = nf
        else:
            sweep_neigh = all_neighbor_idx
            sweep_nfound = all_n_found

        for alpha in [0.3, 0.5, 0.7]:
            X_sw = smooth_expression(X_norm, sweep_neigh, sweep_nfound,
                                      alpha=alpha)
            corr_sw, _ = correlate(X_sw, centroids)
            labels_sw, _, margin_sw = assign_labels(corr_sw, type_names)

            n_changed = (labels_sw != labels_orig).sum()
            pct_changed = 100 * n_changed / len(labels_orig)
            mean_margin = margin_sw.mean()

            # Proportion shift for Sst_20
            sst20_ctrl_orig = []
            sst20_scz_orig = []
            sst20_ctrl_sw = []
            sst20_scz_sw = []
            for j, sid in enumerate(sample_order):
                off = sample_offsets[j]
                end = sample_offsets[j + 1]
                dx = meta.loc[sid, "diagnosis"]
                n_s = end - off
                orig_frac = (labels_orig[off:end] == "Sst_20").sum() / n_s
                sw_frac = (labels_sw[off:end] == "Sst_20").sum() / n_s
                if dx == "Control":
                    sst20_ctrl_orig.append(orig_frac)
                    sst20_ctrl_sw.append(sw_frac)
                else:
                    sst20_scz_orig.append(orig_frac)
                    sst20_scz_sw.append(sw_frac)

            sst20_diff_orig = np.mean(sst20_scz_orig) - np.mean(sst20_ctrl_orig)
            sst20_diff_sw = np.mean(sst20_scz_sw) - np.mean(sst20_ctrl_sw)

            sweep_results.append({
                "K": K, "alpha": alpha,
                "pct_reclassified": pct_changed,
                "mean_margin": mean_margin,
                "sst20_scz_diff_original": sst20_diff_orig,
                "sst20_scz_diff_smoothed": sst20_diff_sw,
            })

            print(f"  K={K}, alpha={alpha}: {pct_changed:.1f}% reclassified, "
                  f"margin={mean_margin:.4f}, "
                  f"Sst_20 SCZ diff: {sst20_diff_orig:+.4f} → {sst20_diff_sw:+.4f}")

    sweep_df = pd.DataFrame(sweep_results)
    sweep_path = os.path.join(OUTPUT_DIR, "sst_neighborhood_sweep.csv")
    sweep_df.to_csv(sweep_path, index=False)
    print(f"\nSaved: {sweep_path}")

    # ── 6. Save full cell-level results ──
    cell_results = pd.DataFrame({
        "sample_id": sample_ids,
        "diagnosis": diagnoses,
        "supertype_original": sup_labels,
        "supertype_recorr_baseline": labels_orig,
        "supertype_smoothed_A": labels_smooth,
        "supertype_voted_B": labels_vote,
        "margin_original": margin_orig,
        "margin_smoothed": margin_smooth,
        "n_sst_neighbors": all_n_found,
    })
    cell_path = os.path.join(OUTPUT_DIR, "sst_neighborhood_cell_results.csv")
    cell_results.to_csv(cell_path, index=False)
    print(f"Saved: {cell_path}")

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.0f}s")


if __name__ == "__main__":
    main()
