#!/usr/bin/env python3
"""
BANKSY Phase 2+3: Domain classification and cortical strip identification.

For each sample:
  1. Run BANKSY (λ=0.8, res=0.3) to get spatially coherent clusters
  2. Classify clusters as Cortical / Vascular / Meningeal / WM
  3. Estimate cortical surface orientation from depth gradient
  4. Identify cortical strips with all L1-L6 layers present
  5. Generate per-sample diagnostic figures

Usage:
    python3 -u code/analysis/banksy_02_domains_and_strips.py [--samples S1 S2 ...]
    python3 -u code/analysis/banksy_02_domains_and_strips.py --all
"""

import os
import sys
import time
import argparse
import warnings
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from scipy.ndimage import gaussian_filter
from scipy.stats import kendalltau

# BANKSY imports
from banksy.initialize_banksy import initialize_banksy
from banksy.embed_banksy import generate_banksy_matrix
from banksy.cluster_methods import run_Leiden_partition
from banksy_utils.umap_pca import pca_umap

# Project imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    H5AD_DIR, LAYER_COLORS, LAYER_ORDER, CORTICAL_LAYERS,
    SUBCLASS_TO_CLASS, CLASS_COLORS, EXCLUDE_SAMPLES,
    SAMPLE_TO_DX,
)

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "modules"))
from depth_model import LAYER_BINS
from spatial_domains import VASCULAR_TYPES, NON_NEURONAL_TYPES

# ── Constants ─────────────────────────────────────────────────────────

OUT_DIR = os.path.expanduser("~/Github/SCZ_Xenium/output/banksy")
os.makedirs(OUT_DIR, exist_ok=True)

# Best BANKSY parameters from Phase 1 pilot
BANKSY_LAMBDA = 0.8
BANKSY_RESOLUTION = 0.3
K_GEOM = 15
PCA_DIMS = [20]

# Domain classification thresholds
VASCULAR_THRESH = 0.50       # lower than old 0.80 — BANKSY clusters are spatially coherent
MENINGEAL_NN_THRESH = 0.50   # non-neuronal fraction for meningeal
MENINGEAL_DEPTH_THRESH = 0.20  # mean depth < 0.20 (meningeal/extra-cortical above pia)
WM_OLIGO_THRESH = 0.40       # oligodendrocyte fraction for WM detection
WM_DEPTH_THRESH = 0.80       # mean depth for WM

# Cortical strip parameters
STRIP_WIDTH_UM = 750         # strip width in micrometers (wider for L1 coverage)
MIN_CELLS_PER_LAYER = 15     # minimum cells per cortical layer in a strip
DEPTH_GRID_BINS = 50         # grid resolution for depth gradient estimation
DEPTH_SMOOTH_SIGMA = 2.0     # Gaussian smoothing for gradient

# Layer definitions for strip scoring
REQUIRED_LAYERS = ["L1", "L2/3", "L4", "L5", "L6"]
LAYER_DEPTH_ORDER = {"L1": 0, "L2/3": 1, "L4": 2, "L5": 3, "L6": 4}

# Layer boundary positions
LAYER_BOUNDARIES = [0.10, 0.30, 0.45, 0.65, 0.85]
LAYER_TICKS = [0.05, 0.20, 0.375, 0.55, 0.75, 0.925]
LAYER_NAMES = ["L1", "L2/3", "L4", "L5", "L6", "WM"]

# Domain colors
DOMAIN_COLORS = {
    "Cortical": "#4CAF50",
    "Vascular": "#F24C99",
    "Meningeal": "#FF9800",
    "WM": "#999999",
}


# ── Data Loading ──────────────────────────────────────────────────────

def load_sample(sample_id):
    """Load a single sample, subset to QC-pass cells."""
    path = os.path.join(H5AD_DIR, f"{sample_id}_annotated.h5ad")
    adata = ad.read_h5ad(path)

    if "hybrid_qc_pass" in adata.obs.columns:
        qc_mask = adata.obs["hybrid_qc_pass"].values.astype(bool)
    elif "qc_pass" in adata.obs.columns:
        qc_mask = adata.obs["qc_pass"].values.astype(bool)
    else:
        qc_mask = np.ones(adata.n_obs, dtype=bool)

    adata = adata[qc_mask].copy()
    return adata


def preprocess_for_banksy(adata):
    """Normalize, log-transform, and z-score for BANKSY."""
    adata_b = adata.copy()
    adata_b.layers["counts"] = adata_b.X.copy()
    sc.pp.normalize_total(adata_b, target_sum=1e4)
    sc.pp.log1p(adata_b)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.scale(adata_b)
    return adata_b


# ── BANKSY Clustering ─────────────────────────────────────────────────

def run_banksy(adata_b, lam=BANKSY_LAMBDA, res=BANKSY_RESOLUTION):
    """Run BANKSY with fixed parameters, return cluster labels."""
    coord_keys = ("x", "y", "spatial")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        banksy_dict = initialize_banksy(
            adata_b, coord_keys, K_GEOM,
            nbr_weight_decay="scaled_gaussian", max_m=1,
            plt_edge_hist=False, plt_nbr_weights=False, plt_agf_angles=False,
        )

        banksy_dict, _ = generate_banksy_matrix(
            adata_b, banksy_dict, [lam], max_m=1
        )

        pca_umap(banksy_dict, pca_dims=PCA_DIMS, add_umap=False)

        results_df, _ = run_Leiden_partition(
            banksy_dict, [res],
            num_nn=50, num_iterations=-1, partition_seed=42,
            match_labels=False,
        )

    # Extract labels
    row = results_df.iloc[0]
    labels = row["labels"]
    if hasattr(labels, "dense"):
        labels = labels.dense
    labels = np.asarray(labels).astype(int)

    return labels


# ── Phase 2: Domain Classification ────────────────────────────────────

def classify_banksy_domains(adata, banksy_labels):
    """Classify BANKSY clusters into domains based on composition and depth.

    Returns
    -------
    domains : np.ndarray of str
        Per-cell domain labels.
    cluster_info : dict
        Per-cluster statistics.
    """
    subclass = adata.obs.get("corr_subclass",
                              adata.obs.get("subclass_label")).values.astype(str)
    pred_depth = pd.to_numeric(adata.obs["predicted_norm_depth"],
                                errors="coerce").values
    coords = adata.obsm["spatial"]

    unique_cl = np.unique(banksy_labels)
    cluster_info = {}

    for cl in unique_cl:
        mask = banksy_labels == cl
        n_cl = mask.sum()
        cl_sub = subclass[mask]
        cl_depth = pred_depth[mask]
        cl_coords = coords[mask]

        # Cell type fractions
        sub_counts = {}
        for s in cl_sub:
            sub_counts[s] = sub_counts.get(s, 0) + 1

        vasc_frac = sum(sub_counts.get(v, 0) for v in VASCULAR_TYPES) / n_cl
        nn_frac = sum(sub_counts.get(v, 0) for v in NON_NEURONAL_TYPES) / n_cl
        oligo_frac = sub_counts.get("Oligodendrocyte", 0) / n_cl

        valid_depth = cl_depth[~np.isnan(cl_depth)]
        mean_depth = float(np.mean(valid_depth)) if len(valid_depth) > 0 else np.nan

        # Top 3 types
        top3 = sorted(sub_counts.items(), key=lambda x: -x[1])[:3]
        top3_str = ", ".join(f"{s}({100*c/n_cl:.0f}%)" for s, c in top3)

        # Neuronal fraction (for cortical vs non-cortical distinction)
        neuronal_frac = 1.0 - nn_frac

        # Classification (checked in priority order)
        if vasc_frac > VASCULAR_THRESH:
            domain = "Vascular"
        elif nn_frac > MENINGEAL_NN_THRESH and (np.isnan(mean_depth) or mean_depth < MENINGEAL_DEPTH_THRESH):
            domain = "Meningeal"
        elif oligo_frac > WM_OLIGO_THRESH and mean_depth > WM_DEPTH_THRESH:
            domain = "WM"
        elif neuronal_frac > 0.20 and not np.isnan(mean_depth) and 0.0 <= mean_depth <= 0.90:
            domain = "Cortical"
        elif mean_depth > WM_DEPTH_THRESH:
            domain = "WM"
        else:
            domain = "Cortical"  # default

        cluster_info[cl] = {
            "n_cells": n_cl,
            "domain": domain,
            "vasc_frac": vasc_frac,
            "nn_frac": nn_frac,
            "oligo_frac": oligo_frac,
            "neuronal_frac": neuronal_frac,
            "mean_depth": mean_depth,
            "top3": top3_str,
        }

    # Print per-cluster summary
    print(f"  {'Cl':>4} | {'N':>6} | {'Domain':<10} | {'Vasc%':>6} | {'NN%':>6} "
          f"| {'Oligo%':>6} | {'Depth':>6} | Top types")
    print(f"  {'-'*95}")
    for cl in sorted(cluster_info.keys()):
        info = cluster_info[cl]
        print(f"  {cl:>4} | {info['n_cells']:>6,} | {info['domain']:<10} "
              f"| {info['vasc_frac']*100:>5.1f}% | {info['nn_frac']*100:>5.1f}% "
              f"| {info['oligo_frac']*100:>5.1f}% | {info['mean_depth']:>6.3f} "
              f"| {info['top3']}")

    # Map to per-cell domains
    domains = np.array([cluster_info[cl]["domain"] for cl in banksy_labels])
    return domains, cluster_info


# ── Phase 3: Cortical Strip Identification ────────────────────────────

def estimate_cortical_orientation(coords, pred_depth, n_bins=DEPTH_GRID_BINS,
                                   sigma=DEPTH_SMOOTH_SIGMA):
    """Estimate the pia-to-WM axis direction from the depth gradient field.

    Returns
    -------
    angle : float
        Median gradient angle in radians (pia → WM direction).
    grad_x, grad_y : 2D arrays
        Gradient field components.
    x_edges, y_edges : 1D arrays
        Bin edges for the grid.
    mean_depth_grid : 2D array
        Smoothed mean depth grid.
    """
    valid = ~np.isnan(pred_depth)
    x, y = coords[valid, 0], coords[valid, 1]
    d = pred_depth[valid]

    # Build depth grid
    depth_sum, x_edges, y_edges = np.histogram2d(x, y, bins=n_bins, weights=d)
    count, _, _ = np.histogram2d(x, y, bins=n_bins)
    mean_depth_grid = np.where(count > 5, depth_sum / count, np.nan)

    # Fill NaN with nearest valid neighbor for smoothing
    from scipy.ndimage import generic_filter
    def fill_nan(values):
        center = values[len(values) // 2]
        if np.isnan(center):
            valid_v = values[~np.isnan(values)]
            return np.mean(valid_v) if len(valid_v) > 0 else np.nan
        return center

    filled = generic_filter(mean_depth_grid, fill_nan, size=5,
                             mode='constant', cval=np.nan)
    smoothed = gaussian_filter(np.nan_to_num(filled, nan=0.5), sigma=sigma)

    # Compute gradient (depth increases from pia toward WM)
    grad_y, grad_x = np.gradient(smoothed)

    # Get median gradient angle (where depth field is reliable)
    valid_mask = count > 10
    if valid_mask.sum() > 0:
        angles = np.arctan2(grad_y[valid_mask], grad_x[valid_mask])
        angle = np.median(angles)
    else:
        angle = np.pi / 2  # default: assume vertical

    return angle, grad_x, grad_y, x_edges, y_edges, smoothed


def identify_cortical_strips(adata, domains, banksy_labels, pred_depth,
                              angle, strip_width=STRIP_WIDTH_UM,
                              min_cells_per_layer=MIN_CELLS_PER_LAYER):
    """Identify cortical strips perpendicular to the depth gradient.

    Tiles the tissue into strips along the depth gradient axis, then scores
    each strip by layer completeness and ordering quality.

    Returns
    -------
    strip_ids : np.ndarray of int
        Strip ID per cell (-1 if not in any cortical strip).
    strip_scores : list of dict
        Quality metrics per strip.
    """
    coords = adata.obsm["spatial"]
    x, y = coords[:, 0], coords[:, 1]

    # Project coordinates onto axis PERPENDICULAR to depth gradient
    # (strips run parallel to cortical surface, perpendicular to depth axis)
    perp_angle = angle + np.pi / 2  # perpendicular to depth gradient
    proj = x * np.cos(perp_angle) + y * np.sin(perp_angle)

    # Tile into strips
    proj_min, proj_max = proj.min(), proj.max()
    n_strips = max(1, int((proj_max - proj_min) / strip_width))
    strip_edges = np.linspace(proj_min, proj_max, n_strips + 1)

    # Assign cells to depth-bin layers for strip scoring
    depth_layers = np.full(len(pred_depth), "Unknown", dtype=object)
    for layer_name, (lo, hi) in LAYER_BINS.items():
        mask = (pred_depth >= lo) & (pred_depth < hi)
        depth_layers[mask] = layer_name

    # Assign cells to strips
    strip_assignment = np.digitize(proj, strip_edges) - 1  # 0-indexed
    strip_assignment = np.clip(strip_assignment, 0, n_strips - 1)

    # Score each strip
    cortical_mask = domains == "Cortical"
    strip_scores = []

    for s in range(n_strips):
        in_strip = (strip_assignment == s) & cortical_mask
        n_total = in_strip.sum()

        if n_total < min_cells_per_layer * 3:
            strip_scores.append({
                "strip_id": s,
                "n_cells": n_total,
                "completeness": 0.0,
                "depth_range": 0.0,
                "order_score": 0.0,
                "purity": 0.0,
                "composite_score": 0.0,
                "layers_present": [],
            })
            continue

        # Count cells per required layer
        strip_layers = depth_layers[in_strip]
        strip_depths = pred_depth[in_strip]
        layer_counts = {}
        for lname in REQUIRED_LAYERS:
            layer_counts[lname] = (strip_layers == lname).sum()

        # Completeness: fraction of required layers with enough cells
        layers_present = [l for l in REQUIRED_LAYERS
                         if layer_counts.get(l, 0) >= min_cells_per_layer]
        completeness = len(layers_present) / len(REQUIRED_LAYERS)

        # Depth coverage (5th to 95th percentile)
        valid_d = strip_depths[~np.isnan(strip_depths)]
        if len(valid_d) > 10:
            depth_range = np.percentile(valid_d, 95) - np.percentile(valid_d, 5)
        else:
            depth_range = 0.0

        # Laminar order quality (Kendall tau)
        if len(layers_present) >= 3:
            expected_order = sorted(layers_present, key=lambda l: LAYER_DEPTH_ORDER[l])
            actual_medians = {l: np.nanmedian(strip_depths[strip_layers == l])
                             for l in layers_present
                             if (strip_layers == l).sum() >= min_cells_per_layer}
            if len(actual_medians) >= 3:
                actual_order = sorted(actual_medians.keys(),
                                     key=lambda l: actual_medians[l])
                expected_ranks = [expected_order.index(l) for l in actual_order]
                tau, _ = kendalltau(range(len(expected_ranks)), expected_ranks)
                order_score = max(0, tau)
            else:
                order_score = 0.0
        else:
            order_score = 0.0

        # Purity: fraction of cells in strip that are cortical
        all_in_strip = strip_assignment == s
        purity = cortical_mask[all_in_strip].sum() / max(1, all_in_strip.sum())

        composite = completeness * order_score * purity

        strip_scores.append({
            "strip_id": s,
            "n_cells": n_total,
            "completeness": completeness,
            "depth_range": depth_range,
            "order_score": order_score,
            "purity": purity,
            "composite_score": composite,
            "layers_present": layers_present,
            "layer_counts": layer_counts,
        })

    # Print strip scores for strips with any cells
    scored = [s for s in strip_scores if s["n_cells"] >= min_cells_per_layer * 3]
    if scored:
        print(f"  {'S':>4} | {'N':>6} | {'Cmpl':>5} | {'Order':>5} | {'Purity':>6} "
              f"| {'DepRng':>6} | {'Layers':<20} | Status")
        print(f"  {'-'*85}")
        for s in scored:
            layers_str = ",".join(s.get("layers_present", []))
            tier = _classify_strip_tier(s)
            status = ""
            if tier == 1:
                status = "COMPLETE"
            elif tier == 2:
                status = "PARTIAL"
            elif s["completeness"] >= 0.80:
                fails = []
                if s["completeness"] < 1.0:
                    fails.append(f"compl={s['completeness']:.2f}")
                if s["order_score"] <= 0.6:
                    fails.append(f"order={s['order_score']:.2f}")
                if s["purity"] <= 0.60:
                    fails.append(f"purity={s['purity']:.2f}")
                status = "near: " + ", ".join(fails) if fails else "near"
            print(f"  {s['strip_id']:>4} | {s['n_cells']:>6,} | {s['completeness']:>5.2f} "
                  f"| {s['order_score']:>5.2f} | {s['purity']*100:>5.1f}% "
                  f"| {s['depth_range']:>6.2f} | {layers_str:<20} | {status}")

    # Two-tier strip selection:
    # Tier 1 (COMPLETE): all L1-L6, good order, cortical-dominated
    # Tier 2 (PARTIAL): 4/5 layers (L1 often sparse), decent order
    complete_ids = set()
    partial_ids = set()
    for score in strip_scores:
        tier = _classify_strip_tier(score)
        if tier == 1:
            complete_ids.add(score["strip_id"])
        elif tier == 2:
            partial_ids.add(score["strip_id"])

    selected_strip_ids = complete_ids | partial_ids

    # Build per-cell strip IDs and tier labels
    strip_ids = np.full(len(pred_depth), -1, dtype=int)
    strip_tiers = np.full(len(pred_depth), "", dtype=object)
    for i, s in enumerate(strip_assignment):
        if s in selected_strip_ids and cortical_mask[i]:
            strip_ids[i] = s
            strip_tiers[i] = "complete" if s in complete_ids else "partial"

    return strip_ids, strip_scores, strip_edges, perp_angle, complete_ids, partial_ids


def _classify_strip_tier(score):
    """Classify a strip score into tier 1 (complete), 2 (partial), or 0 (rejected)."""
    # Tier 1: All L1-L6 present, good laminar order, cortical-dominated
    if (score["completeness"] == 1.0 and
        score["order_score"] > 0.8 and
        score["purity"] > 0.75):
        return 1
    # Tier 2: 4/5 layers (L1 often sparse in tissue), decent order
    if (score["completeness"] >= 0.80 and
        score["order_score"] > 0.6 and
        score["purity"] > 0.60):
        return 2
    return 0


# ── Plotting ──────────────────────────────────────────────────────────

def plot_sample_summary(adata, sample_id, banksy_labels, domains,
                        strip_ids, strip_scores, strip_edges, perp_angle,
                        cluster_info, angle):
    """Generate a 2×3 diagnostic figure for one sample."""
    coords = adata.obsm["spatial"]
    x, y = coords[:, 0], coords[:, 1]
    pred_depth = pd.to_numeric(adata.obs["predicted_norm_depth"],
                                errors="coerce").values

    fig = plt.figure(figsize=(24, 16))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.25)
    dx = SAMPLE_TO_DX.get(sample_id, "?")

    def setup_spatial_ax(ax, title):
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_facecolor("#0a0a0a")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title, fontsize=16, fontweight="bold", color="white",
                     pad=8, backgroundcolor="#333333")

    # Panel 1: BANKSY clusters (colored by depth order)
    ax1 = fig.add_subplot(gs[0, 0])
    unique_cl = np.unique(banksy_labels)
    cl_depth = {cl: np.nanmedian(pred_depth[banksy_labels == cl]) for cl in unique_cl}
    ordered = sorted(unique_cl, key=lambda c: cl_depth[c])
    cmap = matplotlib.colormaps.get_cmap("Spectral_r").resampled(len(ordered))
    for i, cl in enumerate(ordered):
        mask = banksy_labels == cl
        color = cmap(i / max(1, len(ordered) - 1))
        ax1.scatter(x[mask], y[mask], c=[color], s=0.3, alpha=0.5, rasterized=True)
    setup_spatial_ax(ax1, f"BANKSY Clusters ({len(ordered)} clusters)")

    # Panel 2: Domain classification
    ax2 = fig.add_subplot(gs[0, 1])
    for domain, color in DOMAIN_COLORS.items():
        mask = domains == domain
        if mask.sum() > 0:
            ax2.scatter(x[mask], y[mask], c=[color], s=0.3, alpha=0.5,
                        rasterized=True, label=f"{domain} ({mask.sum():,})")
    setup_spatial_ax(ax2, "BANKSY Domain Classification")
    ax2.legend(loc="upper right", fontsize=11, markerscale=10, framealpha=0.8,
               facecolor="#333333", labelcolor="white")

    # Panel 3: Current depth-bin layers (reference)
    ax3 = fig.add_subplot(gs[0, 2])
    layers = adata.obs["layer"].values.astype(str)
    for layer_name, color in LAYER_COLORS.items():
        mask = layers == layer_name
        if mask.sum() > 0:
            ax3.scatter(x[mask], y[mask], c=[color], s=0.3, alpha=0.5,
                        rasterized=True)
    setup_spatial_ax(ax3, "Current Depth-Bin Layers")

    # Panel 4: Cortical strips (selected in green, others dim)
    ax4 = fig.add_subplot(gs[1, 0])
    in_strip = strip_ids >= 0
    not_in_strip = ~in_strip
    if not_in_strip.sum() > 0:
        ax4.scatter(x[not_in_strip], y[not_in_strip], c="#333333", s=0.2,
                    alpha=0.3, rasterized=True)
    if in_strip.sum() > 0:
        # Color by depth within strips
        depth_colors = plt.cm.viridis(np.clip(pred_depth[in_strip], 0, 1))
        ax4.scatter(x[in_strip], y[in_strip], c=depth_colors, s=0.5,
                    alpha=0.7, rasterized=True)

    # Draw strip boundaries
    proj = x * np.cos(perp_angle) + y * np.sin(perp_angle)
    for edge in strip_edges:
        # Convert projection back to line in x-y space
        x_line = np.array([x.min(), x.max()])
        y_line = (edge - x_line * np.cos(perp_angle)) / np.sin(perp_angle) if abs(np.sin(perp_angle)) > 0.01 else np.array([y.min(), y.max()])
        ax4.plot(x_line, y_line, color="#555555", linewidth=0.3, alpha=0.5)

    n_selected = len(set(strip_ids[in_strip]))
    n_in_strip = in_strip.sum()
    setup_spatial_ax(ax4, f"Cortical Strips ({n_selected} selected, {n_in_strip:,} cells)")

    # Panel 5: Strip quality scores
    ax5 = fig.add_subplot(gs[1, 1])
    completeness_vals = [s["completeness"] for s in strip_scores if s["n_cells"] > 0]
    order_vals = [s["order_score"] for s in strip_scores if s["n_cells"] > 0]
    composite_vals = [s["composite_score"] for s in strip_scores if s["n_cells"] > 0]

    bar_x = np.arange(len(completeness_vals))
    width = 0.25
    ax5.bar(bar_x - width, completeness_vals, width, label="Completeness",
            color="#4CAF50", alpha=0.8)
    ax5.bar(bar_x, order_vals, width, label="Order",
            color="#2196F3", alpha=0.8)
    ax5.bar(bar_x + width, composite_vals, width, label="Composite",
            color="#FF9800", alpha=0.8)
    ax5.axhline(0.8, color="red", linestyle="--", linewidth=1, alpha=0.5)
    ax5.set_xlabel("Strip Index", fontsize=14)
    ax5.set_ylabel("Score", fontsize=14)
    ax5.set_title("Strip Quality Scores", fontsize=16, fontweight="bold")
    ax5.legend(fontsize=12)
    ax5.set_ylim(0, 1.05)
    ax5.spines["top"].set_visible(False)
    ax5.spines["right"].set_visible(False)

    # Panel 6: Domain comparison summary (BANKSY vs current pipeline)
    ax6 = fig.add_subplot(gs[1, 2])
    current_domains = adata.obs["spatial_domain"].values.astype(str)
    banksy_domain_counts = {}
    current_domain_counts = {}
    for d in ["Cortical", "Vascular", "Meningeal", "WM", "Extra-cortical"]:
        banksy_domain_counts[d] = (domains == d).sum()
        current_domain_counts[d] = (current_domains == d).sum()

    categories = ["Cortical", "Vascular", "Meningeal", "WM", "Extra-cortical"]
    banksy_vals = [banksy_domain_counts.get(c, 0) for c in categories]
    current_vals = [current_domain_counts.get(c, 0) for c in categories]

    bar_x = np.arange(len(categories))
    ax6.barh(bar_x + 0.2, [v / max(1, sum(banksy_vals)) * 100 for v in banksy_vals],
             0.35, label="BANKSY", color="#4CAF50", alpha=0.8)
    ax6.barh(bar_x - 0.2, [v / max(1, sum(current_vals)) * 100 for v in current_vals],
             0.35, label="Current Pipeline", color="#2196F3", alpha=0.8)
    ax6.set_yticks(bar_x)
    ax6.set_yticklabels(categories, fontsize=13)
    ax6.set_xlabel("% of Cells", fontsize=14)
    ax6.set_title("Domain Classification Comparison", fontsize=16, fontweight="bold")
    ax6.legend(fontsize=12, loc="lower right")
    ax6.spines["top"].set_visible(False)
    ax6.spines["right"].set_visible(False)

    fig.suptitle(f"{sample_id} ({dx}) — BANKSY Domains & Cortical Strips",
                 fontsize=22, fontweight="bold")

    return fig


# ── Per-Sample Processing ─────────────────────────────────────────────

def process_sample(sample_id):
    """Full pipeline for one sample. Returns summary dict."""
    t0 = time.time()
    print(f"\n{'='*60}")
    print(f"Processing {sample_id} ({SAMPLE_TO_DX.get(sample_id, '?')})...")
    print(f"{'='*60}")

    # Load
    adata = load_sample(sample_id)
    print(f"  {adata.n_obs:,} QC-pass cells")

    # Preprocess
    adata_b = preprocess_for_banksy(adata)

    # Run BANKSY
    print(f"  Running BANKSY (λ={BANKSY_LAMBDA}, res={BANKSY_RESOLUTION})...")
    t_banksy = time.time()
    banksy_labels = run_banksy(adata_b, BANKSY_LAMBDA, BANKSY_RESOLUTION)
    n_clusters = len(np.unique(banksy_labels))
    print(f"  BANKSY: {n_clusters} clusters ({time.time()-t_banksy:.0f}s)")

    # Domain classification
    print("  Classifying domains...")
    domains, cluster_info = classify_banksy_domains(adata, banksy_labels)

    domain_counts = {}
    for d in ["Cortical", "Vascular", "Meningeal", "WM"]:
        domain_counts[d] = (domains == d).sum()
    print(f"  Domains: " + ", ".join(f"{d}={n:,}" for d, n in domain_counts.items()))

    # Depth gradient orientation
    pred_depth = pd.to_numeric(adata.obs["predicted_norm_depth"],
                                errors="coerce").values
    print("  Estimating cortical orientation...")
    angle, grad_x, grad_y, x_edges, y_edges, depth_grid = \
        estimate_cortical_orientation(adata.obsm["spatial"], pred_depth)
    print(f"  Depth gradient angle: {np.degrees(angle):.1f}°")

    # Cortical strips
    print("  Identifying cortical strips...")
    strip_ids, strip_scores, strip_edges, perp_angle, complete_ids, partial_ids = \
        identify_cortical_strips(adata, domains, banksy_labels, pred_depth,
                                  angle)
    n_complete = len(complete_ids)
    n_partial = len(partial_ids)
    n_in_strip = (strip_ids >= 0).sum()
    n_cortical = (domains == "Cortical").sum()
    strip_coverage = n_in_strip / max(1, n_cortical) * 100
    print(f"  Strips: {n_complete} complete + {n_partial} partial = "
          f"{n_complete + n_partial} total, {n_in_strip:,} cells "
          f"({strip_coverage:.1f}% of cortical)")

    # Generate diagnostic figure
    print("  Generating diagnostic figure...")
    fig = plot_sample_summary(adata, sample_id, banksy_labels, domains,
                               strip_ids, strip_scores, strip_edges,
                               perp_angle, cluster_info, angle)
    fig_path = os.path.join(OUT_DIR, f"domains_strips_{sample_id}.png")
    fig.savefig(fig_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {fig_path}")

    elapsed = time.time() - t0
    print(f"  Total: {elapsed:.0f}s")

    # Summary
    return {
        "sample_id": sample_id,
        "diagnosis": SAMPLE_TO_DX.get(sample_id, "?"),
        "n_cells": adata.n_obs,
        "n_clusters": n_clusters,
        "n_cortical": domain_counts.get("Cortical", 0),
        "n_vascular": domain_counts.get("Vascular", 0),
        "n_meningeal": domain_counts.get("Meningeal", 0),
        "n_wm": domain_counts.get("WM", 0),
        "pct_cortical": domain_counts.get("Cortical", 0) / adata.n_obs * 100,
        "pct_vascular": domain_counts.get("Vascular", 0) / adata.n_obs * 100,
        "n_complete_strips": n_complete,
        "n_partial_strips": n_partial,
        "n_total_strips": n_complete + n_partial,
        "n_cells_in_strips": n_in_strip,
        "strip_coverage_pct": strip_coverage,
        "gradient_angle_deg": np.degrees(angle),
        "time_sec": elapsed,
    }


# ── Main ──────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", nargs="+", default=None,
                        help="Sample IDs to process")
    parser.add_argument("--all", action="store_true",
                        help="Process all 24 samples")
    args = parser.parse_args()

    if args.all:
        h5ad_files = sorted(
            f.replace("_annotated.h5ad", "")
            for f in os.listdir(H5AD_DIR)
            if f.endswith("_annotated.h5ad")
        )
        sample_ids = h5ad_files
    elif args.samples:
        sample_ids = args.samples
    else:
        # Default: a few representative samples
        sample_ids = ["Br2039", "Br6389", "Br8433", "Br6437"]

    t_start = time.time()
    print(f"BANKSY Phase 2+3: Domains & Cortical Strips")
    print(f"Samples: {sample_ids}")
    print(f"BANKSY: λ={BANKSY_LAMBDA}, res={BANKSY_RESOLUTION}")
    print(f"Strips: width={STRIP_WIDTH_UM}μm, min_cells/layer={MIN_CELLS_PER_LAYER}")

    results = []
    for sample_id in sample_ids:
        try:
            result = process_sample(sample_id)
            results.append(result)
        except Exception as e:
            import traceback
            print(f"\n  ERROR processing {sample_id}: {e}")
            traceback.print_exc()

    # Summary table
    if results:
        df = pd.DataFrame(results)
        print(f"\n{'='*60}")
        print("SUMMARY")
        print(f"{'='*60}")
        print(df.to_string(index=False))

        csv_path = os.path.join(OUT_DIR, "domains_strips_summary.csv")
        df.to_csv(csv_path, index=False)
        print(f"\nSaved: {csv_path}")

    total_time = time.time() - t_start
    print(f"\nTotal time: {total_time:.0f}s ({total_time/60:.1f} min)")


if __name__ == "__main__":
    main()
