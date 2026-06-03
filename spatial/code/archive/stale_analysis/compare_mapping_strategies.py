#!/usr/bin/env python3
"""
Compare MapMyCells mapping strategies on a single Xenium sample (Br5400).

Runs 4 strategies, applies QC filtering based on mapping confidence,
propagates labels through the depth model, and compares depth consistency
with the MERFISH reference.

Strategies:
  1. Hierarchical (current baseline) — re-run with 300-gene mapping
  2. Hierarchical + more markers (n_per_utility=50)
  3. HANN algorithm (full-tree correlation walk per bootstrap iteration)
  4. Flat mapping (flatten=True, no hierarchical tree traversal)

Usage:
    /Users/shreejoy/venv/bin/python3 -u compare_mapping_strategies.py
"""

import os
import sys
import time
import json
import tempfile
import pickle
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

# Pipeline config
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "pipeline"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "modules"))
from pipeline_config import (
    H5AD_DIR, OUTPUT_DIR, PRECOMPUTED_STATS_PATH, GENE_MAPPING_PATH,
    TAXONOMY_LEVELS,
)
from depth_model import (
    load_model, predict_depth, build_neighborhood_features,
    assign_layers_with_ood, LAYER_BINS,
)

BASE_DIR = os.path.expanduser("~/Github/SCZ_Xenium")
MERFISH_PATH = os.path.join(BASE_DIR, "data", "reference",
                             "SEAAD_MTG_MERFISH.2024-12-11.h5ad")
OUT_DIR = os.path.join(BASE_DIR, "output", "plots")
DEPTH_MODEL_PATH = os.path.join(OUTPUT_DIR, "depth_model_normalized.pkl")

SAMPLE_ID = "Br5400"

# Subclass confidence thresholds for QC filtering
CONF_THRESHOLDS = [0.28, 0.50]  # current pipeline + Allen Institute default

# ──────────────────────────────────────────────────────────────────────
# Strategy definitions
# ──────────────────────────────────────────────────────────────────────

STRATEGIES = {
    "S1_hierarchical": {
        "label": "Hierarchical (baseline)",
        "config_overrides": {},
    },
    "S2_hierarchical_50markers": {
        "label": "Hierarchical (50 markers/node)",
        "config_overrides": {
            "query_markers": {"n_per_utility": 50},
        },
    },
    "S3_hann": {
        "label": "HANN (full-tree per iteration)",
        "config_overrides": {
            "type_assignment": {"algorithm": "hann"},
        },
    },
    "S4_flat": {
        "label": "Flat (no hierarchy)",
        "config_overrides": {
            "flatten": True,
        },
    },
}


# ──────────────────────────────────────────────────────────────────────
# Gene mapping
# ──────────────────────────────────────────────────────────────────────

def load_gene_mapping():
    with open(GENE_MAPPING_PATH) as f:
        return json.load(f)


def convert_genes_to_ensembl(adata, gene_mapping):
    mappable = [g for g in adata.var_names if g in gene_mapping]
    adata_sub = adata[:, mappable].copy()
    adata_sub.var_names = [gene_mapping[g] for g in mappable]
    adata_sub.var_names_make_unique()
    return adata_sub


# ──────────────────────────────────────────────────────────────────────
# MapMyCells runner
# ──────────────────────────────────────────────────────────────────────

def _load_taxonomy_tree():
    """Load the taxonomy tree from precomputed stats for HANN parsing."""
    import h5py
    with h5py.File(PRECOMPUTED_STATS_PATH, 'r') as f:
        tree = json.loads(f['taxonomy_tree'][()].decode())
    return tree


def _build_taxonomy_lookups(tree):
    """
    Build lookups from the taxonomy tree:
      - supertype_id → supertype_name
      - supertype_id → subclass_name
      - supertype_id → class_name

    The tree structure has:
      CCN..._CLAS: {class_id: [list of subclass_ids]}
      CCN..._SUBC: {subclass_id: [list of supertype_ids]}
      CCN..._SUPT: {supertype_id: [list of cell barcodes]}
    """
    hierarchy = tree['hierarchy']  # ['CCN...CLAS', 'CCN...SUBC', 'CCN...SUPT']
    nm = tree['name_mapper']

    clas_level = hierarchy[0]
    subc_level = hierarchy[1]
    supt_level = hierarchy[2]

    clas_data = tree[clas_level]  # {class_id: [subclass_ids]}
    subc_data = tree[subc_level]  # {subclass_id: [supertype_ids]}

    supt_to_name = {sid: info['name']
                     for sid, info in nm[supt_level].items()}
    supt_to_subclass = {}
    supt_to_class = {}

    for clas_id, subc_ids in clas_data.items():
        clas_name = nm[clas_level][clas_id]['name']
        for subc_id in subc_ids:
            subc_name = nm[subc_level][subc_id]['name']
            for supt_id in subc_data.get(subc_id, []):
                supt_to_subclass[supt_id] = subc_name
                supt_to_class[supt_id] = clas_name

    return supt_to_name, supt_to_subclass, supt_to_class


def run_mapmycells(query_h5ad_path, output_dir, strategy_key,
                   config_overrides, n_processors=1):
    """
    Run MapMyCells with the given strategy configuration.

    Returns output_path (csv for hierarchical, hdf5 for hann) and algorithm.
    """
    from cell_type_mapper.cli.map_to_on_the_fly_markers import OnTheFlyMapper

    # Determine algorithm
    algorithm = config_overrides.get("type_assignment", {}).get(
        "algorithm", "hierarchical")

    if algorithm == "hann":
        # HANN requires HDF5 output only
        hdf5_path = os.path.join(output_dir, f"{strategy_key}_output.h5")
        config = {
            "query_path": query_h5ad_path,
            "hdf5_result_path": hdf5_path,
            "precomputed_stats": {
                "path": PRECOMPUTED_STATS_PATH,
            },
            "type_assignment": {
                "normalization": "raw",
                "bootstrap_iteration": 100,
                "bootstrap_factor": 0.5,
                "algorithm": "hann",
            },
            "n_processors": n_processors,
            "query_markers": {
                "n_per_utility": 30,
            },
            "reference_markers": {
                "precomputed_path_list": None,
            },
        }
        # Apply non-algorithm overrides
        for key, val in config_overrides.items():
            if key == "type_assignment":
                for k2, v2 in val.items():
                    config["type_assignment"][k2] = v2
            elif isinstance(val, dict) and key in config:
                config[key].update(val)
            else:
                config[key] = val

        print(f"    Config: algorithm=hann, "
              f"n_per_utility={config['query_markers'].get('n_per_utility', 30)}")

        runner = OnTheFlyMapper(args=[], input_data=config)
        runner.run()

        return hdf5_path, "hann"
    else:
        # Hierarchical/flat: CSV + JSON output
        csv_path = os.path.join(output_dir, f"{strategy_key}_output.csv")
        json_path = os.path.join(output_dir, f"{strategy_key}_output.json")

        config = {
            "query_path": query_h5ad_path,
            "extended_result_path": json_path,
            "csv_result_path": csv_path,
            "precomputed_stats": {
                "path": PRECOMPUTED_STATS_PATH,
            },
            "type_assignment": {
                "normalization": "raw",
                "bootstrap_iteration": 100,
                "bootstrap_factor": 0.5,
                "algorithm": "hierarchical",
            },
            "n_processors": n_processors,
            "query_markers": {
                "n_per_utility": 30,
            },
            "reference_markers": {
                "precomputed_path_list": None,
            },
        }

        # Apply strategy-specific overrides
        for key, val in config_overrides.items():
            if isinstance(val, dict) and key in config:
                config[key].update(val)
            else:
                config[key] = val

        print(f"    Config: algorithm={config['type_assignment'].get('algorithm', 'hierarchical')}, "
              f"n_per_utility={config['query_markers'].get('n_per_utility', 30)}, "
              f"flatten={config.get('flatten', False)}")

        runner = OnTheFlyMapper(args=[], input_data=config)
        runner.run()

        return csv_path, "hierarchical"


def parse_mapmycells_output(output_path, algorithm):
    """Parse MapMyCells output (CSV for hierarchical, HDF5 for HANN)."""
    if algorithm == "hann":
        return _parse_hann_hdf5(output_path)
    else:
        return _parse_hierarchical_csv(output_path)


def _parse_hierarchical_csv(csv_path):
    """Parse MapMyCells CSV and extract labels + confidence."""
    df = pd.read_csv(csv_path, comment='#')

    result = pd.DataFrame(index=range(len(df)))

    for level in TAXONOMY_LEVELS:
        name_col = f"{level}_name"
        label_col = f"{level}_label"
        prob_col = f"{level}_bootstrapping_probability"

        if name_col in df.columns:
            result[f"{level}_label"] = df[name_col].values
        elif label_col in df.columns:
            result[f"{level}_label"] = df[label_col].values

        if prob_col in df.columns:
            result[f"{level}_label_confidence"] = df[prob_col].values.astype(np.float32)
        else:
            result[f"{level}_label_confidence"] = np.float32(0.0)

    # Handle older taxonomies: "cluster" → "supertype"
    if "supertype_label" not in result.columns:
        for alt in ["cluster_name", "cluster_label"]:
            if alt in df.columns:
                result["supertype_label"] = df[alt].values
                break
        if "cluster_bootstrapping_probability" in df.columns:
            result["supertype_label_confidence"] = \
                df["cluster_bootstrapping_probability"].values.astype(np.float32)

    return result


def _parse_hann_hdf5(hdf5_path):
    """
    Parse HANN algorithm HDF5 output.

    The HDF5 contains:
      - votes: (n_cells x n_clusters) vote counts across bootstrap iterations
      - correlation: (n_cells x n_clusters) average correlation per cluster
      - cell_identifiers: cell IDs (bytes)
      - cluster_identifiers: cluster IDs (bytes, supertype-level taxonomy IDs)

    We assign each cell to the cluster with the most votes, using correlation
    as a tiebreaker. Then map cluster IDs → human-readable names using the
    taxonomy tree.
    """
    import h5py

    tree = _load_taxonomy_tree()
    supt_to_name, supt_to_subclass, supt_to_class = _build_taxonomy_lookups(tree)

    with h5py.File(hdf5_path, 'r') as f:
        votes = f['votes'][:]          # (n_cells, n_clusters)
        correlation = f['correlation'][:]  # (n_cells, n_clusters)
        cluster_ids = [c.decode() if isinstance(c, bytes) else c
                       for c in f['cluster_identifiers'][:]]

    n_cells = votes.shape[0]
    n_boot = int(votes.sum(axis=1).max())  # total bootstrap iterations

    # For each cell, pick cluster with most votes (correlation as tiebreaker)
    # Add tiny correlation to break ties
    score = votes.astype(float) + correlation * 1e-6
    best_idx = np.argmax(score, axis=1)

    # Confidence = fraction of votes for winning cluster
    max_votes = votes[np.arange(n_cells), best_idx]
    total_votes = votes.sum(axis=1)
    confidence = np.where(total_votes > 0,
                          max_votes / total_votes,
                          0.0).astype(np.float32)

    # Map to human-readable names
    best_cluster_ids = [cluster_ids[i] for i in best_idx]

    result = pd.DataFrame(index=range(n_cells))

    # Supertype (leaf level)
    result["supertype_label"] = [supt_to_name.get(cid, cid)
                                  for cid in best_cluster_ids]
    result["supertype_label_confidence"] = confidence

    # Subclass (inferred from supertype via taxonomy)
    subclass_labels = [supt_to_subclass.get(cid, "Unknown")
                       for cid in best_cluster_ids]
    result["subclass_label"] = subclass_labels

    # Subclass confidence: sum votes for all supertypes in same subclass
    # Vectorized: build cluster → subclass mapping, then group and sum
    cluster_subclasses = [supt_to_subclass.get(cid, "") for cid in cluster_ids]
    unique_subclasses = sorted(set(cluster_subclasses))
    sub_to_idx_local = {s: i for i, s in enumerate(unique_subclasses)}

    # Build subclass vote matrix: (n_cells, n_subclasses)
    n_subclasses = len(unique_subclasses)
    sub_votes = np.zeros((n_cells, n_subclasses), dtype=np.float32)
    for j, sc in enumerate(cluster_subclasses):
        if sc in sub_to_idx_local:
            sub_votes[:, sub_to_idx_local[sc]] += votes[:, j]

    # Get confidence for each cell's assigned subclass
    subclass_confidence = np.zeros(n_cells, dtype=np.float32)
    for i in range(n_cells):
        sc = subclass_labels[i]
        if sc in sub_to_idx_local and total_votes[i] > 0:
            subclass_confidence[i] = sub_votes[i, sub_to_idx_local[sc]] / total_votes[i]
    result["subclass_label_confidence"] = subclass_confidence

    # Class (inferred from supertype via taxonomy)
    class_labels = [supt_to_class.get(cid, "Unknown")
                    for cid in best_cluster_ids]
    result["class_label"] = class_labels

    # Class confidence: vectorized
    cluster_classes = [supt_to_class.get(cid, "") for cid in cluster_ids]
    unique_classes = sorted(set(cluster_classes))
    cls_to_idx_local = {c: i for i, c in enumerate(unique_classes)}
    n_classes = len(unique_classes)
    cls_votes = np.zeros((n_cells, n_classes), dtype=np.float32)
    for j, cc in enumerate(cluster_classes):
        if cc in cls_to_idx_local:
            cls_votes[:, cls_to_idx_local[cc]] += votes[:, j]

    class_confidence = np.zeros(n_cells, dtype=np.float32)
    for i in range(n_cells):
        cc = class_labels[i]
        if cc in cls_to_idx_local and total_votes[i] > 0:
            class_confidence[i] = cls_votes[i, cls_to_idx_local[cc]] / total_votes[i]
    result["class_label_confidence"] = class_confidence

    print(f"    HANN: {n_cells} cells, {len(cluster_ids)} clusters, "
          f"{n_boot} bootstrap iterations")
    print(f"    Mean supertype confidence: {confidence.mean():.3f}")
    print(f"    Mean subclass confidence: {subclass_confidence.mean():.3f}")

    return result


# ──────────────────────────────────────────────────────────────────────
# MERFISH reference data
# ──────────────────────────────────────────────────────────────────────

def load_merfish_reference():
    """Load MERFISH reference: median depth per subclass and supertype."""
    print("Loading MERFISH reference...")
    adata = ad.read_h5ad(MERFISH_PATH, backed="r")

    obs = adata.obs[["Subclass", "Supertype", "Normalized depth from pia",
                      "Donor ID"]].copy()
    obs.columns = ["subclass", "supertype", "depth", "donor"]
    obs["depth"] = obs["depth"].astype(float)
    obs = obs.dropna(subset=["depth"])
    obs["subclass"] = obs["subclass"].astype(str)
    obs["supertype"] = obs["supertype"].astype(str)

    print(f"  MERFISH cells with depth: {len(obs):,}")

    # Per-subclass median depth
    sub_depth = obs.groupby("subclass")["depth"].median().to_dict()
    # Per-supertype median depth
    sup_depth = obs.groupby("supertype")["depth"].median().to_dict()
    # Per-subclass proportions (all cells, not just depth-annotated)
    all_obs = adata.obs[["Subclass", "Supertype"]].copy()
    all_obs.columns = ["subclass", "supertype"]
    all_obs["subclass"] = all_obs["subclass"].astype(str)
    all_obs["supertype"] = all_obs["supertype"].astype(str)
    total = len(all_obs)
    sub_props = (all_obs["subclass"].value_counts() / total).to_dict()
    sup_props = (all_obs["supertype"].value_counts() / total).to_dict()

    # L6b depth distribution for histogram comparison
    l6b_depths = obs[obs["subclass"] == "L6b"]["depth"].values

    return {
        "subclass_depth": sub_depth,
        "supertype_depth": sup_depth,
        "subclass_props": sub_props,
        "supertype_props": sup_props,
        "l6b_depths": l6b_depths,
        "n_total": total,
    }


# ──────────────────────────────────────────────────────────────────────
# Depth prediction
# ──────────────────────────────────────────────────────────────────────

def run_depth_prediction(adata, subclass_labels, model_bundle):
    """
    Run depth prediction using the given subclass labels.

    Creates a temporary copy of adata with the new labels and runs
    the depth model.

    Returns (predicted_depth, ood_scores, layers).
    """
    # Create a modified copy with the new labels
    adata_copy = adata.copy()
    adata_copy.obs["subclass_label"] = subclass_labels

    pred_depth, ood_scores = predict_depth(
        adata_copy, model_bundle,
        subclass_col="subclass_label", compute_ood=True
    )

    layers = assign_layers_with_ood(
        pred_depth, ood_scores,
        model_bundle=model_bundle
    )

    return pred_depth, ood_scores, layers


# ──────────────────────────────────────────────────────────────────────
# Metrics
# ──────────────────────────────────────────────────────────────────────

def compute_metrics(labels_df, pred_depth, layers, merfish_ref,
                    conf_thresh=0.28):
    """
    Compute comparison metrics between Xenium strategy and MERFISH.

    Returns a dict with all metrics.
    """
    # Apply confidence filter
    conf_col = "subclass_label_confidence"
    if conf_col in labels_df.columns:
        conf_mask = labels_df[conf_col].astype(float) >= conf_thresh
    else:
        conf_mask = np.ones(len(labels_df), dtype=bool)

    n_total = len(labels_df)
    n_pass = int(conf_mask.sum())
    n_fail = n_total - n_pass

    # Filtered labels and depths
    sub_labels = labels_df["subclass_label"].values[conf_mask]
    sup_labels = labels_df["supertype_label"].values[conf_mask]
    depths = pred_depth[conf_mask]

    # -- Subclass depth comparison --
    xen_sub_depth = {}
    for s in np.unique(sub_labels):
        mask = sub_labels == s
        if mask.sum() >= 5:
            xen_sub_depth[s] = float(np.median(depths[mask]))

    shared_sub = sorted(set(xen_sub_depth.keys()) &
                         set(merfish_ref["subclass_depth"].keys()))
    if len(shared_sub) >= 3:
        xen_vals = [xen_sub_depth[s] for s in shared_sub]
        mer_vals = [merfish_ref["subclass_depth"][s] for s in shared_sub]
        sub_r, _ = pearsonr(xen_vals, mer_vals)
        sub_rho, _ = spearmanr(xen_vals, mer_vals)
        sub_mad = float(np.mean(np.abs(np.array(xen_vals) - np.array(mer_vals))))
    else:
        sub_r = sub_rho = sub_mad = np.nan

    # -- Supertype depth comparison --
    xen_sup_depth = {}
    for s in np.unique(sup_labels):
        mask = sup_labels == s
        if mask.sum() >= 5:
            xen_sup_depth[s] = float(np.median(depths[mask]))

    shared_sup = sorted(set(xen_sup_depth.keys()) &
                         set(merfish_ref["supertype_depth"].keys()))
    if len(shared_sup) >= 3:
        xen_vals_sup = [xen_sup_depth[s] for s in shared_sup]
        mer_vals_sup = [merfish_ref["supertype_depth"][s] for s in shared_sup]
        sup_r, _ = pearsonr(xen_vals_sup, mer_vals_sup)
        sup_rho, _ = spearmanr(xen_vals_sup, mer_vals_sup)
        sup_mad = float(np.mean(np.abs(np.array(xen_vals_sup) -
                                        np.array(mer_vals_sup))))
    else:
        sup_r = sup_rho = sup_mad = np.nan

    # -- Proportion comparison --
    total_filt = len(sub_labels)
    xen_sub_props = {}
    for s in np.unique(sub_labels):
        xen_sub_props[s] = float((sub_labels == s).sum()) / total_filt

    # L6b proportion
    l6b_prop_xen = xen_sub_props.get("L6b", 0)
    l6b_prop_mer = merfish_ref["subclass_props"].get("L6b", 0)

    # Deep-layer excitatory types
    deep_types = ["L5 IT", "L5 ET", "L5/6 NP", "L6 IT", "L6 IT Car3",
                  "L6 CT", "L6b"]
    upper_types = ["L2/3 IT", "L4 IT"]

    deep_prop_xen = sum(xen_sub_props.get(s, 0) for s in deep_types)
    deep_prop_mer = sum(merfish_ref["subclass_props"].get(s, 0)
                        for s in deep_types)
    upper_prop_xen = sum(xen_sub_props.get(s, 0) for s in upper_types)
    upper_prop_mer = sum(merfish_ref["subclass_props"].get(s, 0)
                         for s in upper_types)

    # L6b depth stats (for L6b-specific comparison)
    l6b_mask = sub_labels == "L6b"
    l6b_depths_xen = depths[l6b_mask] if l6b_mask.sum() > 0 else np.array([])

    return {
        "n_total": n_total,
        "n_pass_conf": n_pass,
        "n_fail_conf": n_fail,
        "pct_fail_conf": 100 * n_fail / n_total if n_total > 0 else 0,
        # Depth correlations
        "sub_r": sub_r,
        "sub_rho": sub_rho,
        "sub_mad": sub_mad,
        "sup_r": sup_r,
        "sup_rho": sup_rho,
        "sup_mad": sup_mad,
        "n_shared_sub": len(shared_sub),
        "n_shared_sup": len(shared_sup),
        # Proportions
        "l6b_prop_xen": l6b_prop_xen,
        "l6b_prop_mer": l6b_prop_mer,
        "l6b_ratio": l6b_prop_xen / l6b_prop_mer if l6b_prop_mer > 0 else np.nan,
        "deep_prop_xen": deep_prop_xen,
        "deep_prop_mer": deep_prop_mer,
        "upper_prop_xen": upper_prop_xen,
        "upper_prop_mer": upper_prop_mer,
        # L6b depth
        "l6b_median_depth_xen": float(np.median(l6b_depths_xen)) if len(l6b_depths_xen) > 0 else np.nan,
        "l6b_median_depth_mer": float(np.median(merfish_ref["l6b_depths"])),
        "l6b_n_xen": int(l6b_mask.sum()),
        # Per-type depth/proportion data for plotting
        "_xen_sub_depth": xen_sub_depth,
        "_xen_sup_depth": xen_sup_depth,
        "_xen_sub_props": xen_sub_props,
        "_shared_sub": shared_sub,
        "_shared_sup": shared_sup,
        "_l6b_depths_xen": l6b_depths_xen,
    }


# ──────────────────────────────────────────────────────────────────────
# Plotting
# ──────────────────────────────────────────────────────────────────────

def plot_comparison(all_metrics, merfish_ref, out_dir, thresh):
    """Generate multi-panel comparison figure for a specific threshold."""
    n_strat = len(all_metrics)
    fig, axes = plt.subplots(n_strat, 3, figsize=(22, 6 * n_strat))
    if n_strat == 1:
        axes = axes[np.newaxis, :]

    colors = ["#e74c3c", "#3498db", "#2ecc71", "#9b59b6"]

    for i, (skey, metrics_by_thresh) in enumerate(all_metrics.items()):
        strategy = STRATEGIES[skey]
        metrics = metrics_by_thresh[thresh]
        color = colors[i % len(colors)]

        # --- Col 0: Subclass depth scatter ---
        ax = axes[i, 0]
        shared = metrics["_shared_sub"]
        xen_d = metrics["_xen_sub_depth"]
        mer_d = merfish_ref["subclass_depth"]

        for s in shared:
            ax.scatter(mer_d[s], xen_d[s], c=color, s=60, alpha=0.8,
                       edgecolors="white", linewidth=0.5, zorder=3)
            ax.annotate(s, (mer_d[s], xen_d[s]), fontsize=7,
                        xytext=(3, 3), textcoords="offset points",
                        color="#333333")

        ax.plot([0, 1], [0, 1], "k--", alpha=0.3, linewidth=1)
        ax.set_xlabel("MERFISH median depth", fontsize=13)
        ax.set_ylabel("Xenium median depth", fontsize=13)
        ax.set_title(f"{strategy['label']}\nSubclass depth (r={metrics['sub_r']:.3f})",
                      fontsize=14, fontweight="bold")
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        ax.set_aspect("equal")

        # --- Col 1: Supertype depth scatter ---
        ax = axes[i, 1]
        shared_sup = metrics["_shared_sup"]
        xen_d_sup = metrics["_xen_sup_depth"]
        mer_d_sup = merfish_ref["supertype_depth"]

        for s in shared_sup:
            ax.scatter(mer_d_sup[s], xen_d_sup[s], c=color, s=30, alpha=0.6,
                       edgecolors="white", linewidth=0.3, zorder=3)

        ax.plot([0, 1], [0, 1], "k--", alpha=0.3, linewidth=1)
        ax.set_xlabel("MERFISH median depth", fontsize=13)
        ax.set_ylabel("Xenium median depth", fontsize=13)
        ax.set_title(f"Supertype depth (r={metrics['sup_r']:.3f}, "
                      f"n={metrics['n_shared_sup']})",
                      fontsize=14, fontweight="bold")
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        ax.set_aspect("equal")

        # --- Col 2: L6b depth distribution ---
        ax = axes[i, 2]
        bins = np.linspace(0, 1, 50)
        if len(merfish_ref["l6b_depths"]) > 0:
            ax.hist(merfish_ref["l6b_depths"], bins=bins, density=True,
                    alpha=0.5, color="#4fc3f7",
                    label=f"MERFISH (n={len(merfish_ref['l6b_depths']):,})",
                    edgecolor="none")
        l6b_xen = metrics["_l6b_depths_xen"]
        if len(l6b_xen) > 0:
            ax.hist(l6b_xen, bins=bins, density=True, alpha=0.5,
                    color=color,
                    label=f"Xenium (n={len(l6b_xen):,})",
                    edgecolor="none")
        ax.set_xlabel("Normalized depth from pia", fontsize=13)
        ax.set_ylabel("Density", fontsize=13)
        ax.set_title(f"L6b depth distribution\n"
                      f"L6b proportion: Xen={metrics['l6b_prop_xen']*100:.1f}% "
                      f"vs MER={metrics['l6b_prop_mer']*100:.1f}% "
                      f"({metrics['l6b_ratio']:.1f}x)",
                      fontsize=14, fontweight="bold")
        ax.legend(fontsize=11, loc="upper left")
        ax.set_xlim(0, 1.05)

        # Layer boundaries
        for b in [0.10, 0.30, 0.45, 0.65, 0.85]:
            ax.axvline(b, color="#888888", linestyle="--",
                       linewidth=0.6, alpha=0.5)

    fig.suptitle(f"MapMyCells Strategy Comparison — {SAMPLE_ID}\n"
                  f"Subclass confidence filter ≥ {thresh}",
                  fontsize=20, fontweight="bold", y=1.01)
    plt.tight_layout()
    out_path = os.path.join(out_dir,
                            f"mapping_strategy_comparison_conf{thresh:.2f}.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


def plot_proportions(all_metrics, merfish_ref, out_dir, thresh):
    """Plot subclass proportion comparison across strategies."""
    # Get union of all subclasses
    all_subclasses = set()
    for metrics_by_thresh in all_metrics.values():
        metrics = metrics_by_thresh[thresh]
        all_subclasses.update(metrics["_xen_sub_props"].keys())
    all_subclasses.update(merfish_ref["subclass_props"].keys())
    all_subclasses = sorted(all_subclasses)

    n_strat = len(all_metrics)
    fig, ax = plt.subplots(figsize=(20, 8))

    x = np.arange(len(all_subclasses))
    width = 0.8 / (n_strat + 1)

    # MERFISH bars
    mer_vals = [merfish_ref["subclass_props"].get(s, 0) * 100
                for s in all_subclasses]
    ax.bar(x - width * n_strat / 2, mer_vals, width, color="#4fc3f7",
           label="MERFISH", edgecolor="white", linewidth=0.3)

    colors = ["#e74c3c", "#3498db", "#2ecc71", "#9b59b6"]
    for i, (skey, metrics_by_thresh) in enumerate(all_metrics.items()):
        metrics = metrics_by_thresh[thresh]
        vals = [metrics["_xen_sub_props"].get(s, 0) * 100
                for s in all_subclasses]
        offset = -width * n_strat / 2 + width * (i + 1)
        ax.bar(x + offset, vals, width, color=colors[i % len(colors)],
               alpha=0.7, label=STRATEGIES[skey]["label"],
               edgecolor="white", linewidth=0.3)

    ax.set_xticks(x)
    ax.set_xticklabels(all_subclasses, rotation=60, ha="right", fontsize=10)
    ax.set_ylabel("Proportion (%)", fontsize=14)
    ax.set_title(f"Subclass Proportions — {SAMPLE_ID} (conf ≥ {thresh})",
                  fontsize=18, fontweight="bold")
    ax.legend(fontsize=11, loc="upper right")

    plt.tight_layout()
    out_path = os.path.join(out_dir,
                            f"mapping_strategy_proportions_conf{thresh:.2f}.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────

def main():
    t_start = time.time()

    # 1. Load sample
    h5ad_path = os.path.join(H5AD_DIR, f"{SAMPLE_ID}_annotated.h5ad")
    print(f"Loading {SAMPLE_ID}...")
    adata = ad.read_h5ad(h5ad_path)
    n_total = adata.shape[0]

    # 2. Get QC mask
    if "qc_pass" in adata.obs.columns:
        qc_mask = adata.obs["qc_pass"].values.astype(bool)
    else:
        qc_mask = np.ones(n_total, dtype=bool)

    n_pass = int(qc_mask.sum())
    print(f"  {n_total:,} total, {n_pass:,} QC pass")

    # 3. Subset to QC-pass cells
    adata_pass = adata[qc_mask].copy()

    # 4. Load gene mapping and convert to Ensembl
    gene_mapping = load_gene_mapping()
    print(f"  Gene mapping: {len(gene_mapping)} entries")

    adata_ensembl = convert_genes_to_ensembl(adata_pass, gene_mapping)
    print(f"  Mapped {adata_ensembl.shape[1]}/{adata_pass.shape[1]} genes to Ensembl IDs")

    # 5. Load MERFISH reference
    merfish_ref = load_merfish_reference()

    # 6. Load depth model
    print("Loading depth model...")
    model_bundle = load_model(DEPTH_MODEL_PATH)
    print(f"  Subclass names in model: {len(model_bundle['subclass_names'])}")

    # 7. Run each strategy
    all_labels = {}
    all_metrics = {}

    with tempfile.TemporaryDirectory() as tmpdir:
        # Save query h5ad once (all strategies use same input)
        query_path = os.path.join(tmpdir, f"{SAMPLE_ID}_qcpass.h5ad")
        adata_ensembl.write_h5ad(query_path)
        print(f"  Saved query h5ad: {adata_ensembl.shape}")

        for skey, strategy in STRATEGIES.items():
            print(f"\n{'='*60}")
            print(f"  Strategy: {strategy['label']} ({skey})")
            print(f"{'='*60}")

            t_s = time.time()

            # Run MapMyCells
            print("  Running MapMyCells...")
            output_path, algorithm = run_mapmycells(
                query_path, tmpdir, skey,
                strategy["config_overrides"],
                n_processors=1,
            )
            t_map = time.time() - t_s
            print(f"  MapMyCells done in {t_map:.0f}s")

            # Parse results
            labels_df = parse_mapmycells_output(output_path, algorithm)
            all_labels[skey] = labels_df

            # Log label summary
            for level in TAXONOMY_LEVELS:
                col = f"{level}_label"
                conf_col = f"{level}_label_confidence"
                if col in labels_df.columns:
                    n_unique = labels_df[col].nunique()
                    mean_conf = labels_df[conf_col].mean() if conf_col in labels_df.columns else 0
                    print(f"    {level}: {n_unique} types, "
                          f"mean confidence={mean_conf:.3f}")

            # Run depth prediction with new labels
            print("  Running depth prediction...")
            t_d = time.time()

            # Build adata for depth prediction (need spatial coords)
            adata_depth = adata_pass.copy()
            # Write strategy labels into adata
            for col in labels_df.columns:
                adata_depth.obs[col] = labels_df[col].values

            pred_depth, ood_scores = predict_depth(
                adata_depth, model_bundle,
                subclass_col="subclass_label", compute_ood=True
            )
            layers = assign_layers_with_ood(
                pred_depth, ood_scores, model_bundle=model_bundle
            )
            t_depth = time.time() - t_d
            print(f"  Depth prediction done in {t_depth:.0f}s")

            # Compute metrics at each confidence threshold
            metrics_by_thresh = {}
            for thresh in CONF_THRESHOLDS:
                m = compute_metrics(
                    labels_df, pred_depth, layers, merfish_ref,
                    conf_thresh=thresh
                )
                m["time_mapping"] = t_map
                m["time_depth"] = t_depth
                m["conf_thresh"] = thresh
                metrics_by_thresh[thresh] = m

            all_metrics[skey] = metrics_by_thresh

            # Print summary for each threshold
            for thresh in CONF_THRESHOLDS:
                m = metrics_by_thresh[thresh]
                print(f"\n  --- Results (conf ≥ {thresh}) ---")
                print(f"  Cells after filter: {m['n_pass_conf']:,} "
                      f"({m['pct_fail_conf']:.1f}% removed)")
                print(f"  Subclass depth: r={m['sub_r']:.3f}, "
                      f"ρ={m['sub_rho']:.3f}, "
                      f"MAD={m['sub_mad']:.4f} "
                      f"(n={m['n_shared_sub']})")
                print(f"  Supertype depth: r={m['sup_r']:.3f}, "
                      f"ρ={m['sup_rho']:.3f}, "
                      f"MAD={m['sup_mad']:.4f} "
                      f"(n={m['n_shared_sup']})")
                print(f"  L6b: prop={m['l6b_prop_xen']*100:.2f}% "
                      f"(MERFISH: {m['l6b_prop_mer']*100:.2f}%, "
                      f"ratio={m['l6b_ratio']:.1f}x)")
                print(f"  L6b median depth: Xen={m['l6b_median_depth_xen']:.3f} "
                      f"vs MER={m['l6b_median_depth_mer']:.3f}")

    # 8. Generate figures (per threshold)
    print(f"\n{'='*60}")
    print("Generating figures...")
    print(f"{'='*60}")

    os.makedirs(OUT_DIR, exist_ok=True)

    for thresh in CONF_THRESHOLDS:
        plot_comparison(all_metrics, merfish_ref, OUT_DIR, thresh)
        plot_proportions(all_metrics, merfish_ref, OUT_DIR, thresh)

    # 9. Summary table (all strategies × all thresholds)
    summary_rows = []
    for skey, metrics_by_thresh in all_metrics.items():
        for thresh, metrics in metrics_by_thresh.items():
            row = {
                "strategy": STRATEGIES[skey]["label"],
                "conf_thresh": thresh,
                "n_pass_conf": metrics["n_pass_conf"],
                "pct_fail_conf": metrics["pct_fail_conf"],
                "sub_r": metrics["sub_r"],
                "sub_rho": metrics["sub_rho"],
                "sub_mad": metrics["sub_mad"],
                "sup_r": metrics["sup_r"],
                "sup_rho": metrics["sup_rho"],
                "sup_mad": metrics["sup_mad"],
                "l6b_prop_xen_pct": metrics["l6b_prop_xen"] * 100,
                "l6b_prop_mer_pct": metrics["l6b_prop_mer"] * 100,
                "l6b_ratio": metrics["l6b_ratio"],
                "l6b_median_depth_xen": metrics["l6b_median_depth_xen"],
                "l6b_median_depth_mer": metrics["l6b_median_depth_mer"],
                "deep_prop_xen_pct": metrics["deep_prop_xen"] * 100,
                "deep_prop_mer_pct": metrics["deep_prop_mer"] * 100,
                "upper_prop_xen_pct": metrics["upper_prop_xen"] * 100,
                "upper_prop_mer_pct": metrics["upper_prop_mer"] * 100,
                "time_mapping_s": metrics["time_mapping"],
                "time_depth_s": metrics["time_depth"],
            }
            summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    csv_path = os.path.join(OUT_DIR, "mapping_strategy_summary.csv")
    summary_df.to_csv(csv_path, index=False)
    print(f"\nSaved: {csv_path}")

    # Print final summary table (per threshold)
    for thresh in CONF_THRESHOLDS:
        print(f"\n{'='*80}")
        print(f"SUMMARY — {SAMPLE_ID} — Confidence ≥ {thresh}")
        print(f"{'='*80}")
        print(f"{'Strategy':<35s} {'Sub r':>6s} {'Sup r':>6s} {'Sub MAD':>8s} "
              f"{'L6b%':>5s} {'L6b×':>5s} {'L6b dep':>7s} {'%removed':>8s}")
        print(f"{'-'*80}")
        for skey, metrics_by_thresh in all_metrics.items():
            m = metrics_by_thresh[thresh]
            print(f"{STRATEGIES[skey]['label']:<35s} "
                  f"{m['sub_r']:>6.3f} "
                  f"{m['sup_r']:>6.3f} "
                  f"{m['sub_mad']:>8.4f} "
                  f"{m['l6b_prop_xen']*100:>5.1f} "
                  f"{m['l6b_ratio']:>5.1f} "
                  f"{m['l6b_median_depth_xen']:>7.3f} "
                  f"{m['pct_fail_conf']:>7.1f}%")
        print(f"{'MERFISH reference':<35s} "
              f"{'':>6s} {'':>6s} {'':>8s} "
              f"{merfish_ref['subclass_props'].get('L6b', 0)*100:>5.1f} "
              f"{'1.0':>5s} "
              f"{np.median(merfish_ref['l6b_depths']):>7.3f}")

    total_time = time.time() - t_start
    print(f"\nTotal time: {total_time:.0f}s ({total_time/60:.1f} min)")


if __name__ == "__main__":
    main()
