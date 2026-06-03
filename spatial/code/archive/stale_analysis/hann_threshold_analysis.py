#!/usr/bin/env python3
"""
Analyze HANN confidence thresholds on Br5400 — breakdown by subclass.

Re-runs HANN mapping on Br5400, then sweeps multiple confidence thresholds
to show:
  1. Per-subclass cell removal rates at each threshold
  2. Confidence distribution per subclass
  3. How threshold choice affects depth consistency with MERFISH

Usage:
    /Users/shreejoy/venv/bin/python3 -u hann_threshold_analysis.py
"""

import os
import sys
import time
import json
import tempfile
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# Pipeline config
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "pipeline"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "modules"))
from pipeline_config import (
    H5AD_DIR, OUTPUT_DIR, PRECOMPUTED_STATS_PATH, GENE_MAPPING_PATH,
    TAXONOMY_LEVELS,
)
from depth_model import load_model, predict_depth, assign_layers_with_ood

BASE_DIR = os.path.expanduser("~/Github/SCZ_Xenium")
MERFISH_PATH = os.path.join(BASE_DIR, "data", "reference",
                             "SEAAD_MTG_MERFISH.2024-12-11.h5ad")
OUT_DIR = os.path.join(BASE_DIR, "output", "plots")
DEPTH_MODEL_PATH = os.path.join(OUTPUT_DIR, "depth_model_normalized.pkl")

SAMPLE_ID = "Br5400"

# Thresholds to sweep
THRESHOLDS = [0.0, 0.10, 0.20, 0.28, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80]


# ──────────────────────────────────────────────────────────────────────
# Re-use functions from compare_mapping_strategies.py
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


def _load_taxonomy_tree():
    import h5py
    with h5py.File(PRECOMPUTED_STATS_PATH, 'r') as f:
        tree = json.loads(f['taxonomy_tree'][()].decode())
    return tree


def _build_taxonomy_lookups(tree):
    hierarchy = tree['hierarchy']
    nm = tree['name_mapper']
    clas_level = hierarchy[0]
    subc_level = hierarchy[1]
    supt_level = hierarchy[2]
    clas_data = tree[clas_level]
    subc_data = tree[subc_level]

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


def run_hann(query_h5ad_path, output_dir):
    """Run HANN mapping, return parsed labels DataFrame."""
    from cell_type_mapper.cli.map_to_on_the_fly_markers import OnTheFlyMapper

    hdf5_path = os.path.join(output_dir, "hann_output.h5")
    config = {
        "query_path": query_h5ad_path,
        "hdf5_result_path": hdf5_path,
        "precomputed_stats": {"path": PRECOMPUTED_STATS_PATH},
        "type_assignment": {
            "normalization": "raw",
            "bootstrap_iteration": 100,
            "bootstrap_factor": 0.5,
            "algorithm": "hann",
        },
        "n_processors": 1,
        "query_markers": {"n_per_utility": 30},
        "reference_markers": {"precomputed_path_list": None},
    }

    runner = OnTheFlyMapper(args=[], input_data=config)
    runner.run()

    return _parse_hann_hdf5(hdf5_path)


def _parse_hann_hdf5(hdf5_path):
    """Parse HANN HDF5, return labels + confidence per cell."""
    import h5py

    tree = _load_taxonomy_tree()
    supt_to_name, supt_to_subclass, supt_to_class = _build_taxonomy_lookups(tree)

    with h5py.File(hdf5_path, 'r') as f:
        votes = f['votes'][:]
        correlation = f['correlation'][:]
        cluster_ids = [c.decode() if isinstance(c, bytes) else c
                       for c in f['cluster_identifiers'][:]]

    n_cells = votes.shape[0]

    # Assign each cell to cluster with most votes
    score = votes.astype(float) + correlation * 1e-6
    best_idx = np.argmax(score, axis=1)

    # Supertype confidence = fraction of votes for winning cluster
    max_votes = votes[np.arange(n_cells), best_idx]
    total_votes = votes.sum(axis=1)
    supt_confidence = np.where(total_votes > 0,
                                max_votes / total_votes, 0.0).astype(np.float32)

    best_cluster_ids = [cluster_ids[i] for i in best_idx]

    result = pd.DataFrame(index=range(n_cells))
    result["supertype_label"] = [supt_to_name.get(cid, cid) for cid in best_cluster_ids]
    result["supertype_label_confidence"] = supt_confidence

    subclass_labels = [supt_to_subclass.get(cid, "Unknown") for cid in best_cluster_ids]
    result["subclass_label"] = subclass_labels

    class_labels = [supt_to_class.get(cid, "Unknown") for cid in best_cluster_ids]
    result["class_label"] = class_labels

    # Subclass confidence: sum votes for all supertypes in same subclass
    cluster_subclasses = [supt_to_subclass.get(cid, "") for cid in cluster_ids]
    unique_subclasses = sorted(set(cluster_subclasses))
    sub_to_idx = {s: i for i, s in enumerate(unique_subclasses)}

    sub_votes = np.zeros((n_cells, len(unique_subclasses)), dtype=np.float32)
    for j, sc in enumerate(cluster_subclasses):
        if sc in sub_to_idx:
            sub_votes[:, sub_to_idx[sc]] += votes[:, j]

    subclass_confidence = np.zeros(n_cells, dtype=np.float32)
    for i in range(n_cells):
        sc = subclass_labels[i]
        if sc in sub_to_idx and total_votes[i] > 0:
            subclass_confidence[i] = sub_votes[i, sub_to_idx[sc]] / total_votes[i]
    result["subclass_label_confidence"] = subclass_confidence

    # Class confidence
    cluster_classes = [supt_to_class.get(cid, "") for cid in cluster_ids]
    unique_classes = sorted(set(cluster_classes))
    cls_to_idx = {c: i for i, c in enumerate(unique_classes)}
    cls_votes = np.zeros((n_cells, len(unique_classes)), dtype=np.float32)
    for j, cc in enumerate(cluster_classes):
        if cc in cls_to_idx:
            cls_votes[:, cls_to_idx[cc]] += votes[:, j]

    class_confidence = np.zeros(n_cells, dtype=np.float32)
    for i in range(n_cells):
        cc = class_labels[i]
        if cc in cls_to_idx and total_votes[i] > 0:
            class_confidence[i] = cls_votes[i, cls_to_idx[cc]] / total_votes[i]
    result["class_label_confidence"] = class_confidence

    return result


# ──────────────────────────────────────────────────────────────────────
# MERFISH reference
# ──────────────────────────────────────────────────────────────────────

def load_merfish_reference():
    print("Loading MERFISH reference...")
    adata = ad.read_h5ad(MERFISH_PATH, backed="r")
    obs = adata.obs[["Subclass", "Supertype", "Normalized depth from pia"]].copy()
    obs.columns = ["subclass", "supertype", "depth"]
    obs["depth"] = obs["depth"].astype(float)
    obs = obs.dropna(subset=["depth"])
    obs["subclass"] = obs["subclass"].astype(str)

    sub_depth = obs.groupby("subclass")["depth"].median().to_dict()

    all_obs = adata.obs[["Subclass"]].copy()
    all_obs.columns = ["subclass"]
    all_obs["subclass"] = all_obs["subclass"].astype(str)
    total = len(all_obs)
    sub_props = (all_obs["subclass"].value_counts() / total).to_dict()

    return {"subclass_depth": sub_depth, "subclass_props": sub_props,
            "l6b_depths": obs[obs["subclass"] == "L6b"]["depth"].values}


# ──────────────────────────────────────────────────────────────────────
# Main analysis
# ──────────────────────────────────────────────────────────────────────

def main():
    t_start = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)

    # 1. Load sample
    h5ad_path = os.path.join(H5AD_DIR, f"{SAMPLE_ID}_annotated.h5ad")
    print(f"Loading {SAMPLE_ID}...")
    adata = ad.read_h5ad(h5ad_path)

    # QC filter
    qc_mask = adata.obs["qc_pass"].values.astype(bool)
    adata_pass = adata[qc_mask].copy()
    n_pass = adata_pass.shape[0]
    print(f"  {n_pass:,} QC-pass cells")

    # Gene mapping
    gene_mapping = load_gene_mapping()
    adata_ensembl = convert_genes_to_ensembl(adata_pass, gene_mapping)
    print(f"  {adata_ensembl.shape[1]} genes mapped to Ensembl")

    # MERFISH reference
    merfish_ref = load_merfish_reference()

    # Depth model
    print("Loading depth model...")
    model_bundle = load_model(DEPTH_MODEL_PATH)

    # 2. Run HANN
    print(f"\n{'='*60}")
    print("Running HANN mapping...")
    print(f"{'='*60}")

    with tempfile.TemporaryDirectory() as tmpdir:
        query_path = os.path.join(tmpdir, f"{SAMPLE_ID}_qcpass.h5ad")
        adata_ensembl.write_h5ad(query_path)
        t_map = time.time()
        labels_df = run_hann(query_path, tmpdir)
        print(f"  HANN done in {time.time() - t_map:.0f}s")

    print(f"  {len(labels_df):,} cells labeled")
    print(f"  {labels_df['subclass_label'].nunique()} subclasses")
    print(f"  Mean subclass confidence: {labels_df['subclass_label_confidence'].mean():.3f}")
    print(f"  Mean supertype confidence: {labels_df['supertype_label_confidence'].mean():.3f}")

    # 3. Run depth prediction (once, using all cells)
    print("\nRunning depth prediction...")
    adata_depth = adata_pass.copy()
    for col in labels_df.columns:
        adata_depth.obs[col] = labels_df[col].values

    pred_depth, ood_scores = predict_depth(
        adata_depth, model_bundle,
        subclass_col="subclass_label", compute_ood=True
    )

    # 4. Per-subclass confidence analysis
    print(f"\n{'='*60}")
    print("Per-subclass confidence analysis")
    print(f"{'='*60}")

    subclasses = sorted(labels_df["subclass_label"].unique())
    sub_conf = labels_df["subclass_label_confidence"].values
    sub_labels = labels_df["subclass_label"].values

    # Table: for each subclass × threshold, show removal %
    rows = []
    for sc in subclasses:
        sc_mask = sub_labels == sc
        n_total_sc = int(sc_mask.sum())
        confs = sub_conf[sc_mask]

        row = {"subclass": sc, "n_total": n_total_sc,
               "median_conf": float(np.median(confs)),
               "mean_conf": float(np.mean(confs)),
               "p10_conf": float(np.percentile(confs, 10)),
               "p25_conf": float(np.percentile(confs, 25)),
               }

        for thresh in THRESHOLDS:
            n_pass = int((confs >= thresh).sum())
            pct_removed = 100 * (1 - n_pass / n_total_sc) if n_total_sc > 0 else 0
            row[f"pct_removed_{thresh:.2f}"] = round(pct_removed, 1)
            row[f"n_pass_{thresh:.2f}"] = n_pass

        rows.append(row)

    thresh_df = pd.DataFrame(rows)

    # Print the table
    print(f"\n{'Subclass':<18s} {'N total':>8s} {'Med conf':>9s} "
          + " ".join(f"{'≥'+str(t):>7s}" for t in [0.28, 0.40, 0.50, 0.60, 0.70]))
    print("-" * 80)
    for _, r in thresh_df.iterrows():
        removal_strs = []
        for t in [0.28, 0.40, 0.50, 0.60, 0.70]:
            pct = r[f"pct_removed_{t:.2f}"]
            removal_strs.append(f"{pct:>6.1f}%")
        print(f"{r['subclass']:<18s} {r['n_total']:>8,d} {r['median_conf']:>9.3f} "
              + " ".join(removal_strs))

    # Overall
    print("-" * 80)
    for t in [0.28, 0.40, 0.50, 0.60, 0.70]:
        n_pass_t = int((sub_conf >= t).sum())
        pct_rem = 100 * (1 - n_pass_t / len(sub_conf))
        print(f"  Overall at ≥{t}: {n_pass_t:,} pass ({pct_rem:.1f}% removed)")

    # Save table
    csv_path = os.path.join(OUT_DIR, "hann_threshold_by_subclass.csv")
    thresh_df.to_csv(csv_path, index=False)
    print(f"\nSaved: {csv_path}")

    # 5. Depth metrics at each threshold
    print(f"\n{'='*60}")
    print("Depth consistency vs MERFISH at each threshold")
    print(f"{'='*60}")

    depth_rows = []
    for thresh in THRESHOLDS:
        conf_mask = sub_conf >= thresh
        n_pass_t = int(conf_mask.sum())
        pct_removed = 100 * (1 - n_pass_t / len(sub_conf))

        # Median depth per subclass (filtered)
        sl = sub_labels[conf_mask]
        dp = pred_depth[conf_mask]

        xen_sub_depth = {}
        xen_sub_props = {}
        for s in np.unique(sl):
            mask = sl == s
            if mask.sum() >= 5:
                xen_sub_depth[s] = float(np.median(dp[mask]))
                xen_sub_props[s] = float(mask.sum()) / len(sl)

        shared = sorted(set(xen_sub_depth.keys()) &
                         set(merfish_ref["subclass_depth"].keys()))

        if len(shared) >= 3:
            xv = [xen_sub_depth[s] for s in shared]
            mv = [merfish_ref["subclass_depth"][s] for s in shared]
            sub_r, _ = pearsonr(xv, mv)
            sub_mad = float(np.mean(np.abs(np.array(xv) - np.array(mv))))
        else:
            sub_r = sub_mad = np.nan

        l6b_prop = xen_sub_props.get("L6b", 0)
        l6b_prop_mer = merfish_ref["subclass_props"].get("L6b", 0)
        l6b_ratio = l6b_prop / l6b_prop_mer if l6b_prop_mer > 0 else np.nan

        l6b_mask = sl == "L6b"
        l6b_depth = float(np.median(dp[l6b_mask])) if l6b_mask.sum() > 0 else np.nan

        depth_rows.append({
            "threshold": thresh,
            "n_pass": n_pass_t,
            "pct_removed": pct_removed,
            "sub_r": sub_r,
            "sub_mad": sub_mad,
            "l6b_prop_pct": l6b_prop * 100,
            "l6b_ratio": l6b_ratio,
            "l6b_median_depth": l6b_depth,
            "n_subclasses": len(shared),
        })

    depth_df = pd.DataFrame(depth_rows)

    print(f"\n{'Thresh':>7s} {'N pass':>8s} {'%removed':>9s} "
          f"{'sub_r':>7s} {'sub_MAD':>8s} {'L6b%':>6s} {'L6b×':>6s} {'L6b dep':>8s}")
    print("-" * 80)
    for _, r in depth_df.iterrows():
        print(f"{r['threshold']:>7.2f} {r['n_pass']:>8,.0f} {r['pct_removed']:>8.1f}% "
              f"{r['sub_r']:>7.3f} {r['sub_mad']:>8.4f} "
              f"{r['l6b_prop_pct']:>5.1f}% {r['l6b_ratio']:>5.1f}x "
              f"{r['l6b_median_depth']:>8.3f}")
    print(f"\n  MERFISH: L6b={merfish_ref['subclass_props'].get('L6b',0)*100:.1f}%, "
          f"L6b depth={np.median(merfish_ref['l6b_depths']):.3f}")

    depth_csv = os.path.join(OUT_DIR, "hann_threshold_depth_metrics.csv")
    depth_df.to_csv(depth_csv, index=False)
    print(f"Saved: {depth_csv}")

    # 6. Figures
    print(f"\n{'='*60}")
    print("Generating figures...")
    print(f"{'='*60}")

    # --- Figure 1: Confidence distribution per subclass ---
    fig, axes = plt.subplots(4, 6, figsize=(28, 18))
    axes_flat = axes.ravel()

    # Sort subclasses by median confidence
    sorted_sc = thresh_df.sort_values("median_conf", ascending=True)

    for idx, (_, row) in enumerate(sorted_sc.iterrows()):
        if idx >= len(axes_flat):
            break
        ax = axes_flat[idx]
        sc = row["subclass"]
        sc_confs = sub_conf[sub_labels == sc]

        ax.hist(sc_confs, bins=50, range=(0, 1), color="#3498db",
                edgecolor="none", alpha=0.7)
        ax.axvline(0.28, color="#e74c3c", linestyle="--", linewidth=1.5,
                   label="0.28")
        ax.axvline(0.50, color="#e67e22", linestyle="--", linewidth=1.5,
                   label="0.50")

        pct_28 = row["pct_removed_0.28"]
        pct_50 = row["pct_removed_0.50"]
        ax.set_title(f"{sc}\nn={row['n_total']:,}, med={row['median_conf']:.2f}\n"
                     f"rem@0.28: {pct_28:.0f}%, @0.50: {pct_50:.0f}%",
                     fontsize=10, fontweight="bold")
        ax.set_xlim(0, 1)
        ax.set_xlabel("Subclass conf", fontsize=8)

    # Hide unused axes
    for idx in range(len(sorted_sc), len(axes_flat)):
        axes_flat[idx].set_visible(False)

    fig.suptitle(f"HANN Subclass Confidence Distributions — {SAMPLE_ID}\n"
                 f"(sorted by median confidence, lowest first)",
                 fontsize=18, fontweight="bold")
    plt.tight_layout()
    fig1_path = os.path.join(OUT_DIR, "hann_confidence_distributions.png")
    fig.savefig(fig1_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {fig1_path}")

    # --- Figure 2: Cell removal by subclass at key thresholds ---
    fig, ax = plt.subplots(figsize=(14, 8))

    # Sort by removal at 0.50
    sorted_sc2 = thresh_df.sort_values("pct_removed_0.50", ascending=True)
    y_pos = np.arange(len(sorted_sc2))

    bar_thresholds = [0.28, 0.40, 0.50, 0.60]
    bar_colors = ["#27ae60", "#f39c12", "#e74c3c", "#8e44ad"]
    width = 0.2

    for j, (t, col) in enumerate(zip(bar_thresholds, bar_colors)):
        vals = sorted_sc2[f"pct_removed_{t:.2f}"].values
        ax.barh(y_pos + j * width - width * 1.5, vals, width,
                color=col, alpha=0.8, label=f"≥{t}")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(sorted_sc2["subclass"], fontsize=11)
    ax.set_xlabel("% Cells Removed", fontsize=14)
    ax.set_title(f"HANN Cell Removal by Subclass at Different Thresholds — {SAMPLE_ID}",
                 fontsize=16, fontweight="bold")
    ax.legend(fontsize=12, title="Threshold", loc="lower right")
    ax.axvline(25, color="#999999", linestyle=":", linewidth=1, alpha=0.5)
    ax.axvline(50, color="#999999", linestyle=":", linewidth=1, alpha=0.5)

    plt.tight_layout()
    fig2_path = os.path.join(OUT_DIR, "hann_removal_by_subclass.png")
    fig.savefig(fig2_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {fig2_path}")

    # --- Figure 3: Depth metrics vs threshold (trade-off curve) ---
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    # Panel 1: sub_r vs threshold
    ax = axes[0]
    ax.plot(depth_df["threshold"], depth_df["sub_r"], "o-", color="#e74c3c",
            markersize=8, linewidth=2)
    ax.set_xlabel("Subclass Confidence Threshold", fontsize=14)
    ax.set_ylabel("Pearson r (Xenium vs MERFISH depth)", fontsize=14)
    ax.set_title("Subclass Depth Correlation", fontsize=16, fontweight="bold")
    ax.set_ylim(0.85, 1.0)
    ax.grid(True, alpha=0.3)

    # Panel 2: L6b ratio vs threshold
    ax = axes[1]
    ax.plot(depth_df["threshold"], depth_df["l6b_ratio"], "s-", color="#2ecc71",
            markersize=8, linewidth=2)
    ax.axhline(1.0, color="#888888", linestyle="--", linewidth=1)
    ax.set_xlabel("Subclass Confidence Threshold", fontsize=14)
    ax.set_ylabel("L6b Proportion Ratio (Xenium/MERFISH)", fontsize=14)
    ax.set_title("L6b Inflation", fontsize=16, fontweight="bold")
    ax.grid(True, alpha=0.3)

    # Panel 3: % removed vs threshold
    ax = axes[2]
    ax.plot(depth_df["threshold"], depth_df["pct_removed"], "D-", color="#3498db",
            markersize=8, linewidth=2)
    ax.set_xlabel("Subclass Confidence Threshold", fontsize=14)
    ax.set_ylabel("% Cells Removed", fontsize=14)
    ax.set_title("Data Loss", fontsize=16, fontweight="bold")
    ax.grid(True, alpha=0.3)

    fig.suptitle(f"HANN Threshold Trade-offs — {SAMPLE_ID}",
                 fontsize=20, fontweight="bold", y=1.02)
    plt.tight_layout()
    fig3_path = os.path.join(OUT_DIR, "hann_threshold_tradeoffs.png")
    fig.savefig(fig3_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {fig3_path}")

    # --- Figure 4: Depth scatter at key thresholds ---
    key_thresholds = [0.28, 0.40, 0.50, 0.60]
    fig, axes = plt.subplots(1, 4, figsize=(24, 6))

    for idx, thresh in enumerate(key_thresholds):
        ax = axes[idx]
        conf_mask = sub_conf >= thresh
        sl = sub_labels[conf_mask]
        dp = pred_depth[conf_mask]

        xen_sub_depth = {}
        for s in np.unique(sl):
            mask = sl == s
            if mask.sum() >= 5:
                xen_sub_depth[s] = float(np.median(dp[mask]))

        shared = sorted(set(xen_sub_depth.keys()) &
                         set(merfish_ref["subclass_depth"].keys()))
        xv = [xen_sub_depth[s] for s in shared]
        mv = [merfish_ref["subclass_depth"][s] for s in shared]
        r_val, _ = pearsonr(xv, mv) if len(shared) >= 3 else (np.nan, None)

        # Color deep-layer types differently
        deep_types = {"L5 IT", "L5 ET", "L5/6 NP", "L6 IT", "L6 IT Car3",
                      "L6 CT", "L6b"}

        for s in shared:
            c = "#e74c3c" if s in deep_types else "#3498db"
            ax.scatter(merfish_ref["subclass_depth"][s], xen_sub_depth[s],
                       c=c, s=60, alpha=0.8, edgecolors="white", linewidth=0.5, zorder=3)
            ax.annotate(s, (merfish_ref["subclass_depth"][s], xen_sub_depth[s]),
                        fontsize=6, xytext=(3, 3), textcoords="offset points",
                        color="#333333")

        ax.plot([0, 1], [0, 1], "k--", alpha=0.3, linewidth=1)
        n_pass_t = int(conf_mask.sum())
        pct_rem = 100 * (1 - n_pass_t / len(sub_conf))
        ax.set_xlabel("MERFISH depth", fontsize=12)
        ax.set_ylabel("Xenium depth", fontsize=12)
        ax.set_title(f"Threshold ≥{thresh}\nr={r_val:.3f}, "
                     f"{pct_rem:.1f}% removed",
                     fontsize=13, fontweight="bold")
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        ax.set_aspect("equal")

    fig.suptitle(f"HANN Subclass Depth vs MERFISH at Different Thresholds — {SAMPLE_ID}\n"
                 f"(red = deep-layer excitatory types)",
                 fontsize=16, fontweight="bold", y=1.04)
    plt.tight_layout()
    fig4_path = os.path.join(OUT_DIR, "hann_depth_scatter_thresholds.png")
    fig.savefig(fig4_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {fig4_path}")

    total_time = time.time() - t_start
    print(f"\nDone in {total_time:.0f}s ({total_time/60:.1f} min)")


if __name__ == "__main__":
    main()
