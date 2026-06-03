"""
Test correlation-based label transfer on a MERFISH section from SEA-AD.

Compares our reimplemented correlation-based approach against the
ground truth labels that SEA-AD assigned using scrattch.mapping.
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from scipy import sparse
from sklearn.metrics import (
    confusion_matrix, classification_report,
    accuracy_score, balanced_accuracy_score
)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from label_transfer_corr import (
    load_reference, correlation_label_transfer, _knn_label_transfer,
    build_centroids, correlate_to_centroids, normalize_log,
    TAXONOMY_LEVELS
)

# ─── Paths ───
_BASE = os.path.expanduser("~/Github/SCZ_Xenium")
REF_PATH = os.path.join(_BASE, "data", "reference", "Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad")
MERFISH_PATH = os.path.join(_BASE, "data", "reference", "SEAAD_MTG_MERFISH.2024-12-11.h5ad")
OUTPUT_DIR = os.path.join(_BASE, "output")


def pick_test_section(merfish, n_min=5000, n_max=50000):
    """Pick a well-populated section for testing."""
    barcodes = merfish.obs["Specimen Barcode"].value_counts()
    candidates = barcodes[(barcodes >= n_min) & (barcodes <= n_max)]
    if len(candidates) == 0:
        candidates = barcodes[barcodes >= n_min]
    chosen = candidates.index[0]
    print(f"  Chosen section: {chosen} ({candidates[chosen]:,} cells)")
    return chosen


def prepare_merfish_as_query(merfish_section):
    """
    Prepare a MERFISH section as if it were an unannotated query.

    Strips labels and keeps data as-is. The MERFISH .X contains
    log2(CPM/volume + 1) normalized data (see 01_map_filtering.R).
    We'll pass this directly and handle normalization differences
    in the label transfer function.
    """
    # Keep only real genes (not blanks)
    gene_mask = ~merfish_section.var_names.str.startswith("Blank")
    section = merfish_section[:, gene_mask].copy()

    # Store the ground truth
    truth = {}
    for col in ["Class", "Subclass", "Supertype"]:
        if col in section.obs.columns:
            truth[col] = section.obs[col].values.copy()

    # .X is already log2(CPM+1) normalized
    X = section.X
    if sparse.issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float32)
    # Replace any inf/nan with 0
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

    # Create a clean query AnnData with the expression as-is
    query = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=section.obs_names),
        var=section.var.copy(),
    )
    # Mark that this is already normalized
    query.uns["already_normalized"] = True

    if "spatial" in section.obsm:
        query.obsm["spatial"] = section.obsm["spatial"]
    elif "X_spatial_raw" in section.obsm:
        query.obsm["spatial"] = section.obsm["X_spatial_raw"]

    return query, truth


def map_taxonomy_levels(ref_levels, merfish_levels):
    """
    Map between reference and MERFISH taxonomy column names.

    Reference uses: class_label, subclass_label, cluster_label
    MERFISH uses: Class, Subclass, Supertype
    """
    return {
        "class_label": "Class",
        "subclass_label": "Subclass",
        "cluster_label": "Supertype",  # closest match
    }


def evaluate_transfer(predicted_obs, truth, level_map):
    """
    Evaluate label transfer accuracy against ground truth.

    Returns dict of metrics per level.
    """
    results = {}
    for pred_col, truth_col in level_map.items():
        if pred_col not in predicted_obs.columns or truth_col not in truth:
            continue

        pred = predicted_obs[pred_col].values.astype(str)
        true = truth[truth_col].astype(str)

        # Compute accuracy
        acc = accuracy_score(true, pred)
        # Balanced accuracy (accounts for class imbalance)
        # Only compute for shared labels
        shared_labels = sorted(set(pred) & set(true))

        # For subclass, we need to handle naming differences
        # Reference uses e.g. "L2/3 IT", MERFISH might use the same
        results[pred_col] = {
            "accuracy": acc,
            "n_cells": len(true),
            "n_pred_types": len(set(pred)),
            "n_true_types": len(set(true)),
            "n_shared_types": len(shared_labels),
        }

        print(f"\n  {pred_col} → {truth_col}:")
        print(f"    Accuracy: {acc:.3f}")
        print(f"    Types: {len(set(pred))} predicted, {len(set(true))} true, "
              f"{len(shared_labels)} shared")

        # Show per-type accuracy for subclass
        if "subclass" in pred_col.lower():
            print(f"\n    Per-type accuracy (subclass):")
            for t in sorted(set(true)):
                mask = true == t
                if mask.sum() > 0:
                    type_acc = (pred[mask] == true[mask]).mean()
                    n = mask.sum()
                    # What did they get mapped to?
                    pred_for_type = pd.Series(pred[mask]).value_counts()
                    top_pred = pred_for_type.index[0]
                    top_frac = pred_for_type.values[0] / n
                    match_str = "✓" if top_pred == t else f"→ {top_pred}"
                    print(f"      {t:20s}: {type_acc:.2f} ({n:5d} cells) {match_str} ({top_frac:.1%})")

    return results


def plot_evaluation(predicted, truth, level_map, spatial_coords, output_path):
    """
    Generate a 2×3 comparison figure: truth vs predicted for each level.
    """
    n_levels = sum(1 for k in level_map if k in predicted.columns and level_map[k] in truth)
    if n_levels == 0:
        print("  No levels to plot!")
        return

    fig, axes = plt.subplots(2, n_levels, figsize=(8 * n_levels, 16))
    if n_levels == 1:
        axes = axes.reshape(2, 1)

    col_idx = 0
    for pred_col, truth_col in level_map.items():
        if pred_col not in predicted.columns or truth_col not in truth:
            continue

        pred_labels = predicted[pred_col].values.astype(str)
        true_labels = truth[truth_col].astype(str)

        # Build consistent color map
        all_types = sorted(set(true_labels) | set(pred_labels))
        cmap = plt.cm.get_cmap("tab20", len(all_types))
        color_map = {t: cmap(i) for i, t in enumerate(all_types)}

        x = spatial_coords[:, 0]
        y = spatial_coords[:, 1]

        for row, (labels, title_prefix) in enumerate([
            (true_labels, "Ground Truth"),
            (pred_labels, "Predicted (Corr)")
        ]):
            ax = axes[row, col_idx]
            colors = np.array([color_map.get(l, (0.5, 0.5, 0.5, 1.0)) for l in labels])

            # Shuffle for better visibility
            idx = np.random.permutation(len(x))
            ax.scatter(x[idx], y[idx], c=colors[idx], s=0.3, alpha=0.7, rasterized=True)
            ax.set_aspect("equal")
            ax.set_facecolor((0.07, 0.07, 0.15))
            level_name = truth_col
            acc = accuracy_score(true_labels, pred_labels)
            ax.set_title(f"{title_prefix}: {level_name}\n(accuracy={acc:.3f})",
                         fontsize=16, fontweight="bold")
            ax.tick_params(labelsize=10)

        col_idx += 1

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Saved: {output_path}")


def plot_confusion_heatmap(predicted, truth, level_map, output_path):
    """Plot confusion matrices for each level."""
    n_levels = sum(1 for k in level_map if k in predicted.columns and level_map[k] in truth)
    if n_levels == 0:
        return

    fig, axes = plt.subplots(1, n_levels, figsize=(10 * n_levels, 10))
    if n_levels == 1:
        axes = [axes]

    col_idx = 0
    for pred_col, truth_col in level_map.items():
        if pred_col not in predicted.columns or truth_col not in truth:
            continue

        pred_labels = predicted[pred_col].values.astype(str)
        true_labels = truth[truth_col].astype(str)

        # Get all types, ordered
        all_types = sorted(set(true_labels) | set(pred_labels))

        # Confusion matrix (rows=true, cols=predicted)
        cm = confusion_matrix(true_labels, pred_labels, labels=all_types)
        # Normalize by row (per true type)
        row_sums = cm.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1
        cm_norm = cm / row_sums

        ax = axes[col_idx]
        im = ax.imshow(cm_norm, cmap="Blues", vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(len(all_types)))
        ax.set_yticks(range(len(all_types)))
        ax.set_xticklabels(all_types, rotation=90, fontsize=8)
        ax.set_yticklabels(all_types, fontsize=8)
        ax.set_xlabel("Predicted", fontsize=14)
        ax.set_ylabel("True (MERFISH)", fontsize=14)
        acc = accuracy_score(true_labels, pred_labels)
        ax.set_title(f"{truth_col} Confusion\n(accuracy={acc:.3f})",
                     fontsize=16, fontweight="bold")

        # Add diagonal accuracy text
        for i in range(len(all_types)):
            for j in range(len(all_types)):
                val = cm_norm[i, j]
                if val > 0.05:
                    color = "white" if val > 0.5 else "black"
                    ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                            fontsize=6, color=color)

        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        col_idx += 1

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Saved: {output_path}")


def main():
    print("=" * 60)
    print("Testing correlation-based label transfer on MERFISH section")
    print("=" * 60)

    # ─── Load reference ───
    print("\n1. Loading snRNAseq reference...")
    ref = load_reference(REF_PATH)

    # ─── Load MERFISH ───
    print("\n2. Loading MERFISH data...")
    t0 = time.time()
    merfish = ad.read_h5ad(MERFISH_PATH)
    print(f"  MERFISH: {merfish.shape[0]:,} cells x {merfish.shape[1]} genes "
          f"({time.time()-t0:.1f}s)")

    # ─── Pick a test section ───
    print("\n3. Selecting test section...")
    test_barcode = pick_test_section(merfish)
    section = merfish[merfish.obs["Specimen Barcode"] == test_barcode].copy()
    print(f"  Section shape: {section.shape}")

    # ─── Check what ground truth labels exist ───
    print("\n4. Ground truth labels in MERFISH:")
    for col in ["Class", "Subclass", "Supertype"]:
        if col in section.obs.columns:
            n_types = section.obs[col].nunique()
            print(f"    {col}: {n_types} types")

    # ─── Prepare as query ───
    print("\n5. Preparing MERFISH section as unannotated query...")
    query, truth = prepare_merfish_as_query(section)
    print(f"  Query: {query.shape[0]:,} cells x {query.shape[1]} genes")

    # Get spatial coords
    spatial_coords = None
    if "spatial" in query.obsm:
        spatial_coords = query.obsm["spatial"]
    elif "X_spatial_raw" in section.obsm:
        spatial_coords = section.obsm["X_spatial_raw"]

    # ─── Run correlation-based label transfer ───
    print("\n6. Running CORRELATION-based label transfer...")
    t0 = time.time()
    pred_corr = correlation_label_transfer(query.copy(), ref)
    t_corr = time.time() - t0
    print(f"  Total time: {t_corr:.1f}s")

    # ─── Run kNN-based label transfer for comparison ───
    print("\n7. Running kNN-based label transfer (for comparison)...")
    t0 = time.time()
    pred_knn = _knn_label_transfer(query.copy(), ref)
    t_knn = time.time() - t0
    print(f"  Total time: {t_knn:.1f}s")

    # ─── Evaluate ───
    level_map = map_taxonomy_levels(None, None)

    print("\n" + "=" * 60)
    print("CORRELATION-BASED RESULTS")
    print("=" * 60)
    results_corr = evaluate_transfer(pred_corr.obs, truth, level_map)

    print("\n" + "=" * 60)
    print("kNN-BASED RESULTS (comparison)")
    print("=" * 60)
    results_knn = evaluate_transfer(pred_knn.obs, truth, level_map)

    # ─── Summary comparison ───
    print("\n" + "=" * 60)
    print("COMPARISON SUMMARY")
    print("=" * 60)
    print(f"{'Level':<20s} {'Corr Acc':>10s} {'kNN Acc':>10s} {'Winner':>10s}")
    print("-" * 52)
    for level in TAXONOMY_LEVELS:
        if level in results_corr and level in results_knn:
            a_corr = results_corr[level]["accuracy"]
            a_knn = results_knn[level]["accuracy"]
            winner = "Corr" if a_corr > a_knn else "kNN" if a_knn > a_corr else "Tie"
            truth_col = level_map[level]
            print(f"{truth_col:<20s} {a_corr:>10.3f} {a_knn:>10.3f} {winner:>10s}")

    # ─── Generate plots ───
    print("\n8. Generating evaluation plots...")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    if spatial_coords is not None:
        plot_evaluation(
            pred_corr.obs, truth, level_map, spatial_coords,
            os.path.join(OUTPUT_DIR, "merfish_label_transfer_spatial.png")
        )

    plot_confusion_heatmap(
        pred_corr.obs, truth, level_map,
        os.path.join(OUTPUT_DIR, "merfish_label_transfer_confusion.png")
    )

    # Also plot kNN confusion for comparison
    plot_confusion_heatmap(
        pred_knn.obs, truth,
        {"subclass_label": "Subclass"},  # just subclass for comparison
        os.path.join(OUTPUT_DIR, "merfish_label_transfer_knn_confusion.png")
    )

    # ─── Correlation confidence distribution ───
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for i, level in enumerate(TAXONOMY_LEVELS):
        conf_col = f"{level}_confidence"
        if conf_col in pred_corr.obs.columns:
            ax = axes[i]
            vals = pred_corr.obs[conf_col].values
            ax.hist(vals, bins=50, color="steelblue", edgecolor="white", alpha=0.8)
            ax.axvline(np.median(vals), color="red", ls="--", lw=2,
                       label=f"median={np.median(vals):.3f}")
            ax.set_xlabel("Pearson Correlation", fontsize=14)
            ax.set_ylabel("# Cells", fontsize=14)
            truth_col = level_map.get(level, level)
            ax.set_title(f"{truth_col} Confidence", fontsize=16, fontweight="bold")
            ax.legend(fontsize=12)
            ax.tick_params(labelsize=11)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "merfish_label_transfer_confidence.png"),
                dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved confidence distribution plot")

    # ─── Save numerical results ───
    summary_rows = []
    for level in TAXONOMY_LEVELS:
        truth_col = level_map.get(level, level)
        if level in results_corr:
            summary_rows.append({
                "level": truth_col,
                "method": "correlation",
                "accuracy": results_corr[level]["accuracy"],
                "n_pred_types": results_corr[level]["n_pred_types"],
                "n_true_types": results_corr[level]["n_true_types"],
                "n_shared": results_corr[level]["n_shared_types"],
                "n_cells": results_corr[level]["n_cells"],
                "time_s": t_corr,
            })
        if level in results_knn:
            summary_rows.append({
                "level": truth_col,
                "method": "knn",
                "accuracy": results_knn[level]["accuracy"],
                "n_pred_types": results_knn[level]["n_pred_types"],
                "n_true_types": results_knn[level]["n_true_types"],
                "n_shared": results_knn[level]["n_shared_types"],
                "n_cells": results_knn[level]["n_cells"],
                "time_s": t_knn,
            })
    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(OUTPUT_DIR, "merfish_label_transfer_comparison.csv")
    summary_df.to_csv(summary_path, index=False)
    print(f"  Saved: {summary_path}")

    print("\nDone!")


if __name__ == "__main__":
    main()
