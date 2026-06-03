#!/usr/bin/env python3
"""
Benchmark our cell typing pipeline on SEA-AD MERFISH data.

Runs the complete pipeline (MapMyCells -> correlation classifier -> doublet QC)
on the 1.89M-cell MERFISH dataset, then benchmarks against ground-truth labels.

Phases:
  A: Prepare MERFISH data (use existing ground-truth labels as HANN-equivalent input)
  B: Run correlation classifier (build centroids from ground-truth exemplars, 2-stage classify)
  C: Use existing depth/layer (avoid circularity — model trained on this data)
  D: Benchmark against MERFISH ground-truth Subclass/Supertype labels
  E: (Future) Benchmark against matched snRNAseq proportions

Design: The MERFISH data already has Allen Institute ground-truth labels at
Subclass/Supertype level. Rather than re-running MapMyCells (which requires
cell_type_mapper and would only use 23/180 MERFISH genes that map to Ensembl),
we use these existing labels as the "HANN equivalent" — treating them as the
first-pass annotation that feeds into the correlation classifier. This:
  1. Uses all 180 MERFISH genes for centroid building (not just 23 mapped to Ensembl)
  2. Tests whether the correlation classifier can reproduce ground-truth labels
  3. Gives us a direct accuracy measure against the Allen Institute's annotation
  4. Mirrors the Xenium pipeline where HANN labels seed the correlation classifier

Key differences from Xenium pipeline:
  - MERFISH has 180 genes vs Xenium 300 (only 23 overlap)
  - Centroids built from MERFISH data itself (all 180 genes)
  - Doublet detection uses adapted markers (GAD1 and SST missing from MERFISH)
  - No cell QC step needed (MERFISH data is already curated)
  - Uses `Section` (69 sections) as the per-sample unit for QC flagging

Usage:
    python3 -u benchmark_merfish_classification.py [--phase A|B|C|D|all]

Output: output/merfish_benchmark/
"""

import os
import sys
import time
import argparse
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import sparse
from sklearn.metrics import (
    confusion_matrix, classification_report,
    f1_score,
)

# ── Paths ──
BASE_DIR = os.path.expanduser("~/Github/SCZ_Xenium")
sys.path.insert(0, os.path.join(BASE_DIR, "code", "pipeline"))
sys.path.insert(0, os.path.join(BASE_DIR, "code", "modules"))
sys.path.insert(0, os.path.join(BASE_DIR, "code", "analysis"))

from pipeline_config import MERFISH_PATH
from config import SUBCLASS_TO_CLASS, CLASS_COLORS

from correlation_classifier import (
    build_subclass_centroids,
    build_supertype_centroids,
    build_flat_centroids,
    run_two_stage_classifier,
    run_flat_classifier,
    flag_low_margin_cells,
    flag_doublet_cells,
    DOUBLET_THRESHOLD_MERFISH,
)

# Output directory
BENCHMARK_DIR = os.path.join(BASE_DIR, "output", "merfish_benchmark")
RECLASSIFIED_PATH = os.path.join(BENCHMARK_DIR, "merfish_reclassified.h5ad")


def ensure_dirs():
    os.makedirs(BENCHMARK_DIR, exist_ok=True)


# ═══════════════════════════════════════════════════════════════════════
# PHASE A: Prepare MERFISH labels for correlation classifier
# ═══════════════════════════════════════════════════════════════════════

def phase_a_prepare(adata):
    """Copy existing MERFISH ground-truth labels into hann_* columns.

    The MERFISH data already has Allen Institute ground-truth labels at
    Class/Subclass/Supertype level. We use these as the "HANN equivalent"
    input to the correlation classifier, just like we use MapMyCells HANN
    labels as the first-pass annotation for Xenium data.

    This approach:
    - Uses all 180 MERFISH genes for centroid building (not just 23 in Ensembl mapping)
    - Tests whether correlation classifier can reproduce ground-truth labels
    - Mirrors the Xenium pipeline (HANN labels -> centroids -> reclassify)

    For confidence, we use the existing 'Supertype confidence' column from MERFISH
    (which represents the Allen Institute's mapping confidence), and compute
    subclass confidence as the max supertype confidence within each subclass.

    Returns adata with hann_* columns added.
    """
    print("\n" + "=" * 70)
    print("PHASE A: Prepare MERFISH Ground-Truth Labels")
    print("=" * 70)

    t0 = time.time()

    # Copy ground-truth labels to hann_* columns
    adata.obs['hann_subclass'] = adata.obs['Subclass'].astype(str).values
    adata.obs['hann_supertype'] = adata.obs['Supertype'].astype(str).values
    adata.obs['hann_class'] = adata.obs['Class'].astype(str).values

    # Use existing Supertype confidence as hann confidence
    if 'Supertype confidence' in adata.obs.columns:
        adata.obs['hann_supertype_confidence'] = \
            adata.obs['Supertype confidence'].astype(np.float32).values
        # Subclass confidence = mean of supertype confidences for that subclass
        # (approximation — in real HANN this would be aggregated votes)
        adata.obs['hann_subclass_confidence'] = \
            adata.obs['Supertype confidence'].astype(np.float32).values
        mean_sup_conf = adata.obs['hann_supertype_confidence'].mean()
        mean_sub_conf = adata.obs['hann_subclass_confidence'].mean()
    else:
        # No confidence column — use 1.0 for all cells
        print("  WARNING: No 'Supertype confidence' column. Using 1.0 for all.")
        adata.obs['hann_supertype_confidence'] = np.float32(1.0)
        adata.obs['hann_subclass_confidence'] = np.float32(1.0)
        mean_sup_conf = 1.0
        mean_sub_conf = 1.0

    # Summary
    n_sub = adata.obs['hann_subclass'].nunique()
    n_sup = adata.obs['hann_supertype'].nunique()
    n_cls = adata.obs['hann_class'].nunique()
    print(f"\n  Ground-truth labels copied to hann_* columns:")
    print(f"    Classes:    {n_cls}")
    print(f"    Subclasses: {n_sub}")
    print(f"    Supertypes: {n_sup}")
    print(f"    Mean confidence — supertype: {mean_sup_conf:.3f}, "
          f"subclass: {mean_sub_conf:.3f}")

    # Class distribution
    print(f"\n  Class distribution:")
    for cls, count in adata.obs['hann_class'].value_counts().items():
        print(f"    {cls:30s}: {count:>9,} ({100*count/adata.n_obs:.1f}%)")

    elapsed = time.time() - t0
    print(f"\n  Phase A complete in {elapsed:.0f}s")

    return adata


# ═══════════════════════════════════════════════════════════════════════
# PHASE B: Correlation Classifier
# ═══════════════════════════════════════════════════════════════════════

def phase_b_correlation_classifier(adata):
    """Run two-stage correlation classifier on MERFISH.

    Mirrors Xenium pipeline step 02b:
    1. Build subclass centroids from top-100 HANN exemplars (180-gene space)
    2. Build supertype centroids within each subclass
    3. Run two-stage Pearson correlation classifier
    4. Low-margin QC (bottom 5% per section)
    5. Doublet detection (adapted for MERFISH panel)

    Returns adata with corr_* and doublet_* columns added.
    """
    print("\n" + "=" * 70)
    print("PHASE B: Correlation Classifier on MERFISH")
    print("=" * 70)

    t0 = time.time()

    # ── Build centroids from HANN exemplars ──
    print("\nBuilding subclass centroids (top-100 HANN exemplars)...")
    sub_centroids, sub_counts, gene_names = build_subclass_centroids(
        adata, top_n=100,
        subclass_col='hann_subclass',
        confidence_col='hann_subclass_confidence',
    )

    print(f"\nBuilding supertype centroids (top-100 exemplars)...")
    sup_centroids, sup_to_sub = build_supertype_centroids(
        adata, top_n=100,
        subclass_col='hann_subclass',
        supertype_col='hann_supertype',
        confidence_col='hann_supertype_confidence',
    )

    # ── Run two-stage classifier ──
    print("\nRunning two-stage classifier on all cells...")
    results = run_two_stage_classifier(
        adata, sub_centroids, sup_centroids, gene_names)

    # Derive corr_class
    results['corr_class'] = results['corr_subclass'].map(
        lambda x: SUBCLASS_TO_CLASS.get(x, 'Unknown'))

    # ── QC: Low margin flagging (per section) ──
    print(f"\nFlagging bottom 5% margin per section...")
    section_ids = adata.obs['Section'].astype(str).values
    margins = results['corr_subclass_margin'].values
    corr_qc_pass, thresholds = flag_low_margin_cells(
        margins, section_ids, percentile=1.0)

    # ── QC: Doublet detection (MERFISH panel) ──
    print(f"\nDetecting spatial doublets (MERFISH panel)...")
    class_labels = results['corr_class'].values
    doublet_suspect, doublet_type, doublet_stats = flag_doublet_cells(
        adata, class_labels, SUBCLASS_TO_CLASS,
        gaba_threshold=DOUBLET_THRESHOLD_MERFISH,
        panel='merfish',
    )

    # Merge: doublets also fail QC
    n_doublet_only = (doublet_suspect & corr_qc_pass).sum()
    corr_qc_pass = corr_qc_pass & (~doublet_suspect)

    # Store results
    for col in results.columns:
        adata.obs[col] = results[col].values
    adata.obs['corr_qc_pass'] = corr_qc_pass
    adata.obs['doublet_suspect'] = doublet_suspect
    adata.obs['doublet_type'] = doublet_type

    # Summary
    print(f"\n  Correlation classifier summary:")
    hann_sub = adata.obs['hann_subclass'].astype(str).values
    corr_sub = adata.obs['corr_subclass'].astype(str).values
    agree = (hann_sub == corr_sub).mean()
    print(f"    HANN vs Corr subclass agreement: {100*agree:.1f}%")
    print(f"    QC flagged: {(~corr_qc_pass).sum():,} / {len(corr_qc_pass):,} "
          f"({100*(~corr_qc_pass).sum()/len(corr_qc_pass):.1f}%)")
    print(f"    Doublets: {doublet_suspect.sum():,} "
          f"({100*doublet_suspect.sum()/len(doublet_suspect):.2f}%)")

    elapsed = time.time() - t0
    print(f"\n  Phase B complete in {elapsed:.0f}s ({elapsed/60:.1f} min)")

    return adata


# ═══════════════════════════════════════════════════════════════════════
# PHASE C: Depth / Layer (use existing)
# ═══════════════════════════════════════════════════════════════════════

def phase_c_depth(adata):
    """Use existing MERFISH depth and layer annotations.

    Our depth model was trained ON this MERFISH data, so re-running it
    would be circular. We use the existing predicted_norm_depth and
    predicted_layer columns.

    Validates that predicted cell types land at biologically plausible depths.
    """
    print("\n" + "=" * 70)
    print("PHASE C: Depth / Layer (using existing annotations)")
    print("=" * 70)

    # Check what depth columns exist
    depth_col = None
    layer_col = None

    for col_candidate in ['predicted_norm_depth', 'Normalized depth from pia']:
        if col_candidate in adata.obs.columns:
            depth_col = col_candidate
            break

    for col_candidate in ['predicted_layer', 'Layer annotation']:
        if col_candidate in adata.obs.columns:
            layer_col = col_candidate
            break

    if depth_col:
        print(f"  Using depth column: {depth_col}")
        depth = adata.obs[depth_col].astype(float)
        print(f"    Non-null: {depth.notna().sum():,} / {len(depth):,}")
        print(f"    Range: [{depth.min():.3f}, {depth.max():.3f}]")
        print(f"    Median: {depth.median():.3f}")
    else:
        print("  WARNING: No depth column found")

    if layer_col:
        print(f"  Using layer column: {layer_col}")
        layers = adata.obs[layer_col].astype(str)
        print(f"    Layer distribution:")
        for layer, count in layers.value_counts().head(10).items():
            print(f"      {layer}: {count:,} ({100*count/len(layers):.1f}%)")
    else:
        print("  WARNING: No layer column found")

    # Validate: median depth of predicted cell types
    if depth_col:
        print(f"\n  Depth validation — median predicted depth by corr_subclass:")
        corr_sub = adata.obs['corr_subclass'].astype(str)
        depth_vals = adata.obs[depth_col].astype(float)

        # Expected depth ordering (superficial to deep)
        expected_order = {
            'L2/3 IT': 0.2, 'L4 IT': 0.35, 'L5 IT': 0.5,
            'L5 ET': 0.55, 'L5/6 NP': 0.6, 'L6 IT': 0.7,
            'L6 CT': 0.75, 'L6b': 0.85,
        }

        for sub in sorted(expected_order.keys()):
            mask = corr_sub == sub
            if mask.sum() > 0:
                med = depth_vals[mask].median()
                expected = expected_order[sub]
                status = "OK" if abs(med - expected) < 0.2 else "CHECK"
                print(f"    {sub:15s}: {med:.3f} (expected ~{expected:.2f}) [{status}]")

    print(f"\n  Phase C complete")
    return adata


# ═══════════════════════════════════════════════════════════════════════
# PHASE D: Benchmark vs MERFISH Ground Truth
# ═══════════════════════════════════════════════════════════════════════

def phase_d_benchmark(adata):
    """Benchmark pipeline predictions against MERFISH ground-truth labels.

    Compares corr_subclass/corr_supertype against existing Subclass/Supertype.
    Generates confusion matrices, F1 scores, and accuracy stratified by
    depth, donor, and confidence.
    """
    print("\n" + "=" * 70)
    print("PHASE D: Benchmark vs MERFISH Ground Truth")
    print("=" * 70)

    t0 = time.time()

    # ── Prepare labels ──
    gt_subclass = adata.obs['Subclass'].astype(str).values
    gt_supertype = adata.obs['Supertype'].astype(str).values
    pred_subclass = adata.obs['corr_subclass'].astype(str).values
    pred_supertype = adata.obs['corr_supertype'].astype(str).values
    corr_qc = adata.obs['corr_qc_pass'].values

    n_total = len(gt_subclass)
    print(f"  Total cells: {n_total:,}")
    print(f"  QC-pass cells: {corr_qc.sum():,} ({100*corr_qc.sum()/n_total:.1f}%)")

    # Use all cells for benchmarking (not just QC-pass) to see full picture
    # But also report QC-pass-only accuracy

    # ── Subclass-level metrics ──
    print(f"\n  === SUBCLASS-LEVEL ACCURACY ===")

    all_subclasses = sorted(set(gt_subclass) | set(pred_subclass))
    n_subclasses = len(all_subclasses)
    print(f"  Ground truth subclasses: {len(set(gt_subclass))}")
    print(f"  Predicted subclasses: {len(set(pred_subclass))}")

    # Overall accuracy
    sub_acc_all = (gt_subclass == pred_subclass).mean()
    sub_acc_qc = (gt_subclass[corr_qc] == pred_subclass[corr_qc]).mean()
    print(f"  Overall subclass accuracy (all cells): {100*sub_acc_all:.1f}%")
    print(f"  Overall subclass accuracy (QC-pass):   {100*sub_acc_qc:.1f}%")

    # Per-subclass metrics
    print(f"\n  Per-subclass F1, precision, recall (all cells):")
    report = classification_report(
        gt_subclass, pred_subclass, output_dict=True, zero_division=0)

    sub_metrics = []
    for sub in sorted(set(gt_subclass)):
        if sub in report:
            m = report[sub]
            n_gt = (gt_subclass == sub).sum()
            n_pred = (pred_subclass == sub).sum()
            sub_metrics.append({
                'subclass': sub,
                'f1': m['f1-score'],
                'precision': m['precision'],
                'recall': m['recall'],
                'n_gt': n_gt,
                'n_pred': n_pred,
                'class': SUBCLASS_TO_CLASS.get(sub, 'Unknown'),
            })
            print(f"    {sub:20s}: F1={m['f1-score']:.3f}  P={m['precision']:.3f}  "
                  f"R={m['recall']:.3f}  n_gt={n_gt:>7,}  n_pred={n_pred:>7,}")

    sub_metrics_df = pd.DataFrame(sub_metrics)

    # Macro averages
    macro_f1 = f1_score(gt_subclass, pred_subclass, average='macro', zero_division=0)
    weighted_f1 = f1_score(gt_subclass, pred_subclass, average='weighted', zero_division=0)
    print(f"\n  Macro F1:    {macro_f1:.3f}")
    print(f"  Weighted F1: {weighted_f1:.3f}")

    # ── Supertype-level metrics ──
    print(f"\n  === SUPERTYPE-LEVEL ACCURACY ===")
    sup_acc_all = (gt_supertype == pred_supertype).mean()
    sup_acc_qc = (gt_supertype[corr_qc] == pred_supertype[corr_qc]).mean()
    print(f"  Overall supertype accuracy (all cells): {100*sup_acc_all:.1f}%")
    print(f"  Overall supertype accuracy (QC-pass):   {100*sup_acc_qc:.1f}%")

    # ── Accuracy by class ──
    print(f"\n  === ACCURACY BY CLASS ===")
    gt_class = np.array([SUBCLASS_TO_CLASS.get(s, 'Unknown') for s in gt_subclass])
    for cls in ['Glutamatergic', 'GABAergic', 'Non-neuronal']:
        mask = gt_class == cls
        if mask.sum() > 0:
            acc = (gt_subclass[mask] == pred_subclass[mask]).mean()
            sup_acc = (gt_supertype[mask] == pred_supertype[mask]).mean()
            print(f"    {cls:20s}: subclass={100*acc:.1f}%  supertype={100*sup_acc:.1f}%  "
                  f"n={mask.sum():,}")

    # ── Accuracy by depth stratum ──
    depth_col = None
    for col in ['predicted_norm_depth', 'Normalized depth from pia']:
        if col in adata.obs.columns:
            depth_col = col
            break

    depth_acc = None
    if depth_col:
        print(f"\n  === ACCURACY BY DEPTH STRATUM ===")
        depth = adata.obs[depth_col].astype(float).values
        strata = {
            'L1 (<0.10)': (0.0, 0.10),
            'L2/3 (0.10-0.40)': (0.10, 0.40),
            'L4 (0.40-0.55)': (0.40, 0.55),
            'L5 (0.55-0.70)': (0.55, 0.70),
            'L6 (0.70-0.90)': (0.70, 0.90),
            'WM (>0.90)': (0.90, 1.5),
        }
        depth_acc_records = []
        for name, (lo, hi) in strata.items():
            mask = (depth >= lo) & (depth < hi) & ~np.isnan(depth)
            if mask.sum() > 0:
                acc = (gt_subclass[mask] == pred_subclass[mask]).mean()
                print(f"    {name:25s}: {100*acc:.1f}%  (n={mask.sum():,})")
                depth_acc_records.append({
                    'stratum': name, 'accuracy': acc, 'n_cells': mask.sum()})
        depth_acc = pd.DataFrame(depth_acc_records)

    # ── Accuracy by donor ──
    print(f"\n  === ACCURACY BY DONOR ===")
    donors = adata.obs['Donor ID'].astype(str).values
    donor_acc_records = []
    for donor in sorted(set(donors)):
        mask = donors == donor
        acc = (gt_subclass[mask] == pred_subclass[mask]).mean()
        n = mask.sum()
        print(f"    {donor}: {100*acc:.1f}% (n={n:,})")
        donor_acc_records.append({'donor': donor, 'accuracy': acc, 'n_cells': n})
    donor_acc = pd.DataFrame(donor_acc_records)

    # ── Accuracy by MERFISH Supertype confidence ──
    print(f"\n  === ACCURACY BY MERFISH SUPERTYPE CONFIDENCE ===")
    if 'Supertype confidence' in adata.obs.columns:
        conf = adata.obs['Supertype confidence'].astype(float).values
        bins = [(0.0, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 0.9), (0.9, 1.01)]
        conf_acc_records = []
        for lo, hi in bins:
            mask = (conf >= lo) & (conf < hi)
            if mask.sum() > 0:
                acc = (gt_subclass[mask] == pred_subclass[mask]).mean()
                print(f"    Conf [{lo:.1f}, {hi:.1f}): {100*acc:.1f}%  (n={mask.sum():,})")
                conf_acc_records.append({
                    'conf_bin': f'[{lo:.1f},{hi:.1f})', 'accuracy': acc,
                    'n_cells': mask.sum()})
        conf_acc = pd.DataFrame(conf_acc_records)
    else:
        print("    No 'Supertype confidence' column found")
        conf_acc = None

    # ── Save results ──
    # Summary CSV
    summary_path = os.path.join(BENCHMARK_DIR, "benchmark_summary.csv")
    sub_metrics_df.to_csv(summary_path, index=False)
    print(f"\n  Saved: {summary_path}")

    if depth_acc is not None:
        depth_path = os.path.join(BENCHMARK_DIR, "accuracy_by_depth.csv")
        depth_acc.to_csv(depth_path, index=False)

    donor_path = os.path.join(BENCHMARK_DIR, "accuracy_by_donor.csv")
    donor_acc.to_csv(donor_path, index=False)

    if conf_acc is not None:
        conf_path = os.path.join(BENCHMARK_DIR, "accuracy_by_confidence.csv")
        conf_acc.to_csv(conf_path, index=False)

    # ── Generate figures ──
    _plot_confusion_matrix(gt_subclass, pred_subclass, adata)
    _plot_f1_bar(sub_metrics_df)
    if depth_acc is not None:
        _plot_accuracy_by_depth(depth_acc)
    if conf_acc is not None:
        _plot_accuracy_by_confidence(conf_acc)

    elapsed = time.time() - t0
    print(f"\n  Phase D complete in {elapsed:.0f}s")

    return {
        'sub_acc_all': sub_acc_all,
        'sub_acc_qc': sub_acc_qc,
        'sup_acc_all': sup_acc_all,
        'sup_acc_qc': sup_acc_qc,
        'macro_f1': macro_f1,
        'weighted_f1': weighted_f1,
        'sub_metrics': sub_metrics_df,
    }


# ═══════════════════════════════════════════════════════════════════════
# Plotting functions
# ═══════════════════════════════════════════════════════════════════════

def _plot_confusion_matrix(gt_subclass, pred_subclass, adata):
    """Plot row-normalized confusion matrix at subclass level."""
    print("  Generating confusion matrix plot...")

    # Use only subclasses present in ground truth
    labels = sorted(set(gt_subclass))
    cm = confusion_matrix(gt_subclass, pred_subclass, labels=labels)

    # Row-normalize
    row_sums = cm.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    cm_norm = cm / row_sums

    fig, ax = plt.subplots(figsize=(16, 14))
    im = ax.imshow(cm_norm, cmap='Blues', vmin=0, vmax=1, aspect='auto')
    plt.colorbar(im, ax=ax, label='Row-normalized proportion', shrink=0.8)

    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=9, ha='center')
    ax.set_yticklabels(labels, fontsize=9)

    ax.set_xlabel("Predicted (our pipeline)", fontsize=14)
    ax.set_ylabel("Ground truth (MERFISH)", fontsize=14)

    # Add diagonal accuracy values
    for i in range(len(labels)):
        val = cm_norm[i, i]
        color = 'white' if val > 0.5 else 'black'
        ax.text(i, i, f'{val:.2f}', ha='center', va='center',
                fontsize=8, color=color, fontweight='bold')

    overall_acc = (gt_subclass == pred_subclass).mean()
    ax.set_title(f"Subclass Confusion Matrix (row-normalized)\n"
                 f"Overall accuracy: {100*overall_acc:.1f}% | "
                 f"{len(gt_subclass):,} cells", fontsize=16)

    plt.tight_layout()
    path = os.path.join(BENCHMARK_DIR, "subclass_confusion_matrix.png")
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved: {path}")


def _plot_f1_bar(sub_metrics_df):
    """Plot per-subclass F1 scores sorted and colored by class."""
    print("  Generating F1 bar plot...")

    df = sub_metrics_df.sort_values('f1', ascending=True).copy()
    colors = [CLASS_COLORS.get(c, '#888888') for c in df['class']]

    fig, ax = plt.subplots(figsize=(10, max(8, len(df) * 0.35)))
    bars = ax.barh(range(len(df)), df['f1'], color=colors, edgecolor='none')

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df['subclass'], fontsize=10)
    ax.set_xlabel("F1 Score", fontsize=14)
    ax.set_xlim(0, 1.05)

    # Add value labels
    for i, (f1, n) in enumerate(zip(df['f1'], df['n_gt'])):
        ax.text(f1 + 0.01, i, f'{f1:.2f} (n={n:,})', va='center', fontsize=9)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=CLASS_COLORS[c], label=c)
                       for c in ['Glutamatergic', 'GABAergic', 'Non-neuronal']]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=12)

    mean_f1 = df['f1'].mean()
    ax.set_title(f"Per-Subclass F1 Score (Macro F1 = {mean_f1:.3f})", fontsize=16)
    ax.axvline(mean_f1, color='gray', linestyle='--', alpha=0.5, label=f'Macro F1')

    plt.tight_layout()
    path = os.path.join(BENCHMARK_DIR, "subclass_f1_bar.png")
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved: {path}")


def _plot_accuracy_by_depth(depth_acc):
    """Bar chart of subclass accuracy per depth stratum."""
    print("  Generating accuracy by depth plot...")

    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(range(len(depth_acc)), depth_acc['accuracy'],
                  color='#4fc3f7', edgecolor='none')

    ax.set_xticks(range(len(depth_acc)))
    ax.set_xticklabels(depth_acc['stratum'], fontsize=11, rotation=30, ha='right')
    ax.set_ylabel("Subclass Accuracy", fontsize=14)
    ax.set_ylim(0, 1.05)

    for i, (acc, n) in enumerate(zip(depth_acc['accuracy'], depth_acc['n_cells'])):
        ax.text(i, acc + 0.01, f'{100*acc:.0f}%\n(n={n:,})',
                ha='center', va='bottom', fontsize=10)

    ax.set_title("Subclass Accuracy by Cortical Depth Stratum", fontsize=16)
    ax.axhline(0.8, color='gray', linestyle='--', alpha=0.3)

    plt.tight_layout()
    path = os.path.join(BENCHMARK_DIR, "accuracy_by_depth.png")
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved: {path}")


def _plot_accuracy_by_confidence(conf_acc):
    """Accuracy vs MERFISH supertype confidence."""
    print("  Generating accuracy by confidence plot...")

    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(range(len(conf_acc)), conf_acc['accuracy'],
                  color='#ef5350', edgecolor='none')

    ax.set_xticks(range(len(conf_acc)))
    ax.set_xticklabels(conf_acc['conf_bin'], fontsize=11)
    ax.set_ylabel("Subclass Accuracy", fontsize=14)
    ax.set_xlabel("MERFISH Supertype Confidence Bin", fontsize=14)
    ax.set_ylim(0, 1.05)

    for i, (acc, n) in enumerate(zip(conf_acc['accuracy'], conf_acc['n_cells'])):
        ax.text(i, acc + 0.01, f'{100*acc:.0f}%\n(n={n:,})',
                ha='center', va='bottom', fontsize=10)

    ax.set_title("Subclass Accuracy by Ground-Truth Confidence", fontsize=16)

    plt.tight_layout()
    path = os.path.join(BENCHMARK_DIR, "accuracy_by_confidence.png")
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved: {path}")


# ═══════════════════════════════════════════════════════════════════════
# PHASE E: Ablation Benchmark
# ═══════════════════════════════════════════════════════════════════════

# 8 ablation configurations
ABLATION_CONFIGS = [
    {
        'name': '1. Baseline (current)',
        'short': 'Baseline',
        'flat': False, 'agg': 'mean', 'norm': 'log_cp10k',
        'blanks': True, 'top_n': 100,
    },
    {
        'name': '2. Flat only',
        'short': 'Flat',
        'flat': True, 'agg': 'mean', 'norm': 'log_cp10k',
        'blanks': True, 'top_n': 100,
    },
    {
        'name': '3. Median only',
        'short': 'Median',
        'flat': False, 'agg': 'median', 'norm': 'log_cp10k',
        'blanks': True, 'top_n': 100,
    },
    {
        'name': '4. log2_cpm only',
        'short': 'log2CPM',
        'flat': False, 'agg': 'mean', 'norm': 'log2_cpm',
        'blanks': True, 'top_n': 100,
    },
    {
        'name': '5. No blanks only',
        'short': 'NoBlanks',
        'flat': False, 'agg': 'mean', 'norm': 'log_cp10k',
        'blanks': False, 'top_n': 100,
    },
    {
        'name': '6. All cells only',
        'short': 'AllCells',
        'flat': False, 'agg': 'mean', 'norm': 'log_cp10k',
        'blanks': True, 'top_n': None,
    },
    {
        'name': '7. Full SEA-AD recipe',
        'short': 'SEA-AD',
        'flat': True, 'agg': 'median', 'norm': 'log2_cpm',
        'blanks': False, 'top_n': None,
    },
    {
        'name': '8. Best hybrid',
        'short': 'Hybrid',
        'flat': True, 'agg': 'median', 'norm': 'log_cp10k',
        'blanks': False, 'top_n': 100,
    },
]


def _run_single_ablation(adata, config, eval_mask=None):
    """Run a single ablation configuration and return metrics.

    Parameters
    ----------
    adata : AnnData
        MERFISH data with hann_subclass/hann_supertype/hann_*_confidence columns
    config : dict
        Ablation configuration
    eval_mask : np.ndarray of bool, optional
        If provided, only evaluate metrics on cells where eval_mask is True.
        Centroids are still built from ALL cells; only evaluation is restricted.

    Returns
    -------
    dict with metrics and per-subclass/supertype results
    """
    name = config['name']
    flat = config['flat']
    agg = config['agg']
    norm = config['norm']
    exclude_blank = not config['blanks']
    top_n = config['top_n']

    print(f"\n  ── {name} ──")
    print(f"     flat={flat}, agg={agg}, norm={norm}, "
          f"exclude_blank={exclude_blank}, top_n={top_n}")

    t0 = time.time()

    if flat:
        # Build flat centroids for all supertypes at once
        centroids, type_to_sub, gene_names, cell_counts = build_flat_centroids(
            adata,
            label_col='hann_supertype',
            confidence_col='hann_supertype_confidence',
            subclass_col='hann_subclass',
            top_n=top_n, agg_func=agg, norm_method=norm,
            exclude_blank=exclude_blank,
        )
        # Run flat classifier
        results = run_flat_classifier(
            adata, centroids, type_to_sub, gene_names, norm_method=norm)
    else:
        # Build hierarchical centroids
        sub_centroids, sub_counts, gene_names = build_subclass_centroids(
            adata, top_n=top_n if top_n is not None else 999999,
            subclass_col='hann_subclass',
            confidence_col='hann_subclass_confidence',
            agg_func=agg, norm_method=norm, exclude_blank=exclude_blank,
        )
        sup_centroids, sup_to_sub = build_supertype_centroids(
            adata, top_n=top_n if top_n is not None else 999999,
            subclass_col='hann_subclass',
            supertype_col='hann_supertype',
            confidence_col='hann_supertype_confidence',
            agg_func=agg, norm_method=norm, exclude_blank=exclude_blank,
        )
        # Run two-stage classifier
        results = run_two_stage_classifier(
            adata, sub_centroids, sup_centroids, gene_names, norm_method=norm)

    elapsed = time.time() - t0

    # ── Compute metrics (optionally on subset) ──
    gt_sub = adata.obs['Subclass'].astype(str).values
    gt_sup = adata.obs['Supertype'].astype(str).values
    pred_sub = results['corr_subclass'].values
    pred_sup = results['corr_supertype'].values

    # Apply evaluation mask if provided
    if eval_mask is not None:
        n_total = len(gt_sub)
        n_eval = eval_mask.sum()
        print(f"     Evaluating on {n_eval:,} / {n_total:,} cells "
              f"({100*n_eval/n_total:.1f}%)")
        gt_sub = gt_sub[eval_mask]
        gt_sup = gt_sup[eval_mask]
        pred_sub = pred_sub[eval_mask]
        pred_sup = pred_sup[eval_mask]

    sub_acc = (gt_sub == pred_sub).mean()
    sup_acc = (gt_sup == pred_sup).mean()
    macro_f1_sub = f1_score(gt_sub, pred_sub, average='macro', zero_division=0)
    weighted_f1_sub = f1_score(gt_sub, pred_sub, average='weighted', zero_division=0)
    macro_f1_sup = f1_score(gt_sup, pred_sup, average='macro', zero_division=0)
    weighted_f1_sup = f1_score(gt_sup, pred_sup, average='weighted', zero_division=0)

    print(f"     Subclass acc: {100*sub_acc:.1f}% | Supertype acc: {100*sup_acc:.1f}%")
    print(f"     Subclass macro F1: {macro_f1_sub:.3f} | weighted F1: {weighted_f1_sub:.3f}")
    print(f"     Supertype macro F1: {macro_f1_sup:.3f} | weighted F1: {weighted_f1_sup:.3f}")
    print(f"     Elapsed: {elapsed:.1f}s")

    # Per-subclass F1
    report = classification_report(gt_sub, pred_sub, output_dict=True, zero_division=0)
    per_sub_f1 = {}
    for sub in sorted(set(gt_sub)):
        if sub in report:
            per_sub_f1[sub] = report[sub]['f1-score']

    # Per-class accuracy
    gt_class = np.array([SUBCLASS_TO_CLASS.get(s, 'Unknown') for s in gt_sub])
    class_acc = {}
    for cls in ['Glutamatergic', 'GABAergic', 'Non-neuronal']:
        mask = gt_class == cls
        if mask.sum() > 0:
            class_acc[cls] = (gt_sub[mask] == pred_sub[mask]).mean()

    # Per-depth accuracy
    depth_col = None
    for col in ['predicted_norm_depth', 'Normalized depth from pia']:
        if col in adata.obs.columns:
            depth_col = col
            break

    depth_acc = {}
    if depth_col:
        depth = adata.obs[depth_col].astype(float).values
        # If eval_mask was applied, subset depth too
        if eval_mask is not None:
            depth = depth[eval_mask]
        strata = {
            'L1': (0.0, 0.10), 'L2/3': (0.10, 0.40), 'L4': (0.40, 0.55),
            'L5': (0.55, 0.70), 'L6': (0.70, 0.90), 'WM': (0.90, 1.5),
        }
        for sname, (lo, hi) in strata.items():
            mask = (depth >= lo) & (depth < hi) & ~np.isnan(depth)
            if mask.sum() > 0:
                depth_acc[sname] = (gt_sub[mask] == pred_sub[mask]).mean()

    return {
        'name': name,
        'short': config['short'],
        'sub_acc': sub_acc,
        'sup_acc': sup_acc,
        'macro_f1_sub': macro_f1_sub,
        'weighted_f1_sub': weighted_f1_sub,
        'macro_f1_sup': macro_f1_sup,
        'weighted_f1_sup': weighted_f1_sup,
        'per_sub_f1': per_sub_f1,
        'class_acc': class_acc,
        'depth_acc': depth_acc,
        'elapsed': elapsed,
        'config': config,
    }


def phase_e_ablation_benchmark(adata, depth_only=False):
    """Run 8 ablation configurations and compare.

    Tests: flat vs hierarchical, mean vs median, log_cp10k vs log2_cpm,
    blank inclusion/exclusion, top-100 vs all-cells, and the full SEA-AD recipe.

    Parameters
    ----------
    adata : AnnData
        MERFISH data with ground-truth labels
    depth_only : bool
        If True, restrict to cells with manually annotated depths (non-NaN
        'Normalized depth from pia'). Centroids are still built from ALL cells,
        but evaluation is only on depth-annotated cells.
    """
    suffix = " [depth-annotated only]" if depth_only else ""
    print("\n" + "=" * 70)
    print(f"PHASE E: Ablation Benchmark (8 configurations){suffix}")
    print("=" * 70)

    t0 = time.time()

    # Ensure HANN labels exist
    if 'hann_subclass' not in adata.obs.columns:
        print("  Setting up HANN labels from ground truth...")
        adata = phase_a_prepare(adata)

    # Build evaluation mask for depth-annotated subset
    eval_mask = None
    if depth_only:
        depth_col = None
        for col in ['Normalized depth from pia', 'predicted_norm_depth']:
            if col in adata.obs.columns:
                depth_col = col
                break
        if depth_col:
            depth_vals = pd.to_numeric(adata.obs[depth_col], errors='coerce')
            eval_mask = depth_vals.notna().values
            n_eval = eval_mask.sum()
            print(f"\n  Depth-annotated subset: {n_eval:,} / {adata.n_obs:,} cells "
                  f"({100*n_eval/adata.n_obs:.1f}%)")
        else:
            print("  WARNING: No depth column found — running on all cells")

    all_results = []
    for i, config in enumerate(ABLATION_CONFIGS):
        print(f"\n{'─'*60}")
        print(f"  Config {i+1}/{len(ABLATION_CONFIGS)}: {config['name']}")
        print(f"{'─'*60}")
        result = _run_single_ablation(adata, config, eval_mask=eval_mask)
        all_results.append(result)

    # ── Build summary table ──
    summary_records = []
    for r in all_results:
        rec = {
            'config': r['name'],
            'short': r['short'],
            'sub_acc': r['sub_acc'],
            'sup_acc': r['sup_acc'],
            'macro_f1_sub': r['macro_f1_sub'],
            'weighted_f1_sub': r['weighted_f1_sub'],
            'macro_f1_sup': r['macro_f1_sup'],
            'weighted_f1_sup': r['weighted_f1_sup'],
            'elapsed_s': r['elapsed'],
            'flat': r['config']['flat'],
            'agg': r['config']['agg'],
            'norm': r['config']['norm'],
            'blanks': r['config']['blanks'],
            'top_n': r['config']['top_n'] if r['config']['top_n'] is not None else 'all',
        }
        # Add class accuracy
        for cls in ['Glutamatergic', 'GABAergic', 'Non-neuronal']:
            rec[f'acc_{cls}'] = r['class_acc'].get(cls, np.nan)
        # Add depth accuracy
        for sname in ['L1', 'L2/3', 'L4', 'L5', 'L6', 'WM']:
            rec[f'depth_{sname}'] = r['depth_acc'].get(sname, np.nan)
        summary_records.append(rec)

    summary_df = pd.DataFrame(summary_records)
    file_suffix = "_depth" if depth_only else ""
    summary_path = os.path.join(BENCHMARK_DIR, f"ablation_summary{file_suffix}.csv")
    summary_df.to_csv(summary_path, index=False)
    print(f"\n  Saved: {summary_path}")

    # ── Print comparison table ──
    label = " [depth-annotated]" if depth_only else ""
    print(f"\n{'='*70}")
    print(f"ABLATION RESULTS SUMMARY{label}")
    print(f"{'='*70}")
    print(f"{'Config':<25s} {'Sub Acc':>8s} {'Sup Acc':>8s} {'Sub mF1':>8s} {'Sup mF1':>8s} {'Time':>6s}")
    print(f"{'─'*65}")
    for r in all_results:
        print(f"{r['short']:<25s} {100*r['sub_acc']:>7.1f}% {100*r['sup_acc']:>7.1f}% "
              f"{r['macro_f1_sub']:>7.3f}  {r['macro_f1_sup']:>7.3f}  {r['elapsed']:>5.0f}s")

    # ── Generate figures ──
    _plot_ablation_comparison(all_results, suffix=file_suffix)

    elapsed = time.time() - t0
    print(f"\n  Phase E complete in {elapsed:.0f}s ({elapsed/60:.1f} min)")

    return all_results


def _plot_ablation_comparison(all_results, suffix=''):
    """Generate multi-panel ablation comparison figure."""
    import matplotlib.gridspec as gridspec

    print("\n  Generating ablation comparison plots...")

    short_names = [r['short'] for r in all_results]
    sub_accs = [r['sub_acc'] for r in all_results]
    sup_accs = [r['sup_acc'] for r in all_results]

    # ── Panel A: Overall accuracy bar chart ──
    fig, axes = plt.subplots(2, 2, figsize=(20, 16))

    # A: Subclass + supertype accuracy
    ax = axes[0, 0]
    x = np.arange(len(short_names))
    w = 0.35
    bars1 = ax.bar(x - w/2, [100*a for a in sub_accs], w,
                    label='Subclass', color='#42a5f5', edgecolor='none')
    bars2 = ax.bar(x + w/2, [100*a for a in sup_accs], w,
                    label='Supertype', color='#ef5350', edgecolor='none')
    ax.set_xticks(x)
    ax.set_xticklabels(short_names, fontsize=11, rotation=30, ha='right')
    ax.set_ylabel("Accuracy (%)", fontsize=14)
    ax.set_title("A. Overall Accuracy by Configuration", fontsize=16)
    ax.legend(fontsize=12)
    ax.set_ylim(0, 105)
    # Add value labels
    for bar_group in [bars1, bars2]:
        for bar in bar_group:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                    f'{height:.1f}', ha='center', va='bottom', fontsize=9)

    # B: Per-subclass F1 heatmap
    ax = axes[0, 1]
    all_subclasses = sorted(set().union(*[r['per_sub_f1'].keys() for r in all_results]))
    f1_matrix = np.zeros((len(all_subclasses), len(all_results)))
    for j, r in enumerate(all_results):
        for i, sub in enumerate(all_subclasses):
            f1_matrix[i, j] = r['per_sub_f1'].get(sub, 0)

    im = ax.imshow(f1_matrix, cmap='RdYlGn', vmin=0, vmax=1, aspect='auto')
    ax.set_xticks(range(len(short_names)))
    ax.set_xticklabels(short_names, fontsize=9, rotation=30, ha='right')
    ax.set_yticks(range(len(all_subclasses)))
    ax.set_yticklabels(all_subclasses, fontsize=8)
    ax.set_title("B. Per-Subclass F1 Score", fontsize=16)
    plt.colorbar(im, ax=ax, label='F1 Score', shrink=0.8)

    # C: Flat vs Hierarchical paired comparison
    ax = axes[1, 0]
    # Compare baseline (#1) vs flat-only (#2)
    baseline_f1 = all_results[0]['per_sub_f1']
    flat_f1 = all_results[1]['per_sub_f1']
    common_subs = sorted(set(baseline_f1.keys()) & set(flat_f1.keys()))
    if common_subs:
        b_vals = [baseline_f1[s] for s in common_subs]
        f_vals = [flat_f1[s] for s in common_subs]
        # Color by class
        colors_list = [CLASS_COLORS.get(SUBCLASS_TO_CLASS.get(s, 'Unknown'), '#888')
                       for s in common_subs]
        ax.scatter(b_vals, f_vals, c=colors_list, s=80, edgecolors='black',
                   linewidth=0.5, zorder=3)
        ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, zorder=1)
        for i, sub in enumerate(common_subs):
            ax.annotate(sub, (b_vals[i], f_vals[i]), fontsize=6,
                       xytext=(3, 3), textcoords='offset points')
    ax.set_xlabel(f"Baseline (hierarchical) F1", fontsize=12)
    ax.set_ylabel(f"Flat classifier F1", fontsize=12)
    ax.set_title("C. Flat vs Hierarchical (per-subclass F1)", fontsize=16)
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, 1.05)
    ax.set_aspect('equal')

    # D: Accuracy by depth stratum (selected configs)
    ax = axes[1, 1]
    strata_order = ['L1', 'L2/3', 'L4', 'L5', 'L6', 'WM']
    # Show baseline, SEA-AD recipe, and best hybrid
    configs_to_show = [0, 6, 7]  # indices into all_results
    colors_depth = ['#42a5f5', '#66bb6a', '#ffa726']
    bar_width = 0.25
    for k, idx in enumerate(configs_to_show):
        r = all_results[idx]
        vals = [100 * r['depth_acc'].get(s, 0) for s in strata_order]
        positions = np.arange(len(strata_order)) + (k - 1) * bar_width
        ax.bar(positions, vals, bar_width, label=r['short'],
               color=colors_depth[k], edgecolor='none')

    ax.set_xticks(range(len(strata_order)))
    ax.set_xticklabels(strata_order, fontsize=12)
    ax.set_ylabel("Subclass Accuracy (%)", fontsize=14)
    ax.set_xlabel("Depth Stratum", fontsize=14)
    ax.set_title("D. Accuracy by Depth (selected configs)", fontsize=16)
    ax.legend(fontsize=11)
    ax.set_ylim(0, 105)

    title_suffix = " [depth-annotated]" if suffix else ""
    plt.suptitle(f"Ablation Benchmark: Replicating SEA-AD Flat Classifier{title_suffix}",
                 fontsize=20, y=1.01)
    plt.tight_layout()
    path = os.path.join(BENCHMARK_DIR, f"ablation_comparison{suffix}.png")
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved: {path}")

    # ── Supplementary: Delta from baseline ──
    fig, ax = plt.subplots(figsize=(12, 6))
    baseline_sub_acc = all_results[0]['sub_acc']
    baseline_sup_acc = all_results[0]['sup_acc']
    deltas_sub = [100 * (r['sub_acc'] - baseline_sub_acc) for r in all_results[1:]]
    deltas_sup = [100 * (r['sup_acc'] - baseline_sup_acc) for r in all_results[1:]]
    names_delta = [r['short'] for r in all_results[1:]]

    x = np.arange(len(names_delta))
    w = 0.35
    ax.bar(x - w/2, deltas_sub, w, label='Δ Subclass', color='#42a5f5')
    ax.bar(x + w/2, deltas_sup, w, label='Δ Supertype', color='#ef5350')
    ax.set_xticks(x)
    ax.set_xticklabels(names_delta, fontsize=12, rotation=30, ha='right')
    ax.set_ylabel("Δ Accuracy vs Baseline (pp)", fontsize=14)
    ax.set_title("Change in Accuracy Relative to Baseline", fontsize=16)
    ax.axhline(0, color='black', linewidth=0.5)
    ax.legend(fontsize=12)

    # Add value labels
    for i in range(len(names_delta)):
        ax.text(x[i] - w/2, deltas_sub[i] + (0.3 if deltas_sub[i] >= 0 else -0.8),
                f'{deltas_sub[i]:+.1f}', ha='center', fontsize=9, color='#1565c0')
        ax.text(x[i] + w/2, deltas_sup[i] + (0.3 if deltas_sup[i] >= 0 else -0.8),
                f'{deltas_sup[i]:+.1f}', ha='center', fontsize=9, color='#c62828')

    plt.tight_layout()
    delta_path = os.path.join(BENCHMARK_DIR, f"ablation_delta_from_baseline{suffix}.png")
    fig.savefig(delta_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"    Saved: {delta_path}")


# ═══════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Benchmark cell typing pipeline on SEA-AD MERFISH data")
    parser.add_argument("--phase", type=str, default="all",
                        choices=["A", "B", "C", "D", "E", "all"],
                        help="Which phase to run (default: all)")
    parser.add_argument("--resume", action="store_true",
                        help="Resume from saved checkpoint (merfish_reclassified.h5ad)")
    parser.add_argument("--depth-only", action="store_true",
                        help="Phase E: evaluate only on cells with manual depth annotations")
    args = parser.parse_args()

    ensure_dirs()

    t_start = time.time()

    # ── Load MERFISH data ──
    if args.resume and os.path.exists(RECLASSIFIED_PATH):
        print(f"Resuming from checkpoint...")
        print(f"  Path: {RECLASSIFIED_PATH}")
        t_load = time.time()
        adata = ad.read_h5ad(RECLASSIFIED_PATH)
        print(f"  Loaded in {time.time()-t_load:.0f}s")
        print(f"  {adata.n_obs:,} cells x {adata.n_vars} genes")
        has_hann = 'hann_subclass' in adata.obs.columns
        has_corr = 'corr_subclass' in adata.obs.columns
        print(f"  Has HANN labels: {has_hann}")
        print(f"  Has corr labels: {has_corr}")
    else:
        print(f"Loading MERFISH h5ad...")
        print(f"  Path: {MERFISH_PATH}")
        t_load = time.time()
        adata = ad.read_h5ad(MERFISH_PATH)
        print(f"  Loaded in {time.time()-t_load:.0f}s")
        print(f"  {adata.n_obs:,} cells x {adata.n_vars} genes")
        print(f"  Donors: {adata.obs['Donor ID'].nunique()}")
        print(f"  Sections: {adata.obs['Section'].nunique()}")
        has_hann = False
        has_corr = False

    # Check for existing ground truth columns
    for col in ['Subclass', 'Supertype', 'Class']:
        if col in adata.obs.columns:
            n_unique = adata.obs[col].nunique()
            print(f"  Ground truth '{col}': {n_unique} unique values")

    # ── Phase A: Prepare ground-truth labels ──
    run_phase = args.phase
    if run_phase in ['A', 'all'] and not has_hann:
        adata = phase_a_prepare(adata)
        has_hann = True
    elif has_hann:
        print(f"\n  Skipping Phase A (HANN labels already present)")

    # ── Phase B: Correlation Classifier ──
    if run_phase in ['B', 'all'] and not has_corr:
        if 'hann_subclass' not in adata.obs.columns:
            print("\n  ERROR: Phase B requires Phase A (HANN labels). "
                  "Run with --phase all first.")
            sys.exit(1)
        adata = phase_b_correlation_classifier(adata)
        # Save checkpoint (this is the expensive step)
        print(f"\n  Saving checkpoint (with corr labels)...")
        t_save = time.time()
        adata.write_h5ad(RECLASSIFIED_PATH)
        print(f"  Saved in {time.time()-t_save:.0f}s: {RECLASSIFIED_PATH}")
        has_corr = True
    elif has_corr:
        print(f"\n  Skipping Phase B (corr labels already present)")

    # ── Phase C: Depth ──
    if run_phase in ['C', 'all']:
        adata = phase_c_depth(adata)

    # ── Phase D: Benchmark ──
    if run_phase in ['D', 'all']:
        if 'corr_subclass' not in adata.obs.columns:
            print("\n  ERROR: Phase D requires Phases A+B. Run with --phase all first.")
            sys.exit(1)
        metrics = phase_d_benchmark(adata)

        # Print final summary
        print("\n" + "=" * 70)
        print("FINAL BENCHMARK SUMMARY")
        print("=" * 70)
        print(f"  Subclass accuracy (all):     {100*metrics['sub_acc_all']:.1f}%")
        print(f"  Subclass accuracy (QC-pass): {100*metrics['sub_acc_qc']:.1f}%")
        print(f"  Supertype accuracy (all):    {100*metrics['sup_acc_all']:.1f}%")
        print(f"  Supertype accuracy (QC-pass):{100*metrics['sup_acc_qc']:.1f}%")
        print(f"  Macro F1:                     {metrics['macro_f1']:.3f}")
        print(f"  Weighted F1:                  {metrics['weighted_f1']:.3f}")

    # ── Phase E: Ablation Benchmark ──
    if run_phase == 'E':
        # Phase E only needs ground truth labels — load fresh if needed
        if 'hann_subclass' not in adata.obs.columns:
            adata = phase_a_prepare(adata)
        ablation_results = phase_e_ablation_benchmark(
            adata, depth_only=args.depth_only)

    elapsed = time.time() - t_start
    print(f"\n{'='*70}")
    print(f"Total elapsed: {elapsed:.0f}s ({elapsed/60:.1f} min)")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
