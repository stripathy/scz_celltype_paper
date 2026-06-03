#!/usr/bin/env python3
"""
Step 2c: Hierarchical Harmony-based label transfer from snRNAseq to Xenium.

Three-stage hierarchical classification:
  Stage 1: Class (3 types) — fresh PCA + Harmony + kNN
  Stage 2: Subclass within each class — fresh PCA + Harmony + kNN per class
  Stage 3: Supertype within each subclass — fresh PCA + Harmony + kNN per subclass

Each stage recomputes PCA on just the relevant cell subset, giving higher
resolution for within-group distinctions. Harmony is re-run at each stage
to correct for modality-specific effects within each population.

New columns added to each h5ad:
  - harmony_subclass, harmony_supertype, harmony_class (str)
  - harmony_subclass_confidence, harmony_supertype_confidence (float32)
  - harmony_class_confidence (float32)

Requires: Step 01 (QC) must have been run. Does NOT require steps 02/02b.
Requires: pip install harmonypy

Usage:
    python3 -u 02c_run_harmony_transfer.py
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import sparse
from sklearn.neighbors import KNeighborsClassifier
import harmonypy as hm

# Pipeline config
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipeline_config import H5AD_DIR, SNRNASEQ_REF_PATH

# Shared constants
MODULES_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "modules")
sys.path.insert(0, MODULES_DIR)
from constants import SUBCLASS_TO_CLASS

# ── Settings ──
N_PCS = 30          # PCA components (fewer than usual given 300 genes)
N_NEIGHBORS = 15    # kNN neighbors for label transfer
THETA = 4           # Harmony theta (higher = more aggressive batch correction)
MIN_CELLS_FOR_HARMONY = 50  # Minimum cells in a group to run Harmony


# ── Reference class mapping ──
# Reference uses full names; map to our simplified names
REF_CLASS_MAP = {
    "Neuronal: Glutamatergic": "Glutamatergic",
    "Neuronal: GABAergic": "GABAergic",
    "Non-neuronal and Non-neural": "Non-neuronal",
}


def _harmony_correct(X_pca, source_labels, theta=THETA):
    """Run Harmony on PCA matrix, return corrected embeddings."""
    # Ensure contiguous C-order array (harmonypy/torch needs positive strides)
    X_pca = np.ascontiguousarray(X_pca)
    obs_df = pd.DataFrame({"source": source_labels})
    ho = hm.run_harmony(X_pca, obs_df, "source", theta=theta, max_iter_harmony=20)
    Z = ho.Z_corr
    # Handle potential torch tensor output
    try:
        import torch
        if isinstance(Z, torch.Tensor):
            Z = Z.cpu().numpy()
    except ImportError:
        pass
    if Z.ndim == 2 and Z.shape[0] != len(source_labels):
        Z = Z.T
    return Z.astype(np.float32)


def _normalize_pca_harmony(X_raw, source_labels, n_pcs=N_PCS, theta=THETA):
    """
    Full pipeline: normalize raw counts → PCA → Harmony.
    Returns corrected embedding matrix (n_cells x n_pcs).
    """
    adata = ad.AnnData(X_raw.copy())
    adata.obs["source"] = source_labels

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)

    actual_pcs = min(n_pcs, min(adata.shape) - 1)
    sc.pp.pca(adata, n_comps=actual_pcs)

    # Only run Harmony if we have cells from both sources
    sources = set(source_labels)
    if len(sources) > 1:
        ref_count = sum(1 for s in source_labels if s == "reference")
        query_count = sum(1 for s in source_labels if s == "xenium")
        if ref_count >= MIN_CELLS_FOR_HARMONY and query_count >= MIN_CELLS_FOR_HARMONY:
            Z = _harmony_correct(adata.obsm["X_pca"], source_labels, theta=theta)
            return Z

    return adata.obsm["X_pca"].astype(np.float32)


def _knn_transfer(X_ref, X_query, labels_ref, label_name, n_neighbors=N_NEIGHBORS):
    """Transfer labels via kNN. Returns (predictions, confidence)."""
    valid = pd.notna(labels_ref) & (labels_ref.astype(str) != "")
    X_ref_v = X_ref[valid]
    labels_v = labels_ref[valid]

    unique_labels = sorted(set(labels_v))
    n_labels = len(unique_labels)

    if n_labels == 0:
        return np.full(X_query.shape[0], "Unknown"), np.zeros(X_query.shape[0])

    # If only 1 unique label, assign it to all
    if n_labels == 1:
        return np.full(X_query.shape[0], unique_labels[0]), np.ones(X_query.shape[0])

    k = min(n_neighbors, X_ref_v.shape[0])

    print(f"    kNN({label_name}): {X_ref_v.shape[0]:,} ref -> {X_query.shape[0]:,} query, "
          f"k={k}, {n_labels} labels")

    knn = KNeighborsClassifier(n_neighbors=k, weights="distance")
    knn.fit(X_ref_v, labels_v)

    predictions = knn.predict(X_query)
    probabilities = knn.predict_proba(X_query)
    confidence = probabilities.max(axis=1).astype(np.float32)

    return predictions, confidence


def load_reference(xenium_genes):
    """Load snRNAseq reference, subset to shared genes."""
    t0 = time.time()
    print(f"Loading snRNAseq reference from {SNRNASEQ_REF_PATH}...")
    ref = ad.read_h5ad(SNRNASEQ_REF_PATH)
    print(f"  Reference: {ref.shape[0]:,} cells x {ref.shape[1]:,} genes")

    # Store labels, mapping class names
    ref_labels = ref.obs[["Class", "Subclass", "Supertype"]].copy()
    ref_labels["Class"] = ref_labels["Class"].map(REF_CLASS_MAP)

    # Subset to shared genes
    shared = sorted(set(ref.var_names) & set(xenium_genes))
    print(f"  Shared genes: {len(shared)} / {len(xenium_genes)} Xenium genes")
    ref = ref[:, shared].copy()

    if sparse.issparse(ref.X):
        ref.X = ref.X.toarray()
    ref.X = ref.X.astype(np.float32)

    print(f"  Loaded in {time.time() - t0:.1f}s")
    return ref, shared, ref_labels


def load_xenium_samples(shared_genes):
    """Load all QC-pass Xenium cells, concatenate."""
    t0 = time.time()
    h5ad_files = sorted(
        f for f in os.listdir(H5AD_DIR) if f.endswith("_annotated.h5ad")
    )

    adatas = []
    sample_obs_list = []
    for fname in h5ad_files:
        sample_id = fname.replace("_annotated.h5ad", "")
        adata = ad.read_h5ad(os.path.join(H5AD_DIR, fname))

        if "qc_pass" in adata.obs.columns:
            adata = adata[adata.obs["qc_pass"] == True].copy()

        # Save original obs for writing back & comparison
        orig_obs = adata.obs[["sample_id"]].copy()
        for col in ["corr_subclass", "subclass_label"]:
            if col in adata.obs.columns:
                orig_obs[col] = adata.obs[col]
        sample_obs_list.append(orig_obs)

        # Subset to shared genes
        adata = adata[:, shared_genes].copy()
        if sparse.issparse(adata.X):
            adata.X = adata.X.toarray()
        adata.X = adata.X.astype(np.float32)

        adatas.append(adata)
        print(f"  {sample_id}: {adata.shape[0]:,} QC-pass cells")

    # Concatenate X matrices (not AnnData, to avoid obs conflicts)
    xenium_X = np.vstack([a.X for a in adatas])
    xenium_orig_obs = pd.concat(sample_obs_list)

    print(f"  Total Xenium: {xenium_X.shape[0]:,} cells x {xenium_X.shape[1]} genes")
    print(f"  Loaded in {time.time() - t0:.1f}s")
    return xenium_X, xenium_orig_obs


def run_hierarchical_transfer(ref_X, ref_labels, xenium_X):
    """
    Three-stage hierarchical Harmony + kNN label transfer.

    Stage 1: Class (Glutamatergic / GABAergic / Non-neuronal)
    Stage 2: Subclass within each class (e.g., L2/3 IT, L4 IT, ... within Glut)
    Stage 3: Supertype within each subclass
    """
    n_ref = ref_X.shape[0]
    n_query = xenium_X.shape[0]

    # Output arrays
    class_pred = np.empty(n_query, dtype=object)
    class_conf = np.zeros(n_query, dtype=np.float32)
    subclass_pred = np.empty(n_query, dtype=object)
    subclass_conf = np.zeros(n_query, dtype=np.float32)
    supertype_pred = np.empty(n_query, dtype=object)
    supertype_conf = np.zeros(n_query, dtype=np.float32)

    # ── Stage 1: Class ──
    print("\n" + "=" * 60)
    print("STAGE 1: Class-level classification")
    print("=" * 60)

    source = np.array(["reference"] * n_ref + ["xenium"] * n_query)
    X_all = np.vstack([ref_X, xenium_X])

    Z = _normalize_pca_harmony(X_all, source)
    Z_ref = Z[:n_ref]
    Z_query = Z[n_ref:]

    class_pred[:], class_conf[:] = _knn_transfer(
        Z_ref, Z_query, ref_labels["Class"].values, "Class"
    )

    # Report
    for cls in sorted(set(class_pred)):
        n = (class_pred == cls).sum()
        print(f"    {cls}: {n:,} ({100*n/n_query:.1f}%)")
    print(f"    Class confidence: median={np.median(class_conf):.3f}")

    # ── Stage 2: Subclass within each class ──
    print("\n" + "=" * 60)
    print("STAGE 2: Subclass within each class (fresh PCA + Harmony per class)")
    print("=" * 60)

    for cls in sorted(set(class_pred)):
        print(f"\n  --- {cls} ---")
        ref_mask = ref_labels["Class"].values == cls
        query_mask = class_pred == cls

        n_ref_cls = ref_mask.sum()
        n_query_cls = query_mask.sum()
        print(f"    Ref: {n_ref_cls:,}, Query: {n_query_cls:,}")

        if n_query_cls == 0:
            continue

        # Fresh raw counts for this class
        X_ref_cls = ref_X[ref_mask]
        X_query_cls = xenium_X[query_mask]
        labels_cls = ref_labels.loc[ref_mask, "Subclass"].values

        # Fresh normalize → PCA → Harmony on just this class
        source_cls = (["reference"] * n_ref_cls + ["xenium"] * n_query_cls)
        X_cls = np.vstack([X_ref_cls, X_query_cls])
        Z_cls = _normalize_pca_harmony(X_cls, source_cls)

        Z_ref_cls = Z_cls[:n_ref_cls]
        Z_query_cls = Z_cls[n_ref_cls:]

        pred, conf = _knn_transfer(
            Z_ref_cls, Z_query_cls, labels_cls, f"Subclass({cls})"
        )

        subclass_pred[query_mask] = pred
        subclass_conf[query_mask] = conf

        # Report top types
        for sc_type in pd.Series(pred).value_counts().head(5).index:
            n = (pred == sc_type).sum()
            print(f"      {sc_type}: {n:,} ({100*n/n_query_cls:.1f}%)")

    # ── Stage 3: Supertype within each subclass ──
    print("\n" + "=" * 60)
    print("STAGE 3: Supertype within each subclass (fresh PCA + Harmony per subclass)")
    print("=" * 60)

    for sc_type in sorted(set(subclass_pred)):
        ref_mask = ref_labels["Subclass"].values == sc_type
        query_mask = subclass_pred == sc_type

        n_ref_sc = ref_mask.sum()
        n_query_sc = query_mask.sum()

        if n_query_sc == 0 or n_ref_sc == 0:
            supertype_pred[query_mask] = sc_type  # fallback
            supertype_conf[query_mask] = 0.0
            continue

        # Check how many unique supertypes exist in the reference
        labels_sc = ref_labels.loc[ref_mask, "Supertype"].values
        n_unique = len(set(labels_sc))

        if n_unique == 1:
            # Only one supertype — assign directly
            supertype_pred[query_mask] = labels_sc[0]
            supertype_conf[query_mask] = 1.0
            continue

        # Fresh normalize → PCA → Harmony on just this subclass
        X_ref_sc = ref_X[ref_mask]
        X_query_sc = xenium_X[query_mask]
        source_sc = (["reference"] * n_ref_sc + ["xenium"] * n_query_sc)
        X_sc = np.vstack([X_ref_sc, X_query_sc])

        n_pcs_sc = min(N_PCS, min(X_sc.shape) - 1, n_unique * 2)
        Z_sc = _normalize_pca_harmony(X_sc, source_sc, n_pcs=max(n_pcs_sc, 5))

        Z_ref_sc = Z_sc[:n_ref_sc]
        Z_query_sc = Z_sc[n_ref_sc:]

        pred, conf = _knn_transfer(
            Z_ref_sc, Z_query_sc, labels_sc, f"Supertype({sc_type})"
        )

        supertype_pred[query_mask] = pred
        supertype_conf[query_mask] = conf

    print(f"\n  Supertype confidence: median={np.median(supertype_conf):.3f}")

    return class_pred, class_conf, subclass_pred, subclass_conf, supertype_pred, supertype_conf


def write_results_to_h5ads(xenium_cell_ids, xenium_predictions, xenium_orig_obs):
    """Write Harmony labels back to per-sample h5ad files."""
    t0 = time.time()
    print("\n=== Writing results to h5ad files ===")

    pred_df = pd.DataFrame(xenium_predictions, index=xenium_cell_ids)
    pred_df["sample_id_orig"] = xenium_orig_obs["sample_id"].values

    for sample_id, group in pred_df.groupby("sample_id_orig"):
        h5ad_path = os.path.join(H5AD_DIR, f"{sample_id}_annotated.h5ad")
        if not os.path.exists(h5ad_path):
            print(f"  WARNING: {h5ad_path} not found, skipping")
            continue

        adata = ad.read_h5ad(h5ad_path)

        for col in xenium_predictions:
            if col in adata.obs.columns:
                # Convert existing Categorical to object to accept new values
                if hasattr(adata.obs[col], "cat"):
                    adata.obs[col] = adata.obs[col].astype(object)
            else:
                if "confidence" in col:
                    adata.obs[col] = np.float32(np.nan)
                else:
                    adata.obs[col] = pd.NA

        cell_ids = group.index
        for col in xenium_predictions:
            adata.obs.loc[cell_ids, col] = group[col].values

        adata.write_h5ad(h5ad_path)
        print(f"  {sample_id}: wrote {group.shape[0]:,} Harmony labels")

    print(f"  Write done in {time.time() - t0:.1f}s")


def main():
    t_start = time.time()
    print("=" * 70)
    print("Step 2c: Hierarchical Harmony-based label transfer")
    print("  Stage 1: Class  →  Stage 2: Subclass  →  Stage 3: Supertype")
    print("  Fresh PCA + Harmony at each stage")
    print("=" * 70)

    # 1. Get Xenium gene names
    h5ad_files = sorted(
        f for f in os.listdir(H5AD_DIR) if f.endswith("_annotated.h5ad")
    )
    first_sample = ad.read_h5ad(
        os.path.join(H5AD_DIR, h5ad_files[0]), backed="r"
    )
    xenium_genes = list(first_sample.var_names)
    first_sample.file.close()
    print(f"Xenium panel: {len(xenium_genes)} genes")

    # 2. Load reference
    ref, shared_genes, ref_labels = load_reference(xenium_genes)

    # 3. Load Xenium samples (raw count matrices + obs metadata)
    xenium_X, xenium_orig_obs = load_xenium_samples(shared_genes)

    # 4. Run hierarchical transfer
    (class_pred, class_conf,
     subclass_pred, subclass_conf,
     supertype_pred, supertype_conf) = run_hierarchical_transfer(
        ref.X, ref_labels, xenium_X
    )

    # 5. Store predictions
    xenium_predictions = {
        "harmony_class": class_pred,
        "harmony_class_confidence": class_conf,
        "harmony_subclass": subclass_pred,
        "harmony_subclass_confidence": subclass_conf,
        "harmony_supertype": supertype_pred,
        "harmony_supertype_confidence": supertype_conf,
    }

    # 6. Write back to h5ads
    write_results_to_h5ads(xenium_orig_obs.index, xenium_predictions, xenium_orig_obs)

    # 7. Summary
    n = len(subclass_pred)
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total Xenium cells labeled: {n:,}")

    for level, pred, conf in [
        ("Class", class_pred, class_conf),
        ("Subclass", subclass_pred, subclass_conf),
        ("Supertype", supertype_pred, supertype_conf),
    ]:
        n_types = len(set(pred))
        print(f"\n  {level}: {n_types} unique types, "
              f"confidence median={np.median(conf):.3f}, mean={np.mean(conf):.3f}")
        for name, count in pd.Series(pred).value_counts().head(10).items():
            print(f"    {name}: {count:,} ({100*count/n:.1f}%)")

    # Agreement with existing labels
    if "corr_subclass" in xenium_orig_obs.columns:
        corr = xenium_orig_obs["corr_subclass"].astype(str).values
        harm = subclass_pred.astype(str)
        valid = (corr != "nan") & (corr != "") & (corr != "<NA>")
        if valid.sum() > 0:
            agree = (corr[valid] == harm[valid]).mean()
            print(f"\n  Harmony vs Corr Classifier (subclass): {100*agree:.1f}% agreement")

    if "subclass_label" in xenium_orig_obs.columns:
        hann = xenium_orig_obs["subclass_label"].astype(str).values
        harm = subclass_pred.astype(str)
        valid = (hann != "nan") & (hann != "") & (hann != "<NA>")
        if valid.sum() > 0:
            agree = (hann[valid] == harm[valid]).mean()
            print(f"  Harmony vs HANN (subclass): {100*agree:.1f}% agreement")

    elapsed = time.time() - t_start
    print(f"\nTotal time: {elapsed:.1f}s ({elapsed/60:.1f}min)")


if __name__ == "__main__":
    main()
