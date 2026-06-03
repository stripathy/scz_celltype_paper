#!/usr/bin/env python3
"""
Create snRNAseq reference subset from the full SEA-AD MTG dataset.

Filters the full SEA-AD snRNAseq dataset (~1.3M cells, 89 donors) to the
5 neurotypical reference donors, producing a ~137K cell reference for
label transfer and downstream analyses.

Input:  data/reference/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad
Output: data/reference/seaad_mtg_snrnaseq_reference.h5ad

Usage:
    python code/pipeline/create_snrnaseq_reference.py
"""

import gc
import time
import anndata as ad
import numpy as np
import pandas as pd
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[2]
FULL_DATASET = PROJECT_ROOT / "data" / "reference" / "SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"
OUTPUT_PATH = PROJECT_ROOT / "data" / "reference" / "seaad_mtg_snrnaseq_reference.h5ad"

# The 5 neurotypical reference donors
REFERENCE_DONORS = [
    "H18.30.001",
    "H18.30.002",
    "H19.30.001",
    "H19.30.002",
    "H200.1023",
]

# Columns to retain in obs (all columns present in the full dataset that
# are needed by the pipeline or useful for QC/analysis)
OBS_COLUMNS_TO_KEEP = [
    "Donor ID",
    "Neurotypical reference",
    "Class",
    "Subclass",
    "Supertype",
    "Supertype (non-expanded)",
    "Class confidence",
    "Subclass confidence",
    "Supertype confidence",
    "Sex",
    "Age at Death",
    "PMI",
    "Brain Region",
    "Number of UMIs",
    "Genes detected",
    "Fraction mitochondrial UMIs",
    "Doublet score",
    "Used in analysis",
    "Overall AD neuropathological Change",
    "Continuous Pseudo-progression Score",
    "sample_id",
    "sample_name",
]


def main():
    t0 = time.time()

    # ------------------------------------------------------------------
    # 1. Load obs only (backed mode) to identify reference cells
    # ------------------------------------------------------------------
    print(f"Loading obs metadata (backed): {FULL_DATASET}")
    adata_backed = ad.read_h5ad(FULL_DATASET, backed="r")
    obs_full = adata_backed.obs
    print(f"  Full dataset shape: {adata_backed.shape[0]:,} cells x {adata_backed.shape[1]:,} genes")
    print(f"  Loaded metadata in {time.time() - t0:.1f}s")

    # ------------------------------------------------------------------
    # 2. Identify neurotypical reference cells
    # ------------------------------------------------------------------
    assert "Neurotypical reference" in obs_full.columns, \
        "Column 'Neurotypical reference' not found in obs"
    assert "Donor ID" in obs_full.columns, \
        "Column 'Donor ID' not found in obs"

    ref_mask = obs_full["Neurotypical reference"] == "True"
    flagged_donors = sorted(obs_full.loc[ref_mask, "Donor ID"].unique())
    print(f"\n  Donors flagged as 'Neurotypical reference': {flagged_donors}")
    assert set(REFERENCE_DONORS) == set(flagged_donors), \
        f"Donor mismatch! Expected {REFERENCE_DONORS}, got {flagged_donors}"

    ref_indices = np.where(ref_mask.values)[0]
    n_ref = len(ref_indices)
    print(f"  Reference cells: {n_ref:,} / {adata_backed.shape[0]:,}")

    # Get obs for reference cells before closing backed object
    obs_ref = obs_full.iloc[ref_indices].copy()
    var_ref = adata_backed.var.copy()

    # Read the X matrix for reference cells only (memory efficient)
    print(f"\n  Reading expression matrix for {n_ref:,} reference cells...")
    t1 = time.time()
    X_ref = adata_backed.X[ref_indices, :]
    # Convert to dense if sparse, then back to sparse for storage
    import scipy.sparse as sp
    if sp.issparse(X_ref):
        X_ref = X_ref.tocsr()
    else:
        X_ref = sp.csr_matrix(X_ref)
    print(f"  Read expression matrix in {time.time() - t1:.1f}s")

    # Close the backed file
    adata_backed.file.close()
    del adata_backed
    gc.collect()

    # Build the filtered AnnData
    adata_ref = ad.AnnData(X=X_ref, obs=obs_ref, var=var_ref)
    print(f"  After filtering: {adata_ref.shape[0]:,} cells x {adata_ref.shape[1]:,} genes")

    # ------------------------------------------------------------------
    # 3. Subset obs columns
    # ------------------------------------------------------------------
    available_cols = [c for c in OBS_COLUMNS_TO_KEEP if c in adata_ref.obs.columns]
    missing_cols = [c for c in OBS_COLUMNS_TO_KEEP if c not in adata_ref.obs.columns]
    if missing_cols:
        print(f"  Warning: requested obs columns not found: {missing_cols}")

    adata_ref.obs = adata_ref.obs[available_cols].copy()

    # Rename 'Donor ID' -> 'donor_id' for pipeline consistency
    adata_ref.obs = adata_ref.obs.rename(columns={"Donor ID": "donor_id"})

    # Remove unused categories from all categorical columns
    for col in adata_ref.obs.columns:
        if hasattr(adata_ref.obs[col], "cat"):
            adata_ref.obs[col] = adata_ref.obs[col].cat.remove_unused_categories()

    # ------------------------------------------------------------------
    # 4. Save
    # ------------------------------------------------------------------
    print(f"\nSaving to: {OUTPUT_PATH}")
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    adata_ref.write_h5ad(OUTPUT_PATH)
    print(f"  Saved in {time.time() - t0:.1f}s total")

    # ------------------------------------------------------------------
    # 5. Summary statistics
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Cells:      {adata_ref.shape[0]:,}")
    print(f"  Genes:      {adata_ref.shape[1]:,}")
    print(f"  Donors:     {adata_ref.obs['donor_id'].nunique()}")
    print(f"  Classes:    {adata_ref.obs['Class'].nunique()}")
    print(f"  Subclasses: {adata_ref.obs['Subclass'].nunique()}")
    print(f"  Supertypes: {adata_ref.obs['Supertype'].nunique()}")
    print()
    print("Donor distribution:")
    print(adata_ref.obs["donor_id"].value_counts().to_string())
    print()
    print("Class distribution:")
    print(adata_ref.obs["Class"].value_counts().to_string())


if __name__ == "__main__":
    main()
