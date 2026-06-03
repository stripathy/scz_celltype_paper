#!/usr/bin/env python
"""
02_magma_enrichment.py — Run MAGMA-style gene property analysis for SCZ.

Loads pre-computed MAGMA gene-level Z-scores (from PGC3 SCZ GWAS), maps ENTREZ
IDs to gene symbols, intersects with SEA-AD specificity matrix, and runs
OLS regression for each cell type.

Input:  MAGMA .genes.out, NCBI gene loc, specificity CSV (from step 01)
Output: results/tables/seaad_magma_scz_enrichment_all_supertypes.csv
        results/figures/magma_scz_enrichment_all_supertypes.png
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
import numpy as np

from scz_celltype_enrichment.config import (
    MAGMA_GENES_OUT, GENE_LOC_FILE, SPECIFICITY_CSV,
    TABLES_DIR, FIGURES_DIR, INTERMEDIATES_DIR,
)
from scz_celltype_enrichment.enrichment.gwas import (
    load_magma_genes, load_entrez_to_symbol, map_and_merge, build_covariates,
)
from scz_celltype_enrichment.enrichment.celltype import run_enrichment_all_types
from scz_celltype_enrichment.utils import classify_supertypes_by_class
from scz_celltype_enrichment.plotting import plot_enrichment_all_supertypes
from scz_celltype_enrichment.utils import ensure_dirs


def main():
    ensure_dirs(TABLES_DIR, FIGURES_DIR, INTERMEDIATES_DIR)

    print("=" * 60)
    print("Step 02: MAGMA-Style Gene Property Analysis")
    print("=" * 60)

    t0 = time.time()

    # Load specificity matrix (from step 01)
    print("\n1. Loading specificity matrix...")
    specificity = pd.read_csv(str(SPECIFICITY_CSV), index_col=0)
    print(f"  Shape: {specificity.shape[0]:,} genes × {specificity.shape[1]} supertypes")

    # Load GWAS gene-level results
    print("\n2. Loading GWAS data...")
    magma_df = load_magma_genes(MAGMA_GENES_OUT)
    entrez_to_symbol = load_entrez_to_symbol(GENE_LOC_FILE)

    # Map and merge
    print("\n3. Mapping and intersecting...")
    spec_matched, gwas_matched = map_and_merge(magma_df, specificity, entrez_to_symbol)

    # Build covariates
    print("\n4. Building covariates...")
    covariates = build_covariates(gwas_matched)
    zstat = gwas_matched["ZSTAT"].values

    # Run enrichment
    print("\n5. Running enrichment analysis...")
    results = run_enrichment_all_types(spec_matched, zstat, covariates)

    # Validate: check class distribution of significant types
    sig_types = results[results["p_fdr"] < 0.05]["supertype"].tolist()
    class_map = classify_supertypes_by_class(sig_types)
    class_counts = pd.Series(class_map).value_counts()
    print(f"\n  FDR-significant types by class:")
    for cls, count in class_counts.items():
        print(f"    {cls}: {count}")

    # Save results
    output_csv = TABLES_DIR / "seaad_magma_scz_enrichment_all_supertypes.csv"
    results.to_csv(output_csv, index=False)
    print(f"\n  Saved: {output_csv}")

    # Save matched data for downstream scripts
    # Save numeric arrays as npz
    np.savez(
        str(INTERMEDIATES_DIR / "gwas_numeric.npz"),
        zstat=zstat,
        covariates=covariates,
    )
    # Save GWAS gene info as CSV (contains strings)
    gwas_matched.to_csv(INTERMEDIATES_DIR / "gwas_matched.csv", index=False)
    print(f"  Saved intermediate data: gwas_numeric.npz + gwas_matched.csv")

    # Generate figure
    print("\n6. Generating figure...")
    plot_enrichment_all_supertypes(
        results, FIGURES_DIR / "magma_scz_enrichment_all_supertypes.png"
    )

    print(f"\nStep 02 complete. Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
