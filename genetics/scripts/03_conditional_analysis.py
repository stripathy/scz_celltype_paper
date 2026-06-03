#!/usr/bin/env python
"""
03_conditional_analysis.py — Forward selection for independent cell types.

Starting from all FDR-significant types, iteratively adds cell types that
remain significant (p < 0.05) when conditioned on all previously selected types.
This identifies independently enriched cell types whose signals are not redundant.

Input:  Enrichment results (from step 02), specificity matrix, GWAS data
Output: results/tables/independent_supertypes_conditional.csv
        results/figures/conditional_analysis_results.png
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
import numpy as np

from scz_celltype_enrichment.config import (
    SPECIFICITY_CSV, TABLES_DIR, FIGURES_DIR, INTERMEDIATES_DIR,
    FDR_THRESHOLD, CONDITIONAL_P_THRESHOLD,
)
from scz_celltype_enrichment.enrichment.celltype import get_significant_types
from scz_celltype_enrichment.enrichment.conditional import forward_selection
from scz_celltype_enrichment.plotting import plot_conditional_results
from scz_celltype_enrichment.utils import ensure_dirs


def main():
    ensure_dirs(TABLES_DIR, FIGURES_DIR)

    print("=" * 60)
    print("Step 03: Conditional Analysis (Forward Selection)")
    print("=" * 60)

    t0 = time.time()

    # Load specificity
    print("\n1. Loading data...")
    specificity = pd.read_csv(str(SPECIFICITY_CSV), index_col=0)

    # Load enrichment results
    enrichment_csv = TABLES_DIR / "seaad_magma_scz_enrichment_all_supertypes.csv"
    enrichment = pd.read_csv(str(enrichment_csv))

    # Load merged GWAS data
    numeric = np.load(str(INTERMEDIATES_DIR / "gwas_numeric.npz"))
    zstat = numeric["zstat"]
    covariates = numeric["covariates"]
    gwas_matched = pd.read_csv(str(INTERMEDIATES_DIR / "gwas_matched.csv"))
    genes = gwas_matched["SYMBOL"].values

    # Align specificity to matched genes
    spec_matched = specificity.loc[genes]
    print(f"  Matched: {len(genes):,} genes, {spec_matched.shape[1]} supertypes")

    # Get FDR-significant types
    sig_types = get_significant_types(enrichment, fdr_threshold=FDR_THRESHOLD)
    print(f"\n2. FDR-significant types: {len(sig_types)}")

    # Forward selection
    print(f"\n3. Running forward selection (p_threshold={CONDITIONAL_P_THRESHOLD})...")
    independent = forward_selection(
        sig_types, enrichment, spec_matched, zstat, covariates,
        p_threshold=CONDITIONAL_P_THRESHOLD,
    )

    # Save
    output_csv = TABLES_DIR / "independent_supertypes_conditional.csv"
    independent.to_csv(output_csv, index=False)
    print(f"\n  Saved: {output_csv}")

    # Generate figure
    print("\n4. Generating figure...")
    plot_conditional_results(
        independent, FIGURES_DIR / "conditional_analysis_results.png"
    )

    print(f"\nStep 03 complete. Total time: {time.time() - t0:.1f}s")
    print(f"  Independent cell types: {len(independent)}")
    for _, row in independent.iterrows():
        print(f"    {row['supertype']:25s} ({row['subclass']:12s}) "
              f"marginal p={row['marginal_p']:.2e}  conditional p={row['conditional_p']:.2e}")


if __name__ == "__main__":
    main()
