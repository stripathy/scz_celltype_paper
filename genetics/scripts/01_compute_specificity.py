#!/usr/bin/env python
"""
01_compute_specificity.py — Compute cell-type specificity from SEA-AD snRNA-seq.

This is the most expensive step (~5-10 min for 8.3 GB h5ad load).
Output is cached; subsequent runs skip computation unless --force is passed.

Input:  nicole_sea_ad_snrnaseq_reference.h5ad (137,303 cells × 36,601 genes)
Output: results/intermediates/seaad_supertype_specificity.csv

Usage:
    python scripts/01_compute_specificity.py          # Use cache if available
    python scripts/01_compute_specificity.py --force   # Recompute from scratch
"""
import sys
import time
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from scz_celltype_enrichment.config import (
    SEAAD_H5AD, SPECIFICITY_CSV, MEAN_EXPRESSION_CSV, INTERMEDIATES_DIR
)
from scz_celltype_enrichment.enrichment.specificity import (
    compute_and_save_specificity,
    compute_and_save_mean_expression,
)
from scz_celltype_enrichment.utils import ensure_dirs


def main():
    ensure_dirs(INTERMEDIATES_DIR)
    force = "--force" in sys.argv

    print("=" * 60)
    print("Step 01: Compute Cell-Type Specificity & Mean Expression")
    print("=" * 60)
    print(f"  Input:  {SEAAD_H5AD}")
    print(f"  Output: {SPECIFICITY_CSV}")
    print(f"          {MEAN_EXPRESSION_CSV}")
    print(f"  Force:  {force}")
    print()

    t0 = time.time()
    specificity = compute_and_save_specificity(SEAAD_H5AD, SPECIFICITY_CSV, force=force)

    print()
    mean_expr = compute_and_save_mean_expression(SEAAD_H5AD, MEAN_EXPRESSION_CSV, force=force)

    print(f"\nStep 01 complete. Total time: {time.time() - t0:.1f}s")
    print(f"  Specificity matrix: {specificity.shape[0]:,} genes × {specificity.shape[1]} supertypes")
    print(f"  Mean expression matrix: {mean_expr.shape[0]:,} genes × {mean_expr.shape[1]} supertypes")


if __name__ == "__main__":
    main()
