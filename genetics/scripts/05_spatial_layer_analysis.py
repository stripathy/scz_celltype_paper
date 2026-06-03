#!/usr/bin/env python
"""
05_spatial_layer_analysis.py — Layer-stratified SST analysis using spatial depth.

Uses MERFISH + Xenium spatial transcriptomics depth data to classify SST
supertypes as upper/middle/deep, then identifies SCZ risk genes that are
preferentially expressed in upper-layer SST interneurons.

Input:  Specificity matrix, GWAS data, depth CSV, enrichment results
Output: results/tables/sst_supertype_layer_info.csv
        results/tables/sst_upper_layer_spatial_gene_analysis.csv
        results/figures/sst_scz_magma_enrichment.png
        results/figures/sst_supertype_scz_pvalues.png
        results/figures/sst_upper_layer_scz_enrichment.png
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
import numpy as np

from scz_celltype_enrichment.config import (
    SPECIFICITY_CSV, DEPTH_SUPERTYPE_CSV, SCZ_GWAS_GENE_SET_CSV,
    TABLES_DIR, FIGURES_DIR, INTERMEDIATES_DIR,
)
from scz_celltype_enrichment.enrichment.gwas import load_scz_gwas_gene_set
from scz_celltype_enrichment.integration.spatial import (
    load_depth_data, get_sst_supertypes, classify_sst_layers,
    compute_upper_layer_gene_scores,
)
from scz_celltype_enrichment.enrichment.gene_drivers import compute_gene_contributions
from scz_celltype_enrichment.plotting import (
    plot_sst_enrichment_panel, plot_sst_pvalues_by_layer,
    plot_upper_layer_enrichment,
)
from scz_celltype_enrichment.utils import ensure_dirs


def main():
    ensure_dirs(TABLES_DIR, FIGURES_DIR)

    print("=" * 60)
    print("Step 05: Spatial Layer-Stratified SST Analysis")
    print("=" * 60)

    t0 = time.time()

    # Load data
    print("\n1. Loading data...")
    specificity = pd.read_csv(str(SPECIFICITY_CSV), index_col=0)

    gwas_df = pd.read_csv(str(INTERMEDIATES_DIR / "gwas_matched.csv"))
    genes = gwas_df["SYMBOL"].values
    spec_matched = specificity.loc[genes]

    enrichment = pd.read_csv(str(TABLES_DIR / "seaad_magma_scz_enrichment_all_supertypes.csv"))
    scz_genes = load_scz_gwas_gene_set(SCZ_GWAS_GENE_SET_CSV)

    # Load and classify depth
    print("\n2. Loading spatial depth data...")
    depth_df = load_depth_data(DEPTH_SUPERTYPE_CSV)

    print("\n3. Classifying SST supertypes by layer...")
    sst_depth = get_sst_supertypes(depth_df)
    sst_layers = classify_sst_layers(sst_depth)

    # Save layer info
    sst_layers.to_csv(TABLES_DIR / "sst_supertype_layer_info.csv", index=False)
    print(f"  Saved: sst_supertype_layer_info.csv")

    # Upper/deep classification
    upper_types = sst_layers[sst_layers["layer_group"] == "Upper"]["supertype"].tolist()
    deep_types = sst_layers[sst_layers["layer_group"] == "Deep"]["supertype"].tolist()
    all_sst_types = sst_layers["supertype"].tolist()

    print(f"\n  Upper SST: {upper_types}")
    print(f"  Deep SST:  {deep_types}")

    # Compute upper-layer gene scores
    print("\n4. Computing upper-layer gene scores...")
    spatial_genes = compute_upper_layer_gene_scores(
        spec_matched, gwas_df, upper_types, deep_types, all_sst_types,
        scz_gene_set=scz_genes,
    )
    spatial_genes.to_csv(
        TABLES_DIR / "sst_upper_layer_spatial_gene_analysis.csv", index=False
    )
    print(f"  Saved: sst_upper_layer_spatial_gene_analysis.csv")

    # Compute gene contributions for SST enrichment figure
    print("\n5. Computing SST gene contributions for figures...")
    sst_enrichment = enrichment[enrichment["supertype"].str.startswith("Sst")]
    sig_sst = sst_enrichment[sst_enrichment["p_fdr"] < 0.05]["supertype"].tolist()

    gene_contributions = {}
    for st in sig_sst:
        gene_contributions[st] = compute_gene_contributions(
            spec_matched, gwas_df, st, scz_gene_set=scz_genes
        )

    # Generate figures
    print("\n6. Generating figures...")
    plot_sst_enrichment_panel(
        enrichment, sst_layers, gene_contributions,
        FIGURES_DIR / "sst_scz_magma_enrichment.png"
    )
    plot_sst_pvalues_by_layer(
        enrichment, sst_layers,
        FIGURES_DIR / "sst_supertype_scz_pvalues.png"
    )
    plot_upper_layer_enrichment(
        spatial_genes, FIGURES_DIR / "sst_upper_layer_scz_enrichment.png"
    )

    print(f"\nStep 05 complete. Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
