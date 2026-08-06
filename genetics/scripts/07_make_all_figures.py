#!/usr/bin/env python
"""
07_make_all_figures.py — Regenerate all figures from saved table outputs.

This script does NO analysis computation — it only loads CSV results
and recreates all publication-quality figures. Useful for tweaking
figure aesthetics without re-running the analysis pipeline.

Input:  All CSV files in results/tables/
Output: All PNG files in results/figures/
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
import numpy as np

from scz_celltype_enrichment.config import (
    SPECIFICITY_CSV, TABLES_DIR, FIGURES_DIR, DEPTH_SUPERTYPE_CSV,
)
from scz_celltype_enrichment.plotting import (
    plot_enrichment_all_supertypes,
    plot_sst_enrichment_panel,
    plot_sst_pvalues_by_layer,
    plot_upper_layer_enrichment,
    plot_shared_vs_unique,
    plot_hcn1_sag_relationship,
    plot_sst_layer_ephys,
)
from scz_celltype_enrichment.enrichment.gene_drivers import compute_gene_contributions
from scz_celltype_enrichment.integration.spatial import load_depth_data, get_sst_supertypes, classify_sst_layers
from scz_celltype_enrichment.utils import ensure_dirs


def main():
    ensure_dirs(FIGURES_DIR)

    print("=" * 60)
    print("Step 07: Regenerate All Figures")
    print("=" * 60)

    t0 = time.time()

    # Load tables
    enrichment = pd.read_csv(str(TABLES_DIR / "seaad_magma_scz_enrichment_all_supertypes.csv"))

    # Figure 1: All supertypes
    print("\n1. All supertypes enrichment...")
    plot_enrichment_all_supertypes(
        enrichment, FIGURES_DIR / "magma_scz_enrichment_all_supertypes.png"
    )

    # Figure 2: Conditional analysis

    # Figure 3: Shared vs unique
    jaccard_path = TABLES_DIR / "jaccard_similarity.csv"
    shared_path = TABLES_DIR / "shared_pan_interneuron_genes.csv"
    if jaccard_path.exists() and shared_path.exists():
        print("\n3. Shared vs unique signal...")
        jaccard = pd.read_csv(str(jaccard_path), index_col=0)
        shared = pd.read_csv(str(shared_path))
        plot_shared_vs_unique(
            jaccard, shared, FIGURES_DIR / "shared_vs_unique_scz_signal.png"
        )

    # Figures 4-6: SST spatial analysis
    layer_path = TABLES_DIR / "sst_supertype_layer_info.csv"
    spatial_path = TABLES_DIR / "sst_upper_layer_spatial_gene_analysis.csv"

    if layer_path.exists():
        print("\n4. SST enrichment panels...")
        sst_layers = pd.read_csv(str(layer_path))

        # For the SST enrichment panel, we need gene contributions
        # Try to load them if they exist
        gene_contributions = {}
        sst_enrichment = enrichment[enrichment["supertype"].str.startswith("Sst")]
        sig_sst = sst_enrichment[sst_enrichment["p_fdr"] < 0.05]["supertype"].tolist()
        for st in sig_sst:
            safe_name = st.replace(" ", "_").replace("/", "_")
            gc_path = TABLES_DIR / f"scz_enrichment_genes_{safe_name}.csv"
            if gc_path.exists():
                gene_contributions[st] = pd.read_csv(str(gc_path))

        plot_sst_enrichment_panel(
            enrichment, sst_layers, gene_contributions,
            FIGURES_DIR / "sst_scz_magma_enrichment.png"
        )

        print("\n5. SST p-values by layer...")
        plot_sst_pvalues_by_layer(
            enrichment, sst_layers,
            FIGURES_DIR / "sst_supertype_scz_pvalues.png"
        )

    if spatial_path.exists():
        print("\n6. Upper-layer enrichment...")
        spatial_genes = pd.read_csv(str(spatial_path))
        plot_upper_layer_enrichment(
            spatial_genes, FIGURES_DIR / "sst_upper_layer_scz_enrichment.png"
        )

    # Figure 7: HCN1 vs SAG
    ephys_path = TABLES_DIR / "sst_supertype_ephys_summary.csv"
    if ephys_path.exists() and SPECIFICITY_CSV.exists():
        print("\n7. HCN1 vs SAG...")
        specificity = pd.read_csv(str(SPECIFICITY_CSV), index_col=0)
        ephys = pd.read_csv(str(ephys_path))
        if "supertype" in ephys.columns:
            ephys_indexed = ephys.set_index("supertype")
        elif "supertype_scANVI" in ephys.columns:
            ephys_indexed = ephys.set_index("supertype_scANVI")
        else:
            ephys_indexed = ephys

        if "mean_sag" in ephys_indexed.columns:
            plot_hcn1_sag_relationship(
                ephys_indexed, specificity,
                FIGURES_DIR / "hcn1_vs_sag.png"
            )

    # Figure 8: SST layer-ephys
    if layer_path.exists() and ephys_path.exists():
        print("\n8. SST layer-ephys relationships...")
        sst_layers = pd.read_csv(str(layer_path))
        ephys_summary = pd.read_csv(str(ephys_path))

        # Merge if needed
        if "mean_sag" not in sst_layers.columns and "mean_sag" in ephys_summary.columns:
            key_col = "supertype" if "supertype" in ephys_summary.columns else "supertype_scANVI"
            sst_layers = sst_layers.merge(
                ephys_summary.rename(columns={key_col: "supertype"})[["supertype", "mean_sag", "mean_tau"]],
                on="supertype", how="left"
            )

        if "mean_sag" in sst_layers.columns:
            sst_with_ephys = sst_layers.dropna(subset=["mean_sag"])
            if len(sst_with_ephys) > 0:
                plot_sst_layer_ephys(
                    sst_with_ephys,
                    FIGURES_DIR / "sst_layer_tau_sag_relationships.png"
                )

    print(f"\nStep 07 complete. Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
