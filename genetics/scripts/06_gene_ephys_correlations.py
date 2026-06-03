#!/usr/bin/env python
"""
06_gene_ephys_correlations.py — Gene-electrophysiology correlation analysis.

Computes correlations between gene expression and electrophysiology features
(sag, tau) from patch-seq recordings. Identifies genes driving HCN channel
activity and membrane dynamics, and flags SCZ GWAS genes among them.

Input:  Patch-seq CSVs + scANVI labels + specificity matrix + SCZ gene set
Output: results/tables/sag_gene_correlations_sst.csv
        results/tables/sag_gene_correlations_seaad_supertypes.csv
        results/tables/tau_gene_correlations_sst.csv
        results/tables/sag_pos_tau_neg_genes.csv
        results/tables/sst_supertype_ephys_summary.csv
        results/figures/hcn1_vs_sag.png
        results/figures/sst_layer_tau_sag_relationships.png
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
import numpy as np

from scz_celltype_enrichment.config import (
    PROJECT_ROOT, SPECIFICITY_CSV, DEPTH_SUPERTYPE_CSV, SCZ_GWAS_GENE_SET_CSV,
    PATCHSEQ_METADATA_CSV, PATCHSEQ_EPHYS_CSV, SCANVI_RESULTS_CSV,
    TABLES_DIR, FIGURES_DIR,
)
from scz_celltype_enrichment.enrichment.gwas import load_scz_gwas_gene_set
from scz_celltype_enrichment.integration.gene_ephys import (
    load_patchseq_data,
    compute_supertype_correlations,
    find_discordant_genes,
    annotate_with_gwas,
)
from scz_celltype_enrichment.integration.spatial import load_depth_data, get_sst_supertypes, classify_sst_layers
from scz_celltype_enrichment.plotting import plot_hcn1_sag_relationship, plot_sst_layer_ephys
from scz_celltype_enrichment.utils import ensure_dirs


def main():
    ensure_dirs(TABLES_DIR, FIGURES_DIR)

    print("=" * 60)
    print("Step 06: Gene-Electrophysiology Correlations")
    print("=" * 60)

    t0 = time.time()

    # Load data
    print("\n1. Loading data...")
    specificity = pd.read_csv(str(SPECIFICITY_CSV), index_col=0)
    scz_genes = load_scz_gwas_gene_set(SCZ_GWAS_GENE_SET_CSV)

    # Load patch-seq data
    print("\n2. Loading patch-seq data...")
    merged_data, scanvi = load_patchseq_data(
        PATCHSEQ_METADATA_CSV, PATCHSEQ_EPHYS_CSV, SCANVI_RESULTS_CSV
    )

    print(f"  Patch-seq cells with ephys: {len(merged_data)}")
    print(f"  scANVI labeled cells: {len(scanvi)}")

    # Load SST supertype ephys summary from patch-seq data
    # This was computed from patch-seq recordings with scANVI cell type labels
    # and contains mean SAG/TAU per SST supertype
    patchseq_ephys_path = PROJECT_ROOT / "sst_supertype_layer_info.csv"
    if patchseq_ephys_path.exists():
        sst_layer_info = pd.read_csv(str(patchseq_ephys_path))
        # Rename supertype column if needed
        if "supertype_scANVI" in sst_layer_info.columns:
            sst_layer_info = sst_layer_info.rename(
                columns={"supertype_scANVI": "supertype"}
            )
        print(f"  Loaded patch-seq SST ephys summary: {len(sst_layer_info)} supertypes")
    else:
        sst_layer_info = None
        print("  NOTE: sst_supertype_layer_info.csv not found in project root.")

    # Supertype-level correlations using SEA-AD specificity
    print("\n3. Computing supertype-level correlations...")

    if sst_layer_info is not None and "mean_sag" in sst_layer_info.columns:
        print("  Using pre-computed SST ephys from patch-seq...")
        ephys_by_supertype = sst_layer_info.set_index("supertype")[["mean_sag", "mean_tau"]]

        # Filter specificity to SST supertypes
        sst_types = [t for t in ephys_by_supertype.index if t in specificity.columns]
        sst_spec = specificity[sst_types]

        # SAG correlations
        sag_corr = compute_supertype_correlations(
            sst_spec, ephys_by_supertype.loc[sst_types, "mean_sag"],
            feature_name="SAG"
        )
        sag_corr.to_csv(TABLES_DIR / "sag_gene_correlations_seaad_supertypes.csv", index=False)
        print(f"  Saved: sag_gene_correlations_seaad_supertypes.csv ({len(sag_corr):,} genes)")

        # TAU correlations
        tau_corr = compute_supertype_correlations(
            sst_spec, ephys_by_supertype.loc[sst_types, "mean_tau"],
            feature_name="TAU"
        )
        tau_corr.to_csv(TABLES_DIR / "tau_gene_correlations_seaad_supertypes.csv", index=False)
        print(f"  Saved: tau_gene_correlations_seaad_supertypes.csv ({len(tau_corr):,} genes)")

        # Discordant genes (sag+/tau-)
        print("\n4. Finding sag+/tau- discordant genes...")
        discordant = find_discordant_genes(sag_corr, tau_corr, sag_direction="+", tau_direction="-")
        discordant = annotate_with_gwas(discordant, scz_genes)
        discordant.to_csv(TABLES_DIR / "sag_pos_tau_neg_genes.csv", index=False)
        print(f"  Saved: sag_pos_tau_neg_genes.csv ({len(discordant)} genes)")

        # Print top discordant genes
        if len(discordant) > 0:
            print(f"\n  Top sag+/tau- genes:")
            for _, row in discordant.head(10).iterrows():
                gwas_flag = " *SCZ*" if row.get("is_scz", False) else ""
                print(f"    {row['gene']:15s}  sag_rho={row['spearman_rho_sag']:+.3f}  "
                      f"tau_rho={row['spearman_rho_tau']:+.3f}  "
                      f"score={row['combined_score']:.3f}{gwas_flag}")

        # Ephys summary
        ephys_summary = sst_layer_info.copy()
        ephys_summary.to_csv(TABLES_DIR / "sst_supertype_ephys_summary.csv", index=False)

        # Generate figures
        print("\n5. Generating figures...")
        plot_hcn1_sag_relationship(
            ephys_by_supertype, specificity,
            FIGURES_DIR / "hcn1_vs_sag.png"
        )

        # Load spatial depth for layer-ephys figure
        if DEPTH_SUPERTYPE_CSV.exists():
            depth_df = load_depth_data(DEPTH_SUPERTYPE_CSV)
            sst_depth = get_sst_supertypes(depth_df)
            sst_with_layers = classify_sst_layers(sst_depth)

            # Merge ephys into depth data
            sst_merged = sst_with_layers.merge(
                ephys_by_supertype.reset_index().rename(columns={"index": "supertype"}),
                on="supertype", how="left"
            )
            sst_merged = sst_merged.dropna(subset=["mean_sag"])

            if len(sst_merged) > 0:
                plot_sst_layer_ephys(
                    sst_merged, FIGURES_DIR / "sst_layer_tau_sag_relationships.png"
                )
    else:
        print("  WARNING: No pre-computed ephys data available.")
        print("  Supertype-level correlations require running step 05 first,")
        print("  or providing SST ephys data manually.")

    print(f"\nStep 06 complete. Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
