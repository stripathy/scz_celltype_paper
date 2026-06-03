#!/usr/bin/env python
"""
04_gene_drivers.py — Identify genes driving SCZ enrichment per cell type.

For each independent cell type, computes gene-level contribution scores
(specificity × GWAS_Z), uniqueness scores, and shared vs unique gene
decomposition across cell types.

Input:  Specificity matrix, GWAS data, independent types (from step 03)
Output: results/tables/scz_enrichment_genes_{type}.csv (per type)
        results/tables/jaccard_similarity.csv
        results/tables/shared_pan_interneuron_genes.csv
        results/figures/shared_vs_unique_scz_signal.png
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
import numpy as np

from scz_celltype_enrichment.config import (
    SPECIFICITY_CSV, SCZ_GWAS_GENE_SET_CSV,
    TABLES_DIR, FIGURES_DIR, INTERMEDIATES_DIR, TOP_GENES_PER_TYPE,
)
from scz_celltype_enrichment.enrichment.gwas import load_scz_gwas_gene_set
from scz_celltype_enrichment.enrichment.gene_drivers import (
    compute_gene_contributions,
    compute_uniqueness_scores,
    get_top_gene_sets,
    compute_jaccard_matrix,
    identify_shared_genes,
    identify_type_specific_genes,
)
from scz_celltype_enrichment.utils import classify_supertypes_by_class
from scz_celltype_enrichment.plotting import plot_shared_vs_unique
from scz_celltype_enrichment.utils import ensure_dirs


def main():
    ensure_dirs(TABLES_DIR, FIGURES_DIR)

    print("=" * 60)
    print("Step 04: Gene Driver Analysis")
    print("=" * 60)

    t0 = time.time()

    # Load data
    print("\n1. Loading data...")
    specificity = pd.read_csv(str(SPECIFICITY_CSV), index_col=0)

    # Load GWAS data from intermediates
    gwas_df = pd.read_csv(str(INTERMEDIATES_DIR / "gwas_matched.csv"))
    genes = gwas_df["SYMBOL"].values
    spec_matched = specificity.loc[genes]

    # Load independent types
    independent_df = pd.read_csv(str(TABLES_DIR / "independent_supertypes_conditional.csv"))
    independent_types = independent_df["supertype"].tolist()
    print(f"  Independent types: {len(independent_types)}")

    # Load SCZ GWAS gene set
    scz_genes = load_scz_gwas_gene_set(SCZ_GWAS_GENE_SET_CSV)

    # Per-type gene contributions
    print("\n2. Computing gene contributions per type...")
    gene_contributions = {}
    for st in independent_types:
        gc = compute_gene_contributions(spec_matched, gwas_df, st, scz_gene_set=scz_genes)
        gene_contributions[st] = gc

        # Save per-type file
        safe_name = st.replace(" ", "_").replace("/", "_")
        output_csv = TABLES_DIR / f"scz_enrichment_genes_{safe_name}.csv"
        gc.to_csv(output_csv, index=False)

        # Print top 5
        top5 = gc.head(5)
        print(f"\n  {st}:")
        for _, row in top5.iterrows():
            gwas_flag = " *SCZ*" if row.get("is_scz_gwas", False) else ""
            print(f"    {row['gene']:15s}  spec={row['specificity']:.4f}  "
                  f"Z={row['gwas_zstat']:.2f}  contrib={row['contribution']:.4f}{gwas_flag}")

    print(f"\n  Saved {len(independent_types)} per-type gene files")

    # Uniqueness scores
    print("\n3. Computing uniqueness scores...")
    uniqueness = compute_uniqueness_scores(
        spec_matched, gwas_df, independent_types, top_n=TOP_GENES_PER_TYPE
    )

    # Top gene sets and Jaccard similarity
    print("\n4. Computing gene set overlap (Jaccard)...")
    gene_sets = get_top_gene_sets(
        spec_matched, gwas_df, independent_types, top_n=TOP_GENES_PER_TYPE
    )
    jaccard = compute_jaccard_matrix(gene_sets)
    jaccard.to_csv(TABLES_DIR / "jaccard_similarity.csv")
    print(f"  Saved: jaccard_similarity.csv")

    # Print mean Jaccard by GABAergic vs others
    class_map = classify_supertypes_by_class(independent_types)
    gaba_types = [t for t in independent_types if class_map.get(t) == "GABAergic"]
    non_gaba_types = [t for t in independent_types if class_map.get(t) != "GABAergic"]
    if len(gaba_types) >= 2:
        gaba_jaccard = jaccard.loc[gaba_types, gaba_types]
        within_gaba = gaba_jaccard.values[np.triu_indices(len(gaba_types), k=1)]
        print(f"  Mean within-GABAergic Jaccard: {within_gaba.mean():.3f}")

    # Shared genes
    print("\n5. Identifying shared genes...")
    shared_genes = identify_shared_genes(gene_sets, min_types=len(independent_types) // 2)
    shared_genes.to_csv(TABLES_DIR / "shared_pan_interneuron_genes.csv", index=False)
    print(f"  Found {len(shared_genes)} genes shared across >= {len(independent_types)//2} types")

    # SST-specific genes
    sst_types = [t for t in independent_types if t.startswith("Sst")]
    non_sst_types = [t for t in independent_types if not t.startswith("Sst")]
    if sst_types and non_sst_types:
        sst_specific = identify_type_specific_genes(gene_sets, sst_types, non_sst_types)
        if len(sst_specific) > 0:
            sst_specific.to_csv(TABLES_DIR / "sst_specific_genes.csv", index=False)
            print(f"  Found {len(sst_specific)} SST-specific genes")

    # Generate figure
    print("\n6. Generating figure...")
    plot_shared_vs_unique(
        jaccard, shared_genes, FIGURES_DIR / "shared_vs_unique_scz_signal.png"
    )

    print(f"\nStep 04 complete. Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
