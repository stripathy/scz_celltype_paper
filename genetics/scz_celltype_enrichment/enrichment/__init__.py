"""
enrichment — Core GWAS enrichment pipeline (Parts 1-4).

Modules:
  specificity.py  — SEA-AD cell-type specificity computation
  gwas.py         — MAGMA gene-level loading, ENTREZ mapping
  celltype.py     — MAGMA-style cell-type enrichment regression
  conditional.py  — Forward selection for independent cell types
  gene_drivers.py — Gene contribution & uniqueness scores
"""
from .specificity import compute_and_save_specificity, compute_and_save_mean_expression
from .gwas import (
    load_magma_genes, load_entrez_to_symbol, map_and_merge,
    build_covariates, load_scz_gwas_gene_set, compute_gene_effect_direction,
)
from .celltype import run_enrichment_all_types, get_significant_types
from .conditional import forward_selection
from .gene_drivers import (
    compute_gene_contributions, compute_uniqueness_scores,
    get_top_gene_sets, compute_jaccard_matrix,
    identify_shared_genes, identify_type_specific_genes,
)
