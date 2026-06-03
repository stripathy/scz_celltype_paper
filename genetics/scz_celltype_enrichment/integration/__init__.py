"""
integration — Downstream integration analyses (Parts 5-8).

Modules:
  spatial.py    — Layer-stratified SST analysis using MERFISH/Xenium depth
  gene_ephys.py — Patch-seq gene-electrophysiology correlations
"""
from .spatial import (
    load_depth_data, get_sst_supertypes, classify_sst_layers,
    compute_upper_layer_gene_scores,
)
from .gene_ephys import (
    load_patchseq_data, compute_supertype_correlations,
    compute_residualized_correlations, find_discordant_genes,
    annotate_with_gwas,
)
