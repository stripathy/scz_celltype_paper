"""
config.py — Central configuration for the SCZ cell-type enrichment pipeline.

All file paths, analysis thresholds, and figure style constants are defined here.
Scripts and modules import from this file rather than hardcoding values.

Path resolution:
  PROJECT_ROOT is auto-detected from this file's location (scz_celltype_enrichment/
  is a direct child of the project root). All other paths are relative to it.
"""
from pathlib import Path

# ============================================================================
# Project root (auto-detected, no hardcoded absolute paths)
# ============================================================================

PROJECT_ROOT = Path(__file__).parent.parent

# ============================================================================
# Input data paths
# ============================================================================

# SEA-AD snRNA-seq reference (137,303 cells x 36,601 genes, symlink to shared_data)
SEAAD_H5AD = PROJECT_ROOT / "data" / "seaad_reference.h5ad"

# Pre-computed MAGMA gene-level results for PGC3 SCZ GWAS (18,449 genes)
MAGMA_GENES_OUT = (
    PROJECT_ROOT
    / "linking_cell_types_to_brain_phenotypes"
    / "Example_results"
    / "PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.no_heading.step2.genes.out"
)

# NCBI gene location file (ENTREZ -> SYMBOL mapping, MHC excluded, 19,175 genes)
GENE_LOC_FILE = (
    PROJECT_ROOT
    / "linking_cell_types_to_brain_phenotypes"
    / "Data"
    / "NCBI37.3.gene.loc.extendedMHCexcluded"
)

# Spatial transcriptomics depth data (MERFISH + Xenium, per supertype)
DEPTH_SUPERTYPE_CSV = PROJECT_ROOT / "data" / "spatial" / "median_depth_supertype.csv"

# SCZ GWAS gene set (genes within ±100kb of genome-wide significant loci)
SCZ_GWAS_GENE_SET_CSV = PROJECT_ROOT / "data" / "gwas" / "scz_gwas_gene_set.csv"

# PGC3 SCZ SNP-level summary statistics (for effect direction)
PGC3_SNP_SUMSTATS = PROJECT_ROOT / "data" / "gwas" / "PGC3_SCZ_wave3.european.autosome.public.v3.vcf (1).tsv.gz"

# Patch-seq data files (copies of 3 CSVs needed for gene-ephys correlations)
PATCHSEQ_DATA_DIR = PROJECT_ROOT / "data" / "patchseq"
PATCHSEQ_METADATA_CSV = PATCHSEQ_DATA_DIR / "LeeDalley_manuscript_metadata_v2.csv"
PATCHSEQ_EPHYS_CSV = PATCHSEQ_DATA_DIR / "LeeDalley_ephys_fx.csv"
SCANVI_RESULTS_CSV = PATCHSEQ_DATA_DIR / "iterative_scANVI_results_patchseq_only.2022-11-22.csv"

# ============================================================================
# Output directories
# ============================================================================

RESULTS_DIR = PROJECT_ROOT / "results"
INTERMEDIATES_DIR = RESULTS_DIR / "intermediates"
TABLES_DIR = RESULTS_DIR / "tables"
FIGURES_DIR = RESULTS_DIR / "figures"

# ============================================================================
# Cached intermediate files
# ============================================================================

SPECIFICITY_CSV = INTERMEDIATES_DIR / "seaad_supertype_specificity.csv"
MEAN_EXPRESSION_CSV = INTERMEDIATES_DIR / "seaad_supertype_mean_expression.csv"
MERGED_GWAS_CSV = INTERMEDIATES_DIR / "gwas_specificity_merged.npz"
GENE_DIRECTION_CSV = INTERMEDIATES_DIR / "gene_gwas_direction.csv"

# ============================================================================
# Analysis parameters
# ============================================================================

# Normalization
NORMALIZE_TARGET_SUM = 1e4          # scanpy normalize_total target
SPECIFICITY_COLUMN_SUM = 1000       # column normalization before row normalization

# Statistical thresholds
FDR_THRESHOLD = 0.05
BONFERRONI_ALPHA = 0.05
CONDITIONAL_P_THRESHOLD = 0.05      # p-value cutoff for forward selection

# Spatial layer classification
UPPER_LAYER_DEPTH_THRESHOLD = 0.35  # normalized depth < this = "upper"
MIDDLE_LAYER_DEPTH_THRESHOLD = 0.55 # depth < this but >= upper = "middle"

# Gene driver analysis
TOP_GENES_PER_TYPE = 200            # number of top contributing genes per cell type

# Cell type groupby column in SEA-AD h5ad
SUPERTYPE_COLUMN = "Supertype"

# ============================================================================
# Figure style constants
# ============================================================================

FIGURE_DPI = 150
TITLE_FONTSIZE = 20
LABEL_FONTSIZE = 16
TICK_FONTSIZE = 13
LEGEND_FONTSIZE = 12
