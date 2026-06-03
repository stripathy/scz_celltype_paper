"""
Shared configuration for the SCZ Xenium processing pipeline.

Centralizes all paths, constants, and settings used across pipeline steps
(00_create_h5ad through 07_export_boundaries). Import from here instead of
hardcoding paths in each script.

This is the pipeline equivalent of code/analysis/config.py.
"""

import os

# ──────────────────────────────────────────────────────────────────────
# Base paths
# ──────────────────────────────────────────────────────────────────────

BASE_DIR = os.path.expanduser("~/Github/SCZ_Xenium")

# Input data
RAW_DIR = os.path.join(BASE_DIR, "data", "raw")
METADATA_PATH = os.path.join(BASE_DIR, "data", "sample_metadata.xlsx")

# Reference data (SEA-AD)
# MERFISH spatial reference — REQUIRED for depth model training (step 04)
MERFISH_PATH = os.path.join(BASE_DIR, "data", "reference",
                            "SEAAD_MTG_MERFISH.2024-12-11.h5ad")
# snRNAseq reference — NOT required for core pipeline; used for validation only
SNRNASEQ_REF_PATH = os.path.join(BASE_DIR, "data", "reference",
                                  "seaad_mtg_snrnaseq_reference.h5ad")
# MapMyCells precomputed stats — REQUIRED for cell type annotation (step 02)
PRECOMPUTED_STATS_PATH = os.path.join(
    BASE_DIR, "data", "reference",
    "precomputed_stats.20231120.sea_ad.MTG.h5"
)

# Reference availability flags
HAS_MERFISH_REF = os.path.exists(MERFISH_PATH)
HAS_SNRNASEQ_REF = os.path.exists(SNRNASEQ_REF_PATH)
HAS_PRECOMPUTED_STATS = os.path.exists(PRECOMPUTED_STATS_PATH)
GENE_MAPPING_PATH = os.path.join(
    BASE_DIR, "data", "reference",
    "gene_symbol_to_ensembl.json"
)
# Comprehensive gene mapping (36K genes, built from snRNAseq reference + precomputed stats)
GENE_MAPPING_COMPREHENSIVE_PATH = os.path.join(
    BASE_DIR, "data", "reference",
    "gene_symbol_to_ensembl_comprehensive.json"
)

# Output directories
OUTPUT_DIR = os.path.join(BASE_DIR, "output")
H5AD_DIR = os.path.join(OUTPUT_DIR, "h5ad")
VIEWER_DIR = os.path.join(OUTPUT_DIR, "viewer")
PLOTS_DIR = os.path.join(OUTPUT_DIR, "plots")

# Model artifacts
DEPTH_MODEL_PATH = os.path.join(OUTPUT_DIR, "depth_model_normalized.pkl")
CENTROID_PATH = os.path.join(H5AD_DIR, "correlation_centroids.pkl")

# ──────────────────────────────────────────────────────────────────────
# Processing settings
# ──────────────────────────────────────────────────────────────────────

N_WORKERS = 4  # for multiprocessing in steps 04/05

# MapMyCells settings (step 02)
MAPMYCELLS_BOOTSTRAP_ITER = 100
MAPMYCELLS_BOOTSTRAP_FACTOR = 0.5
MAPMYCELLS_N_PER_UTILITY = 30

# Taxonomy levels assigned by MapMyCells
TAXONOMY_LEVELS = ["class", "subclass", "supertype"]

# Correlation classifier settings (step 02b)
CORR_CLASSIFIER_TOP_N = 100       # exemplar cells per type for centroid building
CORR_CLASSIFIER_QC_PERCENTILE = 5.0  # bottom % of margin to flag per sample
L6B_MARGIN_THRESHOLD = 0.02           # absolute margin threshold for L6b cells

# Nuclear doublet resolution settings (optional nuclear resolution)
NUCLEAR_CHUNK_SIZE = 500_000      # transcripts per STRtree query batch
NUCLEAR_MIN_UMI = 50              # min nuclear UMI for reliable doublet assessment
NUCLEAR_MIN_CORR = 0.3            # min correlation for high-confidence resolution

# Transcript export directory (step 03 output, used by optional nuclear resolution)
TRANSCRIPT_DIR = os.path.join(VIEWER_DIR, "transcripts")

# ──────────────────────────────────────────────────────────────────────
# Modules path helper
# ──────────────────────────────────────────────────────────────────────

PIPELINE_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.dirname(PIPELINE_DIR)
MODULES_DIR = os.path.join(CODE_DIR, "modules")
