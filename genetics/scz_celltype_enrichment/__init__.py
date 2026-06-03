"""
scz_celltype_enrichment — Cell-type enrichment analysis for schizophrenia GWAS.

This package implements a MAGMA-style gene property analysis pipeline to identify
brain cell types enriched for schizophrenia genetic risk using SEA-AD snRNA-seq
reference data, and links these findings to patch-seq electrophysiology.

Subpackages:
  enrichment/   — Core GWAS enrichment pipeline (Parts 1-4)
  integration/  — Downstream integration analyses (Parts 5-8)

Top-level modules:
  config.py     — Centralized configuration (paths, thresholds, figure styles)
  utils.py      — Shared statistical and I/O utilities
  plotting.py   — Publication-quality figure generation

See scripts/ for executable pipeline scripts.
"""

__version__ = "0.1.0"
