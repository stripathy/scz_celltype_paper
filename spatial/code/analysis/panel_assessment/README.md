# Panel Assessment & Supertype Classification Confidence

## Motivation

Spatial transcriptomics panels vary dramatically in their ability to resolve cell type supertypes within subclasses. Our analysis of the Kwon Xenium dataset (300-gene custom panel) revealed that supertype-level composition results are unreliable when the panel lacks within-subclass discriminating markers — and that this is a universal problem, not specific to SST or L6b supertypes.

We systematically assessed six spatial panels (Kwon 300, Xenium v1, MERFISH 180, MERSCOPE 250, Xenium 5K, MERSCOPE 4K) and found that **gene count does not determine supertype resolution**: MERSCOPE 250 (250 genes) outperforms Xenium 5K (5,001 genes) because its genes were specifically curated for within-subclass discrimination.

This directory contains scripts for:

1. **Cross-platform calibration**: Which genes survive the snRNAseq → spatial transition?
2. **Within-subclass marker identification**: Wilcoxon markers that distinguish sibling supertypes
3. **Panel coverage assessment**: How many discriminating markers does each panel contain?
4. **Classification confidence rating**: HIGH/MEDIUM/LOW ratings for each supertype on each panel
5. **Optimized add-on gene selection**: Greedy selection of spatial-validated genes to maximize supertype coverage
6. **Classification fragility analysis**: Direct evidence of misclassification in the Kwon dataset

## Scripts

### Phase 1: Cross-Platform Calibration

| Script | Purpose |
|--------|---------|
| `compare_snrnaseq_merscope_expression.py` | Pseudobulk correlations between snRNAseq and MERSCOPE 4K by cell type and gene |
| `cross_platform_gene_corr.py` | Compare gene-level correlations across MERSCOPE 4K, SEA-AD MERFISH, and Xenium |
| `characterize_poor_genes.py` | Identify gene properties (biotype, expression, specificity) that predict poor spatial performance |

### Phase 2: Marker Discovery & Panel Coverage

| Script | Purpose |
|--------|---------|
| `merscope_panel_assessment.py` | Panel overlap matrices, supertype marker coverage, detection rate calibration |
| `supertype_markers_panel_overlap.py` | Wilcoxon markers at supertype level; check which are missing from each panel |
| `supertype_classification_confidence.py` | Within-subclass Wilcoxon markers for all 126 supertypes; confidence ratings combining marker coverage + layer specificity |
| `cross_platform_marker_adequacy.py` | Compare within-subclass marker coverage across 6 spatial panels; generates heatmap and bar chart figures |

### Phase 3: Add-On Gene Selection

| Script | Purpose |
|--------|---------|
| `hierarchical_probe_selection.py` | Original two-round hierarchical probe design with dropout-aware scoring (SST/L6b focused) |
| `compute_addon_gap.py` | Greedy gene selection to match target panel coverage; diminishing-returns analysis |

### Phase 4: Classification Fragility Analysis

| Script | Purpose |
|--------|---------|
| `investigate_sst_confusion.py` | Classification margin analysis by diagnosis; per-gene SCZ effects within Sst supertypes |
| `prototype_sst_neighborhood_classifier.py` | Tests spatial neighborhood smoothing and label voting as rescue strategies (negative result: does not improve classification) |

### Phase 5: Visualization

| Script | Purpose |
|--------|---------|
| `visualize_marker_quality.py` | Dotplots, matrixplots, violin plots of selected markers |

## Key Outputs (in `output/marker_analysis/`)

**Primary report:**
- `XENIUM_PANEL_DESIGN_AND_SUPERTYPE_CLASSIFICATION.md` — Comprehensive panel design recommendation with cross-platform evidence and cost analysis

**Add-on gene lists:**
- `v1_addon_100_spatial_validated.csv` — Recommended 100-gene add-on list (all FISH-validated)
- `v1_addon_spatial_only_curve.csv` — Coverage progression at each gene addition

**Within-subclass markers:**
- `within_subclass_markers_all.csv` — Wilcoxon markers for all 126 supertypes (snRNAseq-derived)
- `within_subclass_markers_merfish.csv` — Independent markers computed from MERFISH spatial data

**Cross-platform assessment:**
- `cross_platform_marker_coverage.csv` — Coverage of top-N markers across 6 panels
- `supertype_classification_confidence.csv` — Per-supertype confidence rating

**Supplementary:**
- `addon_gap_analysis.csv` — Per-supertype deficit analysis
- `addon_gene_merfish_validation.csv` — MERFISH validation of add-on gene detectability

## Dependencies

- Modules: `code/modules/panel_utils.py`, `code/modules/reference_utils.py`, `code/modules/gene_properties.py`, `code/modules/pseudobulk.py`, `code/modules/correlation_classifier.py`
- Data: `data/merscope_4k_probe_testing/`, `data/reference/` (SEA-AD snRNAseq + MERFISH h5ad files)
- Panel metadata: `~/Downloads/Xenium_hBrain_v1_metadata.csv`, `~/Downloads/XeniumPrimeHuman5Kpan_tissue_pathways_metadata.csv`
