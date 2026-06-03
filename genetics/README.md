# SCZ Cell-Type GWAS Enrichment Pipeline

Linking schizophrenia (SCZ) genetic risk to specific brain cell types using MAGMA-style gene property analysis on single-cell RNA-seq transcriptomic taxonomies.

## Overview

This pipeline tests whether genes specific to particular brain cell types carry excess schizophrenia GWAS signal. Starting from cell-type specificity scores derived from single-nucleus RNA-seq reference datasets, it runs OLS regression of GWAS Z-scores against specificity for each cell type, then layers on spatial, electrophysiological, chromatin accessibility, and case-control compositional data to build a multi-modal picture of which cell types are most implicated in SCZ.

### Key Findings

- **GABAergic interneurons** carry the vast majority of SCZ genetic enrichment — 41/44 FDR-significant SEA-AD cell types are GABAergic
- **SST interneurons** are the most enriched subclass, with multiple independently significant subtypes
- **Upper-layer SST subtypes** (cortical layers 1–3) show disproportionately strong enrichment (Spearman r = −0.50 between depth and enrichment)
- **SST subtypes most enriched for GWAS signal are also the most depleted in SCZ post-mortem brains** (r = 0.65, p = 3.85e-3; Endresz et al., in prep)
- Extending to a **503-type combined taxonomy** (SEA-AD + Siletti whole-brain atlas) reveals additional enrichment in subcortical MGE interneuron populations — see [`franken_taxonomy/`](franken_taxonomy/) for the canonical type list and RBH dedup mapping that downstream projects can consume

## Data Sources

| Dataset | Description | Reference |
|---------|-------------|-----------|
| **SEA-AD snRNA-seq** | 137,303 nuclei, 36,601 genes, 137 supertypes from middle temporal gyrus | [Gabitto et al. 2024, *Nature Neuroscience*](https://doi.org/10.1038/s41593-024-01774-5) |
| **PGC3 SCZ GWAS** | Gene-level results from 76,755 cases / 243,649 controls (MAGMA) | [Trubetskoy et al. 2022, *Nature*](https://doi.org/10.1038/s41586-022-04434-5) |
| **SEA-AD MERFISH** | Spatial transcriptomics with cortical depth per cell type (0 = pial, 1 = WM) | [Gabitto et al. 2024](https://doi.org/10.1038/s41593-024-01774-5) |
| **Siletti whole-brain atlas** | 461 clusters spanning neocortex, hippocampus, thalamus, amygdala, and more | [Siletti et al. 2023, *Science*](https://doi.org/10.1126/science.add7046) |
| **SEA-AD snATAC-seq** | 516K nuclei, 219K peaks — chromatin accessibility per cell type | [Gabitto et al. 2024](https://doi.org/10.1038/s41593-024-01774-5) |
| **Patch-seq electrophysiology** | Sag, tau, and other intrinsic properties per transcriptomic type | [Lee & Bhatt Dalley et al.](https://portal.brain-map.org/) |
| **Case-control snRNA-seq composition** | 7-cohort meta-analysis of cell-type proportion changes in SCZ (crumblr model) | Endresz et al., in prep |
| **SCZ fine-mapping** | FINEMAP credible sets from PGC3 SCZ GWAS | [Trubetskoy et al. 2022](https://doi.org/10.1038/s41586-022-04434-5) |

## Pipeline

The analysis consists of 13 scripts run sequentially. Scripts 01–07 use the SEA-AD taxonomy; scripts 08–13 extend to the Siletti atlas and multi-modal integration.

```
python scripts/run_all.py          # Run core pipeline (steps 01-07)
python scripts/run_all.py --all    # Run full pipeline (steps 01-13)
```

| Step | Script | Description | Key Output |
|------|--------|-------------|------------|
| 01 | `01_compute_specificity.py` | Cell-type specificity from SEA-AD snRNA-seq (column-normalize, then row-normalize) | Specificity matrix (genes × 137 types) |
| 02 | `02_magma_enrichment.py` | MAGMA-style OLS: `GWAS_Z ~ specificity + log(gene_size) + log(n_snps) + log(N)` | Enrichment p-values per type |
| 03 | `03_conditional_analysis.py` | Forward selection for independently enriched types | 18 independent types |
| 04 | `04_gene_drivers.py` | Gene contribution scores (specificity × GWAS_Z), Jaccard similarity | Driver gene lists per type |
| 05 | `05_spatial_layer_analysis.py` | SST layer classification using MERFISH depth, upper-layer gene scores | Layer info + spatial gene rankings |
| 06 | `06_gene_ephys_correlations.py` | Gene–electrophysiology correlations (sag, tau) from patch-seq | Sag/tau gene correlation tables |
| 07 | `07_make_all_figures.py` | Regenerate all figures from saved CSVs | Publication figures |
| 08 | `08_siletti_enrichment.py` | Enrichment on Siletti 461-cluster atlas (Conti vs reprocessed) | Siletti enrichment tables |
| 09 | `09_metaneighbor_integration.py` | MetaNeighbor reciprocal best hits: SEA-AD ↔ Siletti matching | 95 matched type pairs |
| 10 | `10_combined_taxonomy.py` | Build 503-type RBH combined taxonomy + enrichment | Combined enrichment (503 types) |
| 11 | `11_gene_driver_scatter.py` | Gene driver scatter plots for top enriched types | Scatter figures per type |
| 12 | `12_atac_peak_overlap.py` | ATAC-seq peak overlap with FINEMAP credible sets | PIP-weighted accessibility scores |
| 13 | `13_gwas_vs_composition.py` | Correlate GWAS enrichment with case-control composition changes | Correlation statistics + figure |

## Repository Structure

```
scz_cell_type_enrichment/
├── scz_celltype_enrichment/        # Python package
│   ├── config.py                   # Central configuration (paths, thresholds, figure styles)
│   ├── utils.py                    # OLS regression, FDR correction, utilities
│   ├── enrichment/                 # Core GWAS enrichment pipeline
│   │   ├── specificity.py          #   Cell-type specificity computation
│   │   ├── gwas.py                 #   MAGMA gene-level results loading, ENTREZ mapping
│   │   ├── celltype.py             #   OLS regression per cell type
│   │   ├── conditional.py          #   Forward selection algorithm
│   │   └── gene_drivers.py         #   Gene contribution scores, Jaccard similarity
│   ├── integration/                # Multi-modal integration
│   │   ├── spatial.py              #   Cortical depth classification (MERFISH)
│   │   └── gene_ephys.py           #   Gene–electrophysiology correlations
│   ├── siletti/                    # Siletti atlas processing
│   │   ├── specificity.py          #   Ensembl→ENTREZ mapping, specificity
│   │   ├── clusters.py             #   Cluster metadata
│   │   └── metaneighbor.py         #   HVG selection, AUROC, reciprocal best hits
│   └── plotting/                   # Publication figures
│       ├── enrichment.py           #   Core enrichment bar plots, Manhattan
│       ├── sst_spatial.py          #   SST-focused panels, depth volcano
│       └── gene_ephys.py           #   Gene–ephys scatter plots
├── scripts/                        # Analysis pipeline (01-13)
│   ├── run_all.py                  # Master orchestrator
│   ├── generate_markdown_figures.py # Manuscript figure generation
│   └── archive/                    # Superseded scripts
├── docs/                           # Documentation
│   ├── analysis_approach.md        # Full methodology
│   ├── combined_taxonomy_approach.md # RBH combined taxonomy details
│   └── rnaseq_specificity_enrichment_results.md  # Results summary
├── results/
│   ├── tables/                     # All result CSVs
│   └── figures/                    # All figures including manuscript/
├── app.py                          # Interactive web explorer
└── data/                           # Large data files (not tracked — see below)
```

## Data Setup

The `data/` directory is not tracked in git (140GB+ of reference data). To reproduce the analysis, you need:

| File | Size | Required for |
|------|------|-------------|
| `data/seaad_reference.h5ad` | ~34 GB | Steps 01–07 (SEA-AD snRNA-seq) |
| `linking_cell_types_to_brain_phenotypes/` | 6 MB | Steps 01–02 (MAGMA gene-level results) |
| `data/spatial/median_depth_supertype.csv` | <1 MB | Step 05 (cortical depth) |
| `data/patchseq/*.csv` | <1 MB | Step 06 (electrophysiology) |
| `data/precomputed_stats.siletti.training.h5` | 8.1 GB | Steps 08–10 (Siletti atlas) — [download from Allen S3](https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/mapmycells/WHB-10Xv3/20240831/precomputed_stats.siletti.training.h5) ([browse dir](https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html#mapmycells/WHB-10Xv3/20240831/)) |
| `data/atac/SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad` | 17 GB | Step 12 (ATAC-seq) |
| `data/gwas/PGC3_SCZ_wave3...tsv.gz` | 229 MB | Step 13 (SNP-level sumstats) |

The `linking_cell_types_to_brain_phenotypes/` directory is a separate repository: [Integrative-Mental-Health-Lab/linking_cell_types_to_brain_phenotypes](https://github.com/Integrative-Mental-Health-Lab/linking_cell_types_to_brain_phenotypes). Clone it into the project root.

The case-control composition data (Endresz et al., in prep) is read from a sibling repository (`SCZ_Xenium/data/nicole_scz_snrnaseq_betas/`).

## Interactive Explorer

A web application provides interactive exploration of enrichment results across the 503-type combined taxonomy:

```bash
python app.py
# Open http://localhost:8050
```

## Dependencies

- Python 3.10+
- numpy, pandas, scipy, statsmodels
- scanpy, anndata
- matplotlib, seaborn
- loompy (for Siletti atlas)
- h5py

## Citation

If you use this pipeline, please cite:

- **PGC3 SCZ GWAS:** Trubetskoy et al. (2022). Mapping genomic loci implicates genes and synaptic biology in schizophrenia. *Nature*, 604, 502–508.
- **SEA-AD:** Gabitto et al. (2024). Integrated multimodal cell atlas of Alzheimer's disease. *Nature Neuroscience*, 27, 2366–2383.
- **Siletti atlas:** Siletti et al. (2023). Transcriptomic diversity of cell types across the adult human brain. *Science*, 382, eadd7046.
- **Case-control composition:** Endresz et al. (in prep). 7-cohort snRNA-seq meta-analysis of cell-type composition changes in schizophrenia.
