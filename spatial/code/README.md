# SCZ Xenium Spatial Transcriptomics Pipeline

Analysis pipeline for 24 Xenium spatial transcriptomics samples from schizophrenia (SCZ) and control brains. Uses the SEA-AD MTG taxonomy via MapMyCells for cell type annotation and a MERFISH-trained depth model for cortical depth estimation.

## Pipeline Steps

The pipeline is modular — each step reads/updates per-sample h5ad files:

| Step | Script | Description |
|------|--------|-------------|
| 00 | `pipeline/00_create_h5ad.py` | Load raw Xenium .h5 + cell boundary CSV → h5ad with spatial coords |
| 01 | `pipeline/01_run_qc.py` | Cell-level QC (Kwon et al. approach), adds `qc_pass` column |
| 02 | `pipeline/02_run_mapmycells.py` | MapMyCells hierarchical annotation → `class_label`, `subclass_label`, `supertype_label` |
| 02b | `pipeline/02b_run_correlation_classifier.py` | Two-stage Pearson correlation reclassification + doublet detection |
| 03 | `pipeline/03_export_transcripts.py` | Export per-gene transcript coordinates (for viewer) |
| 04 | `pipeline/04_run_depth_prediction.py` | Retrain MERFISH depth model, predict cortical depth → `predicted_norm_depth` |
| 05 | `pipeline/05_run_spatial_domains.py` | BANKSY spatial domain classification + layer assignment + spatial smoothing → `layer` column |
| 06 | `pipeline/06_export_viewer.py` | Export JSON for interactive HTML viewer |
| 07 | `pipeline/07_export_boundaries.py` | Export cell + nucleus boundary polygons for viewer |

> **Note:** Nuclear doublet resolution is an optional side investigation and is not part of the main pipeline. See `nuclear_resolution/README.md` for details.

## Samples

24 Xenium samples: 12 SCZ, 12 Control (1,339,151 total cells, 300-gene panel)
- **Br2039** (SCZ) is WM-heavy (65%) but included in all analyses — improves snRNAseq concordance

## Modules (`code/modules/`)

### `depth_model.py`
MERFISH-trained cortical depth prediction from K=50 neighborhood composition features.
- GradientBoostingRegressor, test R² ≈ 0.897 on held-out donors
- Predictions NOT clamped to [0,1] (cells outside cortex can be < 0 or > 1)
- `smooth_layers_spatial()`: 3-step spatial layer smoothing (within-domain majority vote, vascular border trim, BANKSY-anchored L1 contiguity)

### `banksy_domains.py`
BANKSY-based spatial domain classification (replaces older K-NN Leiden approach).
- BANKSY clustering (λ=0.8, res=0.3) for spatially coherent domains
- Classifies: Cortical, Vascular (>50% Endo+VLMC), WM (>40% Oligo + deep)
- L1 border detection: shallow non-neuronal clusters correctly identified as L1 cortex
- Used by pipeline step 05

### `spatial_domains.py` (legacy)
Original K-NN composition → PCA → Leiden domain classifier. Superseded by
`banksy_domains.py`. Retained for `VASCULAR_TYPES` and `NON_NEURONAL_TYPES` constants.

### `analysis.py`
Depth-stratified SCZ vs Control comparison with MERFISH validation.
- Cell type fractions per sample per depth stratum
- Mann-Whitney U tests with FDR correction
- Validation against MERFISH reference proportions

### `loading.py`
Data loading for 10x Xenium .h5 files with cell boundary CSV centroids.

### `metadata.py`
Subject metadata loading (handles non-standard Excel XML format).

### `plotting.py`
Spatial visualization with rasterized color-blended images (dark backgrounds, 20μm bins).

## Key Output Files

```
output/
  h5ad/                           # Per-sample annotated h5ad files
    {sample_id}_annotated.h5ad    # .obs: sample_id, class_label, subclass_label,
                                  #        supertype_label, predicted_norm_depth,
                                  #        qc_pass, layer, layer_unsmoothed
  all_samples_annotated.h5ad      # Combined dataset (1.34M cells, ~1.3 GB)
  viewer/                         # Interactive HTML viewer
    index.html                    # Multi-sample spatial viewer
    xenium_viewer_standalone.html # Self-contained standalone version
    {sample_id}.json              # Per-sample cell data (compact)
    index.json                    # Global metadata + color palettes
  depth_model_normalized.pkl      # Trained depth prediction model
  qc_summary.csv                  # Per-sample QC statistics
  spatial_domain_summary.csv      # Per-sample domain classification & layer counts
```

## Cell Type Taxonomy

Uses the SEA-AD MTG hierarchy (via MapMyCells):
- **Class** (3): Glutamatergic, GABAergic, Non-neuronal
- **Subclass** (24): L2/3 IT, L4 IT, ..., Sst, Pvalb, ..., Astrocyte, Oligodendrocyte, ...
- **Supertype** (127): Finest level, named `{Subclass}_{N}` (e.g., `Sst_5`, `L2/3 IT_3`, `Astro_1`)

## Depth Strata for Analyses

| Stratum | Depth Range | Description |
|---------|------------|-------------|
| L1 | < 0.10 | Molecular layer |
| L2/3 | 0.10 - 0.40 | Upper cortical layers |
| L4 | 0.40 - 0.55 | Granular layer |
| L5 | 0.55 - 0.70 | Deep output layer |
| L6 | 0.70 - 0.90 | Deep cortical layers |
| WM | > 0.90 | White matter |

## Dependencies

```
anndata, scanpy, numpy, scipy, scikit-learn,
matplotlib, pandas, statsmodels, h5py, openpyxl,
adjustText, cell_type_mapper
# See environment.yml for full list
```

## Data Requirements

- Xenium data: `GSM*-cell_feature_matrix.h5` + `GSM*-cell_boundaries.csv.gz`
- SEA-AD MERFISH: `SEAAD_MTG_MERFISH.2024-12-11.h5ad` (depth model training)
- SEA-AD snRNAseq: `seaad_mtg_snrnaseq_reference.h5ad` (correlation classifier reference)
- MapMyCells precomputed stats: `precomputed_stats.20231120.sea_ad.MTG.h5`
- Subject metadata: `data/sample_metadata.xlsx`

## Nuclear Resolution (`code/nuclear_resolution/`)

Optional side investigation into using nuclear-only transcript counts to arbitrate doublet calls. Empirically shown to have negligible impact on downstream compositional analysis — the simplified QC pipeline (spatial QC + 5th-percentile margin filter + doublet exclusion) produces equivalent results. See `nuclear_resolution/README.md` for details.

## Archive

Legacy code is preserved in `code/archive/` for reference:
- `stale_analysis/` — Archived diagnostic and analysis scripts (harmony transfer, edgepython DE, nsforest markers, calibration, crumblr builders, proportion comparisons). Removed dependencies: harmonypy, edgepython, fitz, nsforest, markdown
- `label_transfer.py` — Old kNN-based label transfer (superseded by MapMyCells)
- `layers.py` — Old density-based layer segmentation (superseded by depth model)
- `legacy_runners/` — Old monolithic pipeline runners
- `ood_methods/` — Exploratory OOD detection approaches (superseded by BANKSY domains)
- `spatial_domain_exploration/` — Experimental spatial domain scripts
- `banksy_exploration/` — BANKSY parameter tuning, validation, and batch runner (logic now in `modules/banksy_domains.py`)
- `curved_strips/` — Curved cortex strip identification pipeline (experimental)
