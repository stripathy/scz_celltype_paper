# spatial/ — Xenium spatial transcriptomics of SCZ DLPFC

Cell typing, cortical-depth inference and laminar segmentation for 24 Xenium
DLPFC sections (12 SCZ, 12 control; Kwon et al. 2026, LIBD), plus the
compositional and differential-expression analyses built on them.

This component is described in the paper as **Supplementary Methods SM1**, and
it produces **Supplementary Figs. S2 and S3**. It also supplies the data behind
the Xenium panels of Figures 1, 2 and 3 — but of those, only Figure 2 is
rendered in this repo (see *What is not here*).

## What this component produces

| Output | Feeds |
|---|---|
| `output/de/de_results_{subclass,supertype}.csv` | Figure 2 (cross-platform DE) and Supplementary Fig. S8 |
| `output/crumblr/crumblr_{input,results}_*.csv` | Figure 3 (Xenium compositional replication) and Supplementary Fig. S6 |
| `output/depth_platform/supertype_depth_platform_summary.csv` | Supplementary Fig. S3, and the S8 strata definition |
| `output/depth_proportions/proposed_layer_boundaries.csv` | layer gutters in Supplementary Fig. S8 |
| `manuscript/figures/supplementary/S02_xenium_celltype_annotation` | Supplementary Fig. S2 |
| `manuscript/figures/supplementary/S03_xenium_merfish_concordance` | Supplementary Fig. S3 |
| the annotated h5ads (`corr_subclass`, `predicted_norm_depth`, `layer`) | Figure 1e–f |

## Pipeline

Run from `spatial/`. Steps 00–05 need the raw Xenium data and the SEA-AD
reference; everything downstream runs from their outputs.

| Step | Script | What |
|---|---|---|
| 00 | `code/pipeline/00_create_h5ad.py` | Raw Xenium `.h5` + boundaries → per-sample h5ad |
| 01 | `code/pipeline/01_run_qc.py` | QC flags (`qc_pass`); nothing is ever removed, only flagged |
| 02 | `code/pipeline/02_run_mapmycells.py` | MapMyCells/HANN labels against the SEA-AD MTG taxonomy |
| 02b | `code/pipeline/02b_run_correlation_classifier.py` | Two-stage correlation classifier → `corr_subclass`, `corr_supertype`, margin QC, doublet flags |
| 04 | `code/pipeline/04_run_depth_prediction.py` | Neighbourhood-composition depth model (`K = 50`, three smallest donors held out) → `predicted_norm_depth` |
| 05 | `code/pipeline/05_run_spatial_domains.py` | BANKSY spatial domains + depth→layer binning + smoothing |

`code/pipeline/create_snrnaseq_reference.py` builds the five-donor SEA-AD snRNA-seq
reference used by the panel-resolvability benchmark.

## Disease analyses

```bash
python3 code/analysis/build_de_input.py       # pseudobulk, subclass + supertype
Rscript  code/analysis/run_de.R subclass      # edgeR QLF -> output/de/
Rscript  code/analysis/run_de.R supertype
python3 code/analysis/build_crumblr_input.py  # per-donor counts, 4 strata
Rscript  code/analysis/run_crumblr.R          # crumblr + dream -> output/crumblr/
python3 code/analysis/derive_layer_boundaries.py
```

Cell definition is identical across both: cortical, `qc_pass & corr_qc_pass`,
24 donors. Model is `~ diagnosis + sex + age + pmi`.

## Supplementary figures

```bash
# S2 — annotation agreement, marker dot plots, panel resolvability
python3 code/analysis/integrate_lieber_annotations.py
python3 code/analysis/plot_subclass_marker_dotplot.py
python3 code/analysis/plot_sst_supertype_depth_dotplot.py
python3 code/analysis/resolvability_report.py
Rscript  code/analysis/plot_markers_resolvability_combined.R

# S3 — Xenium vs MERFISH concordance
python3 code/analysis/build_celltyping_validation_data.py
python3 code/analysis/build_supertype_depth_platform_data.py
Rscript  code/analysis/plot_xenium_merfish_composite.R
```

Both write straight into `manuscript/figures/supplementary/`. Their small input
CSVs are committed (force-added past the `output/` ignore rule), so the two
renderers run from a clean clone with no external data. **Run them from this
directory** (`cd spatial && Rscript code/analysis/<script>.R`) — unlike the
renderers in `genetics/`, they resolve `source()` and their inputs relative to
the working directory, not to their own location. They also need the `ragg`
package, without which `ggsave` fails after the figure is already built.

`_xenium_exclusions.R` defines the types too rare to test, dropped from every
Xenium figure.

## Modules

### `depth_model.py`
MERFISH-trained cortical depth prediction from K=50 neighborhood composition features.
- GradientBoostingRegressor; R² ≈ 0.897 on the three held-out donors (the deployed split). The GroupKFold out-of-fold R² quoted in SM1 is a different, stricter estimate — see the validation table below.
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

### `loading.py`
Data loading for 10x Xenium .h5 files with cell boundary CSV centroids.

### `metadata.py`
Subject metadata loading (handles non-standard Excel XML format).

## Validation (the SM1 numbers)

`code/analysis/validation/` holds the checks the Supplementary Methods quote.

| Script | Number it backs |
|---|---|
| `validate_depth_model_cv.py` | GroupKFold out-of-fold R² = 0.889; subclass-only anti-circularity baseline R² = 0.32 |
| `validate_depth_marker_anchor.py` | MERFISH-free check: eight laminar markers ordered correctly by predicted depth |
| `benchmark_merfish_classification.py` | 84.9% subclass accuracy against MERFISH ground truth |
| `plot_merfish_xenium_benchmark.py` | correlation classifier r = 0.80 vs Harmony r = 0.73 on MERFISH proportions |
| `02c_run_harmony_transfer.py`, `compare_harmony_vs_corr.py` | why the centroid classifier replaced Harmony+kNN: Sst inflated to 12.1% vs ~2.5% expected, VLMC called Sst 56% of the time, 69% agreement |


## What is not here

- **Figures 1 and 3 themselves.** Figure 1a and Figure 3 are rendered by
  [`snrnaseq/Final_figures/`](../snrnaseq/Final_figures/README.md), on the
  cluster rather than into `manuscript/figures/`; Figure 1b–f is assembled by
  hand. This component supplies the annotated h5ads and the crumblr tables they
  draw on. Note that Fig. 3's Xenium panels read a cluster-side copy of those
  results, not `output/crumblr/` directly — and panels d and e read a Seurat
  conversion of the Xenium object, `Xenium_SCZ_R.rds`, whose provenance is not
  recorded anywhere (issue 3 in [`../KNOWN_ISSUES.md`](../KNOWN_ISSUES.md); it
  matters because a superseded Xenium object has MGE and CGE transposed). The
  per-supertype **mean** depth table behind Fig. 3f is likewise not produced
  here; the committed depth-platform summary carries medians.
- **The interactive cell browser.** The viewer front end, its Playwright suite
  and the three export steps that feed it (`03`, `06`, `07`) live in the
  upstream `SCZ_Xenium` repo. They are not part of the paper.
- **Exploratory work**: density and depth-stratified composition, presentation
  figures, MERSCOPE, probe-panel design, nuclear-doublet resolution, and
  `K = 100` low-CPS depth-model variants. None was wired into the pipeline.

`plot_supertype_depth_casecontrol.R` and its data builder are kept even though
the figure is not in the paper: it is the null control showing that the
compositional result is not a depth artifact. It writes to
`manuscript/figures/not_in_current_version/`.

## Key design decisions

See [`methods_writeup.md`](methods_writeup.md) for the full account. In brief:
two-stage cell typing because Harmony+kNN inflated interneuron proportions;
depth from spatial neighbours rather than expression, so depth is not
circular with cell type; BANKSY for spatial domains; flag-based QC that never
removes cells; and external validation against the SEA-AD MERFISH atlas at
every step.

Object schema: [`all_samples_annotated_guide.md`](all_samples_annotated_guide.md).
Cross-platform numbers: [`cross_platform_concordance.md`](cross_platform_concordance.md).
