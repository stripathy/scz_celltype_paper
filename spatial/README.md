# SCZ Xenium Spatial Transcriptomics Analysis

## What & Why

Schizophrenia involves widespread but subtle alterations in the cellular composition of the cerebral cortex. Understanding which cell types are affected — and where in the cortical laminae those changes occur — requires measuring cell type identity and spatial position at single-cell resolution across large tissue sections. This project asks: **are cell type proportions altered in schizophrenia prefrontal cortex, and if so, which cell types and cortical layers are most affected?**

To answer this, we reanalyzed data from [Kwon et al. (2026)](https://doi.org/10.64898/2026.02.16.706214), which used Xenium spatial transcriptomics to profile 24 human dorsolateral prefrontal cortex (DLPFC) sections (12 schizophrenia, 12 control), yielding approximately 1.3 million cells measured across a 300-gene panel. Xenium provides single-cell spatial resolution over whole tissue sections, enabling both cell type classification and cortical depth modeling in the same dataset. However, 300 genes is far too few for de novo clustering to reliably separate cell types — many types that are easily distinguishable with 20,000+ genes become ambiguous at this scale. Instead, we transfer cell type labels from the Allen Institute's SEA-AD MTG taxonomy, a well-characterized reference built from single-nucleus RNA-seq of the middle temporal gyrus.

The pipeline uses a two-stage classification approach: initial labels come from MapMyCells (Allen Institute's tool implementing HANN — Hierarchical Approximate Nearest Neighbor — a bootstrapped hierarchical mapping algorithm), followed by a self-referencing Pearson correlation classifier that improves accuracy by building centroids directly from the Xenium data itself. Disease comparisons use crumblr (compositional regression using linear mixed models on centered log-ratio transformed proportions), which properly accounts for the compositional nature of cell type proportions and the repeated-measures structure of the experimental design.

**Key documents:**
- **[Cross-Platform Validation](cross_platform_concordance.md)** — Xenium vs MERFISH measurement validation (subclass r = 0.84, depth r = 0.96)
- **[SCZ Compositional Findings](scz_compositional_findings.md)** — Disease effects, snRNA-seq concordance, depth-stratified results, confidence tiers
- **[Methods: Cell Typing, Depth, and Validation](methods_writeup.md)** — Classifier, depth model, BANKSY domains, MERFISH validation
- **[Data Download Instructions](data/README.md)** — How to obtain all input datasets

## Datasets & References

This pipeline is an **analysis pipeline only** — it does not generate any of the primary datasets. All data were produced by other groups and are publicly available.

### Xenium SCZ Data (primary dataset)

- **Source:** [Kwon et al. (2026)](https://doi.org/10.64898/2026.02.16.706214) — *Mapping spatially organized molecular and genetic signatures of schizophrenia across multiple scales in human prefrontal cortex*
- **GEO Accession:** [GSE307404](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307404)
- **Description:** 24 DLPFC sections (12 SCZ, 12 control) profiled with a 300-gene Xenium panel, yielding ~1.34 million cells total

### SEA-AD Reference Datasets (Allen Institute)

Cell type annotation and cortical depth modeling use reference data from the [Seattle Alzheimer's Disease Brain Cell Atlas](https://portal.brain-map.org/) (SEA-AD), produced by the Allen Institute for Brain Science:

- **SEA-AD MERFISH** (`SEAAD_MTG_MERFISH.2024-12-11.h5ad`) — Middle temporal gyrus MERFISH spatial reference (1.9M cells, 180 genes, 27 donors, 69 sections). Used for training the cortical depth model and for cross-modal proportion validation. [Download (3.1 GB)](https://sea-ad-spatial-transcriptomics.s3.us-west-2.amazonaws.com/middle-temporal-gyrus/all_donors-h5ad/SEAAD_MTG_MERFISH.2024-12-11.h5ad)
- **SEA-AD snRNAseq** (`seaad_mtg_snrnaseq_reference.h5ad`) — SEA-AD MTG single-nucleus RNA-seq dataset, subset to the 5 neurotypical reference donors (137K cells, 36K genes). Used for ground-truth validation of doublet detection and as a full-transcriptome reference. Recreated from the full dataset via `code/pipeline/create_snrnaseq_reference.py`; original data: [Download (33.8 GB)](https://sea-ad-single-cell-profiling.s3.us-west-2.amazonaws.com/MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad)
- **MapMyCells precomputed stats** (`precomputed_stats.20231120.sea_ad.MTG.h5`) — Precomputed taxonomy statistics for hierarchical cell type mapping. [Download (251 MB)](https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/mapmycells/SEAAD/20240831/precomputed_stats.20231120.sea_ad.MTG.h5)
- **Cell type mapper:** [AllenInstitute/cell_type_mapper](https://github.com/AllenInstitute/cell_type_mapper) (requires Python 3.10+)

**References:**
- Gabitto et al. (2024) *Integrated multimodal cell atlas of Alzheimer's disease.* Nature Neuroscience. [doi:10.1038/s41593-024-01774-5](https://doi.org/10.1038/s41593-024-01774-5)
- [Allen Brain Cell Atlas portal](https://portal.brain-map.org/)

See **[data/README.md](data/README.md)** for complete download instructions.

## Pipeline Overview

The pipeline consists of 8 numbered steps in `code/pipeline/`, run sequentially. All paths and settings are centralized in `code/pipeline/pipeline_config.py`.

```
Raw .h5 + boundaries --> [00] --> initial h5ad
                         [01] --> + QC columns (qc_pass)
                         [02] --> + MapMyCells labels (HANN)
                        [02b] --> + correlation classifier labels + doublet flags (corr_qc_pass)
                         [03] --> transcript coordinates (viewer only)
                         [04] --> + depth predictions
                         [05] --> + spatial domains + layers
                         [06] --> viewer JSON + HTML (viewer only)
                         [07] --> cell + nucleus boundary polygons (viewer only)
```

| Step | Script | What it does | Why |
|------|--------|--------------|-----|
| **00** | `00_create_h5ad.py` | Creates initial AnnData h5ad files from raw Xenium `.h5` + cell boundary files. | Standardizes the raw data into a single format for all downstream processing. |
| **01** | `01_run_qc.py` | Flags cells with aberrant control probes, extreme UMI counts, or too few genes; adds `qc_pass` column. No cells are removed. | Identifies unreliable cells so downstream analyses can exclude them without losing data. |
| **02** | `02_run_mapmycells.py` | Assigns hierarchical cell type labels (class, subclass, supertype) via MapMyCells using 100-iteration bootstrapping against the SEA-AD MTG taxonomy. | Provides initial cell type identities and confidence scores as the foundation for refined classification. |
| **02b** | `02b_run_correlation_classifier.py` | Reclassifies cells using a two-stage Pearson correlation classifier built from the top-100 highest-confidence HANN exemplars per type; flags doublet suspects and low-confidence cells. | Corrects HANN misclassifications (e.g., L6b in upper layers) and provides per-cell confidence margins for QC. |
| **03** | `03_export_transcripts.py` | Exports per-gene transcript molecule coordinates from raw Xenium `transcripts.zarr` files. *Viewer only.* | Enables on-demand gene visualization in the interactive spatial viewer. |
| **04** | `04_run_depth_prediction.py` | Trains a GradientBoostingRegressor on MERFISH (Multiplexed Error-Robust Fluorescence In Situ Hybridization — a complementary spatial transcriptomics platform from the SEA-AD atlas) reference data, then predicts normalized cortical depth (0=pia, 1=white matter) for each Xenium cell. | Enables layer-specific analyses; the model uses K=50 neighborhood cell type composition rather than expression, which is more robust at 300 genes. |
| **05** | `05_run_spatial_domains.py` | Classifies tissue domains using BANKSY (a spatial domain segmentation method that augments gene expression with spatial neighbor expression to produce spatially coherent clusters); assigns and spatially smooths cortical layers. | Identifies Cortical, Vascular, and White Matter domains; provides spatially coherent layer assignments that respect tissue geometry. |
| **06** | `06_export_viewer.py` | Exports compact JSON files with coordinates, cell type labels, depth, and layers for the interactive viewer. *Viewer only.* | Generates the standalone HTML viewer for visual exploration. |
| **07** | `07_export_boundaries.py` | Exports cell and nucleus boundary polygons. *Viewer only.* | Enables boundary overlay visualization in the spatial viewer. |

**Output:** `output/h5ad/{sample}_annotated.h5ad` (24 files), `output/all_samples_annotated.h5ad` (merged), `output/viewer/xenium_viewer_standalone.html`

## Key Design Decisions

### Two-stage cell typing: MapMyCells followed by correlation classifier

The initial HANN labels from MapMyCells provide a reasonable starting point, but the bootstrapped hierarchical mapping can misclassify cells — for example, placing L6b cells in upper cortical layers. To address this, step 02b builds Pearson correlation centroids from the top-100 highest-confidence HANN exemplars per cell type, using only Xenium expression data. Because these centroids are constructed from Xenium itself, they capture platform-specific expression characteristics without requiring cross-platform normalization (which we found introduces systematic artifacts — e.g., Harmony misclassified non-neuronal types into GABAergic categories and achieved only 69% agreement with our final classifier). The correlation classifier achieves Pearson r = 0.80 against independent MERFISH proportions, compared to r = 0.73 for Harmony-based approaches. See the [methods writeup](methods_writeup.md) for full benchmarking.

### Depth from spatial neighbors, not expression

Cortical depth is predicted using K=50 spatial neighborhood cell type composition as features, not the 300-gene expression profile directly. At 300 genes, expression-based depth prediction is noisy and overfits to a few marker genes. Neighborhood composition is more robust because it leverages the known laminar organization of cell types — a cell surrounded by L4 IT neurons is almost certainly in layer 4, regardless of its own expression noise. The model (GradientBoostingRegressor trained on SEA-AD MERFISH data) achieves R-squared = 0.90 and MAE = 0.031 on held-out donors.

### BANKSY spatial domains

Tissue domain classification uses BANKSY (lambda=0.8, resolution=0.3), which augments gene expression with spatial neighbor expression to produce spatially coherent clusters. This replaced an earlier K-NN Leiden approach that misclassified L1 border cells as "Extra-cortical" and lacked white matter detection. BANKSY correctly identifies L1 border regions (shallow, non-neuronal-dominated clusters) as Cortical with a `banksy_is_l1` flag, validated against SEA-AD MERFISH showing L1 has ~81% non-neuronal composition.

### Correlation margin QC: bottom 5% per sample

Rather than applying a hard correlation threshold (which disproportionately removes rare cell types with inherently lower classifier confidence), the pipeline flags the bottom 5% of cells by subclass correlation margin *within each sample*. This ensures QC pressure is distributed evenly across cell types and samples. Combined with spatial doublet detection (validated at 0.098% false-positive rate against snRNAseq), this produces the `corr_qc_pass` column used as the default QC gate.

### L6b special treatment

L6b cells receive a stricter absolute margin threshold (0.02) in addition to the per-sample percentile filter. This is an ad hoc but empirically effective fix: without it, too many L6b-labeled cells appeared in upper cortical layers where L6b should not exist. The strict threshold removes these misclassified cells at the cost of slightly reduced L6b counts, which is preferable to systematic spatial mislocalization.

### Flag-based QC, never cell removal

The pipeline never removes cells from the h5ad files. Instead, it adds boolean QC flag columns (`qc_pass`, `corr_qc_pass`, `doublet_suspect`) so downstream analyses can choose their own filtering stringency. All data is preserved and QC decisions are transparent and reversible.

### External validation at every step

Cell type proportions are benchmarked against SEA-AD MERFISH (an independent spatial dataset with 341K cortical cells), depth distributions are validated against manually annotated cortical layers (r=0.92), and doublet detection thresholds are calibrated against snRNAseq false-positive rates (0.098%). No result is accepted on internal consistency alone.

## Cell Type Taxonomy

The pipeline uses the [SEA-AD MTG taxonomy](https://portal.brain-map.org/), organized into three hierarchical levels:

- **3 classes:** Neuronal: Glutamatergic, Neuronal: GABAergic, Non-neuronal and Non-neural
- **24 subclasses:** L2/3 IT, L4 IT, L5 IT, L5 ET, L5/6 NP, L6 IT, L6 IT Car3, L6 CT, L6b, Sst, Sst Chodl, Pvalb, Vip, Lamp5, Lamp5 Lhx6, Sncg, Pax6, Chandelier, Astrocyte, Oligodendrocyte, OPC, Microglia-PVM, Endothelial, VLMC
- **127+ supertypes:** Fine-grained subtypes within each subclass (e.g., Sst_1 through Sst_25)

## Depth Strata

Layers are assigned through a combined approach: BANKSY spatial domain classification (step 05) identifies Cortical, Vascular, and White Matter regions; the MERFISH-trained depth model (step 04) predicts normalized cortical depth for cortical cells; and a 3-step spatial smoothing pipeline refines layer boundaries.

| Layer | Depth Range | Description |
|-------|-------------|-------------|
| L1 | < 0.12 | Molecular layer (identified via BANKSY `banksy_is_l1` + depth) |
| L2/3 | 0.12 -- 0.47 | Superficial layers |
| L4 | 0.47 -- 0.54 | Granular layer (L4 IT median depth ~0.51 in MERFISH) |
| L5 | 0.54 -- 0.71 | Deep pyramidal layer |
| L6 | 0.71 -- 0.93 | Deep layer (includes L6 CT, L6 IT, L6b) |
| WM | > 0.93 | White matter |
| Vascular | -- | Vascular domain (from BANKSY, post-smoothing) |

Final layer proportions: L1 (5.9%), L2/3 (30.7%), L4 (6.3%), L5 (17.3%), L6 (16.8%), WM (16.4%), Vascular (6.6%).

Layer boundaries are derived from pairwise excitatory neuron marker crossovers in SEA-AD MERFISH, validated against Xenium predicted depths (see `code/analysis/derive_layer_boundaries.py`).

See the **[methods writeup](methods_writeup.md)** for full details, validation figures, and design decisions.

## Directory Structure

```
SCZ_Xenium/
├── code/
│   ├── modules/              # Core library modules
│   ├── pipeline/             # Numbered pipeline steps (00-07, run sequentially)
│   ├── analysis/             # Downstream statistical analyses & plotting
│   ├── nuclear_resolution/   # Nuclear doublet resolution (optional side investigation)
│   └── archive/              # Legacy/exploratory scripts (reference only)
├── data/
│   ├── raw/              # 72 raw Xenium files (not in repo; see data/README.md)
│   └── reference/        # SEA-AD reference datasets (not in repo; see data/README.md)
├── output/
│   ├── h5ad/             # Per-sample annotated h5ad files (generated)
│   ├── viewer/           # Interactive spatial viewer
│   ├── crumblr/          # Compositional regression results
│   ├── presentation/     # Presentation-ready figures
│   └── ...               # CSVs, model files (generated)
├── environment.yml      # Python/conda dependencies
└── requirements.txt     # Python dependencies (pip)
```

## h5ad Schema

Each `{sample}_annotated.h5ad` contains:
- `.X`: sparse count matrix (cells x 300 genes)
- `.obsm['spatial']`: (x, y) spatial coordinates

### `.obs` columns

| Column | Type | Description |
|--------|------|-------------|
| `sample_id` | category | Sample identifier |
| **QC (step 01)** | | |
| `qc_pass` | bool | Overall QC pass flag |
| `total_counts` | int64 | Total UMI counts per cell |
| `n_genes` | int64 | Number of genes detected |
| `neg_probe_sum` | int64 | Negative control probe sum |
| `neg_codeword_sum` | int64 | Negative codeword sum |
| `unassigned_sum` | int64 | Unassigned transcript sum |
| `fail_neg_probe` | bool | Failed negative probe QC |
| `fail_neg_codeword` | bool | Failed negative codeword QC |
| `fail_unassigned` | bool | Failed unassigned QC |
| `fail_n_genes_low` | bool | Too few genes detected |
| `fail_total_counts_low` | bool | Too few UMIs |
| `fail_total_counts_high` | bool | Too many UMIs |
| **HANN labels (step 02)** | | |
| `class_label` | category | Broad class (3 types) |
| `subclass_label` | category | Subclass (24 types) |
| `supertype_label` | category | Supertype (127+ types) |
| `class_label_confidence` | float32 | HANN bootstrapping confidence |
| `subclass_label_confidence` | float32 | HANN bootstrapping confidence |
| `supertype_label_confidence` | float32 | HANN bootstrapping confidence |
| **Correlation classifier (step 02b)** | | |
| `corr_subclass` | category | Reclassified subclass |
| `corr_supertype` | category | Reclassified supertype |
| `corr_class` | category | Derived from corr_subclass |
| `corr_subclass_corr` | float32 | Best subclass Pearson correlation |
| `corr_subclass_margin` | float32 | Margin between top-2 subclass correlations |
| `corr_supertype_corr` | float32 | Best supertype Pearson correlation |
| `corr_qc_pass` | bool | Passes both margin QC and doublet QC |
| `doublet_suspect` | bool | Suspected spatial doublet |
| `doublet_type` | str | '' / 'Glut+GABA' / 'GABA+GABA' |
| **Depth & layers (steps 04-05)** | | |
| `predicted_norm_depth` | float64 | Predicted normalized cortical depth (0=pia, 1=WM) |
| `banksy_cluster` | int | BANKSY cluster ID (-1 for QC-failed cells) |
| `banksy_domain` | str | Cortical / Vascular / WM (BANKSY-based domain) |
| `banksy_is_l1` | bool | True if cell in L1 border cluster (shallow, non-neuronal) |
| `spatial_domain` | category | Cortical / Vascular (backward-compatible; WM mapped to Cortical) |
| `layer` | category | L1 / L2/3 / L4 / L5 / L6 / WM / Vascular (final spatially-smoothed layer assignment) |

**Notes:**
- `layer` is the final spatially-smoothed layer assignment and the recommended column for all downstream analyses.
- The correlation classifier columns (`corr_*`, `doublet_*`) are present only in samples processed through step 02b. Analysis scripts prefer `corr_subclass`/`corr_supertype` when available, falling back to HANN labels otherwise.
- The nuclear doublet columns (see Nuclear Doublet Resolution section below) are present only in samples where the optional resolution step was previously run. These columns are retained for reference but `corr_qc_pass` is the default QC gate for all analyses.

## Getting Started

### Data

See **[data/README.md](data/README.md)** for complete download instructions covering:
- Raw Xenium data (72 files across 24 samples)
- SEA-AD MERFISH reference (`SEAAD_MTG_MERFISH.2024-12-11.h5ad`)
- SEA-AD snRNAseq reference (`seaad_mtg_snrnaseq_reference.h5ad`)
- MapMyCells precomputed stats (`precomputed_stats.20231120.sea_ad.MTG.h5`)

### Environment

```bash
conda env create -f environment.yml
pip install pybanksy                    # Required for step 05 (BANKSY spatial domains)
pip install "cell_type_mapper @ git+https://github.com/AllenInstitute/cell_type_mapper"  # Python 3.10+, required for step 02
```

### Running the pipeline

```bash
python3 -u code/pipeline/00_create_h5ad.py
python3 -u code/pipeline/01_run_qc.py
python3 -u code/pipeline/02_run_mapmycells.py
python3 -u code/pipeline/02b_run_correlation_classifier.py
python3 -u code/pipeline/03_export_transcripts.py          # viewer only
python3 -u code/pipeline/04_run_depth_prediction.py
python3 -u code/pipeline/05_run_spatial_domains.py
python3 -u code/pipeline/06_export_viewer.py                # viewer only
python3 -u code/pipeline/07_export_boundaries.py            # viewer only

# Open the interactive viewer
open output/viewer/xenium_viewer_standalone.html
```

**Note:** Steps 04-05 use multiprocessing (4 workers by default, configurable in `pipeline_config.py`). The full pipeline processes ~1.34M cells across 24 samples. Step 04 (depth model training) is the most time-intensive (~80 min). Total pipeline runtime is approximately 2-3 hours.

## Web Viewer & Deployment

The interactive spatial viewer can be run locally or deployed to Netlify.

### Local development

```bash
# Serve from the viewer directory (uses fetch-based data loading)
python3 -m http.server 8080 --directory output/viewer

# Or open the standalone HTML directly (all data embedded, no server needed)
open output/viewer/xenium_viewer_standalone.html
```

### Netlify deployment

The viewer is deployed to Netlify from the `output/deploy/` directory. This directory is gitignored — it is regenerated from `output/viewer/` before each deploy.

**Prerequisites:** Install the [Netlify CLI](https://docs.netlify.com/cli/get-started/) (`npm install -g netlify-cli`) and authenticate (`netlify login`). The Netlify site ID is stored in `output/deploy/.netlify/state.json`.

**Steps to deploy:**

```bash
# 1. Regenerate viewer data from h5ad files (if pipeline outputs changed)
python3 -u code/pipeline/06_export_viewer.py    # exports JSON + standalone HTML
python3 -u code/pipeline/07_export_boundaries.py  # exports cell/nucleus boundary polygons

# 2. Sync viewer files to the deploy directory
rsync -av --delete \
  --exclude='.netlify' \
  --exclude='xenium_viewer_standalone.html' \
  output/viewer/ output/deploy/

# 3. Trim transcript data for deploy size
#    - Remove mitochondrial transcripts (MT-*, MTRNR*) — large files, not useful for viewer
#    - Optionally remove transcript directories for samples you don't need molecule overlay for
#    Full transcripts for all 24 samples is ~50 GB; 5-6 samples is ~3-4 GB

# Remove MT transcripts from all samples
find output/deploy/transcripts/ \( -name "MT-*.json" -o -name "MTRNR*.json" \) -delete

# Update gene_index.json files to exclude MT genes
python3 -c "
import json, glob
for idx_path in glob.glob('output/deploy/transcripts/*/gene_index.json'):
    with open(idx_path) as f:
        d = json.load(f)
    d['genes'] = [g for g in d['genes']
                  if not g['gene'].startswith('MT-') and not g['gene'].startswith('MTRNR')]
    d['n_genes'] = len(d['genes'])
    with open(idx_path, 'w') as f:
        json.dump(d, f)
"

# (Optional) Keep transcripts for only a subset of samples to reduce deploy size
# The viewer gracefully hides the Transcripts panel for samples without data
KEEP_TRANSCRIPTS="Br8667 Br6432 Br2039 Br5746 Br5973"
for dir in output/deploy/transcripts/*/; do
    sample=$(basename "$dir")
    echo "$KEEP_TRANSCRIPTS" | grep -qw "$sample" || rm -rf "$dir"
done

# 4. Verify deploy size (Netlify free tier limit: ~15 GB; aim for <10 GB)
du -sh output/deploy/

# 5. Deploy to Netlify
cd output/deploy && netlify deploy --prod --dir=.
```

**What gets deployed:**
- `index.html` — the viewer app (vanilla HTML/CSS/JS, no build step)
- `index.json` — global metadata (sample list, color palettes, taxonomy)
- `{sample_id}.json` (×24) — per-sample cell data (coordinates, labels, depth, layer, QC)
- `boundaries/{sample_id}.json` + `_nucleus.json` (×24) — cell/nucleus boundary polygons
- `transcripts/{sample_id}/{gene}.json` — per-gene molecule coordinates (subset of samples)

**Size breakdown:** ~98 MB (cell JSONs + HTML) + ~740 MB (boundaries) + ~600 MB per sample with transcripts. Without any transcripts, the deploy is under 1 GB.

**Notes:**
- All 24 samples are always viewable (cell dots, all color modes, layer overlay, etc.) regardless of which samples have transcript data
- The standalone HTML (`xenium_viewer_standalone.html`, ~26 MB) is excluded from deploy — it's for offline/local use only
- If layer boundaries change, re-run `update_layer_boundaries.py` → `06_export_viewer.py` → rsync → deploy

## Nuclear Doublet Resolution (optional)

An optional nuclear doublet resolution module lives in `code/nuclear_resolution/`. It uses nuclear-only count matrices to arbitrate doublet calls, producing a `hybrid_qc_pass` column. Empirical testing showed this has negligible impact on downstream compositional analyses (see `docs/pipeline_qc_audit.md`), so it is not part of the standard pipeline. The `corr_qc_pass` column from step 02b is the default QC gate.

See **[code/nuclear_resolution/README.md](code/nuclear_resolution/README.md)** for details on running this step and interpreting results.

### Nuclear doublet columns (present only if the optional step was run)

| Column | Type | Description |
|--------|------|-------------|
| `nuclear_total_counts` | int | Total UMI counts from nuclear-only transcripts |
| `nuclear_n_genes` | int | Number of genes detected in nuclear transcripts |
| `nuclear_fraction` | float32 | Fraction of UMIs in nucleus vs whole cell |
| `nuclear_doublet_suspect` | bool | Nuclear-level doublet detected |
| `nuclear_doublet_type` | str | '' / 'Glut+GABA' / 'GABA+GABA' (nuclear evidence) |
| `nuclear_doublet_status` | str | 'clean' / 'resolved' / 'persistent' / 'nuclear_only' / 'insufficient' |
| `hybrid_qc_pass` | bool | Hybrid QC pass (nuclear-informed; not used by default) |

## Archive (`code/archive/`)

Contains legacy and exploratory scripts from earlier iterations of the analysis:

| Directory | Contents |
|-----------|----------|
| `stale_analysis/` | Archived analysis scripts: diagnostic comparisons, edgepython DE, Harmony transfer, nsforest markers, calibration scripts |
| `one_time_utils/` | One-time data migration scripts (rename columns, annotate MERFISH depth) |
| `legacy_runners/` | Old pipeline runners (pre-numbered-step architecture) |
| `ood_methods/` | OOD method exploration scripts (superseded by BANKSY domain classification) |
| `spatial_domain_exploration/` | Early spatial domain clustering experiments |
| `banksy_exploration/` | BANKSY parameter tuning (pilot), early domain+strip prototype, validation, batch runner. Logic extracted to `modules/banksy_domains.py` |
| `curved_strips/` | Curved cortex strip identification: pia curve fitting, strip boundary detection, MERFISH adapter, gallery visualization. Experimental approach for selecting cortical strips with complete L1-L6 laminar structure |
| *(root)* | Archived modules: `label_transfer.py`, `layers.py`, old pipeline scripts |

These are kept for reference but are not part of the active pipeline.
