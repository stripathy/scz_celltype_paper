# Data Directory

This directory contains input data for the SCZ Xenium spatial transcriptomics pipeline.
Large files are excluded from git — see download instructions below.

## Directory Structure

```
data/
├── raw/              # Raw Xenium output files (from GEO)
├── reference/        # Allen Brain Cell Atlas reference datasets
└── sample_metadata.xlsx   # Subject metadata (age, sex, diagnosis, PMI) — tracked
```

## What You Need

Not all data is required for every use case. Here's a quick guide:

| Dataset | Size | Core pipeline | Supp Figs S2 / S3 |
|---------|------|:---:|:---:|
| Raw Xenium cell matrices + boundaries | 809 MB | **Required** | **Required** |
| SEA-AD MERFISH reference | 3.1 GB | **Required** (step 04) | **Required** (S3) |
| MapMyCells precomputed stats | 251 MB | **Required** (step 02) | — |
| SEA-AD snRNAseq reference | 33.8 GB | — | **Required** (S2 c, e) |
| Gene symbol mappings | <1 MB | **Required** (step 02) | — |

**Minimum for the core pipeline (steps 00–02b, 04–05):** ~4.1 GB — cell
matrices, MERFISH, MapMyCells stats and the gene mappings.

**Everything, including both supplementary figures:** ~38 GB.

Neither supplementary figure needs any of this to *render*: their input CSVs are
committed, so `plot_markers_resolvability_combined.R` and
`plot_xenium_merfish_composite.R` run from a clean clone. The downloads below are
only needed to rebuild those inputs from the raw data.

The pipeline and analysis scripts check for reference file availability at runtime. Scripts that use snRNAseq or MERFISH references will raise a clear error with download instructions if the file is missing, so you can run whatever you have data for.

---

## Step 1: Raw Xenium Data (GEO: GSE307404)

**Source:** [Kwon et al. (2026)](https://doi.org/10.64898/2026.02.16.706214) — *Mapping spatially organized molecular and genetic signatures of schizophrenia across multiple scales in human prefrontal cortex*

**GEO Accession:** [GSE307404](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307404)

**Note:** The GEO data may still be under embargo while the manuscript is in review. Check the GEO page for availability.

Each sample has 4 files:
- `cell_feature_matrix.h5` — Gene expression counts (541 features: 300 genes + controls) **[Required]**
- `cell_boundaries.csv.gz` — Cell boundary polygons **[Required]**
- `nucleus_boundaries.csv.gz` — Nucleus boundary polygons **[Optional]**
- `transcripts.zarr.zip` — Molecule-level transcript coordinates **[Not needed]** — these feed the interactive cell browser, which lives in the upstream `SCZ_Xenium` repo

```bash
mkdir -p data/raw
cd data/raw

# Cell matrices and boundaries (~809 MB) — everything this repo needs:
wget -r -np -nd -A "*.h5,*.csv.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE307nnn/GSE307404/suppl/
```

The pipeline discovers samples by globbing `data/raw/*-cell_feature_matrix.h5`.

---

## Step 2: SEA-AD MERFISH Reference (Required)

**Source:** [Gabitto et al. (2024)](https://doi.org/10.1038/s41593-024-01774-5) — Seattle Alzheimer's Disease Brain Cell Atlas.

**Used for:** Training the cortical depth model (step 04). Also used by validation/comparison plots if available.

**Note on depth bins:** The pipeline uses depth bins derived from SEA-AD MERFISH manual annotations: L1 <0.10, L2/3 0.10-0.40, L4 0.40-0.55, L5 0.55-0.70, L6 0.70-0.90, WM >0.90.

```bash
mkdir -p data/reference

# Download MERFISH reference (~3.1 GB)
wget -O data/reference/SEAAD_MTG_MERFISH.2024-12-11.h5ad \
  "https://sea-ad-spatial-transcriptomics.s3.us-west-2.amazonaws.com/middle-temporal-gyrus/all_donors-h5ad/SEAAD_MTG_MERFISH.2024-12-11.h5ad"

# Or using AWS CLI (faster, supports resume):
aws s3 cp \
  s3://sea-ad-spatial-transcriptomics/middle-temporal-gyrus/all_donors-h5ad/SEAAD_MTG_MERFISH.2024-12-11.h5ad \
  data/reference/ --no-sign-request
```

---

## Step 3: MapMyCells Precomputed Stats (Required)

**Source:** [Allen Brain Cell Atlas](https://portal.brain-map.org/) — Precomputed statistics for hierarchical cell type mapping against the SEA-AD MTG taxonomy.

**Used for:** Cell type annotation via MapMyCells (step 02).

```bash
# Download precomputed stats (~251 MB)
wget -O data/reference/precomputed_stats.20231120.sea_ad.MTG.h5 \
  "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/mapmycells/SEAAD/20240831/precomputed_stats.20231120.sea_ad.MTG.h5"
```

**Note:** Step 02 also requires the `cell_type_mapper` Python package (Python 3.10+):

```bash
pip install "cell_type_mapper @ git+https://github.com/AllenInstitute/cell_type_mapper"
```

---

## Step 4: SEA-AD snRNAseq Reference

**Source:** [Gabitto et al. (2024)](https://doi.org/10.1038/s41593-024-01774-5) — SEA-AD MTG single-nucleus RNA-seq dataset.

**Used for:** the panel-vs-transcriptome resolvability benchmark behind Supplementary Fig. S2 panels c and e — leave-one-donor-out classification on the full transcriptome versus the 300-gene panel. **Not required for the core pipeline.** Scripts check for this file at runtime and raise a clear error if it is missing.

```bash
# Download the full SEA-AD MTG snRNAseq dataset (~33.8 GB)
wget -O data/reference/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad \
  "https://sea-ad-single-cell-profiling.s3.us-west-2.amazonaws.com/MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"

# Or using AWS CLI (faster, supports resume):
aws s3 cp \
  s3://sea-ad-single-cell-profiling/MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad \
  data/reference/ --no-sign-request
```

**Subsetting:** only the 5 neurotypical reference donors are used: `H18.30.001`, `H18.30.002`, `H19.30.001`, `H19.30.002`, `H200.1023` (137,303 cells x 36,601 genes). `code/pipeline/create_snrnaseq_reference.py` does the subsetting, writing `data/reference/seaad_mtg_snrnaseq_reference.h5ad`.

**Provenance note:** This is the full SEA-AD snRNAseq dataset subset to 5 neurotypical reference donors. The Allen Institute provides a version of this dataset, but it lacks the complete set of SEA-AD supertypes. This subset retains all supertype annotations needed for proportion validation.

---

## Other Data Files

### Sample Metadata

`sample_metadata.xlsx` is included in the repository. Contains donor demographics (age, sex, PMI, RIN) and diagnosis (SCZ vs. Control) for all 24 samples.

---

## Summary

| Dataset | Size | Required for | Source |
|---------|------|-------------|--------|
| Xenium cell matrices + boundaries | 809 MB | Core pipeline (steps 00-05) | [GEO GSE307404](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307404) |
| SEA-AD MERFISH | 3.1 GB | Depth model (step 04) | [Allen Brain Cell Atlas](https://sea-ad-spatial-transcriptomics.s3.us-west-2.amazonaws.com/middle-temporal-gyrus/all_donors-h5ad/SEAAD_MTG_MERFISH.2024-12-11.h5ad) |
| MapMyCells stats | 251 MB | Cell type annotation (step 02) | [Allen Brain Cell Atlas](https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/mapmycells/SEAAD/20240831/precomputed_stats.20231120.sea_ad.MTG.h5) |
| SEA-AD snRNAseq | 33.8 GB | Resolvability benchmark, Supp Fig S2 c/e | [Allen Brain Cell Atlas](https://sea-ad-single-cell-profiling.s3.us-west-2.amazonaws.com/MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad) |
| **Minimum (core pipeline)** | **~4.1 GB** | | |
| **Everything** | **~38 GB** | | |
