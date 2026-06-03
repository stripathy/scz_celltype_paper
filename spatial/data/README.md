# Data Directory

This directory contains input data for the SCZ Xenium spatial transcriptomics pipeline.
Large files are excluded from git — see download instructions below.

## Directory Structure

```
data/
├── raw/                          # Raw Xenium output files (from GEO)
├── reference/                    # Allen Brain Cell Atlas reference datasets
├── nicole_scz_snrnaseq_betas/    # snRNAseq SCZ differential expression betas (Kwon et al.)
├── merscope_4k_probe_testing/    # MERSCOPE 4K probe benchmark data
├── sample_metadata.xlsx          # Subject metadata (age, sex, diagnosis, PMI)
└── xenium_probes.xlsx            # Xenium probe panel info
```

## What You Need

Not all data is required for every use case. Here's a quick guide:

| Dataset | Size | Core pipeline | Analysis scripts | Validation plots |
|---------|------|:---:|:---:|:---:|
| Raw Xenium cell matrices + boundaries | 809 MB | **Required** | — | — |
| Raw Xenium transcript coordinates (.zarr.zip) | ~33 GB | Optional (step 03 only) | — | — |
| SEA-AD MERFISH reference | 3.1 GB | **Required** (step 04) | — | Used if available |
| MapMyCells precomputed stats | 251 MB | **Required** (step 02) | — | — |
| SEA-AD snRNAseq reference | 33.8 GB | Not needed | — | Used if available |
| Gene symbol mappings | <1 MB | **Required** (step 02) | — | — |

**Minimum for core pipeline (steps 00-02b, 04-05):** ~4.1 GB (cell matrices + MERFISH + MapMyCells stats + gene mappings)

**With transcript viewer (+ steps 03, 06-07):** ~37 GB (adds zarr files)

**With all validation plots:** ~71 GB (adds snRNAseq reference)

The pipeline and analysis scripts check for reference file availability at runtime. Scripts that use snRNAseq or MERFISH references will raise a clear error with download instructions if the file is missing, so you can run whatever you have data for.

---

## Step 1: Raw Xenium Data (GEO: GSE307404)

**Source:** [Kwon et al. (2026)](https://doi.org/10.64898/2026.02.16.706214) — *Mapping spatially organized molecular and genetic signatures of schizophrenia across multiple scales in human prefrontal cortex*

**GEO Accession:** [GSE307404](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307404)

**Note:** The GEO data may still be under embargo while the manuscript is in review. Check the GEO page for availability.

Each sample has 4 files:
- `cell_feature_matrix.h5` — Gene expression counts (541 features: 300 genes + controls) **[Required]**
- `cell_boundaries.csv.gz` — Cell boundary polygons **[Required]**
- `nucleus_boundaries.csv.gz` — Nucleus boundary polygons **[Required for step 07]**
- `transcripts.zarr.zip` — Molecule-level transcript coordinates **[Optional — only needed for interactive viewer (step 03) and nuclear doublet resolution (see `code/nuclear_resolution/`)]**

```bash
mkdir -p data/raw
cd data/raw

# Download cell matrices and boundaries only (~809 MB):
wget -r -np -nd -A "*.h5,*.csv.gz" \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE307nnn/GSE307404/suppl/

# Also download transcript coordinates if you want the interactive viewer (~33 GB):
wget -r -np -nd -A "*.zarr.zip" \
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

## Step 4: SEA-AD snRNAseq Reference (Optional)

**Source:** [Gabitto et al. (2024)](https://doi.org/10.1038/s41593-024-01774-5) — SEA-AD MTG single-nucleus RNA-seq dataset.

**Used for:** Validation plots only — ground-truth cell type proportions and cross-modal concordance with Xenium. **Not required for the core pipeline.** Analysis scripts check for this file at runtime and will skip or raise a clear error if it's missing.

```bash
# Download the full SEA-AD MTG snRNAseq dataset (~33.8 GB)
wget -O data/reference/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad \
  "https://sea-ad-single-cell-profiling.s3.us-west-2.amazonaws.com/MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"

# Or using AWS CLI (faster, supports resume):
aws s3 cp \
  s3://sea-ad-single-cell-profiling/MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad \
  data/reference/ --no-sign-request
```

**Subsetting:** The pipeline uses only the 5 neurotypical reference donors: `H18.30.001`, `H18.30.002`, `H19.30.001`, `H19.30.002`, `H200.1023` (137,303 cells x 36,601 genes). Subsetting is performed in code, producing `data/reference/seaad_mtg_snrnaseq_reference.h5ad`.

**Provenance note:** This is the full SEA-AD snRNAseq dataset subset to 5 neurotypical reference donors. The Allen Institute provides a version of this dataset, but it lacks the complete set of SEA-AD supertypes. This subset retains all supertype annotations needed for proportion validation.

---

## Other Data Files

### Sample Metadata

`sample_metadata.xlsx` is included in the repository. Contains donor demographics (age, sex, PMI, RIN) and diagnosis (SCZ vs. Control) for all 24 samples.

### snRNAseq SCZ Betas

`nicole_scz_snrnaseq_betas/scz_coefs.xlsx` contains differential expression coefficients from Kwon et al. snRNAseq meta-analysis, used for cross-modal concordance analysis with Xenium crumblr results.

### MERSCOPE Data

`merscope_4k_probe_testing/` contains MERSCOPE spatial transcriptomics data for cross-platform benchmarking.

---

## Summary

| Dataset | Size | Required for | Source |
|---------|------|-------------|--------|
| Xenium cell matrices + boundaries | 809 MB | Core pipeline (steps 00-05) | [GEO GSE307404](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307404) |
| Xenium transcript coordinates | ~33 GB | Viewer (step 03) | [GEO GSE307404](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307404) |
| SEA-AD MERFISH | 3.1 GB | Depth model (step 04) | [Allen Brain Cell Atlas](https://sea-ad-spatial-transcriptomics.s3.us-west-2.amazonaws.com/middle-temporal-gyrus/all_donors-h5ad/SEAAD_MTG_MERFISH.2024-12-11.h5ad) |
| MapMyCells stats | 251 MB | Cell type annotation (step 02) | [Allen Brain Cell Atlas](https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/mapmycells/SEAAD/20240831/precomputed_stats.20231120.sea_ad.MTG.h5) |
| SEA-AD snRNAseq | 33.8 GB | Validation plots only (optional) | [Allen Brain Cell Atlas](https://sea-ad-single-cell-profiling.s3.us-west-2.amazonaws.com/MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad) |
| **Minimum required** | **~4.1 GB** | | |
| **Full pipeline + viewer** | **~37 GB** | | |
| **Everything** | **~71 GB** | | |
