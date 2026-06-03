# Combined SEA-AD + Siletti Taxonomy: Methodology

This document describes the approach used to build a combined cell-type taxonomy
from the SEA-AD snRNA-seq reference (137 supertypes) and the Siletti et al. (2023)
whole-brain atlas (461 clusters), for the purpose of GWAS cell-type enrichment
analysis of schizophrenia.

> **Looking for the canonical type list?** See
> [`franken_taxonomy/`](../franken_taxonomy/) at the repo root —
> `franken_types.csv` (598 candidate types with provenance) and
> `franken_rbh.csv` (the 95 reciprocal-best-hit pairs) are the artifacts
> that downstream projects should consume.

## Overview

The goal is to create a single taxonomy that:
1. Retains SEA-AD's fine-grained neocortical resolution (25 Sst subtypes, 13 Pvalb subtypes, etc.)
2. Adds cell types from Siletti that are genuinely novel (not already captured by SEA-AD)
3. Avoids double-counting: if a Siletti cluster is transcriptomically equivalent to a SEA-AD type, only the SEA-AD type is kept

The approach uses **MetaNeighbor reciprocal best hits** to determine which Siletti clusters
are redundant with SEA-AD types vs. which represent unique cell populations.

## Part 1: Siletti Data Processing

### Where the Siletti reference data lives

The Allen Brain Cell Atlas (ABCA) hosts the WHB-10Xv3 atlas on S3. Two
artifacts are used downstream:

- **MapMyCells precomputed stats** — `precomputed_stats.siletti.training.h5`
  (8.1 GB), used by steps 08–10 here and by the `RSC_Xenium` cell-typing
  pipeline. [Download from Allen S3](https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/mapmycells/WHB-10Xv3/20240831/precomputed_stats.siletti.training.h5)
  ([browse dir](https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html#mapmycells/WHB-10Xv3/20240831/)).
- **WHB-10Xv3 metadata** — `cell_metadata.csv`, `taxonomy/*.csv`, `gene.csv`.
  Downloadable via the [`abc_atlas_access`](https://github.com/AllenInstitute/abc_atlas_access)
  Python package.

Citation: Siletti et al. 2023, *Science*
([doi:10.1126/science.add7046](https://doi.org/10.1126/science.add7046)).
Taxonomy ID: `CCN202210140` (supercluster / cluster / subcluster).

### Data versions

Two versions of the Siletti specificity matrix are used:

| Version | Source | Genes | Clusters | Gene IDs |
|---------|--------|-------|----------|----------|
| **Conti/Duncan** | Original 2022 loom, precomputed for Duncan et al. 2025 | 17,097 | 461 | ENTREZ |
| **Reprocessed** | Updated Allen Institute release, computed from `siletti_cluster_level_stats.npz` | 17,670 | 461 | ENTREZ |

### Specificity computation (reprocessed)

Starting from the cluster-level mean expression (`siletti_cluster_level_stats.npz`,
461 clusters × 59,357 Ensembl genes):

1. **Gene ID mapping**: Ensembl → Symbol (from loom file `Gene` attribute) → ENTREZ
   (from NCBI gene loc file). This chains two lookups:
   - Loom `Accession` (versioned Ensembl) → `Gene` (symbol): 59,436 mappings
   - NCBI `SYMBOL` → `ENTREZ`: 19,175 mappings
   - Combined: 17,801 Ensembl → ENTREZ mappings

2. **Deduplication**: When multiple Ensembl IDs map to the same ENTREZ ID, keep
   the one with the highest total expression across clusters.

3. **Normalization**: Same as SEA-AD pipeline —
   - Column-normalize each cluster to sum to 1,000 (equalize library sizes)
   - Row-normalize each gene to sum to 1.0 (convert to specificity fractions)

4. **Output**: 17,670 ENTREZ genes × 461 clusters, rows sum to 1.0

### Why ENTREZ IDs?

The Conti/Duncan specificity matrix uses ENTREZ gene IDs, and MAGMA's gene-level
results also use ENTREZ IDs. By keeping everything in ENTREZ space, we avoid the
gene symbol mapping step that loses ~3% of genes. This gives a fairer comparison
between Conti and reprocessed, and slightly more statistical power.

### Conti vs. reprocessed validation

Running identical enrichment pipelines on both versions (all 461 clusters, ENTREZ IDs):

| Metric | Value |
|--------|-------|
| Beta correlation (Spearman) | 0.904 |
| Significance correlation (Spearman) | 0.923 |
| Both FDR-significant | 101 clusters |
| Conti-only significant | 24 |
| Reprocessed-only significant | 14 |

The two versions are highly concordant. The reprocessed data is used for the combined
taxonomy since it reflects the most current Allen Institute release.

## Part 2: MetaNeighbor Integration

### Purpose

MetaNeighbor (Crow et al. 2018) identifies which cell types from two different
datasets represent the same biological entity. We use it to determine which of the
461 Siletti clusters are transcriptomically equivalent to one of the 137 SEA-AD types.

### Method (adapted for mean expression profiles)

Since we work with mean expression per cell type (not single cells), we implement
the core MetaNeighbor algorithm on aggregated profiles:

#### Step 1: Build combined mean expression

Combine SEA-AD (137 types, log-normalized mean expression) with Siletti (461 clusters,
similarly normalized from raw mean expression). Intersect to shared gene symbols.

Result: 35,321 shared genes × 598 cell types.

#### Step 2: Select highly variable genes (HVGs)

Compute coefficient of variation (CV = std/mean) across all 598 types for each gene.
Filter to expressed genes (mean > 0.1), then take the top 5,000 by CV.

5,000 HVGs were chosen (rather than the default 3,000) because the two datasets span
very different cell types (neocortical + whole-brain), and more genes help distinguish
subtypes that differ in subtle expression programs.

#### Step 3: Compute cross-dataset Spearman correlation

Using only HVGs, rank-transform each cell type's expression profile, then compute
Pearson correlation on ranks (= Spearman) between all 137 SEA-AD types and all 461
Siletti clusters.

Result: 137 × 461 correlation matrix, range [-0.06, 0.93].

#### Step 4: Compute neighbor-voting AUROC

For each SEA-AD type *i* and Siletti cluster *j*, the AUROC measures:

> Among all 461 Siletti clusters, what fraction have *lower* correlation with
> SEA-AD type *i* than cluster *j* does?

AUROC = 1.0 means cluster *j* is the most similar to type *i* among all Siletti clusters.
AUROC = 0.5 means random (no similarity). Ties receive 0.5 credit.

Result: 137 × 461 AUROC matrix.

#### Step 5: Identify reciprocal best hits

A reciprocal best hit (RBH) occurs when:
- SEA-AD type A's best Siletti match (highest AUROC) is cluster B
- Siletti cluster B's best SEA-AD match (highest AUROC in the reverse direction) is type A

This bidirectional criterion is stringent: it ensures that both datasets "agree" that
these two types represent the same cell population.

### Results

| Metric | Value |
|--------|-------|
| Reciprocal best hits | 95 / 137 SEA-AD types |
| All with mean AUROC > 0.9 | 95 (100%) |
| Non-reciprocal SEA-AD types | 42 |
| Siletti clusters claimed by ≥1 SEA-AD type | 95 |
| Novel Siletti clusters (not claimed) | 366 |

**Why are 42 SEA-AD types non-reciprocal?**

These are SEA-AD subtypes that map to the same Siletti cluster as another SEA-AD type.
For example, multiple Pvalb subtypes (Pvalb_3, Pvalb_5) both map best to MGE_259, but
MGE_259's best SEA-AD match is only one of them (Pvalb_14). This reflects SEA-AD's finer
resolution — it distinguishes subtypes that Siletti groups into one cluster.

## Part 3: RBH Combined Taxonomy

### Construction logic

```
For each SEA-AD type (137):
    Keep it (always included)

For each Siletti cluster (461):
    If it is a reciprocal best hit for any SEA-AD type:
        Exclude it (already represented by the SEA-AD type)
    Else:
        Include it as a novel type
```

Result: **137 SEA-AD types + 366 novel Siletti clusters = 503 total types**

### Specificity computation

The combined specificity matrix is computed from scratch using the full mean expression:
1. Select the 503 columns (137 SEA-AD + 366 novel Siletti) from the combined mean expression
2. Column-normalize each type to sum to 1,000
3. Row-normalize each gene to sum to 1.0
4. Intersect with MAGMA gene-level results for enrichment testing

### Enrichment results

| Metric | RBH Combined (503) | Old 20%-Cortical (485) |
|--------|-------------------|----------------------|
| Total types | 503 | 485 |
| SEA-AD types | 137 | 137 |
| Siletti clusters | 366 | 348 |
| FDR-significant | 120 | 108 |
| Bonferroni-significant | 54 | ~50 |

The RBH taxonomy finds 12 more FDR-significant types, mostly non-neocortical Siletti
clusters (hippocampal, amygdalar, thalamic) that the old 20%-cortical filter excluded.

### Comparison with old approach

The old approach kept Siletti clusters where ≥20% of subclusters were neocortical.
This had two problems:
1. **Double-counting**: Siletti clusters matching SEA-AD types were kept alongside them
2. **Missing types**: Non-neocortical clusters with genuine enrichment were excluded

The RBH approach solves both: it uses transcriptomic evidence (MetaNeighbor) for
deduplication, and includes all non-redundant Siletti clusters regardless of anatomical
location.

For SEA-AD types present in both taxonomies, enrichment betas are virtually identical
(Spearman r = 0.992).

## Part 4: Data Flow

```
INPUT DATA
├── data/seaad_reference.h5ad                    → SEA-AD mean expression + specificity
├── data/adult_human_20221007.loom               → Ensembl→Symbol gene mapping
├── data/conti_specificity_matrix.txt            → Conti/Duncan Siletti specificity (ENTREZ)
├── results/intermediates/siletti_cluster_level_stats.npz → Reprocessed Siletti cluster means
└── linking_.../PGC3_SCZ...genes.out             → MAGMA gene-level GWAS results

STEP 08: SILETTI ENRICHMENT
├── Compute reprocessed specificity (Ensembl→ENTREZ, 461 clusters)
│   → results/intermediates/reprocessed_siletti_specificity_entrez.csv
├── Run enrichment on Conti specificity
│   → results/tables/conti_siletti_enrichment.csv
├── Run enrichment on reprocessed specificity
│   → results/tables/reprocessed_siletti_461_enrichment.csv
└── Compare versions
    → results/tables/conti_vs_reprocessed_461_comparison.csv

STEP 09: METANEIGHBOR INTEGRATION
├── Build full mean expression (137 SEA-AD + 461 Siletti, gene symbols)
│   → results/intermediates/full_seaad_siletti_461_mean_expression.csv
├── Select 5000 HVGs, compute Spearman correlation, compute AUROC
│   → results/intermediates/metaneighbor_full461_auroc_matrix.csv
├── Find reciprocal best hits (95 RBH)
│   → results/tables/metaneighbor_full461_reciprocal_best_hits.csv
│   → results/tables/metaneighbor_full461_all_matches.csv
└── Compute match lookups for app
    → results/intermediates/rbh_seaad_to_siletti_metaneighbor.csv
    → results/intermediates/rbh_siletti_to_seaad_metaneighbor.csv

STEP 10: COMBINED TAXONOMY
├── Build RBH taxonomy (137 SEA-AD + 366 novel Siletti = 503)
├── Compute combined specificity from mean expression
├── Run enrichment
│   → results/tables/rbh_combined_enrichment.csv
│   → results/tables/rbh_combined_enrichment_app.csv
├── Compare with old taxonomy
│   → results/tables/rbh_vs_old_combined_comparison.csv
├── Generate figures
│   → results/figures/rbh_combined_taxonomy_manhattan.png
│   → results/figures/rbh_combined_taxonomy_top_enrichments.png
│   → results/figures/duncan_manhattan_replication.png
└── Prepare webapp data
    → results/intermediates/rbh_combined_type_order_info.csv

STEP 11: GENE DRIVER SCATTER PLOTS
├── For each cell type: specificity (x) vs GWAS -log10(p) (y)
├── Driver genes = top 10% specificity AND GWAS FDR < 0.05
│   → results/figures/gene_drivers/gene_driver_scatter_{TYPE}.png
│   → results/tables/gene_drivers/gene_drivers_{TYPE}.csv
└── Summary across types
    → results/tables/gene_drivers/driver_gene_summary.csv

WEBAPP (app.py)
├── Main view: Manhattan plot of RBH combined enrichment (503 types)
│   Loads rbh_combined_enrichment_app.csv + MetaNeighbor matches + cluster metadata
└── Gene driver view: interactive scatter per cell type (/drivers/{cell_type})
    Computed on-demand from specificity + MAGMA data
```

## Part 5: Gene Driver Scatter Plots

For each cell type, we identify **driver genes** — genes that are both highly specific
to the cell type AND carry strong SCZ GWAS signal. Following the approach recommended
by Laramie Duncan and colleagues:

- **X-axis**: Cell-type specificity score (fraction of total expression in this type)
- **Y-axis**: -log₁₀(GWAS gene-level p-value) from MAGMA

A gene is classified as a "driver" if it is in the **top 10% of specificity** for that
cell type AND passes **GWAS FDR < 0.05** (3,128 genes out of ~17,000 tested).

Three threshold lines are shown:
- **Bonferroni** (purple): p < 0.05/n_genes ≈ 2.9e-6 (-log₁₀ ≈ 5.5)
- **FDR 0.05** (red): the threshold used to define driver genes
- **Nominal** (orange): p < 0.05

Driver gene counts per cell type range from ~250-325, reflecting both the type's
specificity profile and the GWAS signal distribution.

The webapp provides interactive versions at `/drivers/{cell_type}` with hover tooltips,
a filterable/sortable gene table, and summary statistics.

## Part 6: Case-Control Composition Correlation

Correlates GWAS cell-type enrichment from the RBH combined taxonomy with
case-control compositional changes from a 7-cohort snRNA-seq meta-analysis
of SCZ (Endresz et al., in prep; crumblr model). Now implemented as `scripts/13_gwas_vs_composition.py`.

**Key finding:** Across all 109 matched cell types, there is no significant correlation
between GWAS enrichment and composition change (r = -0.10, p = 0.28). However,
**within SST interneurons specifically** (n = 18), there is a significant positive correlation:
- |Composition change| vs GWAS significance: **r = 0.645, p = 3.85e-3**
- Significance vs significance: **r = 0.701, p = 1.20e-3**

SST subtypes with the strongest GWAS enrichment (Sst_2, Sst_25, Sst_3) are also the ones most
decreased in proportion in SCZ brains — convergence of genetic predisposition and
observed cellular pathology on the same interneuron subtypes.

See `docs/rnaseq_specificity_enrichment_results.md` for detailed tables and interpretation.

**Outputs:**
- `results/tables/gwas_vs_casecontrol_composition.csv`
- `results/figures/gwas_vs_casecontrol_composition.png`

## Key Files

### Scripts
| Script | Description |
|--------|-------------|
| `scripts/08_siletti_enrichment.py` | Siletti enrichment (Conti + reprocessed, 461 clusters) |
| `scripts/09_metaneighbor_integration.py` | MetaNeighbor cross-dataset matching (full 461) |
| `scripts/10_combined_taxonomy.py` | RBH combined taxonomy + enrichment + webapp data |
| `scripts/11_gene_driver_scatter.py` | Gene driver scatter plots (static PNG + CSV) |
| `app.py` | Interactive webapp (enrichment explorer + gene driver views) |

### Package modules
| Module | Description |
|--------|-------------|
| `scz_celltype_enrichment/siletti/specificity.py` | Siletti specificity computation (Ensembl→ENTREZ) |
| `scz_celltype_enrichment/siletti/clusters.py` | Cluster ID mapping and metadata |
| `scz_celltype_enrichment/siletti/metaneighbor.py` | HVG selection, AUROC, reciprocal best hits |

### Key outputs
| File | Description |
|------|-------------|
| `results/tables/rbh_combined_enrichment.csv` | Final combined taxonomy enrichment (503 types) |
| `results/tables/metaneighbor_full461_reciprocal_best_hits.csv` | 95 SEA-AD ↔ Siletti matches |
| `results/tables/metaneighbor_full461_all_matches.csv` | All 137 SEA-AD with top 3 Siletti matches |
| `results/tables/conti_vs_reprocessed_461_comparison.csv` | Siletti data version comparison |
| `results/intermediates/metaneighbor_full461_auroc_matrix.csv` | Full AUROC matrix (137 × 461) |
| `results/intermediates/reprocessed_siletti_specificity_entrez.csv` | Reprocessed specificity (ENTREZ) |
