# SCZ Cell-Type Enrichment Analysis: Approach & Methods

## Pipeline Overview

The full pipeline consists of 13 scripts in `scripts/`, run via `run_all.py`:

| Script | Description |
|--------|-------------|
| **01** `compute_specificity.py` | Compute cell-type specificity from SEA-AD snRNA-seq |
| **02** `magma_enrichment.py` | MAGMA-style gene property analysis (SCZ GWAS × specificity) |
| **03** `conditional_analysis.py` | Forward selection for independently enriched types |
| **04** `gene_drivers.py` | Identify genes driving enrichment per cell type |
| **05** `spatial_layer_analysis.py` | SST layer stratification via spatial transcriptomics |
| **06** `gene_ephys_correlations.py` | Gene-electrophysiology correlations (patch-seq) |
| **07** `make_all_figures.py` | Regenerate all publication figures |
| **08** `siletti_enrichment.py` | Siletti atlas enrichment (Conti vs reprocessed) |
| **09** `metaneighbor_integration.py` | MetaNeighbor cross-dataset matching (SEA-AD ↔ Siletti) |
| **10** `combined_taxonomy.py` | Build RBH combined taxonomy (503 types) |
| **11** `gene_driver_scatter.py` | Gene driver scatter plots (specificity × GWAS p-value) |
| **12** `atac_peak_overlap.py` | ATAC-seq peak overlap with FINEMAP credible sets |
| **13** `gwas_vs_composition.py` | GWAS enrichment vs case-control composition |

Scripts 01–07 form the **core SEA-AD pipeline** (`run_all.py`).
Scripts 08–13 add the **Siletti integration, ATAC, and composition analyses** (`run_all.py --all`).
See `docs/combined_taxonomy_approach.md` for methodology on scripts 08–11.
See `docs/rnaseq_specificity_enrichment_results.md` for RNAseq enrichment results and composition intersection.

## Overview

This analysis links schizophrenia (SCZ) GWAS genetic risk to specific brain cell types — and ultimately to specific genes and cortical layers — by integrating six data modalities:

1. **GWAS summary statistics** (PGC3 SCZ, 53K cases / 77K controls)
2. **snRNA-seq reference** (SEA-AD, 137K cells across 137 supertypes)
3. **snATAC-seq chromatin accessibility** (SEA-AD, 516K nuclei × 219K peaks)
4. **Spatial transcriptomics** (MERFISH + Xenium cortical depth estimates)
5. **Patch-seq electrophysiology** (SST interneuron functional properties)
6. **Tissue-specific eQTL models** (GTEx v8 MASHR, brain cortex)

The central question: **Which cell types carry the most SCZ genetic risk signal, and within those types, which genes drive the signal, where are they located in cortex, and do they influence risk or protection?**

---

## Part 1: Cell-Type Specificity from snRNA-seq

### Goal
Quantify how specifically each gene is expressed in each of 137 brain cell types (supertypes).

### Data
- SEA-AD snRNA-seq reference: 137,303 cells, 36,601 genes, 137 supertypes (Allen Institute)

### Method
1. **Normalize**: scanpy `normalize_total(target_sum=1e4)` then `log1p`
2. **Mean expression per supertype**: Average log-normalized expression across all cells of each type
3. **Column-normalize**: Equalize total expression across supertypes (prevents types with larger libraries from dominating)
4. **Row-normalize to specificity**: For each gene, divide by its total expression across all types. Result: a fraction between 0 and 1, where 0.5 means half the gene's total expression is in that one type

**Output**: A 36,601 genes x 137 supertypes specificity matrix. Each row sums to 1.

We also retain the unnormalized mean expression matrix for comparison analyses.

---

## Part 2: MAGMA-Style Gene Property Enrichment

### Goal
Test whether genes that are specifically expressed in each cell type carry more SCZ GWAS signal than expected.

### Data
- Specificity matrix (Part 1)
- MAGMA gene-level results: for each of ~18,000 genes, a Z-statistic summarizing the total GWAS association signal across all SNPs within and near the gene
- Gene covariates: gene size, number of SNPs, sample size

### Method
For each cell type, run an OLS regression:

```
MAGMA_Z ~ specificity_in_type + log(gene_size) + log(n_snps) + log(N)
```

The t-statistic on the specificity coefficient tells us: **after accounting for confounds, do genes more specific to this type have stronger SCZ association?**

### Multiple testing
- Bonferroni correction (p < 0.05/137)
- FDR (Benjamini-Hochberg, q < 0.05)

### Result
A ranked list of all 137 cell types by SCZ enrichment significance. GABAergic interneurons (especially SST subtypes) emerge as the most enriched class.

---

## Part 3: Conditional Analysis (Forward Selection)

### Goal
Identify which enriched cell types carry *independent* signal (vs. correlated signal from shared gene expression).

### Problem
Many SST subtypes are enriched, but they share overlapping gene expression programs. If we condition on one SST type, do the others remain significant?

### Method
Forward stepwise selection:
1. Start with the most significant type
2. Add it to the "accepted" set
3. Re-test all remaining types, conditioning on the accepted set (include accepted types' specificities as additional covariates in the regression)
4. Among types still significant at p < 0.05, accept the most significant
5. Repeat until no types remain significant

### Result
A smaller set of independently enriched cell types — each carrying non-redundant SCZ signal.

---

## Part 4: Gene Driver Analysis

### Goal
Identify which genes are driving the enrichment in each independent cell type.

### Method
For each gene in each enriched type, compute:
- **Contribution** = specificity x |MAGMA_Z| — genes that are both highly specific AND strongly associated with SCZ rank highest
- **Uniqueness** = contribution weighted by how exclusively a gene contributes to *this* type vs. all enriched types

### Cross-type comparisons
- **Jaccard similarity**: Pairwise overlap of top-200 gene sets across cell types. High Jaccard = shared enrichment program; low = type-specific
- **Pan-interneuron genes**: Genes appearing in the top sets of >= 50% of enriched types (shared machinery)
- **Type-specific genes**: Genes with high contribution in one type but not others

---

## Part 5: Spatial Layer Analysis of SST Subtypes

### Goal
Determine whether SCZ enrichment in SST interneurons is uniform across cortical layers or concentrated in specific layers.

### Data
- Spatial transcriptomics (MERFISH + Xenium): continuous cortical depth estimate (0 = pial surface, 1 = white matter) for each of the 18 SST supertypes

### Method
1. Classify SST supertypes by layer position:
   - Upper (depth < 0.35): layers 1-3
   - Middle (0.35-0.55): layer 4-5a
   - Deep (>= 0.55): layers 5b-6
2. Compare SCZ enrichment between layer groups
3. Compute **upper-layer gene scores**: for each gene, what fraction of its SST expression is in upper-layer types, weighted by GWAS signal

### Finding
Upper-layer SST types tend to carry the strongest SCZ enrichment signal.

---

## Part 6: Gene-Depth Volcano Plots

### Goal
For every gene in the genome, ask: is this gene preferentially expressed in upper-layer vs. deep-layer SST interneurons, and does it carry SCZ risk or protective signal?

### Method
1. For each gene, compute **Spearman correlation** between its expression (or specificity) across 18 SST supertypes and their continuous cortical depth
   - Negative rho = upper-layer enriched (near pia)
   - Positive rho = deep-layer enriched (near white matter)
2. FDR-correct across all genes
3. Annotate each gene with:
   - Whether it's in the SCZ GWAS gene set (genome-wide significant locus)
   - Its MAGMA Z-statistic (strength of association)
   - Its **direction of effect** (risk vs. protective, see Part 7)

### Visualization
Volcano plot with:
- X-axis: Spearman rho (gene vs. depth)
- Y-axis: -log10(p-value)
- Point color: red (risk) vs. blue (protective)
- Point size: scales with |MAGMA Z| (GWAS signal strength)
- Black edges: SCZ GWAS loci genes

Two versions generated: one using specificity scores (row-normalized), one using absolute mean expression. The two are highly concordant (Spearman rho = 0.986).

---

## Part 7: Determining Risk vs. Protective Direction

A central question: for a gene associated with SCZ, does the genetic variation *increase* or *decrease* disease risk? Three approaches were considered, in order of sophistication:

### Approach A: Lead-SNP BETA (simplest, used as fallback)

For each gene, find the most significant SNP within gene boundaries from the PGC3 summary statistics. The sign of its BETA (log odds ratio) indicates direction:
- BETA > 0: the tested allele increases SCZ risk
- BETA < 0: the tested allele is protective

**Limitations**: The lead SNP may not be causal; genes can harbor both risk and protective variants; the sign depends on arbitrary allele labeling. This is a rough heuristic.

### Approach B: Signed MAGMA Z (intermediate)

Combine the direction from the lead SNP with the magnitude from MAGMA:

```
signed_MAGMA_Z = sign(lead_SNP_BETA) x |MAGMA_Z|
```

This gives a continuous metric where large positive values = strong risk signal, large negative = strong protective signal. But it still inherits the lead-SNP limitations, and MAGMA's Z-statistic is fundamentally non-directional (it tests association strength, not direction).

### Approach C: S-PrediXcan (most rigorous, used as primary)

S-PrediXcan asks: **"If your genetics cause higher expression of this gene in brain cortex, does that make you more or less likely to get SCZ?"**

#### How it works

**Step 1** (done once, from GTEx reference): Train a model predicting each gene's expression from nearby SNPs:
```
Predicted expression of Gene X = w1*SNP1 + w2*SNP2 + ...
```
These weights (w) capture the cis-eQTL architecture — which SNPs affect expression, and by how much. Models come from GTEx v8 brain tissue (MASHR fine-mapped models, ~1-5 SNPs per gene).

**Step 2** (using GWAS): Combine eQTL weights with GWAS z-scores without needing individual data:
```
z_gene = (w' * z_gwas) / sqrt(w' * Sigma * w)
```
where Sigma is the LD covariance matrix between the SNPs.

This is mathematically equivalent to: "predict each person's expression from their genotype, then test whether predicted expression correlates with SCZ status."

#### Interpretation
- **z > 0**: Genetically predicted *upregulation* increases SCZ risk
- **z < 0**: Genetically predicted *upregulation* is protective (i.e., *downregulation* increases risk)

#### Why this is better
S-PrediXcan explicitly models the **expression-mediated** causal pathway (SNPs -> expression -> disease), rather than just testing proximity (SNPs near gene -> disease). The direction is biologically interpretable: it tells you whether turning a gene *up* or *down* matters for disease.

#### Implementation
We ran S-PrediXcan for two GTEx brain tissues:
- **Brain_Cortex**: 12,086 genes tested, 5,480 FDR-significant
- **Brain_Frontal_Cortex_BA9**: 11,794 genes tested, 5,389 FDR-significant
- Cross-tissue correlation: Spearman rho = 0.69

The GTEx v8 MASHR prediction models were downloaded from Zenodo. We wrote a custom vectorized implementation because the official MetaXcan software had compatibility issues.

#### Validation
C4A (complement component 4A) shows z = +26 (upregulation -> risk), consistent with the landmark finding (Sekar et al. 2016) that increased C4A copy number / expression increases SCZ risk through excessive synaptic pruning.

#### Limitation
Not all genes have cortical eQTL models. Key genes like HCN1, TCF4, SST, and BCL11B are absent — their expression variation may be driven by trans-regulation, cell-type-specific eQTLs diluted in bulk tissue, or low expression heritability. For these genes, we fall back to the lead-SNP BETA approach.

### Integration in the volcano plot

Direction is assigned with a priority system:
1. **S-PrediXcan z-score** (10,697 genes) — expression-mediated, well-calibrated
2. **Lead-SNP BETA direction** (6,334 additional genes) — heuristic fallback
3. **No direction data** (15,419 genes)

---

## Part 8: Gene-Electrophysiology Integration

### Goal
Connect the genetic and transcriptomic findings to functional neurophysiology using patch-seq data.

### Data
- Patch-seq recordings from SST interneurons: simultaneous RNA-seq + electrophysiology
- Key features: SAG (voltage sag, driven by HCN channels) and TAU (membrane time constant)

### Method
For each gene, compute Spearman correlation between its SEA-AD specificity across SST supertypes and the mean electrophysiology feature of those supertypes (from patch-seq).

### Key finding
HCN1 (a known SCZ risk gene and upper-layer SST marker) is strongly correlated with SAG — providing a direct link from genetics (GWAS hit) through transcription (HCN1 expression in upper-layer SST) to electrophysiology (voltage sag properties).

---

## Summary of Data Flow

```
SEA-AD snRNA-seq (137K cells)
    |
    v
[1] Specificity matrix (genes x 137 types)
    |
    +---> [2] MAGMA enrichment (which types carry SCZ signal?)
    |         |
    |         +---> [3] Conditional analysis (which are independent?)
    |         |         |
    |         |         +---> [4] Gene drivers (which genes drive it?)
    |         |
    |         +---> [5] Spatial layer analysis (which layers?)
    |
    +---> [6] Gene-depth volcano (upper vs deep layer genes)
    |         |
    |         +---> Colored by risk/protective direction from:
    |                   [7a] S-PrediXcan (GTEx eQTL x PGC3 GWAS)
    |                   [7b] Lead-SNP BETA (PGC3 SNP-level fallback)
    |
    +---> [8] Gene-ephys correlations (patch-seq functional link)
```

---

## Part 9: ATAC-seq Peak Overlap with Fine-Mapped Variants

### Goal
Test whether SCZ risk variants preferentially fall in open chromatin of specific cell types, providing an orthogonal (chromatin-based) measure of cell-type enrichment.

### Data
- **SEA-AD snATAC-seq** (516,389 nuclei × 218,882 peaks, 139 supertypes)
- **FINEMAP credible sets** (20,591 SNPs across 249 PGC3 SCZ loci with posterior inclusion probabilities)

### Method
1. **Peak-variant overlap**: For each fine-mapped variant, find ATAC peaks that contain its genomic position
2. **Per-cell-type accessibility**: For each overlapping peak, compute mean accessibility per supertype from the sparse peak × cell matrix
3. **PIP-weighted scoring**: Score = Σ(PIP × mean_accessibility) across all overlapping peaks. This weights each variant by its probability of being causal.

### Key Result
Expression-based enrichment (MAGMA) and chromatin-based enrichment (ATAC) identify **different cell type classes**:
- **MAGMA (expression)**: GABAergic interneurons (Sst, Pvalb, Lamp5) dominate
- **ATAC (chromatin)**: Glutamatergic neurons (L5 IT, L2/3 IT, L5 ET) dominate

The two approaches show **near-zero correlation** (Spearman r ≈ 0.09), suggesting SCZ risk variants sit in open chromatin of excitatory neurons but affect genes specifically expressed in inhibitory neurons — consistent with a regulatory model where variants in excitatory neuron enhancers control genes whose products act in interneurons.

### Output
- `results/tables/atac_pip_weighted_scores.csv` — PIP-weighted ATAC scores per cell type
- `results/intermediates/atac_finemap_peak_overlaps.csv` — All variant-peak overlaps
- `results/figures/atac_vs_gwas_manhattan.png` — Side-by-side comparison
- `results/figures/atac_vs_gwas_enrichment.png` — Scatter of ATAC vs GWAS scores

---

## Part 10: GWAS vs Case-Control Composition

### Goal
Test whether cell types with stronger SCZ GWAS enrichment also show altered proportions in SCZ brains.

### Data
- **GWAS enrichment**: RBH combined taxonomy (503 types) from script 10
- **Case-control composition**: Endresz et al. (in prep) 7-cohort snRNA-seq meta-analysis (crumblr model, 109 matched cell types)

### Method
Spearman correlation between GWAS enrichment metrics and composition change metrics, across all matched types and within SST interneurons specifically.

### Key Result
- **All cell types** (n=109): No significant correlation (r = -0.10, p = 0.28)
- **SST interneurons** (n=18): Strong positive correlation between |composition change| and GWAS -log10(p) (r = 0.645, p = 3.85e-03), and between significance measures (r = 0.701, p = 1.20e-03)

SST subtypes most enriched for SCZ genetic risk (Sst_2, Sst_25, Sst_3) are also depleted in SCZ brains — convergence of genetic predisposition and observed cellular pathology.

### Output
- `results/tables/gwas_vs_casecontrol_composition.csv`
- `results/figures/gwas_vs_casecontrol_composition.png`

See `docs/rnaseq_specificity_enrichment_results.md` for detailed tables and interpretation.

---

## Key Software & Dependencies

- **Python**: scanpy/anndata ecosystem for single-cell data
- **MAGMA**: Gene-level association from GWAS (pre-computed)
- **GTEx v8 MASHR models**: eQTL prediction weights (from Zenodo)
- **S-PrediXcan**: Custom implementation (official MetaXcan had compatibility issues)
- **adjustText**: Non-overlapping gene labels on volcano plots
- **statsmodels**: FDR correction, OLS regression
- **scipy**: Spearman correlations, normal distribution
