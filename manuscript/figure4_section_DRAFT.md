# Figure 4 — results section and methods (rewrite, 2026-08-07)

Every statistic below was extracted programmatically from the rendered panel inputs in
`magma_pgc3/r_panels_bigdeli/`. Panel letters follow the nine-panel figure. Items in
`[brackets]` still need a decision or a citation.

---

## Results

### Schizophrenia genetic risk converges on HCN1-expressing upper-layer Sst supertypes that are also depleted in Alzheimer's disease

Lastly, we asked whether the depletion of these upper-layer Sst supertypes is
etiologically relevant to SCZ, or instead a downstream consequence of the illness or its
treatment. Common-variant risk offers a potential avenue for addressing this: because such
risk is germline and therefore precedes the disorder, the cell types it implicates are more
plausibly part of the causal pathway than a downstream consequence [(CITE)]. Prior work has
implicated Sst interneurons in particular. For example, mapping SCZ risk across 461
whole-brain cell types (Siletti et al. 2023) showed that MGE-derived interneuron
populations carry the strongest genetic association of any brain cell type (Duncan et al.
2025); however, because this was assessed in a neurotypical atlas, it could not establish
whether these same populations are also the ones most changed in disease. We therefore
asked whether the Sst supertypes we find compositionally depleted in case/control data are
also those most enriched for SCZ genetic risk. To this end we applied the enrichment
pipeline of Duncan et al. without modification to any analytical step, changing only the
expression data and the cell-type taxonomy so that both matched the rest of our study
(Methods; Supplementary Methods SM3; Supplementary Table T4).

Across the 16 Sst supertypes, those carrying the greatest SCZ common-variant risk —
including Sst_23, Sst_3, Sst_2 and Sst_20 — were generally the ones most compositionally
depleted in SCZ relative to controls (Fig. 4a; Spearman ρ = 0.66, p = 0.0069). This
relationship was graded by cortical depth: SCZ enrichment was strongest among the most
superficial Sst supertypes and declined with increasing soma depth (ρ = −0.82,
p = 1.0 × 10⁻⁴, n = 16; [Supplementary Fig. SX]).

To identify the genes driving this enrichment we decomposed the genetic signal for Sst_3,
the most strongly enriched Sst supertype with reconstructed morphologies available, into
the individual genes that are both specifically expressed in Sst_3 and genetically
associated with SCZ (Fig. 4b; Methods). Of 345 such driver genes, this analysis nominated
HCN1, a subunit of the hyperpolarization-activated current I_h, which is specifically
expressed in Sst_3 (93.5th percentile of gene specificity), carries a strong gene-level SCZ
association (MAGMA −log₁₀p = 10.9), and lies under a fine-mapped SCZ locus that nominates
HCN1 as the likely causal gene (Fig. 4c; 4-variant credible set, cumulative posterior
inclusion probability 0.985; lead variant rs10035564, PIP = 0.56).

Consistent with a role for I_h in these cells, patch-seq data from an independent set of
150 human neocortical Sst interneurons (Lee, Dalley et al. 2023; Chartrand et al. 2023)
showed that voltage sag, the electrophysiological signature of HCN channel function, was
greater in Sst supertypes with higher HCN1 expression (Fig. 4d; ρ = 0.65, p = 0.0064;
exemplar voltage traces in Fig. 4f).

We next asked what other genes distinguish these cells, and identified 580 genes
differentially expressed between depleted and not-depleted Sst supertypes (FDR < 0.05,
|log₂FC| > 0.25; Fig. 4g). HCN1 was itself among them, expressed more highly in the
depleted supertypes (log₂FC = +0.43, FDR = 1.9 × 10⁻¹³⁸). Also among these genes was
CALB1, encoding the calcium-binding protein calbindin (log₂FC = +0.55,
FDR = 4.8 × 10⁻¹⁸³; Fig. 4h), a previously established marker of double-bouquet cells
(DBCs), an upper cortical layer morphological type suggested to be primate-specialized as
it has no clear rodent counterpart [(CITE)]. Consistent with this, the depleted supertypes
include cells with canonical DBC-like morphologies, and the reconstructed exemplars of
depleted and not-depleted supertypes differ systematically in both laminar position and
voltage sag (Fig. 4e,f).

Finally, we asked whether the Sst supertypes depleted in schizophrenia are also depleted in
Alzheimer's disease, where Sst subtypes have likewise been reported to be lost [(CITE)]. We
turned to the SEA-AD resource, which profiled the DLPFC of 80 aged donors spanning a
continuum of AD neuropathology using snRNAseq and assigned each donor a continuous
pseudo-progression score (CPS) indexing degree of AD-related neuropathological severity.
Because SEA-AD employs the same supertype taxonomy used throughout this study, we could
relate each supertype's SCZ case-control composition change directly to the same
supertype's compositional depletion along the AD CPS trajectory. Among Sst supertypes, the
SCZ- and AD-associated depletion in frontal cortex were strongly concordant, agreeing in
direction for 94% of supertypes, and the Sst supertypes we found depleted in schizophrenia
were the same ones that declined most as Alzheimer's disease advances (Spearman ρ = 0.87,
p = 1.3 × 10⁻⁵; n = 16 supertypes; Fig. 4i). These upper-layer Sst supertypes therefore
likely constitute an intrinsically vulnerable population, depleted across disorders with
distinct etiologies and neuropathology.

---

## Methods — SCZ GWAS enrichment and patch-seq integration

**GWAS summary statistics.** Cell-type enrichment used the European-ancestry
autosomal meta-analysis of Bigdeli et al. 2026 (GRCh38; maximum per-variant effective
N = 114,827). The larger PGC3 analysis (Trubetskoy et al. 2022) was retained as a
robustness comparison (Supplementary Methods SM3). We note that the Bigdeli
European meta-analysis is not simply a better-powered version of PGC3: it carries 1.57×
PGC3's effective sample size but only 1.07× the polygenic signal (mean gene χ² − 1), because
the additional cohorts are biobank- and EHR-ascertained rather than clinically ascertained.
We therefore treat the two as complementary rather than ordered.

**Enrichment pipeline.** We applied the pipeline of Duncan et al. 2025
(`github.com/Integrative-Mental-Health-Lab/linking_cell_types_to_brain_phenotypes`)
without modification to any analytical step, changing only the expression data and the
cell-type taxonomy. Specificity was computed by their recipe: per-cell ln(1 + x) on raw
counts averaged within cell type; restriction to genes with unique names and unique
ENTREZ↔ENSEMBL mappings (`org.Hs.eg.db`); removal of unexpressed genes and of the extended
MHC (chr6:25–34 Mb); scaling of each cell type to a common total; and per-gene
normalization across types, with the ENTREZ join performed last. Enrichment was then tested
with MAGMA v1.10 using their exact invocations — SNP-to-gene annotation with a 35 kb
upstream / 10 kb downstream window against `NCBI37.3.gene.loc`; gene analysis with the
SNP-wise mean model against the 1000 Genomes phase 3 European LD panel; and gene-property
analysis with `--model direction=greater --gene-covar`. Because MAGMA's own gene-property
analysis is used rather than a reimplementation, the analysis inherits its defaults
unchanged: truncation of gene Z-scores at 3 SD below and 6 SD above the mean, truncation of
specificity values at 5 SD from the mean, conditioning on gene size, gene density, sample
size and inverse mean minor allele count together with their logarithms, and a
multivariate-normal error model accounting for LD between genes. Of 19,427 gene locations,
19,364 received at least one SNP and 18,534 contained valid SNPs in the genotype data;
15,981 genes were common to the gene-level results and the specificity matrix and entered
the gene-property regression. P-values were corrected across all tested cell types by both
Bonferroni and Benjamini–Hochberg FDR.

Applying our implementation to the Siletti 461-cluster taxonomy of Duncan et al.
reproduced their published enrichment estimates at Pearson r = 0.99, with 17 of 461 cell
types differing in FDR significance. The residual is attributable to a newer `org.Hs.eg.db`
release changing which genes have unique ENTREZ↔ENSEMBL mappings; exact reproduction would
require the 2020 annotation.

**Expression reference and taxonomy.** Specificity was computed over a combined 500-type
taxonomy comprising 125 SEA-AD supertypes and 375 non-redundant Siletti et al. 2023
whole-brain clusters, across 16,666 genes. SEA-AD expression came from the dorsolateral
prefrontal cortex (A9) release, restricted to its three neurotypical reference donors
(90,579 nuclei), so that cell-type specificity is measured in the same cortical region as
every disease measurement in this study; three of the five SEA-AD MTG reference donors are
the same individuals, making this close to a within-donor comparison. Redundancy between
the two references was removed by reciprocal-best-hit matching on per-type mean expression
over the 5,000 most variable shared genes: 95 Siletti clusters were matched to a SEA-AD
supertype (median Spearman r = 0.90) and removed in favor of the more finely resolved
SEA-AD annotation. [DECIDE: whether to report the eight merges below r = 0.8 — all
microglial or vascular — and whether to retain them as separate types.]

**Fine-mapping.** Credible sets at the HCN1 locus are the SuSiE-R European-ancestry results
of Bigdeli et al. 2026 (their Supplementary Table 13), ancestry-matched to the summary
statistics used for enrichment. Their tables additionally report this locus under five other
method and ancestry combinations, with posterior inclusion probabilities for the lead
variant ranging from 0.001 to 0.841 (Supplementary Methods SM3).

**Differential expression between depleted and not-depleted Sst supertypes.** Supertypes
were assigned to depleted and not-depleted groups from the compositional analysis
(Fig. 3). Counts were normalized to 10,000 per nucleus and log1p-transformed, and genes
tested by cell-level Wilcoxon rank-sum with Benjamini–Hochberg correction, with
Seurat-style avg_log₂FC computed on the un-logged CP10K scale. Genes were called
differentially expressed at FDR < 0.05 and |log₂FC| > 0.25. Because the DLPFC reference
comprises three donors, a donor-level paired design is underpowered for genome-wide
correction; donor-paired tests are reported for named genes (CALB1 p = 0.0030,
HCN1 p = 0.0119) but were not used to define the gene set.

**Patch-seq electrophysiology.** [Unchanged from the current draft — 1,155 cells across
221 donors; scANVI and Harmony-corrected kNN assignment; sag ratio from NWB sweeps with the
Allen IPFX library. Note that these recordings are from temporal cortex (MTG), whereas
specificity is measured in frontal cortex.]

---

## Open items in this draft

1. `[CITE]` markers: germline-risk rationale, DBC primate-specialization, AD Sst loss.
2. `[Supplementary Fig. SX]` for the depth-gradient statistic (ρ = −0.82, p = 1.0 × 10⁻⁴) —
   the panel is built but not yet numbered.
3. Whether to report the eight sub-threshold RBH merges in the methods.
4. Verify the Duncan et al. characterization ("MGE-derived interneuron populations carry
   the strongest genetic association") against their results text — our reproduction of
   their 461-cluster analysis gives MGE/CGE clusters as the top hits, but the original
   wording said "cortical Sst subtypes."
5. Supplementary Methods SM2 must be corrected: the "mean AUROC ≥ 0.9" quality bar is
   vacuous (the statistic is 1.0 by construction for every best hit), and the metric is a
   percentile-rank-of-correlation rather than MetaNeighbor neighbor-voting AUROC, so the
   Crow et al. 2018 citation overstates it.
