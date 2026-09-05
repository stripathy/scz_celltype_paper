> ## ⚠️ STALE — the live manuscript is the Google Doc (`1cO5ZStbp9b3cb6vb9s6H2YfbprYGikXBvFHQR7gROjo`)
>
> Since 2026-08-20 every number, figure reference and decision lives in the Doc; this
> file is the 2026-07-31 assembled draft, committed for provenance only. Do not edit it
> or quote from it. The dated `figure5_*_2026-09-01.md`, `figure3_tail_paragraph_*.md`
> and `*_PLAN.md` files in this folder are decision records that were pasted into the
> Doc; treat them the same way.

# Resolving cellular pathology in schizophrenia through single-cell meta-analysis

*(Working title.)*

## Abstract

Schizophrenia (SCZ) is associated with cortical GABAergic dysfunction, but it has
remained unresolved whether inhibitory interneurons are lost or instead persist in
an altered molecular state. Using fine-grained, transcriptomically defined
cell-type taxonomies (SEA-AD supertypes) and cross-dataset meta-analysis, we
integrated cell-type–specific differential expression and cellular composition
across seven post-mortem single-nucleus RNA-seq datasets (298 controls, 171 SCZ
cases), spatial replication by Xenium (12 SCZ, 12 controls), and SCZ
common-variant genetics anchored to human patch-seq electrophysiology and
morphology. We resolve the debate in favour of both phenomena acting on different
cells: a subclass-wide, cell-intrinsic reduction of *SST* mRNA within Sst neurons
coexists with a numerically restricted depletion of upper-layer Sst supertypes —
most robustly Sst_25 `[TODO: Sst_25 effect size and 95% CI from the composition
meta-analysis]` — alongside an increase in deep-layer L6b excitatory neurons. Both
compositional changes replicate spatially. SCZ common-variant risk is enriched in
the same upper-layer Sst supertypes, and genetic enrichment tracks compositional
depletion. Expression of the I_h channel gene *HCN1* tracks these cells' sag
physiology, and the vulnerable, *CALB1*-expressing population corresponds to
primate-specialized double-bouquet interneurons. The vulnerable upper-layer Sst
supertypes overlap with those lost earliest in Alzheimer's disease, indicating a
transdiagnostic axis of upper-layer cortical inhibitory vulnerability.

## Introduction

Schizophrenia (SCZ) is a chronic neuropsychiatric disorder affecting approximately
0.5–1% of adults worldwide. It is characterized by positive symptoms such as
hallucinations and delusions, negative symptoms such as blunted affect and
anhedonia, and pervasive cognitive impairment, and it typically requires lifelong
treatment [ref]. Available therapies — principally second-generation antipsychotics
targeting dopaminergic signaling — leave a substantial fraction of patients
inadequately treated: roughly 10–30% show frank treatment resistance and a further
30–60% achieve only partial responses or cannot tolerate side effects [ref].
Progress toward mechanistically novel therapies has been slowed by an incomplete
picture of the specific cortical cell types and molecular changes that underlie the
disorder.

Unlike Alzheimer's or Parkinson's disease, schizophrenia does not involve gross
neurodegeneration, yet it is consistently associated with dysfunction of cortical
circuits — particularly in the prefrontal cortex (PFC) — and with impairments in
γ-aminobutyric acid (GABA)-mediated inhibition [ref]. GABAergic interneurons
provide the inhibitory control that shapes cortical computation, and their
disruption is among the most reproducible neurobiological findings in
schizophrenia, with the somatostatin (SST) and parvalbumin (PVALB) interneuron
classes most frequently implicated [ref]. How these interneurons are altered,
however, has remained contested. Some studies report a reduced *number* of
GABAergic interneurons (Toker et al., 2017; Batiuk et al., 2022; Kiss et al.,
2025), whereas others find altered interneuron *gene expression* — notably lowered
SST and GABA-synthetic transcripts — without a detectable loss of cells (Dienel et
al., 2022). These observations motivate three competing possibilities: that
interneurons are fewer in number, that interneuron number is preserved but their
molecular state is altered, or that both occur together, potentially in different
cells.

Resolving this question requires the ability to distinguish closely related
inhibitory cell types, because a change confined to a specific subset of neurons
can be masked or diluted when cells are analyzed only at a coarse level. Recent
single-cell and single-nucleus transcriptomic atlases now provide taxonomies of
human neocortex at unprecedented granularity: the Seattle Alzheimer's Disease Brain
Cell Atlas (SEA-AD) expands the conventional ~24 cortical subclasses (e.g., Sst)
into 137 finer "supertypes" (e.g., Sst_25), many of which are spatially restricted
to particular cortical layers (Gabitto et al., 2024). Complementary patch-seq
datasets link these transcriptomic types to their electrophysiological,
morphological, and laminar properties (Lee, Dalley et al., 2023), so that a
transcriptomic label can be tied to a concrete cellular phenotype. In parallel, the
rapid growth of publicly available post-mortem SCZ snRNA-seq datasets creates an
opportunity to test cell-type–specific hypotheses at scale and to mitigate the
small sample sizes and dataset heterogeneity that have limited individual studies.

Here we combine fine-grained cell-type taxonomies with cross-dataset meta-analysis
to ask (i) whether granular cell taxonomies can reconcile the competing hypotheses
of how inhibitory cells change in schizophrenia, and (ii) whether meta-analysis
across many post-mortem datasets can overcome the heterogeneity that has
historically limited reproducibility. We integrate cell-type–specific differential
expression and cellular composition across seven snRNA-seq datasets, replicate the
compositional and expression findings in an independent Xenium spatial
transcriptomics dataset, and connect the implicated cell types to schizophrenia
common-variant genetic risk and to their patch-seq-defined physiology and
morphology.

## Results

### A unified cell-type taxonomy resolves cells across seven snRNA-seq datasets and validated Xenium sections

To resolve how cortical cell types change in schizophrenia (SCZ), we assembled seven post-mortem snRNA-seq datasets of frontal cortex — 469 donors in total (298 neurotypical controls and 171 SCZ cases; approximately 2.3 million nuclei `[TODO: confirm total nuclei from Nicole's pipeline]`) — and mapped every nucleus onto a single reference taxonomy, the SEA-AD middle temporal gyrus (MTG) supertypes, which refine the conventional 24 subclasses (e.g., Sst) into 137 finer supertypes (e.g., Sst_25; Gabitto et al., 2024). After reference-based label transfer and integration, nuclei grouped by cell type rather than by dataset of origin or diagnosis (Fig. 1a–d). Critically, this harmonization removed technical and dataset structure while preserving biological identity, providing the common cell-type framework required for a cell-type-resolved meta-analysis across datasets.

To test whether the snRNA-seq findings replicate in an orthogonal, spatially resolved modality, we added a Xenium spatial transcriptomics cohort of 24 DLPFC sections (12 SCZ, 12 control; LIBD), profiling 1.34 million cells across a 300-gene panel. Reference-based label transfer from the SEA-AD MTG reference assigned each cell to one of 24 subclasses and, within subclass, to one of 137 supertypes (Fig. 1e), and we benchmarked these labels against the independent SEA-AD MERFISH atlas. The annotations were biologically accurate: per-donor subclass proportions agreed at Pearson *r* = 0.85, subclass cortical depth at *r* = 0.96, and label transfer reproduced MERFISH subclass calls at 84.9% accuracy. To place cells in their laminar context, we further partitioned each section into its cortical laminae (L1–L6) plus white-matter and vascular domains (Fig. 1f). Together these panels establish accurate, laminarly organized cell-type annotation at subclass resolution in the Xenium data — the foundation for the compositional replication in Figure 3. Supertype-resolution annotation, agreement with the original authors' annotations, and the panel's within-subclass resolution limits are documented in Supplementary Fig. S4.

### *SST* mRNA is down-regulated in somatostatin interneurons and replicates in Xenium

To test whether SCZ alters expression of the canonical inhibitory markers within the cell types that express them, we performed a cell-type–specific differential-expression (DE) meta-analysis across the seven frontal-cortex snRNA-seq datasets and sought independent replication in the Xenium cohort. Within Sst interneurons, *SST* mRNA was significantly reduced in SCZ (meta-analytic log₂ fold change = −0.46, *P* = 6.3 × 10⁻⁴, FDR = 0.049; Fig. 2a,b). This reduction was directionally consistent across all seven datasets and reached nominal significance (*P* < 0.05) in three: Fröhlich (log₂FC = −0.57, *P* = 0.016), HBCC (−0.61, *P* = 1.1 × 10⁻³) and McLean (−0.85, *P* = 0.026) (Fig. 2b). Notably, Xenium recovered the same direction and magnitude (log₂FC = −0.32, *P* = 0.052, panel-wide FDR = 0.30; 12 vs 12 donors), narrowly short of nominal significance (Fig. 2b–d). Beyond the marker itself, Sst cells down-regulated the mitochondrial protease *AFG3L2* (log₂FC = −0.11, FDR = 0.030) and *NAT16* (−0.81, FDR = 4.6 × 10⁻⁵) and up-regulated *SLC9A9* (+0.25, FDR = 4.3 × 10⁻³), *STAC* (+0.54, FDR = 0.042) and *SMAD1* (+0.22, FDR = 6.6 × 10⁻⁵) (Fig. 2a).

The second canonical interneuron marker, *PVALB*, showed the complementary pattern: undetectable in snRNA-seq, reduced in Xenium. *PVALB* was not differentially expressed within Pvalb interneurons in the snRNA-seq meta-analysis (log₂FC = −0.06, *P* = 0.54, FDR = 0.86, n.s.) and was not significant in any individual dataset (Fig. 2e,f). Conversely, in Xenium *PVALB* was significantly decreased within Pvalb cells (log₂FC = −0.22, *P* = 0.044, panel-wide FDR = 0.39; Fig. 2f,g) — consistent with a partly cytoplasmic signal that Xenium captures but nucleus-restricted snRNA-seq does not, and, for a type-defining marker, with mild label-transfer circularity (Discussion). Other DE genes within Pvalb cells included down-regulated *ANXA2* (−0.39, FDR = 6.8 × 10⁻³), *NAT16* (−0.49, FDR = 0.021), *VGF* (−0.60, FDR = 0.036) and *CIRBP* (−0.18, FDR = 6.8 × 10⁻³), and up-regulated *SCN3A* (+0.17, FDR = 7.5 × 10⁻⁴), *SMAD1* (+0.19, FDR = 7.5 × 10⁻⁴) and *TCAF2* (+0.40, FDR = 0.017) (Fig. 2e).

Beyond the two marker genes, the meta-analysis identified 6,903 DE gene × subclass associations transcriptome-wide at FDR < 0.10 (3,586 at FDR < 0.05), spanning the 22 of 23 tested subclasses with at least one DE gene (Fig. 2i), with the signal dominated by down-regulation (3,916 down vs 2,987 up at FDR < 0.10). This burden tracked cell-type abundance rather than selective vulnerability: astrocytes carried the most DE genes (943; 619 down-regulated) and were also the most abundant subclass, and across subclasses the DE-gene count was strongly correlated with mean per-donor cell-type proportion (Spearman ρ = 0.83, *P* = 2.1 × 10⁻⁶, n = 22 subclasses; Fig. 2i inset). Critically, raw DE tallies therefore largely index the number of nuclei sampled per subclass, and effect sizes — not gene counts — carry the biology.

To assess whether these effects reproduce across platforms, we compared snRNA-seq meta-analytic and Xenium DE effects across the meta-DE gene × subclass pairs testable on both. Although the 300-gene Xenium panel restricts cross-validation to a subset of meta-DE genes, the two platforms agreed closely where they overlapped. Across all 166 gene × subclass pairs that were meta-DE (FDR < 0.10) and testable on the Xenium panel, snRNA-seq meta-analytic and Xenium log₂ fold changes were positively correlated (Pearson *r* = 0.73; ordinary-least-squares slope = 0.81) and agreed in direction for 72% of pairs (120/166; sign-test *P* = 8.4 × 10⁻⁹; Fig. 2j). Concordant examples spanned every cell class — for example down-regulated *VGF* in Chandelier cells (meta log₂FC = −1.04, FDR = 0.002; Xenium −0.77, *P* = 0.027) and up-regulated *ATP2B4* in Sst cells (meta +0.21, FDR = 0.016; Xenium +0.31, *P* = 1.6 × 10⁻³). Importantly, two orthogonal technologies thus recapitulate the same cell-type–resolved transcriptional signature of SCZ. Within-subclass (supertype-level) DE was underpowered for most cell types, reflecting the smaller number of nuclei per donor at that resolution (Supplementary Fig. S8).

### Upper-layer SST supertypes are depleted and deep-layer L6b excitatory supertypes more abundant in schizophrenia

Reduced *SST* mRNA need not imply fewer cells; to test whether cell-type abundances themselves shift in SCZ, we performed a stratified crumblr compositional meta-analysis across the seven frontal-cortex snRNA-seq datasets (298 control / 171 SCZ donors; Endresz et al., this study). This meta-analysis identified a coordinated shift within the GABAergic and deep excitatory compartments. Several somatostatin (Sst) supertypes were significantly depleted in SCZ — Sst_2 (β = −0.25, FDR = 9 × 10⁻⁴), Sst_22 (β = −0.27, FDR = 0.022), Sst_25 (β = −0.26, FDR = 0.037) and Sst_3 (β = −0.19, FDR = 0.046) — while Sst_20 was depleted only at nominal significance (β = −0.17, *P* = 0.019, FDR = 0.19). Conversely, layer-6b (L6b) excitatory supertypes were more abundant: L6b_1 (β = +0.48, FDR = 2 × 10⁻⁴) and L6b_4 (β = +0.26, FDR = 0.022), with L6b_2 at trend level (β = +0.17, *P* = 0.055). These effects are individually subtle — reaching independent significance in only one dataset (Fröhlich; FDR < 0.10; Fig. 3b,c) — but reproducible across datasets, the regime where meta-analysis outperforms underpowered single datasets, motivating orthogonal replication in a spatial platform.

To test these compositional signatures on an orthogonal platform, we applied the same stratified crumblr analysis to the 24-section Xenium DLPFC cohort (12 SCZ / 12 control). At the supertype level, Xenium recovered the same directional pattern: the three Sst supertypes best powered on the panel were concordantly depleted — Sst_25 (logFC = −0.66, FDR = 0.028), Sst_3 (logFC = −0.56, FDR = 0.028) and Sst_2 (logFC = −0.31, *P* = 0.028, FDR = 0.15) — while the L6b supertypes were concordantly more abundant — L6b_4 (logFC = +1.05, FDR = 0.028), L6b_2 (logFC = +0.70, FDR = 0.028) and L6b_1 (logFC = +0.80, FDR = 0.056). Aggregating the five vulnerable Sst supertypes (Sst_2, Sst_20, Sst_22, Sst_25, Sst_3) into a single per-donor proportion, this aggregate was significantly depleted in SCZ (control 1.65% vs SCZ 1.30% of cortical cells; Welch's *t*-test *P* = 0.016, n = 12/12), with a concordant, non-significant reduction in density (Methods). Aggregating all L6b cells at the subclass level (excluding four sections with insufficient deep-cortex coverage; Methods), L6b cells were significantly more abundant in SCZ by both proportion (0.84% vs 1.22%; *P* = 0.009, n = 10/10) and density (5.16 vs 7.79 cells/mm²; *P* = 0.005). The Sst_25 depletion is visible in representative sections closest to each group's median Sst_25 proportion (control 0.199% vs SCZ 0.083%; Fig. 3d).

Two supertypes did not replicate in Xenium — Sst_22 (logFC = −0.04, *P* = 0.72) and Sst_20 (logFC = +0.10, *P* = 0.45) — precisely the two the 300-gene panel resolves worst (versus Sst_25, the best-resolved; Methods). The aggregate Sst depletion is therefore robust and cross-platform, while individual supertype allocation is reliable only for better-resolved types.

To quantify cross-platform agreement directly, we compared the SCZ compositional effect size for each cell type between the snRNA-seq meta-analysis and Xenium. Agreement was strong at the subclass level: the 17 neuronal subclasses correlated at Pearson *r* = 0.70 (*P* = 1.9 × 10⁻³) — deep excitatory types (L6b, L6 CT, L6 IT) more abundant and inhibitory types (Sst, Pvalb, Vip) less abundant on both platforms — whereas the six non-neuronal subclasses did not agree (*r* = −0.20, *P* = 0.71). At supertype resolution the neuronal effects correlated more moderately (*r* = 0.50, *P* = 4.2 × 10⁻⁸, n = 106 neuronal supertypes); the attenuation to supertype level reflects the panel's limited within-subclass markers (Methods). Taken together, the spatial data independently confirm the two headline compositional signatures of the snRNA-seq meta-analysis — vulnerable-Sst loss and L6b gain — in situ.

To test whether this SST depletion extends beyond DLPFC to a second cortical region, we re-analyzed FISH-based SST interneuron counts from the subgenual anterior cingulate cortex (sgACC; Arbabi et al. 2025; 15 control / 11 SCZ after laminar quality control) under the same upper-layer prediction. SST density was lower in SCZ (mixed-model β = −0.60, *P* = 0.071), with a larger deficit in upper (L2/3, β = −0.76) than deep (L5/6, β = −0.43) layers — individually subtle but directionally consistent with the DLPFC result, extending the upper-layer SST depletion to a second cortical region (Supplementary Fig. S14; Methods).

### Schizophrenia genetic risk converges on the vulnerable upper-layer SST cells

Recent work mapped schizophrenia common-variant risk onto specific brain cell types — most strongly a somatostatin (SST) interneuron subtype — using MAGMA gene-property analysis of the Siletti whole-brain single-nucleus atlas (461 cell types), and began to localize the associated SST subtypes to specific cortical layers (Duncan et al. 2025). To place this genetic signal in the cell-type taxonomy used throughout the present study, we ran the same analysis (MAGMA gene-property regression of the PGC3 SCZ GWAS on cell-type expression specificity; Trubetskoy et al. 2022; Methods) on our combined Franken-503 taxonomy (137 SEA-AD cortical supertypes + 366 non-redundant Siletti whole-brain clusters), which retains the whole-brain reference needed for well-calibrated specificity while adding SEA-AD's cortical resolution and enables a direct, supertype-by-supertype comparison with our composition and spatial results. Of the 503 types, 120 were FDR-significant (54 Bonferroni), overwhelmingly interneurons: 68 cortical SEA-AD supertypes (56 GABAergic) and 52 subcortical Siletti clusters. Among cortical supertypes the strongest signals were GABAergic interneurons — parvalbumin, Lamp5 and somatostatin — including the two SST supertypes central to this study, Sst_2 (β = 38.4, *p* = 3.8 × 10⁻¹¹) and Sst_25 (β = 33.7, *p* = 3.2 × 10⁻⁸). Brain-wide, the top signals were subcortical MGE-derived interneurons, reproducing the whole-brain associations of Duncan et al. (Supplementary Table T4; Supplementary Methods SM2).

To test whether genetic enrichment was spatially structured within the SST subclass, we positioned the 18 SST supertypes along cortical depth. Enrichment rose toward the pia in a superficial-to-deep gradient (Spearman *r* = −0.50, *p* = 0.034, n = 18): all 7 upper-layer SST types (depth < 0.35) were FDR-significant with mean −log₁₀*p* = 6.6, nearly double the 9 deep-layer types (mean 3.7), and 17 of 18 SST types were FDR-significant with 10 passing Bonferroni. Notably, these enriched upper-layer types were also the ones preferentially *depleted* in the 7-dataset case–control snRNA-seq composition meta-analysis (Endresz et al., in prep). Within SST (n = 18), the magnitude of compositional change tracked genetic enrichment (|composition β| vs GWAS −log₁₀*p*, Spearman *r* = 0.645, *p* = 3.85 × 10⁻³), and the most genetically enriched SST types were the most reduced in SCZ (signed composition β vs enrichment β, *r* = −0.494, *p* = 0.037). Because no such relationship held across all 109 matched types (*r* = −0.10, *p* = 0.28), this genetics↔depletion convergence is SST-specific, with Sst_2 (composition β = −0.255, FDR = 9.2 × 10⁻⁴) and Sst_25 (composition β = −0.259, FDR = 0.037) the clearest dual-hit types — top-ranked for genetic risk and significantly depleted in tissue.

To nominate candidate gene drivers of the enriched Sst_25 population, we decomposed its enrichment into its constituent driver genes (Methods), yielding a set of candidate drivers `[TODO: Sst_25 driver-gene count and HCN1's driver status pending regeneration of the Sst_25 driver panel — only the Sst_2 driver table (269 genes) is committed and HCN1 is absent from it]`. We focus on the hyperpolarization-activated channel gene *HCN1*, which lies under a credibly fine-mapped SCZ locus: the PGC3 FINEMAP 95% credible set indexed by lead variant rs10035564 comprises 9 variants carrying cumulative posterior probability 0.95, with the lead variant alone at PIP = 0.525 and 4 of the 9 variants intronic to *HCN1*. Importantly, the credible set overlapped SST open chromatin more than astrocyte chromatin (PIP-weighted snATAC accessibility 0.045 vs 0.006), consistent with an interneuron-acting regulatory effect.

To connect *HCN1* to a measurable physiological phenotype, we drew on patch-seq recordings of voltage sag — a hyperpolarization-activated (I_h/HCN-dependent) membrane property — in 150 human cortical SST interneurons spanning 16 SST supertypes (Lee, Dalley et al. 2023; Chartrand et al. 2023). Across supertypes, *HCN1* expression specificity was positively correlated with mean sag (Spearman ρ = 0.62, *p* = 0.010, n = 16 supertypes), and upper-layer SST supertypes showed higher mean sag (0.35) than deep-layer types (0.22). This places the genetically enriched, *HCN1*-high upper-layer SST populations at the high-sag end of the physiological range.

Two reconstructed patch-seq cells illustrate the extremes of this range: an upper-layer Sst_25 cell (specimen 819770858, layer 2, soma 405.6 µm from pia) with high sag, versus a deep Sst_5 cell (specimen 758996755, layer 4) with near-absent sag `[TODO: single-cell sag values 0.564/0.00196 are methods-draft only — verify from the patch-seq source or report supertype means]`. Their supertype means bracket the population (Sst_25 mean sag = 0.456, n = 43 cells; Sst_5 mean sag = 0.167, n = 10 cells). Together, panels a–g trace a coherent chain from SCZ common-variant risk, through *HCN1* expression in upper-layer SST, to a measurable electrophysiological phenotype (sag) in the same cells that are genetically enriched and compositionally depleted in schizophrenia.

We repeated the enrichment under the larger Bigdeli et al. 2026 SCZ GWAS, with largely concordant results (per-cell-type ρ = 0.95 vs PGC3; Supplementary Fig. S13; Methods) `[TODO: finalize primary GWAS — PGC3 vs Bigdeli 2026]`.

## Discussion

By combining fine-grained, transcriptomically defined cell-type taxonomies with cross-dataset meta-analysis, we reconcile two long-competing accounts of GABAergic pathology in schizophrenia (SCZ): a loss of inhibitory interneurons versus an altered molecular state of interneurons that remain present. Across four largely orthogonal lines of evidence — cell-type–specific differential expression (DE) from a seven-dataset single-nucleus RNA-seq meta-analysis, compositional meta-analysis of the same datasets, spatial replication with Xenium, and common-variant genetics anchored to patch-seq electrophysiology — both accounts hold, but for partly different cells: a subclass-wide reduction in somatostatin (*SST*) mRNA within Sst neurons co-occurs with a numerical depletion confined to specific, upper-layer Sst supertypes, set against a coordinated increase in deep-layer L6b excitatory neurons. That these signals converge across dissociation-based and spatial platforms, align with SCZ common-variant risk, and are echoed by an independent histological re-analysis argues they reflect genuine features of the disease rather than any single platform's artifact — and that the resolution they afford was largely hidden from coarser, single-dataset analyses.

The central contribution of this work is to show that the field's two hypotheses are not in conflict. The altered-state view — lowered SST and GABA-synthetic gene expression in interneurons that remain present — rests on two decades of RNAscope and microarray work, principally from the Lewis group (Dienel et al., 2022, 2023), whereas the cell-loss view — a reduced number of interneurons — is supported by an independent set of studies (Toker et al., 2017; Batiuk et al., 2022; Kiss et al., 2025). Our data indicate that both are correct but describe different phenomena in partly different cells. The reduction of *SST* transcript within Sst neurons is subclass-wide and cell-intrinsic (Figure 2), matching the altered-state literature; the abundance loss, by contrast, is not subclass-wide but restricted to a subset of Sst supertypes concentrated in the upper cortical layers, most robustly Sst_25 (Figure 3), while other Sst supertypes remain numerically stable. Granular taxonomy is what dissolves the apparent contradiction: at the subclass level a pan-Sst molecular down-shift and a subtype-restricted loss of cells are confounded, whereas at supertype resolution they cleanly separate.

That said, a deeper ambiguity remains that our data cannot fully resolve. Because SCZ alters the very transcripts used to assign cell identity, a supertype that appears depleted may reflect either bona fide loss of those neurons or a drastic erosion of their transcriptomic identity to the point that they can no longer be confidently classified — and both scenarios reduce the count of confidently labeled Sst_25 cells. Cross-sectional post-mortem sampling likewise cannot order the two changes, leaving open whether the molecular down-regulation precedes, follows, or drives the compositional loss. We mitigate the identity concern through cross-platform replication and depth- and marker-based validation of the spatial labels, but a definitive resolution will require orthogonal, expression-independent labels (protein, morphology, or lineage). We therefore frame the abundance findings as identifying which cells are most affected, and where, while remaining deliberately agnostic about loss versus identity erosion as the underlying mechanism.

Three independent observations converge to identify the vulnerable population as a specific, primate-specialized cell type. First, the depletion is laminar and molecularly coherent: the affected upper-layer Sst supertypes preferentially express calbindin (*CALB1*) `[TODO: quantify CALB1 specificity, vulnerable vs deeper Sst]`, and in primates these *CALB1*-expressing cells include the double-bouquet cell — a morphological specialization whose tight, vertically descending axon provides intracolumnar inhibition and is largely absent in rodents (Yañez et al., 2005; Raghanti et al., 2010; Lee, Dalley et al., 2023). This lineage connects our result to a classical neuropathological observation made more than twenty years ago — a reduced density of calbindin-immunoreactive interneurons in SCZ prefrontal cortex (Beasley et al., 2002; ~20% reduction, *P* ≈ 0.035 `[TODO: confirm from Beasley et al., 2002]`) — which our data now localize to specific, upper-layer Sst transcriptomic supertypes. Second, SCZ common-variant heritability is enriched in these same upper-layer Sst supertypes, and across supertypes the degree of genetic enrichment tracks the degree of compositional depletion (Figure 4); because this alignment is difficult to explain if the abundance changes were purely secondary to chronic illness or treatment, it nominates upper-layer Sst loss as at least partly upstream — though the correlation holds under PGC3 (Trubetskoy et al., 2022) but attenuates under the more recent Bigdeli et al. (2026) statistics (Figure 4; Methods). Third, patch-seq links this population to a specific intrinsic physiology: per-supertype expression of *HCN1*, which encodes a subunit of the hyperpolarization-activated current I_h, tracks the electrophysiological sag potential across Sst supertypes, and *HCN1* is additionally a candidate cell-type–specific SCZ risk gene for these cells `[TODO: HCN1-as-Sst_25-driver pending regeneration of the Sst_25 driver panel; the HCN1-expression-vs-sag correlation, ρ = 0.62, p = 0.010, is verified]`. Together, these observations raise the possibility that the I_h/HCN1 physiology characteristic of upper-layer Sst cells is mechanistically tied to their vulnerability, and they nominate this axis as a concrete, druggable target for follow-up.

Our second compositional finding — a coordinated increase in deep-layer L6b excitatory neurons, replicated spatially (Figure 3) — recovers one of the more reproducible neuropathological observations in schizophrenia. An excess of deep-layer and white-matter (interstitial) neurons has been reported repeatedly across cohorts and methods (reviewed in Kubo et al., 2019, 2020) and is most often interpreted as a neurodevelopmental signature: aberrant migration or failed developmental apoptosis of subplate-lineage neurons. Because L6b/subplate neurons are positioned to coordinate activity between cortex and higher-order thalamus (Zolnik et al., 2026), a relative excess of this population may reflect an early developmental perturbation retained into the adult cortex, complementary to the inhibitory loss rather than a compensatory response to it. Whether the inhibitory loss and the excitatory excess share a developmental origin, or arise independently, remains an open question.

The vulnerability we describe may extend beyond schizophrenia and beyond a single cortical region. The upper-layer Sst supertypes depleted in SCZ overlap substantially with those lost earliest in Alzheimer's disease (Gabitto et al., 2024) `[TODO: confirm the overlap count (talk states "five of eight") and list the intersection]`, and an independent RNAscope re-analysis extends the SST deficit from dorsolateral prefrontal cortex to the subgenual anterior cingulate while resolving a graded, diagnosis-specific pattern — largest in schizophrenia, intermediate in bipolar disorder, and absent in major depression (β = −0.60, −0.38 and +0.08, respectively; Supplementary Fig. S14). Taken together, the recurrence of the same primate-specialized, upper-layer inhibitory population across a neurodegenerative and a psychotic disorder, and the concentration of the deficit on the psychosis spectrum, suggest a shared axis of upper-layer cortical inhibitory fragility rather than strictly disease-specific pathology — with the important caveats that the time-courses and drivers almost certainly differ (early, progression-linked loss in AD; a likely earlier, more static profile in SCZ), that the overlap is at the level of which cells are affected rather than the mechanism, and that the histological effects are trend-level. More broadly, our results argue that pairing fine-grained cell-type taxonomies with cross-dataset meta-analysis is a generalizable strategy for resolving cell-type-specific neuropathology in psychiatric disorders — where individual post-mortem datasets are small and heterogeneous — and for nominating specific, physiologically defined cell populations, here the I_h/HCN1-high upper-layer Sst cell, for mechanistic and therapeutic study.

Limitations. Several caveats bound these conclusions. First, cell typing under transcriptional dysregulation: as noted above, we cannot fully exclude that depleted supertypes persist with altered expression, so the abundance results should be read as layer- and subtype-localized vulnerability rather than as proof of cell death. Second, meta-analysis is tuned for consistency: random-effects pooling privileges effects reproducible across datasets and down-weights large but dataset-specific effects, so the most consistent change is not guaranteed to be the most biologically consequential one. Third, the order of events is unresolved: cross-sectional post-mortem sampling cannot order molecular down-regulation and compositional change, nor fully separate primary disease processes from illness chronicity, medication, or agonal state, although we controlled for age, sex, and PMI and restricted to donors under 70 years to limit survivorship bias. Fourth, supertype resolvability in the spatial data is limited: the 300-gene Xenium panel carries few within-subclass markers for many Sst supertypes, so subclass-level spatial conclusions are firmer than individual supertype allocations, which should be treated as corroborative (Figure 3; Methods). Finally, our datasets are predominantly dorsolateral prefrontal cortex; reassuringly, the independent sgACC RNAscope re-analysis provides preliminary cross-region support, but establishing that the same cell-type–specific changes generalize broadly will require sampling additional cortical areas.

Conclusion. Fine-grained cell-type taxonomies and cross-dataset meta-analysis together resolve a decades-old debate over the nature of GABAergic pathology in schizophrenia: a subclass-wide reduction of *SST* mRNA coexists with a spatially and genetically specific depletion of upper-layer, *CALB1*-expressing, primate-specialized Sst neurons, set against an increase in deep-layer L6b excitatory neurons. The convergence of genetics and patch-seq physiology on this population nominates it — and its *HCN1*/I_h physiology — as a priority target for mechanistic dissection and therapeutic development.

## Methods

### Data sources and cohorts

We compiled seven post-mortem snRNA-seq datasets spanning multiple institutional
cohorts and brain banks: the Harvard Brain Tissue Resource Center at McLean Hospital
(McLean); the Mount Sinai NIH NeuroBioBank as profiled by Ruzicka et al. (MSSM 1) and
the PsychAD Consortium (MSSM 2); the NIMH-IRP Human Brain Collection Core (HBCC); the
brainSCOPE/PsychENCODE Resource multiome data (Multiome); and the open-access Batiuk
and Fröhlich studies. Studies were required to (i) profile prefrontal cortex tissue
via snRNA-seq on the 10x Genomics Chromium v3 or v3.1 chemistry, and (ii) classify
donors as SCZ cases or neurotypical controls. Included datasets predominantly sample
dorsolateral prefrontal cortex (DLPFC), with one dataset (Fröhlich) sampling
orbitofrontal cortex (OFC). To limit survivorship bias (cf. Kiss et al., 2025), only
cases and controls younger than 70 years at death were included. The meta-analysis
comprises 298 controls and 171 SCZ cases across the seven datasets. Per-dataset sex and
post-mortem interval (PMI) distributions are given in Supplementary Table T1.

| Dataset | Data type | CON | SCZ | Mean age (SD) CON | Mean age (SD) SCZ | Region | Source |
|---|---|---|---|---|---|---|---|
| Batiuk | snRNA-seq | 10 | 5 | 57.3 (5.4) | 61.4 (4.4) | DLPFC | Batiuk et al., Sci Adv 2022 |
| HBCC | snRNA-seq | 79 | 51 | 40.3 (13.0) | 48.9 (10.4) | DLPFC | Lee et al. (PsychAD), preprint 2024 |
| MSSM 2 | snRNA-seq | 140 | 45 | 55.8 (11.6) | 57.4 (9.3) | DLPFC | Lee et al. (PsychAD), preprint 2024 |
| McLean | snRNA-seq | 18 | 14 | 55.4 (9.8) | 50.8 (11.4) | DLPFC | Ruzicka et al., Science 2024 |
| MSSM 1 | snRNA-seq | 17 | 18 | 51.5 (14.4) | 59.8 (8.3) | DLPFC | Ruzicka et al., Science 2024 |
| Multiome | snRNA-seq | 5 | 6 | 39.2 (13.7) | 46.3 (13.0) | DLPFC | Emani et al., Science 2024 |
| Fröhlich | snRNA-seq | 29 | 32 | 54.2 (10.2) | 51.7 (9.5) | OFC | Fröhlich et al., Nat Neurosci 2024 |
| **Meta-analysis** | snRNA-seq | **298** | **171** | 51.0 (13.6) | 53.2 (10.6) | DLPFC/OFC | — |
| Xenium (LIBD) | Xenium | 12 | 12 | 46.6 (9.7) | 49.2 (8.1) | DLPFC | Kwon et al., preprint 2026 |
| Arbabi (sgACC) | RNAscope | 15 | 11 | 47.5 (11.4) | 43.8 (9.4) | sgACC | Arbabi et al., Mol Psychiatry 2025 |

To test whether snRNA-seq findings replicate in an orthogonal, spatially resolved
modality, we analyzed a Xenium dataset of post-mortem human DLPFC from the Lieber
Institute for Brain Development (LIBD) Human Brain and Tissue Repository, comprising
12 SCZ cases and 12 neurotypical controls (Kwon et al., preprint 2026); panel design,
cell segmentation, and cell-typing are described under *Xenium spatial
transcriptomics*.

All datasets were mapped to a single unified taxonomy defined by the Seattle
Alzheimer's Disease Brain Cell Atlas (SEA-AD), which defines transcriptomic
"supertypes" — a resolution finer than the conventional class/subclass hierarchy
(e.g., Sst → Sst_25) — from neurotypical middle temporal gyrus (MTG) reference donors
(Gabitto et al., 2024); we mapped to the supertype level (137 SEA-AD MTG supertypes).
To attach functional, morphological, and laminar properties to these transcriptomic
types, we used BRAIN Initiative Cell Census Network (BICCN) human patch-seq recordings
aligned to the same taxonomy (401 donors, 2,602 cells; Lee, Dalley et al., 2023)
`[TODO: confirm BICCN patch-seq donor/cell counts against source]`.

### snRNA-seq processing and cell-type harmonization

For each of the seven datasets we obtained the author-provided, quality-controlled
filtered count matrices. Because each study's native cell-type taxonomy differed, we
harmonized all donors to a single reference by Seurat reference-based label transfer,
projecting every nucleus onto the SEA-AD MTG supertype taxonomy (137 supertypes,
defined from 5 neurotypical MTG reference donors; Gabitto et al., 2024) `[TODO:
confirm Seurat label-transfer parameters against the upstream pipeline]`. This
assigned each nucleus a class, subclass (e.g., Sst) and supertype (e.g., Sst_25) label
used consistently across all datasets and downstream analyses.

### Compositional analysis (crumblr)

Cell-type composition changes were modelled with crumblr, which applies a
centred-log-ratio transform to per-donor cell-type counts and fits precision-weighted
linear models (donors contributing more nuclei are weighted more heavily). Neuronal
and non-neuronal compartments were analysed separately to account for technical
sampling variation (e.g., grey-versus-white-matter content) and so that proportions
reflect within-compartment composition. Within each dataset we regressed cell-type
composition on diagnosis with demographic covariates:

> composition ~ diagnosis + age at death + sex + PMI

PMI was unavailable for the Multiome dataset and was omitted from that dataset's model.
Per-dataset estimates and standard errors were combined across the seven datasets by
random-effects meta-analysis (metafor::rma, restricted maximum likelihood), producing
weighted, heterogeneity-adjusted pooled estimates; nominal p-values were adjusted by
the Benjamini–Hochberg (BH) FDR procedure.

### Cell-type–specific differential expression

To quantify how SCZ alters gene expression within cell types, we tested differential
expression (DE) separately in each of the seven snRNA-seq datasets at SEA-AD subclass
resolution. For each dataset, counts were aggregated into per-donor, per-cell-type
pseudobulk profiles; to ensure stable estimates we retained only donors contributing
≥ 500 cells to a given cell type and, within each cell type, tested only genes
expressed in ≥ 80% of retained donors. Pseudobulk profiles were modelled with
limma-voom and edgeR, with SCZ diagnosis as the predictor of interest and age at
death, sex, and PMI as covariates (PMI omitted for the Multiome dataset). Per-dataset
log₂ fold changes and standard errors were combined for each gene × cell-type pair by
random-effects meta-analysis (metafor::rma, REML), and nominal p-values were
BH-FDR-adjusted across all tested gene × cell-type pairs. This 7-dataset meta-analytic
table (23 subclasses, 16,368 genes) is the canonical DE result; Xenium was never
pooled into it, preserving a strict snRNA-seq discovery / Xenium replication
separation. For display, per-gene forest diamonds (Fig. 2b,f) recompute the pooled
estimate with the DerSimonian–Laird estimator, with significance annotations from the
canonical REML FDR.

For spatial replication, the Xenium cohort (Kwon et al. 2026; GSE307404; 24 DLPFC
sections, 12/12; SEA-AD label transfer) was summed into per-section, per-subclass
pseudobulk and tested with an edgeR quasi-likelihood F-test adjusting for sex, age at
death, and PMI — the same covariates as the snRNA-seq model — retaining, within each
cell type, genes detected in ≥ 80% of the donors contributing a pseudobulk. Results
were FDR-corrected within cell type across the 300-gene panel; because this panel-wide
FDR is computed over far fewer genes, it is not directly comparable to the genome-wide
snRNA-seq FDR. Two snRNA-seq conventions were deliberately not transferred, because
Xenium is a single 24-section cohort rather than seven pooled datasets: the ≥ 500
cells-per-donor threshold (which would remove 8 of 23 subclasses, including L6b, for
which only 1 of 24 donors reaches 500 cells; we required ≥ 10 cells per donor per cell
type), and BH correction across all gene × cell-type pairs (unattainable on a 300-gene
panel, where the smallest observed *P* = 3.6 × 10⁻⁵ across 6,792 pairs cannot reach
FDR < 0.10).
Cross-platform concordance intersected the meta-DE gene × subclass pairs (meta
FDR < 0.10) with the pairs testable in Xenium (166 pairs); we correlated the two
platforms' log₂ fold changes (Pearson) and assessed directional agreement with a
binomial sign test (120/166 concordant; two-sided *P* = 8.4 × 10⁻⁹). DE burden per subclass
(up/down counts at FDR < 0.10 and 0.05) was correlated with each subclass's mean
per-donor Xenium proportion (Spearman) to test whether burden reflects statistical
power rather than selective vulnerability. To confirm the SST and PVALB reductions
were not normalisation artefacts, per-cell marker counts within the Xenium Sst and
Pvalb populations were modelled with negative-binomial mixed models (donor random
intercept) under four normalisations (Supplementary Fig. S7). These per-cell models showed the Sst SST reduction was robust (0.70× raw, P = 0.002; 0.76× library-normalised, P = 0.011), whereas the Pvalb PVALB reduction was modest and borderline (0.86× raw, P = 0.065; 0.87× library-normalised, P = 0.060, n.s.).

### Xenium spatial transcriptomics

We reanalysed the Xenium DLPFC dataset of Kwon et al. (2026; 24 sections, 12 SCZ / 12 control; 1.34 million cells; 300-gene panel). Because 300 genes are too few for de novo clustering, cells were assigned to SEA-AD subclasses and supertypes by a self-referencing two-stage correlation classifier and to a continuous cortical depth (0 = pia, 1 = white matter) by a neighbourhood-composition model; both were validated against the independent SEA-AD MERFISH atlas (subclass proportions *r* = 0.85, depth *r* = 0.96; 84.9% classifier accuracy) and the panel's supertype-resolution limits characterized by a panel-vs-transcriptome F1 benchmark (Supplementary Methods SM1; Supplementary Figs S2–S4). Disease analyses used one shared cell definition (cortical, quality-passing cells): composition via stratified crumblr (neuronal and non-neuronal compartments; covarying age, sex, PMI), cell densities (cells/mm²) as a constraint-free complement (Supplementary Fig. S9), and DE via pseudobulk edgeR (~ diagnosis + sex + age + PMI), mirroring the snRNA-seq meta-analysis. For the L6b aggregation, four sections with < 3% L6 cells (Br2039, Br5973, Br2719, Br5314) were excluded. The full cross-platform concordance matrix (composition and density; subclass and supertype) is Supplementary Table T7.

### SCZ GWAS enrichment and patch-seq integration

GWAS summary statistics. SCZ heritability enrichment used PGC3-SCZ summary statistics
(Trubetskoy et al. 2022; European-ancestry autosomal meta-analysis VCF v3,
GRCh37/hg19) and the same study's FINEMAP posterior inclusion probabilities (PIPs) and
95% credible sets `[TODO: state exact PGC3 EUR effective N — README lists 76,755 cases
/ 243,649 controls; docs list EUR-only ~53K/77K]`. A larger Bigdeli et al. 2026 SCZ
meta-analysis (MAGMA v1.10, dense g1000_eur LD panel) was used for a robustness re-run
(below).

Combined Franken-503 taxonomy. Enrichment was tested across a combined 503-type taxonomy that merges the 137 SEA-AD MTG supertypes with 366 non-redundant Siletti et al. 2023 whole-brain clusters (reciprocal-best-hit deduplication of the 461 Siletti clusters; Supplementary Methods SM2). The full 503-type enrichment table is Supplementary Table T4 (SEA-AD-137 standalone, Supplementary Table T4b; Siletti 461-cluster, Supplementary Table T5).

Per-cell-type enrichment. SCZ heritability enrichment per cell type was tested by MAGMA-style gene-property regression (de Leeuw et al. 2015) of the gene-level SCZ Z-statistic on cell-type expression specificity, with MHC exclusion and Bonferroni/BH-FDR correction (Supplementary Methods SM3) `[TODO: document upstream MAGMA gene-analysis parameters (gene window, LD panel, SNP-to-gene aggregation) for the PGC3 .genes.out]`. In the SEA-AD-only 137-type test, 44 types were FDR-significant (17 Bonferroni; 41 GABAergic, 3 glutamatergic); in the Franken-503 test, 120 (54 Bonferroni). Panel 4a plots −log₁₀*p* for the 137 SEA-AD supertypes using their Franken-503 values, following the cell-type Manhattan convention of Duncan et al. 2025; independently enriched (conditionally significant) types are shown in Supplementary Fig. S11.

SST depth stratification and GWAS–composition convergence. Continuous cortical depth
for each of the 18 SST supertypes (0 = pia) was the per-supertype average of SEA-AD
MERFISH and Xenium median depths; enrichment vs depth was tested by Spearman
correlation (*r* = −0.50, *p* = 0.034, n = 18), with upper-layer types (depth < 0.35,
n = 7) at mean −log₁₀*p* = 6.6 versus 3.7 for deep types (n = 9). Franken-503 SST
enrichment was correlated with case–control compositional change from the 7-dataset
snRNA-seq meta-analysis (Endresz et al., in prep; crumblr β and FDR matched to 109
SEA-AD supertypes): within SST (n = 18), |composition β| vs GWAS −log₁₀*p* Spearman
*r* = 0.645 (*p* = 3.85 × 10⁻³) and signed composition β vs enrichment *r* = −0.494
(*p* = 0.037), with no such relationship across all 109 types (*r* = −0.10). In panels
4b,e, supertypes depleted in the composition meta-analysis (FDR < 0.20, β < 0) are
bold-outlined.

Gene drivers and HCN1 locus. For the featured supertype Sst_25, candidate driver
genes were defined as genes jointly in the top decile of Sst_25 specificity and passing
MAGMA gene-level FDR < 0.05 `[TODO: Sst_25 driver-gene count and HCN1's inclusion
pending regeneration of the Sst_25 driver table; only the Sst_2 table (269 genes) is
committed and HCN1 is absent from it]`. Around *HCN1*, SCZ GWAS variants were extracted
in an asymmetric window and lifted hg19→hg38 for display against the *HCN1* track; the
PGC3 FINEMAP 95% credible set (lead rs10035564) comprised 9 variants (cumulative
PIP = 0.95; lead PIP = 0.525; 4 intronic to *HCN1*), and the credible set's
PIP-weighted snATAC accessibility was higher in SST (0.045) than astrocytes (0.006).

Patch-seq electrophysiology. A harmonized human cortical patch-seq dataset (1,155 cells across 221 donors; Lee, Dalley et al. 2023; Chartrand et al. 2023) was assigned to SEA-AD supertypes by scANVI and Harmony-corrected kNN, and sag ratio extracted from NWB sweeps with the Allen IPFX library (Supplementary Methods SM4; scANVI-vs-kNN agreement, Supplementary Table T6). The 16 SST supertypes with ≥ 1 assigned cell (150 cells) are shown in panel 4e; per-supertype *HCN1* specificity correlated with mean sag (Spearman ρ = 0.62, *p* = 0.010) `[TODO: confirm Fig. 4e axis — committed value is HCN1 specificity-vs-sag; methods draft describes expression-vs-sag]`. Two exemplar cells span the range (Sst_25 specimen 819770858, L2; Sst_5 specimen 758996755, L4; supertype mean sag 0.456 vs 0.167) `[TODO: single-cell exemplar sag values (0.564/0.00196) and the Sst_25 donor ID are methods-draft only — verify from the patch-seq source]`; biocytin-recovered SWC morphologies were drawn soma-aligned to pia (100 µm scale bar).

Robustness to GWAS version. The GABAergic/SST enrichment headline reproduced under the larger Bigdeli et al. 2026 SCZ meta-analysis (per-cell-type enrichment Spearman ρ = 0.95 vs PGC3), but the two SST sub-claims (depth gradient; genetics↔depletion convergence) attenuated to non-significance; because Bigdeli contains PGC3 this is a dilution robustness check, not an independent replication (Supplementary Methods SM3; Supplementary Fig. S13). The PGC3-vs-Bigdeli choice is unresolved (see Outstanding items).

### RNAscope re-analysis (sgACC)

To test whether the upper-layer SST deficit is detectable at the cellular level in an
independent region and modality, we re-analysed RNAscope fluorescent in-situ
hybridization counts of SST interneurons in subgenual anterior cingulate cortex (sgACC)
from Arbabi et al. (2025), who quantified SST and VIP densities across 20 stereological
counting frames per section assigned to upper (L2/3) or deep (L5/6) cortical layers.
Because VIP interneurons are concentrated superficially, we used their L2/3-over-L5/6
enrichment as an internal laminar-quality control and excluded subjects lacking it
(one-sided *t*-test *P* ≥ 0.05), retaining 54 of 68 subjects (15 control, 11
schizophrenia, 13 major depression, 15 bipolar). SST density was modelled with a subject
random-intercept mixed model (SST ~ diagnosis + age + sex + PMI + (1|subject)) fit on all
per-frame observations, both aggregate and layer-stratified. Representative RNAscope
images were lipofuscin-suppressed for display only; all quantification used the original
SlideBook cell counts (Supplementary Methods SM5).

### Statistics and code availability

Multiple-testing correction used the Benjamini–Hochberg FDR procedure throughout, with
both nominal and adjusted values reported; nominally significant results were treated
as leads rather than confirmed findings. Random-effects meta-analyses (composition and
DE) used metafor::rma with REML; per-gene forest diamonds were recomputed with the
DerSimonian–Laird estimator. Correlations are Pearson or Spearman as stated, with n
reported alongside each. Xenium panel-wide FDR is computed over the 300-gene panel and
is not directly comparable to the genome-wide snRNA-seq FDR. Data extraction used
Python (pandas, numpy, scipy, statsmodels, scanpy/anndata, pynwb); figures were
rendered in R (ggplot2, cowplot, ggrepel, scales). Cell-type colours follow the Allen
SEA-AD supertype palette.

Analysis code spans the transcriptomic DE/composition meta-analysis, the spatial
(Xenium) pipeline, and the genetics pipeline (Franken taxonomy, MAGMA enrichment, gene
drivers, fine-mapping, and patch-seq integration; source repository
github.com/stripathy/scz_cell_type_enrichment). Patch-seq NWB/SWC data are on DANDI
(dandisets 000636, 000630, 000228, 000337) and the Brain Image Library; PGC3 SCZ
summary statistics and FINEMAP results are from the PGC (https://pgc.unc.edu/); SEA-AD
and Siletti references are from the Allen Brain Cell Atlas.

## Figure legends

### Figure 1 | Study design, cross-dataset snRNA-seq harmonization, and Xenium spatial annotation

**(a)** Study-design and cross-dataset harmonization schematic. **(b–d)** UMAP of the
integrated snRNA-seq nuclei coloured by **(b)** dataset, **(c)** diagnosis and **(d)**
subclass, showing that nuclei group by cell type rather than by dataset or diagnosis.
**(e)** Representative Xenium DLPFC sections (one control, one SCZ), each cell plotted
at its spatial position and coloured by SEA-AD subclass (24-colour taxonomy;
correlation-classifier `corr_subclass` labels, cortical `corr_qc_pass` cells). **(f)**
The same sections coloured by inferred cortical layer (L1–L6, white matter, vascular),
assigned by binning MERFISH-trained neighbourhood-composition depth predictions and
spatially smoothing; dashed lines mark layer boundaries. Scale bar, `[TODO: scale-bar
length]`. Cross-platform validation (subclass proportion *r* = 0.85 log₁₀; subclass
depth *r* = 0.96 vs SEA-AD MERFISH) is shown in Supplementary Fig. S2.

### Figure 2 | Cross-platform schizophrenia differential expression for the two canonical interneuron markers, in the context of transcriptome-wide DE burden

**(a–d)** SST in Sst cells. (a) Volcano of SCZ log₂ fold change versus −log₁₀ *P* (raw
meta-analytic *p*) for genes tested in Sst cells; points coloured by direction × FDR
tier (up orange, down blue; dark FDR < 0.05, light FDR < 0.10; NS grey); dashed line,
the *p* at the FDR < 0.10 boundary; selected genes labelled (down: *NAT16*, *AFG3L2*,
*DRD3*, *SST*; up: *SMAD1*, *SLC9A9*, *STAC*, *KCTD4*); data-driven axes. (b) Forest
plot for *SST*: seven frontal-cortex snRNA-seq datasets (grey squares), pooled
random-effects meta-analytic estimate (black diamond; log₂FC = −0.46, FDR = 0.049) and
the independent Xenium spatial estimate (green triangle; log₂FC = −0.32, *P* = 0.052);
whiskers, 95% CI. (c) Per-donor library-normalised expression (CP1K) of *SST* in Sst
cells, control versus SCZ (12 vs 12); *p*, Xenium edgeR quasi-likelihood test
(*P* = 0.052). (d) Representative Xenium Sst cells at the pooled group-median transcript
density per diagnosis (control Br6432, 35 molecules; SCZ Br5973, 28); grey outline,
cell boundary; dashed, nucleus; red dots, *SST* molecules; scale bar, 5 µm. **(e–h)**
PVALB in Pvalb cells, same four views. (e) Volcano (down: *ANXA2*, *NAT16*, *VGF*,
*CIRBP*; up: *SMAD1*, *SCN3A*, *TCAF2*, *FGF10*; *PVALB* as the n.s. marker reference).
(f) Forest for *PVALB* (meta log₂FC = −0.06, FDR = 0.86, n.s.; Xenium −0.22,
*P* = 0.044). (g) CP1K of *PVALB* in Pvalb cells (edgeR *P* = 0.044). (h) Exemplar
Pvalb cells (control Br6432, 10 molecules; SCZ Br5973, 9). **(i)** Up- (orange) and
down-regulated (blue) DE-gene counts per subclass (22 of 23 tested with ≥ 1 DE gene),
FDR < 0.10 (light) with the FDR < 0.05 subset overlaid (dark); counts annotated (6,903
associations at FDR < 0.10; 3,586 at FDR < 0.05). Inset: DE-gene count versus mean
per-donor cell-type proportion (Xenium, log₁₀), Spearman ρ = 0.83 (*P* = 2.1 × 10⁻⁶,
n = 22). **(j)** snRNA-seq meta-analytic versus Xenium log₂ fold change for all gene ×
cell-type pairs with meta FDR < 0.10 (n = 166); coloured by cell class, sized by meta
FDR (larger, FDR < 0.05, n = 109; smaller, 0.05–0.10, n = 57); ten pairs labelled;
dashed line, identity. Pearson *r* = 0.73; 72% sign-concordant (120/166). Significance:
\*FDR < 0.10, \*\*FDR < 0.05, \*\*\*FDR < 0.01; •, nominal *P* < 0.05; n.s., not
significant. Xenium FDR is panel-wide, not comparable to genome-wide snRNA-seq FDR.

### Figure 3 | Cell-type composition in schizophrenia: snRNA-seq meta-analysis and Xenium spatial replication

**(a)** Meta-analysed SCZ compositional effect (crumblr CLR β ± SE) per supertype
across seven snRNA-seq datasets (298 control / 171 SCZ), covarying age, sex, PMI; bars
coloured by significance (FDR < 0.01 \*\*\*, < 0.05 \*\*, < 0.10 \*, 0.10–0.20 ·).
Vulnerable Sst supertypes (Sst_2, Sst_22, Sst_25, Sst_3) reduced and L6b_1/L6b_4
increased. **(b)** Per-donor Sst_25 proportion by dataset (control blue, SCZ red);
unadjusted *P*. **(c)** Forest plot of the Sst_25 effect: per-dataset (black) and
random-effects pooled (red). **(d)** Representative Xenium sections (one control, one
SCZ, closest to group-median Sst_25 proportion), cells coloured by Sst supertype
vulnerability (red = vulnerable, blue = non-vulnerable Sst, grey = other),
cortical-layer boundaries dashed; right, per-layer density. **(e)** Per-donor Sst_25
proportion in Xenium (n = 12/12); unadjusted *P* from crumblr covarying age, sex.
**(f)** snRNA-seq meta-analytic β (x) versus Xenium logFC (y) for neuronal supertypes
(GABAergic purple, glutamatergic green); dashed lines, origin; Pearson *r* = 0.50
(n = 106 neuronal supertypes; all 120 supertypes *r* = 0.51). Deep-layer excitatory
supertypes upper-right, Sst supertypes lower-left on both platforms. snRNA-seq β and
Xenium logFC are stratified crumblr estimates (neuronal/non-neuronal analysed
separately). Subclass concordance: neuronal *r* = 0.70 (n = 17); density-based neuronal
supertype *r* = 0.55.

### Figure 4 | SCZ genetic risk converges on vulnerable upper-layer SST interneurons and an HCN1-linked physiological phenotype

**(a)** SCZ heritability enrichment (MAGMA gene-property regression of PGC3 gene-level
Z on cell-type specificity) across the Franken-503 taxonomy (137 SEA-AD MTG supertypes
+ 366 Siletti clusters); the 137 SEA-AD supertypes shown as −log₁₀*p* grouped by
subclass. Overwhelmingly GABAergic (68 of the 137 shown supertypes FDR-significant, 56 GABAergic); Pvalb_3 the strongest cortical signal (β = 47.2, *p* = 3.9 × 10⁻¹²). Dashed lines: Bonferroni, FDR
thresholds. **(b)** GWAS enrichment vs cortical depth for 18 SST supertypes (0 = pia);
upper-layer types carry stronger signal (Spearman *r* = −0.50, *p* = 0.034).
Composition-depleted types (FDR < 0.20, β < 0) bold-outlined; most-enriched SST types
are most depleted (|composition β| vs enrichment, *r* = 0.645, *p* = 3.85 × 10⁻³).
**(c)** Gene-driver scatter for Sst_25 (specificity vs SCZ −log₁₀*p*); *HCN1*
highlighted `[TODO: Sst_25 driver list pending regeneration]`. **(d)** *HCN1* locus
zoom: PGC3 FINEMAP 95% credible set (lead rs10035564, PIP = 0.525; 9 variants,
cumulative PIP = 0.95) coloured by PIP over the *HCN1* gene track. **(e)** Per-supertype
*HCN1* expression vs patch-seq mean sag across 16 SST supertypes (Spearman ρ = 0.62,
*p* = 0.010); point size ∝ n cells; depleted types outlined. **(f)** Morphological
reconstructions: high-sag upper-layer Sst_25 (specimen 819770858, L2) and low-sag deep
Sst_5 (specimen 758996755, L4), soma-aligned to pia; 100 µm scale bar. **(g)** Matched
hyperpolarizing voltage responses; the Sst_25 cell shows pronounced I_h/HCN-dependent
sag and the Sst_5 cell essentially none `[TODO: single-cell sag values 0.564/0.00196
are methods-draft only; supertype means 0.456/0.167 verified]`. SEA-AD supertype
palette throughout.

## Outstanding items

Every inline `[TODO]` left in the manuscript, grouped by section.

**Abstract**
- Sst_25 effect size and 95% CI from the composition meta-analysis (hero number).

**Results — Figure 1**
- Confirm total nuclei analysed (~2.3 million) from Nicole's pipeline.

**Results — Figure 4**
- Sst_25 driver-gene count and HCN1's driver status — pending regeneration of the Sst_25 driver panel (only the Sst_2 table, 269 genes, is committed; HCN1 absent).
- Single-cell exemplar sag values (0.564 / 0.00196) are methods-draft only — verify from the patch-seq source or report the verified supertype means.
- Finalize GWAS version — PGC3 (primary) vs Bigdeli 2026 (SST sub-claims attenuate).

**Discussion** *(AD / CALB1 / double-bouquet material moved here from the former Figure 5)*
- HCN1-as-Sst_25-driver pending regeneration of the Sst_25 driver panel (the HCN1-expression-vs-sag correlation, ρ = 0.62, p = 0.010, is verified).
- Quantify *CALB1* specificity/expression in vulnerable vs deeper Sst supertypes (SEA-AD reference matrix) — used in the "Why upper-layer Sst cells?" paragraph.
- Confirm the Beasley et al., 2002 effect size and *P* (~20% calbindin+ reduction, *P* ≈ 0.035).
- Confirm the SCZ↔early-AD vulnerable-Sst overlap count ("five of eight") and list the intersection (Gabitto 2024 + composition meta-analysis).

**Methods — Data sources and cohorts**
- Confirm BICCN patch-seq donor/cell counts (401 donors / 2,602 cells) against source.
- Confirm the Arbabi (sgACC RNAscope) cohort ages (47.5/43.8) against the histology demographics table (`histology/data/pTable`).

**Methods — snRNA-seq processing and cell-type harmonization**
- Confirm Seurat label-transfer parameters against the upstream pipeline.

**Methods — SCZ GWAS enrichment and patch-seq integration**
- State the exact PGC3 EUR effective N (README 76,755 / 243,649 vs docs EUR-only ~53K/77K).
- Document upstream MAGMA gene-analysis parameters (gene window, LD panel, SNP-to-gene aggregation).
- Sst_25 driver-gene count and HCN1 inclusion — pending regeneration of the Sst_25 driver table.
- Reconcile Fig. 4e axis — committed value is HCN1 specificity-vs-sag (ρ = 0.62); methods draft describes expression-vs-sag.
- Single-cell exemplar sag values (0.564 / 0.00196) and the Sst_25 donor ID — methods-draft only, verify from the patch-seq source.

**Figure legends**
- Figure 1: scale-bar length; Supplementary figure number for the cell-typing validation.
- Figure 4c: Sst_25 driver list pending regeneration.
- Figure 4g: single-cell sag values (0.564 / 0.00196) unverified; supertype means (0.456 / 0.167) verified.

**Hand-offs to collaborators**
- **Send Nicole the Figure 2 inputs and the Xenium DE results for her updated Figure 2.**
  The Xenium DE was rerun on 2026-07-31 and the numbers she is working from are
  now out of date, so this is a send-before-she-redraws item, not a courtesy copy.
  - Files: the six snapshots in `transcriptomic/data/figure_inputs/` (plus
    `MANIFEST.tsv`, which records exactly which upstream file each came from),
    and both Xenium DE tables — `de_results_subclass.csv` and
    `de_results_supertype.csv` in `SCZ_Xenium/output/de/` (also reachable via
    `spatial/output/de/`).
  - What changed: rebuilt on the reinstated 2026-04-01 Xenium object, with PMI
    added to the model (`~ diagnosis + sex + age + PMI`, matching the snRNA-seq
    methods) and the gene filter switched to "detected in ≥ 80% of retained
    donors". Subclass hits went 94 → 208 at FDR < 0.10, almost entirely from the
    PMI covariate rather than the data change.
  - Downstream consequence to flag for her: cross-platform concordance moved from
    76% (126/166) to 72% (120/166), sign-test *P* 7.2 × 10⁻¹² → 8.4 × 10⁻⁹, OLS
    slope 0.79 → 0.81; Pearson *r* stays 0.73. Any panel or number she carries
    over from the previous DE needs regenerating, not just re-labelling.
