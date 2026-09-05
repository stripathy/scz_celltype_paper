# Supplementary Information — DRAFT (Option B preview)

> ## ⚠️ SUPERSEDED — do not cite or edit as source of truth (2026-08-06)
>
> The live Supplementary Methods now live in the Google Doc
> `Draft_Manuscript_SCZ_project_1` (`1cO5ZStbp9b3cb6vb9s6H2YfbprYGikXBvFHQR7gROjo`),
> which is substantially rewritten and **more accurate** than this file. Known
> divergences where **the Google Doc is correct and this file is wrong**:
>
> - **SM1 depth model.** This file describes a K = 100 / 6-lowest-CPS-donor model at
>   held-out R² = 0.907. The deployed pipeline is **K = 50 with the 3 smallest donors
>   held out** (train R² = 0.93, test R² = 0.89, MAE 0.050/0.069, r 0.96/0.95) — see
>   `spatial/code/pipeline/04_run_depth_prediction.py:112` and
>   `spatial/code/modules/depth_model.py:190`. The low-CPS/K = 100 variant is an
>   *uncommitted exploration* (`validate_depth_lowcps_cv.py`, `compare_depth_K.py`),
>   never wired into the pipeline.
> - SM1 is ~1,240 words here vs ~2,300 in the Doc; the Doc adds the full QC section,
>   the L6b margin threshold, marker-detection validation, and Sst_25 GAD1/GAD2 evidence.
>
> Retained here only as a record of material the Doc **dropped** and that may be worth
> reinstating: the anti-circularity package (GroupKFold R² = 0.889; subclass-only
> baseline R² = 0.32; marker-anchor validation), the 84.9% MERFISH classifier accuracy
> (still cited in the Doc's main Methods but no longer supported in its SM1), and the
> full Lamp5 Lhx6 argument.

> **Supplementary Information (Option B) — main-Methods rewiring APPLIED 2026-07-21.**
> The detailed Xenium and GWAS procedures now live here as Supplementary Methods SM1–SM4 (plus SM5 for RNAscope display processing);
> the main Methods in `manuscript.md` point to them (see "Main-Methods rewiring" at the
> bottom for the exact pointer text now in the manuscript).
>
> Markers: **[BUILD]** = figure/table must be generated · **[RETRIEVE]** = lives in a
> sibling repo · **[ASSEMBLE]** = compile from existing data · **[RECONCILE]** = a
> number to lock first (see `README_STATUS.md`). Figure/table numbers are provisional.

## Contents
- **Supplementary Methods** SM1–SM5
- **Supplementary Figures** S1–S14 (+ Option-C spares)
- **Supplementary Tables** T1–T7
- **Main-Methods rewiring (preview)** — what shrinks in `manuscript.md`

---

# Supplementary Methods

## SM1 — Xenium cell typing and spatial inference

*The novel Xenium label-transfer and depth-inference pipeline, relocated in full from the
main Methods. Supports Figure 1e–f and Figure 3; anchors Supplementary Figs S2–S4.*

**Cohort.** We reanalysed the Xenium dataset of Kwon et al. (2026; GEO GSE307404): 24
post-mortem human DLPFC sections from the LIBD Human Brain and Tissue Repository (12 SCZ,
12 control; matched for age), profiled with a 300-gene panel and yielding 1,339,151
segmented cells. All cell-type, depth, and layer inferences below are analytical
additions to this primary dataset.

**Cell-type classification.** Because 300 genes (~1% of the transcriptome) are too few
for de novo clustering to resolve the 24 subclasses and 137 supertypes of the SEA-AD MTG
reference (Gabitto et al., 2024), we transferred labels with a self-referencing two-stage
correlation classifier operating entirely within the Xenium feature space. Per-type
centroids were built from the top-100 highest-confidence Allen MapMyCells (HANN) exemplars
per type (counts-per-10k, log1p); stage 1 assigned each cell to a subclass by z-scored
Pearson correlation to subclass centroids (recording the best-minus-second-best margin),
and stage 2 correlated the cell only against supertype centroids within its assigned
subclass, preventing biologically impossible cross-subclass supertype calls. Spatial QC
(negative-probe/codeword and count filters) plus a per-sample 5th-percentile
classifier-margin filter and a spatial-doublet exclusion (cells co-expressing incompatible
GABAergic/glutamatergic marker sets) retained 1,221,519 cells; restricting to cortical
tissue left 742,103 cells across the 24 donors as the disease-analysis set. QC was applied
as boolean flags rather than by deleting cells. A self-referencing centroid approach was
chosen because Harmony + kNN transfer inflated Sst to 12.1% of cells (versus ~2.5%
expected) and agreed with the correlation classifier on only 69% of subclass calls
(Supplementary Fig. S2).

**Cortical depth model.** Depth was predicted from local neighbourhood cell-type
composition (the subclass fractions among each cell's K nearest same-section neighbours,
concatenated with a one-hot encoding of the cell's own subclass) rather than from its
sparse expression, using a GradientBoostingRegressor (n_estimators = 300, max_depth = 5,
learning_rate = 0.1). To keep AD neurodegeneration from confounding platform calibration,
the model was trained on SEA-AD MERFISH cells from the lowest-pathology donors (ranked by
Continuous Pseudo-progression Score, CPS) carrying manual "normalized depth from pia"
annotations (0 = pia, 1 = white matter). The deployed pipeline trains on the 6 lowest-CPS
donors at K = 100 with a post-hoc k = 30 within-section spatial smoothing; on 3 held-out
low-CPS donors it reached R² = 0.907 **[RECONCILE: 0.907 deployed vs 0.889 committed
K=50/24-donor variant]** (MAE = 0.059, Pearson r = 0.957). Two circularities were addressed
explicitly: a donor-grouped GroupKFold CV predicted all 27 MERFISH donors out-of-fold at
R² = 0.889 (n = 368,795 cells), and a subclass-identity-only baseline reached only
R² = 0.32 versus R² = 0.89 for the full neighbourhood model, showing the depth signal
exceeds the label (Supplementary Fig. S3). As an independent, MERFISH-free anchor,
predicted Xenium depth ordered eight canonical panel layer markers — never used as depth
targets — into their correct laminar sequence (upper-layer CUX2, LAMP5, CALB1 declining;
deep NR4A2, NXPH4 rising; Supplementary Fig. S3). Continuous depth was binned into laminae
(L1 < 0.12, L2/3 0.12–0.47, L4 0.47–0.54, L5 0.54–0.71, L6 0.71–0.93, WM > 0.93) from
SEA-AD MERFISH excitatory-marker crossovers.

**Spatial-domain segmentation.** BANKSY spatial clustering (Singhal et al., 2024; λ = 0.8,
Leiden resolution 0.3, k_geom = 15, 20 PCs) grouped cells into spatially coherent domains
classed as Cortical, Vascular, or White Matter by composition and mean depth; raw
depth-bin layers were then refined by within-domain majority-vote smoothing. The
segmentation partitions each section into its cortical laminae (L1, L2/3, L4, L5, L6) plus
white-matter and vascular domains.

**Cross-platform validation.** Against the independent SEA-AD MERFISH atlas (341,595
cortical cells, 27 neurotypical donors, 180-gene panel, manual layer annotations), Xenium
per-donor subclass proportions agreed at Pearson r = 0.85 (log₁₀; ρ = 0.88, n = 23) and
median cortical depth at r = 0.96 (subclass, n = 23); supertype proportions were markedly
weaker (r = 0.43, n = 131), reflecting the panel's limited within-subclass marker content.
Applying the correlation classifier directly to MERFISH (which carries its own ground-truth
labels) gave 84.9% subclass accuracy, with expected weak points at rare or marker-poor
types (Sst 69%, Lamp5 55%; Supplementary Fig. S2).

**Panel resolvability.** To separate "the panel lacks the genes" from "these types are
intrinsically hard to distinguish," we benchmarked classification F1 from the 300 panel
genes versus the full transcriptome, holding classifier (nearest-centroid Pearson) and
validation (leave-one-donor-out CV on the SEA-AD neurotypical snRNA-seq reference) fixed
and varying only the gene set. Subclass identity was panel-sufficient (median F1 0.96
panel vs 0.99 transcriptome; 15/23 subclasses at ceiling), whereas Sst supertype identity
was only partly resolved (median F1 0.62 vs 0.88) — the high transcriptome ceiling shows
the subtypes are genuinely separable, so the gap is panel coverage, not biological
continuity. The 300-gene panel carries 0–1 within-subclass discriminating markers for most
Sst supertypes: Sst_25 was the best-resolved (F1 0.81), whereas Sst_22 (F1 0.45) and
Sst_20 (F1 0.44) were the worst-resolved — accounting for the non-replication of the
Sst_22 and Sst_20 composition effects in Xenium. Per-supertype Sst allocation in Xenium is
therefore corroborative rather than definitive (Supplementary Fig. S4).

**Cell types too rare to analyse in Xenium.** We required a cell type to be present in at
least half of the donors, which excluded one subclass and four supertypes from our downstream
analyses. The subclass, Lamp5 Lhx6, is reported in our cell-type assignments but was recovered
too sparsely in cortex to support per-donor estimates: of the 2,598 cells assigned to it that
passed quality control, only 18 lay inside a cortical spatial domain whereas the majority,
2,549 (98.1%), lay in white matter annotated domains, and those remaining cortical cells were
present in only 11 of 24 donors. Its single supertype, Lamp5_Lhx6_1, is excluded for the same
reason. Of the 109 neuronal supertypes carried through the snRNA-seq meta-analysis this leaves
106 (Fig. 3e), the others being Pvalb_3, present in 2 donors, and Pvalb_14, which received no
Xenium cells at all. Both Pvalb supertypes are poorly resolved by the 300-gene panel (F1 0.42
and 0.49, versus 0.69 and 0.82 using the full transcriptome), so we expect their cells were
absorbed into neighbouring Pvalb supertypes rather than missing from the tissue. Among
non-neuronal supertypes the criterion excluded Astro_1 and OPC_2, leaving 27 of 29. All
exclusions were applied before testing.

## SM2 — Combined Franken-503 taxonomy construction

*Relocated from main Methods; full detail in `genetics/docs/combined_taxonomy_approach.md`.
Supports Figure 4a; anchors Supplementary Tables T4–T5.*

To test enrichment across cortical and whole-brain cell types, we merged the 137 SEA-AD
MTG supertypes with the Siletti et al. 2023 whole-brain reference (461 clusters). Redundant
types were removed by cross-dataset reciprocal-best-hit (RBH) deduplication using a
MetaNeighbor-style AUROC (Crow et al. 2018) on per-type mean expression (top 5,000
highly-variable genes over 35,321 shared symbols): 95 Siletti clusters were RBH-matched to
a SEA-AD supertype (all mean AUROC ≥ 0.9) and removed, retaining all 137 SEA-AD supertypes
plus 366 novel Siletti clusters — 503 non-redundant types.

## SM3 — SCZ GWAS enrichment specification and version robustness

*Relocated from main Methods; full detail in `genetics/docs/analysis_approach.md`.
Supports Figure 4a–b; anchors Supplementary Fig. S13.*

**Specification.** For each gene, cell-type specificity was computed by column-normalizing
each type's mean expression to a common library size, then row-normalizing each gene across
types (rows sum to 1). Per-cell-type enrichment was tested by MAGMA-style gene-property
regression (de Leeuw et al. 2015) of the gene-level SCZ Z-statistic on specificity,
controlling for log gene size, log SNP count, and log N (one-sided). Genes in the MHC
(chr6:28–34 Mb, hg19) were excluded; p-values were Bonferroni- and BH-FDR-corrected across
tested types. **[BUILD/RETRIEVE: upstream MAGMA gene-analysis parameters — gene window, LD
panel, SNP-to-gene aggregation — see `README_STATUS.md`.]**

**Version robustness.** Because SCZ GWAS is actively updating, we re-ran the enrichment
under the larger Bigdeli et al. 2026 SCZ meta-analysis. The GABAergic/SST enrichment
headline reproduced (per-cell-type enrichment Spearman ρ = 0.95 vs PGC3; gene-Z ρ = 0.855;
SEA-AD FDR-significant types 44 → ~39, still ~90% GABAergic), but the two SST *sub-claims*
attenuated to non-significance (depth gradient r = −0.50 → −0.32, p = 0.20; SST
genetics↔depletion convergence r = 0.645 → 0.387, p = 0.11). Because Bigdeli *contains*
PGC3, this is a dilution robustness check, not an independent replication (Supplementary
Fig. S13). **[DECISION: PGC3 vs Bigdeli as primary — see `README_STATUS.md`.]**

## SM4 — Patch-seq label transfer and electrophysiology

*Relocated from main Methods; anchors Supplementary Table T6.*

A harmonized human cortical patch-seq dataset was consolidated from Lee, Dalley et al. 2023
(779 GABAergic interneurons across human MTG L1–6) and Chartrand et al. 2023 (419 human L1
interneurons) via the `human_int_patch_seq` pipeline, giving 1,155 unique cells across 221
donors after de-duplication (expression harmonized to log₂(FPKM + 1)). Each cell was
assigned a SEA-AD subclass/supertype by a published scANVI run (Xu et al. 2021) and by
Harmony-corrected kNN (15 nearest reference neighbours); the 16 SST supertypes with ≥ 1
assigned cell (150 cells) are shown in Figure 4e. Intrinsic features were extracted from
NWB sweeps with the Allen IPFX library; sag ratio was computed on the largest hyperpolarizing
subthreshold sweep as (V_peak − V_steady)/(V_peak − V_baseline). Agreement between the two
label-transfer methods is reported in Supplementary Table T6 **[RETRIEVE: scANVI-vs-kNN
confusion from `human_int_patch_seq` `compare_knn_vs_scanvi`]**.

## SM5 — RNAscope representative imaging & lipofuscin suppression

*Supports the RNAscope re-analysis (Figure 3; Supplementary Fig. S14). Full pipeline in `histology/results/microscopy/README.md`.*

Quantitative RNAscope analyses used the original SlideBook mask-based SST cell counts from Arbabi et al. (2025). For representative image figures only, RGB composite TIFFs (488 nm SST; 568 nm VIP; 405 nm DAPI) were cropped to the 800×800-px stereological counting frame and processed to suppress lipofuscin autofluorescence by spectral subtraction: lipofuscin was estimated per pixel as the minimum of the SST and VIP channel intensities and subtracted from each signal channel (multiplicative factor 1.2, cross-channel correction 0.15, intensity floor 0.03), applied uniformly to all displayed images. Cell markers overlaid on the raw composites indicate example labelled cells and are illustrative, not exhaustive counts.

---

# Supplementary Figures

**Figure S1 | snRNA-seq integration QC.** UMAP of the integrated 469-donor snRNA-seq nuclei
coloured by dataset, diagnosis, and subclass, showing nuclei group by cell type rather than
by dataset or diagnosis. Supports Fig 1a–d. *(Nicole's; confirm whether this is already main
Fig 1a–d or a separate supplement.)*

**Figure S2 | Xenium cell-type annotation: agreement with the original annotations, panel
marker genes, and panel resolvability.** (a) Reannotated SEA-AD subclass calls against the
dataset authors' (Kwon et al.) independent annotations; each cell is the percentage of a
reannotated subclass assigned to a given author cell type (columns sum to 100), values
≥ 10% labelled, ordered to place the dominant correspondence on the diagonal. (b, d) Marker
dot-plots at subclass level (b; inhibitory → excitatory → non-neuronal) and Sst-supertype
level (d; ordered pia → white matter); dot size, fraction of cells expressing; fill, scaled
mean expression. (c, e) Classification F1 on the 300-gene Xenium panel (blue) versus the
full transcriptome (black), nearest-centroid Pearson classifier under leave-one-donor-out
cross-validation on the SEA-AD neurotypical snRNAseq reference: subclass (c; median 0.96
vs 0.99, 15 of 23 at ceiling) and Sst supertype (e; median 0.62 vs 0.88; Sst_25 = 0.81
highest, Sst_22 = 0.45 and Sst_20 = 0.44 lowest). Rows in b/c and d/e share the y-axis
labels shown at left. Lamp5 Lhx6 is omitted throughout (SM1). Supports Fig 1e; SM1.
READY (canonical): `spatial/supplemental_figures/markers_resolvability_combined`.

**Figure S3 | Xenium vs SEA-AD MERFISH concordance: cell-type proportions and cortical
depth.** (a) Per-donor subclass proportions, Xenium vs MERFISH (Pearson r = 0.85, log₁₀;
ρ = 0.88, n = 23). (b) Neuronal supertype proportions (r = 0.54, log₁₀; ρ = 0.49, n = 106),
illustrating the panel's within-subclass limit. (c) Subclass median cortical depth
(r = 0.96; ρ = 0.95, n = 23). (d, e) Per-supertype cortical depth distributions (0 = pia,
1 = WM) in MERFISH (purple; manual annotation) and Xenium (orange; model prediction), for
the 42 glutamatergic (d) and 64 GABAergic (e) supertypes of b; grouped by subclass and
ordered pia → WM by MERFISH median; violins show the inner 95% of cells, white point =
median; dashed lines, laminar boundaries. Median depths agree at r = 0.97 across these
supertypes and at r = 0.94 within subclass, though Xenium distributions are systematically
narrower (median IQR 0.120 vs 0.145). Supports Fig 1e–f; SM1. READY (canonical):
`spatial/supplemental_figures/xenium_merfish_composite`.

**Figure S4 | Xenium cortical-depth model validation.** (a) Held-out-donor GroupKFold CV of
predicted vs annotated depth (R² = 0.889, n = 368,795 cells). (b) Anti-circularity control:
subclass-identity-only baseline R² = 0.32 vs full neighbourhood model R² = 0.89. (c)
MERFISH-free anchor: eight canonical panel layer markers ordered correctly by predicted
depth (CUX2/LAMP5/CALB1 decline; NR4A2/NXPH4 rise). Supports Fig 1f; SM1. READY:
`spatial/output/depth_validation/`. **[RECONCILE deployed 0.907 vs 0.889 before caption.]**

*(S2 supersedes the former S4 panel markers, S5 resolvability and S6 Kwon agreement
figures; S3 merges the former S2 scatters with the standalone per-supertype depth figure,
dropping its redundant supertype median-depth scatter. Three S-numbers are freed by the two
merges — renumber S5–S14 before submission.)*

**Figure S7 | Per-cell normalisation robustness of the SST and PVALB reductions.** Per-cell
negative-binomial SCZ/control fold change under four normalisations: SST reduction robust
(0.70× raw, P = 0.002; 0.76× library-normalised, P = 0.011), PVALB borderline (0.86× raw,
P = 0.065; 0.87× library-normalised, P = 0.060, n.s.). Supports Fig 2b–c, 2f–g. READY:
`transcriptomic/results/figures/S_percell_metrics`. *(= the existing "Supplementary Fig.
S_percell" callout.)*

**Figure S8 | Supertype-level differential-expression power. [BUILD]**
> ⚠️ **REMINDER — figure to be generated.** Per-supertype nuclei counts vs DE-gene yield
> (or a power curve), showing within-subclass supertype-level DE is underpowered because of
> fewer nuclei per donor. Referenced in main text (Fig 2 Results) as "Supplementary Fig.
> Sx"; no committed figure exists yet.

**Figure S9 | Density-based cross-platform concordance.** snRNA-seq compositional β vs Xenium
cell-density logFC at (a) supertype (neuronal r = 0.55) and (b) subclass (neuronal r = 0.68)
resolution; density corroborates the composition-based concordance (Fig 3f) while avoiding
the compositional zero-sum constraint. Supports Fig 3. READY:
`spatial/output/density_analysis/`.

**Figure S10 | Supertype cortical depth is stable across diagnosis (control).** No supertype's
median depth differs by diagnosis at FDR < 0.05; the strongest shift, Pax6_4, is only marginal
(β = +0.10, p = 9.1 × 10⁻⁴, FDR = 0.094 over 104 supertypes), showing the compositional changes
are not a depth-assignment artifact. Supports Fig 3. READY (canonical):
`manuscript/figures/supplementary/S10_supertype_depth_by_diagnosis`.

**Figure S11 | Conditional (forward-selection) enrichment.** Marginal vs conditional
significance identifying 18 independently enriched SEA-AD supertypes (3 Sst, 4 Pvalb, 3 Pax6,
2 Lamp5, 3 glutamatergic), showing the interneuron enrichment is subtype-specific and
non-redundant. Supports Fig 4a. READY: `genetics/results/figures/conditional_analysis_results`.

**Figure S12 | Gene-driver scatters across top enriched interneuron types.** Specificity ×
SCZ −log₁₀p driver-gene scatter for Sst_2, Pvalb_3, Pvalb_6, Lamp5_5, Pax6_4. Supports Fig 4c.
READY: `genetics/results/figures/gene_drivers/`. *(Main Fig 4c requires the Sst_25 scatter —
see Reminders.)*

**Figure S13 | GWAS-version robustness (PGC3 vs Bigdeli 2026).** (a) Gene-Z concordance
(ρ = 0.855) and per-cell-type enrichment concordance (ρ = 0.95) — the GABAergic/SST
enrichment headline is concordant across GWAS versions. (b) The two version-fragile SST
sub-claims under each GWAS: the depth gradient (*r* = −0.50 → −0.32, *p* = 0.20) and the
genetics↔depletion convergence (*r* = 0.645 → 0.387, *p* = 0.11) attenuate to
non-significance under Bigdeli. Supports Fig 4a–b; SM3. READY: `genetics/data/gwas/magma_bigdeli/`.

**Figure S14 | Independent RNAscope corroboration of the SST reduction in subgenual cingulate cortex.** Re-analysis of sgACC RNAscope FISH (Arbabi et al. 2025; VIP-QC-filtered N = 54; 15 control, 11 SCZ). (a) Aggregate SST density (cells/mm²), control vs SCZ (mixed-model *P* = 0.071). (b–c) Layer-stratified SST density — L2/3 (*P* = 0.126) and L5/6 (*P* = 0.180) — the deficit is numerically larger superficially. (d) SST-density diagnosis coefficients vs control (density units): MDD B = +0.7 (*P* = 0.813), bipolar B = −3.4 (*P* = 0.214), schizophrenia B = −5.4 (*P* = 0.071) — a graded psychosis-spectrum pattern (the equivalent per-frame count-model SCZ coefficient is β = −0.60; Methods). (e) Representative RNAscope micrographs of L2/3 control vs SCZ, lipofuscin-suppressed for display (Supplementary Methods SM5); markers indicate example SST cells and are illustrative. Supports Figure 3 — a cross-region, cross-modality corroboration. READY: `histology/results/fig_main_result` + `histology/results/microscopy/`.

**Option-C spares (hold unless reviewers push):** *S15* depth-stratified non-neuronal changes
(L2/3 oligodendrocyte depletion, endothelial increase — Xenium-only, not cross-platform
validated); *S16* ATAC-vs-expression divergence (chromatin → excitatory, expression →
inhibitory); extra per-gene forest exemplars (BDNF, FKBP5, CX3CR1…).

---

# Supplementary Tables

- **T1 | Dataset demographics** — per-dataset sex and post-mortem interval. **[ASSEMBLE]** (= the "Supplementary Table SX" sex/PMI callout).
- **T2 | Cell-type–specific DE meta-analysis** — all gene × subclass estimates, SE, p, FDR (7-dataset). READY: `transcriptomic/data/figure_inputs/DE_genes_all_cells_scz.csv`.
- **T3 | Compositional results** — per-supertype crumblr β/SE/FDR, both platforms. READY: `spatial/output/crumblr/`.
- **T4 | Franken-503 SCZ enrichment** — all 503 types (β, p, FDR, Bonferroni); 120 of 503 FDR-significant (54 Bonferroni). The cortical GABAergic/SST pattern matches the SEA-AD-only analysis; brain-wide, the strongest signals are subcortical MGE-derived interneuron clusters (MGE_260, β = 61.3, *p* = 4.4 × 10⁻¹⁶; MGE_263; Misc_132), consistent with the whole-brain (Siletti) analysis of Duncan et al. (2025) (T5). READY: `genetics/results/tables/rbh_combined_enrichment.csv` (= existing "Supplementary Table SX"). **T4b** = SEA-AD-137 standalone cortical-only enrichment (companion/robustness; 44/137 FDR-significant, 41 GABAergic, top signal Lamp5_5). Panel 4a and the main analysis use the Franken-503 values (T4), not this cortical-only version — the cortical-only reference re-ranks the supertypes, underscoring why Franken-503 is primary.
- **T5 | Siletti whole-brain enrichment** — 461-cluster enrichment (top brain-wide = MGE_260, MGE_263, Misc_132). READY: `genetics/results/tables/reprocessed_siletti_461_enrichment.csv`.
- **T6 | Patch-seq label-transfer agreement** — scANVI vs Harmony-kNN confusion/agreement. **[RETRIEVE]** from `human_int_patch_seq`.
- **T7 | Cross-platform concordance matrix** — subclass/supertype × composition/density r, p, n. READY: `spatial/output/{crumblr,density_analysis}/`.

---

# Reminders — figures/tables still to produce before submission
- **MAIN Fig 4c: Sst_25 gene-driver scatter with HCN1 labelled — [BUILD].** (HCN1 absent from the committed Sst_2 driver list; main text carries the `[TODO]`.)
- **Supp Fig S8: supertype-level DE power — [BUILD].**
- **Supp Table T6: patch-seq scANVI-vs-kNN agreement — [RETRIEVE].**
- **Supp Table T1: dataset sex/PMI — [ASSEMBLE].**
- **Nicole's supplements — [PENDING].** Her own supplemental figures/tables (from the snRNA-seq DE/composition analyses) will be inserted into this S/T numbering at a later time; leave slots and renumber on merge.

# Reconcile before finalizing any caption (see `README_STATUS.md`)
- Depth-model R²: deployed 0.907 vs committed 0.889 — pick the canonical variant (S3).
- Fig 3f concordance r = 0.50 vs 0.45 (QC gate).
- L2/3-oligo interaction FDR under the redeployed depth model (if S15 is used).
- Delete stale `transcriptomic/results/tables/08_meta_vs_xenium_pairs.csv` (r = 0.70/74%) so it can't leak into T2/T7.

---

# Main-Methods rewiring (APPLIED 2026-07-21)

These main-Methods subsections were shrunk to the pointers below; the removed detail now
lives in SM1–SM4 above.

**`### Xenium spatial transcriptomics` → becomes:**

> We reanalysed the Xenium DLPFC dataset of Kwon et al. (2026; 24 sections, 12 SCZ / 12
> control; 1.34 million cells; 300-gene panel). Because 300 genes are too few for de novo
> clustering, cells were assigned to SEA-AD subclasses and supertypes by a self-referencing
> two-stage correlation classifier and to a continuous cortical depth (0 = pia, 1 = white
> matter) by a neighbourhood-composition model; both were validated against the independent
> SEA-AD MERFISH atlas (subclass proportions r = 0.85, depth r = 0.96; 84.9% classifier
> accuracy) and the panel's supertype-resolution limits characterized by a
> panel-vs-transcriptome F1 benchmark (Supplementary Methods SM1; Supplementary Figs S2–S4).
> Disease analyses used one shared cell definition (cortical, quality-passing cells):
> composition via stratified crumblr (neuronal and non-neuronal compartments; covarying age,
> sex), densities (cells/mm²) as a constraint-free complement, and DE via pseudobulk edgeR
> (~ diagnosis + sex + age), mirroring the snRNA-seq meta-analysis. For the L6b aggregation,
> four sections with < 3% L6 cells (Br2039, Br5973, Br2719, Br5314) were excluded.

**GWAS Methods → keep the enrichment model + results; move construction/specifics/robustness
to SM2–SM4** with pointers, e.g. "…across a combined Franken-503 taxonomy (Supplementary
Methods SM2)…"; "…MAGMA-style gene-property regression (Supplementary Methods SM3)…";
"…patch-seq cells assigned to SEA-AD supertypes by scANVI and Harmony-kNN (Supplementary
Methods SM4)…"; and the robustness paragraph → "…(Supplementary Methods SM3; Supplementary
Fig. S13)."
