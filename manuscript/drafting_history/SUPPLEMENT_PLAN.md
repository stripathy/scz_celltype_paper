# Supplement & Supplementary-Methods plan

Planning doc (Claude, 2026-07-21) — options for where trimmed/moved detail should live:
main Methods vs **Supplementary Figure** vs **Supplementary Table** vs **Supplementary
Methods**. Grounded in an inventory of committed assets across `spatial/`,
`transcriptomic/`, and `genetics/`. Nothing here edits the manuscript; it's a menu to
choose from. Companion to [`README_STATUS.md`](README_STATUS.md).

## Guiding principle (the destination rule)

- **Main text (Results)** = the *claim* + its headline number.
- **Main Methods** = a *concise pointer* — what was done in 2–4 sentences.
- **Supplementary Methods** = the *full novel procedure* (algorithm, parameters, QC,
  anti-circularity controls) that a methods-scrutinizing reviewer needs but that would
  bloat the main story.
- **Supplementary Figure** = a *visual* validation/robustness result that supports a
  main claim but isn't itself a headline (a scatter, dumbbell, CV curve, heatmap).
- **Supplementary Table** = *values/lists* (full result tables, per-type stats,
  agreement matrices).

Litmus test used below: if the moved detail has (or deserves) a **figure**, promote it
to a numbered Supp Fig and *cite it at the point of the claim* rather than describing it
in Methods prose. If it's a **procedure**, it goes to Supp Methods. If it's **values**,
Supp Table.

---

## Proposed Supplementary Methods sections

Four self-contained blocks. Each shrinks a large main-Methods passage to a pointer.

### SM1 — Xenium cell typing & spatial inference (the big one)
The novel, ~2-page block. **Main Methods keeps ~3 sentences**; the rest moves here.

> **Main-Methods residue (keep):** "Xenium cells were assigned to SEA-AD subclasses and
> supertypes by a two-pass correlation classifier, and to a continuous cortical depth
> (0 = pia, 1 = white matter) by a neighbourhood-composition model; both were validated
> against the independent SEA-AD MERFISH atlas and the panel's supertype-resolution
> limits characterized by a panel-vs-transcriptome benchmark (Supplementary Methods SM1;
> Supplementary Figs S2–S4)."

Moves to SM1: two-stage centroid classifier (MapMyCells exemplars, z-scored Pearson,
subclass→within-subclass supertype, margin); QC/doublet flags and the Harmony-kNN
failure that motivated the centroid approach (Sst inflated 12.1% vs ~2.5%; 69% agreement);
the GradientBoosting depth model (neighbourhood composition + one-hot self-subclass;
low-CPS training; K, smoothing); the two anti-circularity controls (GroupKFold out-of-fold;
subclass-only baseline R²=0.32 vs 0.89); the MERFISH-free marker-anchor check; BANKSY
segmentation + layer binning; the panel-resolvability benchmark (information ceiling).

### SM2 — Combined Franken-503 taxonomy construction
RBH / MetaNeighbor-style AUROC dedup merging SEA-AD 137 + Siletti 461 → 503. Fully
drafted already in `genetics/docs/combined_taxonomy_approach.md`. Main Methods keeps one
sentence + pointer.

### SM3 — GWAS enrichment specification & version robustness
MAGMA gene-property regression specifics (specificity normalization, covariates, MHC
exclusion), **and** the PGC3-vs-Bigdeli robustness with the explicit "Bigdeli *contains*
PGC3, so this is a dilution robustness check, not an independent replication" caveat.
Anchors Supp Fig S13.

### SM4 — Patch-seq label transfer & feature extraction
scANVI + Harmony-corrected kNN assignment to SEA-AD supertypes; IPFX sag extraction.
Anchors Supp Table T6 (scANVI-vs-kNN agreement — needs retrieval).

---

## Per-figure supplement recommendations

Status legend: **READY** = rendered on disk (regenerable from committed scripts/CSVs);
**BUILD** = must be generated; **RETRIEVE** = lives in a sibling repo.

### Figure 1 — Xenium annotation & layers

| Moved/detailed item | → Destination | Asset | Status |
|---|---|---|---|
| Classifier accuracy 84.9% + Xenium↔MERFISH proportions r=0.85 / depth r=0.96 | **Supp Fig S2** | `spatial/output/celltyping_supplement/celltyping_validation_supplement`; `output/merfish_benchmark/subclass_confusion_matrix` | READY (promote to canonical) |
| Depth-model validation: GroupKFold R²=0.889, anti-circularity baseline R²=0.32, marker-anchor ordering | **Supp Fig S3** | `output/depth_validation/depth_validation_supplement` + `anchor_marker_gradients` (+ `lowcps/depth_learning_curve`) | READY |
| Agreement with Lieber/Kwon original annotations (a) + panel marker dotplots (b, d) + panel resolvability F1 (c, e) | **Supp Fig S4** | `spatial/supplemental_figures/markers_resolvability_combined` | READY (canonical; merges the former S4/S5/S6) |
| Two-stage classifier, depth model, BANKSY, QC/doublet procedure | **SM1** | `spatial/methods_writeup.md` | text |
| Harmony-vs-classifier justification | **SM1** (fig optional) | `output/presentation/subclass_proportions_by_method` | READY |

### Figure 2 — Cross-platform DE

| Item | → Destination | Asset | Status |
|---|---|---|---|
| Per-cell normalisation robustness (SST/PVALB across 4 norms) | **Supp Fig S7** (= existing `S_percell` placeholder) | `transcriptomic/results/figures/S_percell_metrics` | READY |
| Supertype-level DE underpowered (the `Supplementary Fig. Sx` placeholder) | **Supp Fig S8** | — | **BUILD** (per-supertype nuclei vs DE yield / power) |
| Full DE meta table (all gene×subclass) | **Supp Table T2** | `data/figure_inputs/DE_genes_all_cells_scz.csv` | READY |
| Extra per-gene forests (BDNF, FKBP5, CX3CR1…) | optional Supp Fig | `results/figures/07_forest_composite` | READY |
| Concordance FDR<0.05 sensitivity | optional (fold into S-concordance) | `results/figures/08_*_fdr05` | READY |
| DE-vs-proportion labelled (script 16) | **skip** — redundant with Fig 2i inset | `de_vs_proportion_subclass` | READY |

### Figure 3 — Composition

| Item | → Destination | Asset | Status |
|---|---|---|---|
| Density-based concordance (corroborates composition r=0.50) | **Supp Fig S9** | `spatial/output/density_analysis/snrnaseq_vs_density_{supertype,subclass}` (r=0.55 / 0.68) | READY |
| Supertype depth stable by diagnosis (null control: composition ≠ depth artifact) | **Supp Fig S10** | `spatial/supplemental_figures/supertype_depth_casecontrol` | READY (canonical #4) |
| Full concordance value matrix (subclass/supertype × composition/density) | **Supp Table T7** | `output/crumblr/*`, `output/density_analysis/*` | READY |
| Full compositional results (crumblr, both platforms) | **Supp Table T3** | crumblr result CSVs | READY |
| Panel resolvability (shared) | cite **Supp Fig S4** (panels c, e) | — | READY |
| Depth-stratified non-neuronal (L2/3 oligo depletion, endothelial) | optional Supp Fig (Option C) | `output/depth_proportions/oligo_L23_detail`, depth profiles | READY (Xenium-only — flag) |

### Figure 4 — GWAS + patch-seq

| Item | → Destination | Asset | Status |
|---|---|---|---|
| Full Franken-503 enrichment (503 types) | **Supp Table T4** (= existing `Supplementary Table SX`) | `genetics/results/tables/rbh_combined_enrichment.csv` | READY |
| SEA-AD-137 standalone enrichment (source of panel 4a) | **Supp Table T4b** | `seaad_magma_scz_enrichment_all_supertypes.csv` | READY |
| Conditional/forward-selection (18 independent types) | **Supp Fig S11 + Table** | `results/figures/conditional_analysis_results`; `independent_supertypes_conditional.csv` | READY |
| Siletti whole-brain (MGE) enrichment (461) | **Supp Table T5** (+ optional Manhattan) | `reprocessed_siletti_461_enrichment.csv`; `rbh_combined_taxonomy_manhattan` | READY |
| Gene-driver scatters across top enriched interneuron types | **Supp Fig S12** | `results/figures/gene_drivers/gene_driver_scatter_{Sst_2,Pvalb_3,Pvalb_6,Lamp5_5,Pax6_4}` | READY (5 committed) |
| **Sst_25 driver scatter (main Fig 4c)** | **MAIN FIG** | — | **BUILD** (HCN1 absent from Sst_2 list) |
| GWAS-version robustness (PGC3 vs Bigdeli) | **Supp Fig S13 + SM3** | `data/gwas/magma_bigdeli/validate_*` | READY |
| Franken-503 RBH dedup procedure | **SM2** | `docs/combined_taxonomy_approach.md` | text |
| MAGMA / specificity procedure | **SM3** | `docs/analysis_approach.md` | text |
| Patch-seq scANVI-vs-kNN agreement | **Supp Table T6 + SM4** | — | **RETRIEVE** (`human_int_patch_seq`) |
| ATAC-vs-expression divergence (chromatin→excitatory vs expression→inhibitory) | optional Supp Fig (Option C) | `results/figures/atac_vs_gwas_manhattan` | READY |

---

## Draft numbered manifest (Option B — recommended)

**Supplementary Figures:** S1 snRNA-seq integration QC (Nicole; may be main Fig 1) · S2
Xenium↔MERFISH cell-typing validation · S3 depth-model validation (CV + anchor) · S4
cell-type annotation (author agreement + marker dotplots + panel resolvability; merges the
former S4/S5/S6, leaving S5/S6 vacant pending a renumber) ·
S7 per-cell normalisation robustness · **S8 supertype DE power [BUILD]** · S9 density
concordance · S10 supertype-depth null control · S11 conditional enrichment · S12
gene-driver scatters · S13 PGC3-vs-Bigdeli robustness.

**Supplementary Tables:** T1 dataset demographics/sex/PMI [assemble] · T2 full DE
meta-analysis · T3 full compositional results · T4 Franken-503 enrichment (+T4b SEA-AD-137) ·
T5 Siletti 461 enrichment · **T6 patch-seq scANVI-vs-kNN agreement [RETRIEVE]** · T7
concordance value matrix.

**Supplementary Methods:** SM1 Xenium cell typing & spatial inference · SM2 Franken-503
construction · SM3 GWAS enrichment & robustness · SM4 patch-seq label transfer.

**Figure 3 add-on (Option 4, applied 2026-07-21):** **S14** RNAscope corroboration of the
SST reduction in sgACC (Arbabi 2025 re-analysis; VIP-filtered stratified model) —
`histology/results/`; shipped as a Fig 3 sentence + Discussion passage + dataset-table row
+ Supplementary Fig. S14 + Supplementary Methods SM5.

**Option C add-ons (if scope broadens / reviewers push):** S15 depth-stratified
non-neuronal (L2/3 oligo) · S16 ATAC-vs-expression divergence · extra forest exemplars.

**GSEA pathway story** (OxPhos-in-inhibitory / cholesterol-in-excitatory): recommend
**separate companion paper** — too large and orthogonal to be a supplement. If a hook is
wanted, one Supp Fig (theme heatmap) + one Supp Table, cross-referenced from Discussion,
after clearing its two hygiene flags.

---

## Must-build vs ready

- **Blocks MAIN figures:** Sst_25 gene-driver scatter (Fig 4c); supertype DE-power fig (Fig 2, S8).
- **Supplement-only build/retrieve:** patch-seq scANVI-vs-kNN table (T6); dataset sex/PMI table (T1).
- **Everything else on the S/T list is READY** (rendered; regenerable from committed scripts+CSVs).

## Reconcile before writing any supplement caption
- Depth-model R²: committed 0.889 (K=50/24-donor) vs deployed 0.907 (low-CPS) — pick the canonical variant (3 renders exist).
- Fig 3f concordance r=0.50 vs 0.45 (QC-gate).
- L2/3-oligo interaction FDR under the redeployed depth model (may weaken to ~0.09).
- Delete stale `transcriptomic/results/tables/08_meta_vs_xenium_pairs.csv` (r=0.70/74%) so it can't leak into T-tables.
