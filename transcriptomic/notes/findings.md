# Findings — annotated take-stock

> **Update 2026-05-28**: per-cell-type volcanos now exist for all 23 SEA-AD subclasses (`results/figures/01_volcano_<celltype>.png`). Full GSEA pipeline re-ran from scratch — numerical results identical, confirming reproducibility. Systematic literature review per cell type now lives at [`notes/literature_per_celltype.md`](literature_per_celltype.md); see Section 6 below for additional findings beyond the OxPhos/cholesterol stories that emerged from that review.

Working notes from the analysis session that built this pipeline. Intended as a hand-off / memory aid: what we found, what we ruled out, what to revisit. Numbers are pulled from `results/tables/` — see source files for the underlying tests.

---

## 1. The four robust transcriptional themes

After GSEA with FDR per cell type **and** multi-testing correction (Fisher BH + Stouffer signed-Z, across 23 cell types × ~12K pathways), these survive:

### 1.1 Synaptic / vesicle machinery — pan-cell-type DOWN
- Top Stouffer hits: `GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING` (Z = -15.2, padj = 2e-48), `GOCC_SYNAPTIC_MEMBRANE` (Z = -14), `GOCC_PRESYNAPSE` (Z = -12).
- Strongest in **Astro, Oligo, OPC, Vip, Chandelier** — surprisingly *more* significant in glia than in neurons. Worth thinking about: are glial cells losing their *synaptic-support* gene programs (Astro tripartite synapse), or is this a marker-gene confound that needs explicit testing?
- Leading-edge drivers include classic SCZ synaptic genes: SNAP25, SYT1, DNM1, STX1A, GRM4, ELFN1, CNTNAP1, NRGN, GRIN2A, GABBR2, CACNB2, CPLX1, ERBB4.

### 1.2 OxPhos / mitochondrial — inhibitory-selective DOWN
- HALLMARK_OXIDATIVE_PHOSPHORYLATION: Sst (NES = -2.14, padj = 3e-5), Vip (-1.80, padj = 0.02), Lamp5 (-1.72, padj = 0.06), Pvalb (-1.65, padj = 0.06).
- *All* mitochondrial subprocesses are hit in Sst (Complex I/III/IV/V, TCA, mito translation, mito import, inner membrane). Not isolated to one complex — wholesale mitochondrial program failure.
- **Mechanism dissection** (`exploratory/mechanism_oxphos_tests.R`) ruled out:
  - PGC-1α suppression (axis NES is positive in Sst; biogenesis TFs are mostly flat)
  - Glycolytic compensation (ρ(OxPhos, Glycolysis) = +0.44 — they move together, not oppositely)
  - Activity-driven (ρ(OxPhos, IEGs) = -0.15, n.s.; IEGs slightly UP in Sst)
- **Sharper hypothesis** from gene-level pull: **mitochondrial protein import + ETC complex assembly + quality control failure**. Key genes individually down in Sst:
  - **AFG3L2** (mito AAA protease) — z = -3.74, padj = 0.030
  - **TOMM40** (outer-membrane translocase) — z = -3.26, padj = 0.069
  - **UQCC2** (Complex III assembly) — z = -4.01, padj = 0.019
  - **PINK1** (mitophagy initiator) — z = -2.33
  - 57% of 30 assembly/chaperone genes negative; PPRC1, POLRMT, PPARGC1B subtly down (PGC-1 family at the gene level, even though the whole-axis GSEA misses it because non-PGC1 TFs dilute the signal)
- Curious cross-disease resonances: AFG3L2 (SCA28), TOMM40 (AD GWAS), PINK1 (early-onset PD).

### 1.3 Cholesterol / sterol biosynthesis — excitatory + macroglia DOWN
- Strongest hits: REACTOME_CHOLESTEROL_BIOSYNTHESIS in L5 IT, L5_6 NP, L6 IT, Astro, Oligo (all padj < 0.01).
- Coordinated suppression of the entire mevalonate pathway: ACAT2 → HMGCS1 → **HMGCR** (statin target) → MVK → PMVK → MVD → IDI1 → FDPS → FDFT1 → LSS → ... → DHCR7/DHCR24, capped by master regulator **SREBF2** (also a PGC3 GWAS hit).
- **OPCs go the opposite direction** (NES > 0) — possible compensatory ramp-up for differentiation/myelin biogenesis.
- Broader lipid program (sphingolipid, phospholipid metabolism) also down in macroglia — not purely cholesterol-specific in glia.
- INSIG1 (negative-feedback regulator of SREBP processing) also down → consistent with master TF being below threshold.

### 1.4 Growth / stress / chaperones — UP (glial-biased)
- HALLMARK_MYC_TARGETS_V1, HALLMARK_MTORC1_SIGNALING, REACTOME_HSF1_ACTIVATION, GOBP_CHAPERONE_MEDIATED_PROTEIN_FOLDING all UP in 5+ cell types, particularly Oligo, OPC, Micro-PVM, Astro, plus L6b (an outlier excitatory).
- Reads as compensatory / stress response. Recurrent leading-edge driver: HSPD1 (HSP60, mitochondrial chaperone, PGC3 GWAS hit).

---

## 2. Signals that are real but heterogeneous

### Translation / ribosome (bidirectional)
- Survives Fisher BH (padj = 6e-30) but **fails Stouffer** because direction varies by cell type: DOWN in Sst, Pvalb, Pax6; UP in L6b, Astro, Oligo, OPC, Micro-PVM.
- Interesting biology — not a uniform "translation collapse" but an opposing reorganization. Worth a focused look.

### Axon guidance / SLIT-ROBO
- Survives Fisher but not Stouffer. Strongest signal in Sst (`REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS` NES = -2.05, padj = 1e-4).
- Heterogeneous direction across other cell types.

### BMP / TGF-β signaling — UP (inhibitory + Astro)
- `REACTOME_SIGNALING_BY_BMP` significantly UP in Vip (padj = 0.011), Sncg (0.016), Astro (0.047); trend in Lamp5, Pvalb, Sst, Endo.
- Leading-edge: SMAD1, SMAD5, SMAD9 (R-SMADs), ACVR2A, ACVR2B, BMPR1A (receptors). Clean BMP receptor-SMAD signature.
- Survives Fisher BH (padj = 1e-19) and Stouffer (Z = +9.8).

---

## 3. GWAS convergence — supports the narrative, doesn't drive it

PGC3 SCZ GWAS genes (484 genes from Trubetskoy 2022, no-MHC version) tested as a custom gene set via GSEA:

- Significantly DOWN-shifted in 12 of 23 cell types. Strongest: **Oligo** (NES = -1.71, padj = 2e-4), L2_3 IT (-1.51, padj = 0.005), OPC (-1.48, 0.005), Vip (-1.43, 0.005), Astro (-1.30, 0.027).
- With-MHC and without-MHC versions agree — MHC region not driving it.
- Specific GWAS genes that drive the major themes (from leading-edge × GWAS intersection):
  - **SREBF2** → cholesterol biosynthesis (5 cell types: Astro, L5 IT, L5_6 NP, L6 IT, Oligo)
  - **GOT2** → OxPhos (3 inhibitory: Lamp5, Pvalb, Sst)
  - **NRGN, GRIN2A, GABBR2, CACNB2, CPLX1, ERBB4** → synaptic transmission
  - **HSPD1** → MYC targets (6 cell types)
- **Caveat we explicitly chose not to resolve**: marker-gene confound. Are PGC3 genes enriched for cell-type-specific markers, making them down-shift purely because cell identity is dampened in disease? Test scaffolded but not run — see "Open questions" below.
- **Reframing**: the dominant transcriptional findings (OxPhos in Sst, cholesterol in glia/exc, synaptic broadly) stand on their own statistical and biological merits. GWAS overlap is supporting evidence, not the foundation.

---

## 4. Open questions / next steps

### Mechanistic
- **Mitochondrial machinery failure in Sst**: covariance of TOMM40 + AFG3L2 + PINK1 + UQCC2 across per-study estimates — do they share an upstream regulator? Requires per-study data, not pooled meta-analytic.
- **OPC cholesterol UP-regulation**: is this coupled to differentiation markers? Test OPC cholesterol NES vs OPC differentiation gene set NES.
- **SREBF2 axis as bottleneck**: per-sample correlation of SREBF2 expression vs downstream enzyme expression would test whether one TF is the choke point.
- **Translation bidirectionality**: which specific ribosomal/translation factors are up in glia vs down in inhibitory? Identify the divergent subset.

### Methodological
- **Marker-gene confound test for GWAS overlap** (option (c) from the session — scaffolded but not run). Use SEA-AD specificity matrix in `~/Github/scz_cell_type_enrichment/data/conti_specificity_matrix.txt`. Tests: (1) are PGC3 genes unusually cell-type-specific, (2) does GSEA enrichment survive controlling for specificity, (3) do non-SCZ GWAS sets (AD, BIP, height) show similar broad-down pattern (= generic effect) or only SCZ (= disease-specific).
- **Within-class subtype resolution**: aggregate to subclass loses heterogeneity. Sst > Vip > Lamp5 > Pvalb ordering for OxPhos collapse is interesting (Pvalb has highest metabolic demand but isn't most affected — argues against simple "fast-spiking vulnerability" hypothesis). Subtype-level analysis might sharpen.
- **PGC-1α-axis-specific GSEA**: the diluted 23-gene set hid the signal. Test with a minimal 5-gene set (PPARGC1A, PPARGC1B, PPRC1, POLRMT, ESRRA).

### Cross-modal
- **Compare to PGC3 prioritized genes** (fine-mapped credible sets, not just lead-SNP nearest-gene) — likely cleaner GWAS overlap.
- **PRS-stratified expression**: do high-PRS controls and patients converge on the same DE patterns?
- **Layer/depth analysis**: in Astro (broadest signal) and Sst (strongest OxPhos), does the DE pattern vary by cortical depth? Spatial transcriptomics from SEA-AD MERFISH could test.

### Out of scope but flagged
- **Cell composition confound**: meta-analysis effect sizes are within-cell-type but cell composition changes between SCZ and controls could distort baseline. The Endresz et al. cohort composition data in `~/Github/scz_cell_type_enrichment/results/intermediates/` could be cross-referenced.
- **Antipsychotic medication effects**: not addressable from these data without patient-level metadata. Worth a literature pull on antipsychotic-induced transcriptomic signatures vs. our findings (especially for the metabolic/cholesterol signals — clozapine has well-documented metabolic effects).

---

## 6. Additional themes from per-cell-type literature review (2026-05-28)

The cell-type-by-cell-type review surfaced findings beyond the OxPhos/cholesterol stories. Worth elevating:

### Strongest single-gene anchors (manuscript-ready)
- **NTNG1 ↓ in L6 IT Car3** (estimate = -0.54, padj = 8.3e-12, *strongest single-gene hit in the entire meta-analysis*). Directly replicates [Aoki-Suzuki 2005 Biol Psych](https://pubmed.ncbi.nlm.nih.gov/15705354/) and [Eastwood & Harrison 2008 Neuropsychopharm](https://www.nature.com/articles/1301457). **Novel cell-type substrate** for a long-known SCZ DEG → strong pipeline validation.
- **PDE4B ↓ in L5 ET** — replicates Millar/Fatemi DISC1-PDE4 work; localizes a classic SCZ gene to a specific cell type.
- **ABCG2 ↓ in Endo** — directly replicates [Cai 2019 Mol Psych](https://www.nature.com/articles/s41380-018-0235-x) high-inflammation SCZ subgroup BBB transporter loss.

### New cell-type stories surfaced

1. **Astrocyte interferon-gamma / alpha response UP** (NES = 2.7-2.9, padj < 1e-9) — partial-to-novel; magnitude is unusually strong. SP110 (IFN-inducible PML body protein) is the novel gene anchor. Connects to [Warre-Cornish 2020 iPSC IFN-gamma reproduces SCZ phenotypes](https://www.science.org/doi/10.1126/sciadv.aay9506) and Bast 2025 preprint mitochondrial signatures. Astrocyte inflammation has emerged as a parallel theme to the existing synaptic-down + cholesterol-down signal.

2. **Microglia translation UP** (REACTOME_TRANSLATION NES = 3.17, padj = 1e-14 — *strongest GSEA hit in the dataset*). Partial-to-novel for SCZ specifically. Consistent with activated-microglia signature; partially contradicts [Snijders 2021](https://onlinelibrary.wiley.com/doi/full/10.1002%2Fglia.23962) (homeostatic-marker loss without activation). Meta-analytic power may resolve activation that smaller cohorts missed.

3. **L6 IT Car3 synaptic membrane DOWN** — claustrum-projecting cell type with deep significance (padj = 2e-8) and the NTNG1 anchor above. Notable because Ji 2022 reported decreased L5/6_IT_CAR3 proportion in SCZ patients — proportion *and* per-cell transcriptome both affected.

4. **L6 CT chaperone / HSF1 UP** — partial / novel-for-L6CT extension of Arion 2007 cortical chaperone upregulation. Biologically plausible: long-range corticothalamic projection neurons have high metabolic demand. Worth pairing with the broader chaperone-up signature in L2_3 IT.

5. **L2_3 IT cell-cycle-mitotic UP** — paradoxical because CDC25B is DOWN. **Novel in postmortem SCZ snRNA-seq**. Possible cell-cycle re-entry / stress signature; warrants follow-up.

6. **OPC vs Oligo cholesterol directional divergence** confirmed: Oligo DOWN, OPC UP. Combined with OPC mTORC1 UP + HK2 UP, suggests **compensatory proliferative-stress state in OPCs** as mature Oligos lose function. Consistent with Mauney 2015 OLIG2+ density drop without NG2+ loss (impaired differentiation).

7. **L6b innate immune + chaperone-attenuation UP, glycosyl hydrolase DOWN** — anchored to the [Hoerder-Suabedissen 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3587197/) subplate-SCZ-risk story. Adult human L6b SCZ literature was previously sparse; this work adds substantial cell-type-specific signature.

### Themes that contradict prior literature (worth noting upfront)

- **mTORC1 UP in Chandelier + OPC** vs bulk literature mTOR DOWN ([Chadha 2020](https://www.nature.com/articles/s41386-020-0614-2), [Ibarra-Lecue 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7105616/)). Reconcilable as cell-type heterogeneity: minority cell types ramp mTORC1 while pyramidal-dominated bulk goes down.
- **Microglia activation signal** vs Snijders 2021 homeostatic loss.
- **Cholesterol DOWN in L5/L6 IT neurons** vs clozapine-induced UP in iPSC neurons (Kruse 2021) — argues for disease-not-drug origin.

### Themes likely to be artifacts (flag in any manuscript)

- **Male-gamete-generation / sex-differentiation UP** in multiple inhibitory cells (Sst, Pvalb, Sncg, Pax6) and excitatory cells (L5 IT, L6 IT Car3). Almost certainly shared cilium / dynein / RNA-binding gene overlap with sperm flagella — MSigDB cross-tissue artifact. Recommend re-running GSEA after removing testis-restricted leading-edge genes before reporting as biology.

## 5. Where things live

- **Final GSEA results table**: `results/tables/gsea_gsea_all_celltypes.csv` (135K rows, one row per pathway × cell type × NES + p)
- **Multi-test corrected pathway summary**: `results/tables/gsea_pathway_multitest_summary.csv`
- **Leading-edge genes (long format)**: `results/tables/gsea_leading_edge_genes_long.csv`
- **Leading-edge × GWAS overlap**: `results/tables/gsea_leading_edge_x_gwas_overlap.csv`
- **Curated theme heatmap**: `results/figures/gsea_GSEA_themes_multitest.png`
- **OxPhos story figure**: `results/figures/story_oxphos_figure_oxphos_inhibitory.png`
- **Cholesterol story figure**: `results/figures/story_cholesterol_figure_cholesterol_glia_exc.png`
- **Mechanism testing figure**: `results/figures/story_oxphos_mech_figure_oxphos_mechanism_tests.png`
- **GSEA-vs-ORA explainer**: `results/figures/06_explainer_Sst_OXPHOS.png`
