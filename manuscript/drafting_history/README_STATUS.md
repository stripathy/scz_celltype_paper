> ## ⚠️ STALE — the live manuscript is the Google Doc (`1cO5ZStbp9b3cb6vb9s6H2YfbprYGikXBvFHQR7gROjo`)
>
> Since 2026-08-20 every number, figure reference and decision lives in the Doc; this
> file is the 2026-07-31 assembled draft, committed for provenance only. Do not edit it
> or quote from it. The dated `figure5_*_2026-09-01.md`, `figure3_tail_paragraph_*.md`
> and `*_PLAN.md` files in this folder are decision records that were pasted into the
> Doc; treat them the same way.

# Manuscript — status & open decisions

**Started 2026-07-21 (Claude Code)**, assembled from the repo components, the group
talk (`martinowich_group_talk_june_2026.pptx`), and the two Google Docs (main draft
`1cO5Z…`; GWAS methods `1IaqyjU5…`). Rule throughout: every statistic is extracted from
a source file (never memory/plots); unverified values are visible `[TODO: …]` markers.

**Consolidated to four files (2026-07-22):** the many per-section drafts were stitched
into `manuscript.md` and then removed (archived outside the repo). `manuscript.md` is
now the single source of truth.

## Assembled draft

**[`manuscript.md`](manuscript.md)** is the single polished, submission-style draft —
all sections stitched in reading order, scaffolding stripped, Nicole's rough sections
(Fig 1a–d, Fig 3a–c, snRNA-seq + crumblr methods) redrafted, and the resolved decisions
below applied. It carries inline `[TODO]`s consolidated in its own Outstanding-items
appendix; the per-section drafts that were stitched into it (with their provenance
tables) have been removed and archived outside the repo.

**Supplement (Option B, applied 2026-07-21):** the novel Xenium label-transfer/depth
pipeline and the GWAS taxonomy/specificity/robustness detail now live in
[`supplementary.md`](supplementary.md) as **Supplementary Methods SM1–SM5**, with a
numbered Supplementary Figure (S1–S14) / Table (T1–T7) manifest; the main Methods point
to them. Planning rationale + full asset inventory in [`SUPPLEMENT_PLAN.md`](SUPPLEMENT_PLAN.md).
Two supplement items still to make: **Supp Fig S8** (supertype DE power) and **main Fig 4c**
Sst_25 driver scatter — both `[BUILD]`; **Supp Table T6** (scANVI-vs-kNN) `[RETRIEVE]`;
**Supp Table T1** (sex/PMI) `[ASSEMBLE]`. **Nicole's own supplemental figures/tables**
(from her snRNA-seq DE/composition analyses) will be inserted into the S/T numbering
later — reserve slots and renumber on merge.

**RNAscope corroboration (Option 4, applied 2026-07-21):** the `histology/` RNAscope
re-analysis (Arbabi 2025, sgACC) is now integrated — a Fig 3 sentence + a Discussion
passage + a dataset-table row + **Supplementary Fig. S14** + **Supplementary Methods SM5**
(VIP-filtered stratified model; SST β = −0.60, *P* = 0.071).

**Figure 4 polished figure ported (2026-07-22):** the cleaned-up R-rendered Fig 4
(then `scz_sst_hcn1_multipanel_R_v3`, 6 panels A–F; **now `scz_sst_hcn1_figure4`,
10 panels a–j** — the intermediate 6-panel, 7-panel and compact variants were deleted
2026-08-04, and the figure is assembled by a single `build_figure4()` call in
`fig4_assemble.R`) + its supp enrichment landscape were
ported from the original `scz_cell_type_enrichment` repo into
`genetics/` (script `scz_sst_hcn1_story.R` + extractor `export_for_R.py` + committed
`r_panels/` CSVs; full lineage verified — regenerated CSVs byte-identical, render pixel-matched).
Provenance + regenerate steps in `genetics/scripts/figures/PORT_NOTES.md`. This figure
**supersedes** the old matplotlib `fig4_sst_rbh_depth.png`. **Downstream TODO:** the Fig 4
Results + legend still describe the old panels — reconcile to the A–F layout (A ρ=0.56,
C ρ=0.69, D ρ=0.65; the Sst_25 gene-driver panel is demoted to a separate supp, likely
mooting the "Fig 4c driver pending regeneration" item).

## Directory contents (4 files)

| File | What it is |
|---|---|
| [`manuscript.md`](manuscript.md) | The single assembled draft — Abstract, Introduction, Results (Figs 1–4), Discussion, Methods, Figure legends, and an Outstanding-items `[TODO]` appendix. **Source of truth.** |
| [`supplementary.md`](supplementary.md) | Supplementary Information — Supplementary Methods SM1–SM5, Supplementary Figures S1–S14, Supplementary Tables T1–T7, and the applied Main-Methods rewiring record. |
| [`README_STATUS.md`](README_STATUS.md) | This file — status + the prioritised open-decision list below. |
| [`SUPPLEMENT_PLAN.md`](SUPPLEMENT_PLAN.md) | Supplement asset inventory + build/reconcile rationale (which repo asset backs each Supp Fig/Table). |

**Still Nicole's** (in her own draft, not re-written here): Introduction detail, Fig 1a–d
UMAP panels, Fig 3a–c composition, and the snRNA-seq processing/crumblr/DE methods. Her
own supplemental figures/tables slot into the S/T numbering later.

---

## Decisions & reconciliations (prioritised)

### P1 — scientific decisions that change claims

1. **GWAS version: PGC3 vs Bigdeli 2026.** The GABAergic/SST **enrichment headline
   (Fig 4a) is robust** (per-cell-type ρ = 0.95 vs PGC3). But two SST **sub-claims
   attenuate to n.s. under Bigdeli**: the upper-layer depth gradient (Fig 4b,
   r = −0.50 → −0.32, p = 0.20) and the genetics↔depletion convergence
   (r = 0.645 → 0.387, p = 0.11). **Decide before submission**; if PGC3 is kept,
   state the Bigdeli caveat (already written into Fig 4 Results + Discussion).
2. **HCN1 as a Sst_25 driver gene (Fig 4c) is not verifiable.** No committed Sst_25
   driver table exists and HCN1 is *absent* from the committed Sst_2 driver list.
   **Panel C must be regenerated** to substantiate the "481 drivers / HCN1" text, or
   the claim softened. The HCN1-*expression*-vs-sag correlation (ρ = 0.62, p = 0.010)
   **is** verified and carries the physiology link on its own.
3. **Sst_20 is a trend, not FDR-significant** (nominal P = 0.019, FDR = 0.19). Fix
   the "significant" Sst list wherever it appears (already corrected in Fig 3 text
   and Fig 5 biology; check the draft's Fig 3a wording and slide-derived lists).

### P2 — numbers / pipeline to lock

4. **Fig 3f concordance R.** Use **r = 0.50** (neuronal supertype composition,
   n = 106, P = 4.2 × 10⁻⁸) — the draft's **0.45 is a non-default permissive-QC
   value and should be corrected**. Confirm the final panel's QC gate + metric
   (composition vs density; neuronal-only vs all 120) and lock one number.
5. ~~**Depth model version.**~~ **RESOLVED 2026-08-06 — this item was backwards.**
   The **K = 50 / 3-smallest-donors-held-out** model *is* the deployed one
   (`spatial/code/pipeline/04_run_depth_prediction.py:112` calls
   `train_depth_model(merfish, K=50)` with `test_donors=None`, which
   `spatial/code/modules/depth_model.py:204` documents as holding out the 3 donors with
   fewest depth-annotated cells). Train R² = 0.93, test R² = 0.89, MAE 0.050/0.069,
   r 0.96/0.95 — as stated in `spatial/methods_writeup.md §3.3`, `spatial/README.md:114`,
   and the **Google Doc SM1**, which is correct. The **K = 100 / 6 low-CPS donor /
   R² = 0.907** variant described in `manuscript/supplementary.md` is an *uncommitted
   exploration* (`validate_depth_lowcps_cv.py`, `compare_depth_K.py`,
   `validate_depth_bottomthird_pipeline.py`) that was never wired into the pipeline.
   No committed prose needs changing; `supplementary.md` is the file that is wrong.
6. **Exemplar single-cell sags + Sst_25 donor ID** (Fig 4f,g: 0.564 / 0.00196) are
   **methods-draft only, not in any committed table**. Verify from patch-seq source
   (DANDI 000636) or present the verified supertype means (Sst_25 0.456, Sst_5 0.167).
7. **Fig 4e axis:** committed value is HCN1 **specificity**-vs-sag (ρ = 0.62); the
   methods draft says **expression**-vs-sag. Reconcile the axis label + statistic.

### P3 — nomenclature / consistency

8. **Dataset naming:** the 17/18 Ruzicka dataset is "MSSM 1" (draft) vs "MtSinai3" (talk).
9. **Prevalence:** 0.5–1% (talk) vs 1–1.5% (draft).
10. **Supertype count:** use **137** (SEA-AD MTG, the taxonomy actually mapped to)
    consistently; the talk's "127" is the neocortical-only count.
11. **"Six cortical layers":** the pipeline yields L1, L2/3, L4, L5, L6 (L2+L3 merged)
    + WM + vascular. Reword to avoid implying separate L2 and L3.
12. **Total cells:** snRNA-seq ≈ 2.3 M (fills the draft's `[x cells]`); Xenium = 1.34 M.
    Keep the two modalities' counts distinct.
13. **PGC3 N:** README says 76,755 / 243,649; `docs/analysis_approach.md` says EUR-only
    ~53 K / 77 K. State the exact EUR effective N from Trubetskoy 2022.

### P4 — data hygiene / housekeeping

14. **Stale table:** `transcriptomic/results/tables/08_meta_vs_xenium_pairs.csv` (Jun 2)
    gives r = 0.70 / 74% for DE concordance, disagreeing with the canonical figure
    inputs (r = 0.73 / 76%, n = 166). Remove/regenerate it to avoid a wrong number.
15. **External CSVs:** some validation tables (84.9% classifier accuracy, Harmony
    label-transfer) live in the sibling `~/Github/SCZ_Xenium` repo, **not committed
    here** — commit or document for reproducibility.
16. **Fig 1 scale bar + representative-donor IDs** are `[BLANK]`.
17. **Fig 4 panel-letter mapping:** the genetics↔depletion convergence is currently
    rendered as bold outlines on panels b/e, not a standalone panel — renumber if it
    gets its own panel.

### P5 — scope decisions

18. **GSEA pathway story** (`transcriptomic/`: OxPhos-in-inhibitory, cholesterol-in-
    excitatory/glia) — built, in neither slides nor draft. Decide: supplement, figure,
    or separate paper.
19. **Fig 5 placement** — RESOLVED (2026-07-21): moved into the Discussion; the paper's figures are SCZ-only (Figs 1–4).
20. **Reference list:** formalize all citations (draft uses URLs / "cite X" / `[ref]`).

---

## Key locked numbers (quick reference — all source-verified)

- **DE (Fig 2):** SST-in-Sst meta log₂FC −0.46 (FDR 0.049); Xenium −0.32 (P 0.052).
  PVALB-in-Pvalb: snRNA n.s. (−0.06, P 0.54); Xenium −0.22 (P 0.044). DE-burden vs
  abundance ρ = 0.83 (P 2.1e−6). Cross-platform DE concordance r = 0.73, 76% (n 166).
- **Composition (Fig 3):** Sst_25 snRNA β −0.26 (FDR 0.037) / Xenium logFC −0.66
  (FDR 0.028); L6b_4 +0.26 (0.022) / +1.05 (0.028). Aggregated vulnerable-Sst
  proportion P 0.016; L6b proportion P 0.009 / density P 0.005. Concordance:
  subclass neuronal r 0.70; supertype neuronal r 0.50.
- **Validation:** Xenium vs MERFISH subclass proportion r 0.85, depth r 0.96;
  classifier 84.9%; depth model held-out R² 0.907; Sst_25 F1 0.81.
- **Genetics (Fig 4, PGC3):** 44/137 SEA-AD FDR-sig (41 GABAergic); Pvalb_3
  p 3.9e−12; SST depth r −0.50 (p 0.034); genetics↔depletion r 0.645 (p 3.9e−3);
  HCN1 credible set 9 vars, lead rs10035564 PIP 0.525; HCN1 spec-vs-sag ρ 0.62 (p 0.010).

## Not attempted (Nicole's, or out of scope for this pass)
- Introduction detail, Fig 1a–d UMAP panels, Fig 3a–c composition panels, snRNA-seq
  processing/crumblr/DE methods — Nicole's sections.
- A single concatenated clean manuscript file — the section files above are the
  source of truth; can be assembled on request.
