# Figure 5 finalization plan (Sst strata DE/GSEA)

**Status 2026-08-29.** Mock-up v7 (two rows, pure SCZ: a definition + layer gutter,
b burden, c gene scatter with module + motile-cilium coloring, d pathway blocks,
e exemplar genes) is built by `transcriptomic/scripts/18_sst_strata_supplement_mockup.R`
from `transcriptomic/results/sst_strata_gsea/`. Everything SCZ-quantitative in it rests
on the **IVW shortcut**: per-gene inverse-variance combination of Nicole's per-supertype
meta-DE (`shared/snrnaseq_de/nicole_scz_snrnaseq_betas/sst_meta_de/SST_*_meta.csv`).
The IVW is fine for ranking but anticonservative for inference (supertypes share
donors), and the pooled-subclass column mixes methods (proper RE meta) with the strata
columns (IVW). Finalization = replace the IVW with real stratum-level pseudobulk
meta-DE, then re-verify every number and rerun the checks below.

---

## A. What we need from Nicole

### A1. Preferred ask: raw pseudobulk counts (enables everything below)

Per cohort (all 7), one matrix + one covariate table:

1. **Per-donor × per-supertype pseudobulk UMI counts** for the 16 Sst supertypes
   (Sst_1, 2, 3, 4, 5, 7, 9, 10, 11, 12, 13, 19, 20, 22, 23, 25; Sst Chodl excluded),
   genes × (donor × supertype) raw sums — the same aggregation that fed her
   supertype DE. Include per-group **n_cells**.
2. **Donor covariate table** per cohort: donor id, diagnosis, sex, age, PMI, plus any
   cohort-specific covariates her models used (batch, brain bank, ancestry PCs?).
3. **Label provenance**: which taxonomy/label version produced these supertype calls —
   must be the same labels the crumblr composition used, or the strata definition and
   the DE disagree about who is in which stratum.

With these we aggregate supertypes → strata per donor per cohort ourselves, which is
what makes the sensitivity analyses in §C possible (re-definitions, leave-one-out,
downsampling) without going back to her repeatedly. It also gives the per-donor
exemplar boxplots (panel-c/e companions and the g/h-style raw-data panels), replacing
the Jens placeholder.

### A2. Fallback ask: she runs the stratum DE in her pipeline

If handing out pseudobulks is impractical, we send her the strata mapping (3-line CSV:
supertype → depleted / intermediate / non_depleted; also an `all_sst` pooled group) and
ask for exactly the `sst_meta_de/`-style outputs at stratum grain:
`STRATUM_{depleted,intermediate,non_depleted,all}_meta.csv` with the same columns
(Gene, estimate, se, zval, pval, padj). Her March `Unaffected_DE.csv` /
`Vulnerable_DE_AD.csv` / `Vulnerable_DE_SCZAD.csv` files prove this mode exists —
we are asking for the same run with the SCZ-native grouping. Downside: each
sensitivity in §C becomes another ask.

### A3. Framework details to confirm with her (so our text and reruns match)

For the methods section and to keep the stratum run identical to the paper's
subclass/supertype DE, confirm:

- **Per-cohort model**: dream / limma-voom? Exact formula (diagnosis + sex + age +
  PMI + …?); the observational unit (pseudobulk per donor — no repeated measures?);
  normalization (TMM?) and gene-filtering rule (filterByExpr? min counts / min donors).
- **Meta-analysis**: metafor random-effects across the 7 cohorts (the subclass file
  carries k, tau2, I2 — confirm estimator, e.g., REML vs DL) and how genes present in
  only some cohorts are handled (min k?).
- **FDR universe**: padj computed within cell type across genes (assumed) — confirm.
- **QC gates**: which nuclei enter the pseudobulks (her pipeline's QC + doublet rules),
  and whether MSSM1/2 donor overlap handling (the 12 shared donors) affects strata DE.

---

## B. The swap-in (what gets rerun when the data lands)

1. Rebuild stratum signatures: replace the IVW step in
   `transcriptomic/scripts/17_sst_strata_gsea.R` with the pseudobulk meta-DE
   (per-cohort stratum DE → metafor), or ingest Nicole's stratum meta CSVs directly.
   Keep the IVW output long enough to run the **IVW-vs-pseudobulk concordance check**
   (scatter of z's; NES correlation of the GSEA) — this validates every mock-up-era
   conclusion and is one supplement-ready panel if a reviewer asks about the method.
2. Rerun fgsea (house settings) → all bracketed numbers in
   `manuscript/figure5_section_DRAFT.md` and the v7 figure update:
   burden counts ([112]/[57]/[0]/[50]), NES values and stars in panel d, gene z's in
   panels c/e, the shared/graded categories, the motile-cilium up-in-both values.
3. Now the subclass reference column and the strata columns share one method — remove
   the mixed-methods caveat from the legend.
4. Rerun the Xenium cross-platform checks against the new snRNA side
   (stratum concordance r's; SST/VGF; module panel-coverage numbers are unchanged).
5. Rebuild the exemplar boxplots (`26_concordance_boxplots.R` pattern) with per-donor
   7-cohort pseudobulks. **Flagged verify**: RPL36 and NDUFS8 should now separate
   CTRL vs SCZ in the depleted stratum (meta z −3.8 / −3.3). The Jens placeholder
   could not show this (its depleted labels are the unreliable ones) — expected, not
   a failure. If the pseudobulks ALSO fail to show it, that is a real finding to chase.

---

## C. Sensitivity analyses (agreed across this session)

| # | Analysis | Purpose | Method | Destination |
|---|---|---|---|---|
| S1 | Strata definition: FDR<0.05 (4 depleted, Sst_20 → intermediate) vs FDR<0.20 (current 5) | Sst_20 is borderline (padj 0.19); co-author gate | Rerun stratum DE+GSEA under both; report headline-module NES both ways | supplement sentence or reviewer reserve |
| S2 | Cell-count-matched downsampling | Non-depleted has ~40–50% fewer cells everywhere (Xenium 6.2k vs 10.6k; Jens 5.9k vs 11.2k; A9 16.3k vs 27.7k; meta SE proxy 0.0595 vs 0.046); its zero-burden is partly power | Downsample depleted/intermediate cells (or donors' cells) to non-depleted counts before pseudobulking; re-run GSEA; note the two power-safe facts (intermediate best-powered yet half the burden; translation sign-flip is directional) | supplement or legend sentence |
| S3 | Leave-one-cohort-out | Mirrors the paper's S6 robustness idiom | Drop each cohort from the stratum meta; track headline-module NES + burden gradient | supplement panel |
| S4 | Age dependence (Kiss 2026) | SCZ Sst depletion is age-dependent (<~70y); strata state effects may be too | Stratum DE with age interaction, or split younger/older donors; check gradient + modules | supplement or internal check |
| S5 | Leave-Sst_22-out | Sst_22's DE profile is noise-like (near-zero correlation with every other supertype) | Rerun depleted stratum without Sst_22 | internal check; footnote if it matters |
| S6 | Independent-cohort replication (DONE, refresh) | Jens (83 donors; scripts 20a–c): burden gradient, shared synaptic/SST/VGF, OxPhos replicate; translation localization unresolved (depleted-label porosity; signal sits in intermediate/affected) | Refresh comparison against the pseudobulk side; optionally add GSE158516 (already brisc-mapped: `full_mapped_gse.parquet`) as a second independent cohort | companion supplement |
| S7 | Re-annotation consistency (S7-analog, optional) | Strata rest on labels; abundance already has the DE-gene-excluded re-annotation control | If reviewers press: stratum DE using re-annotated labels (needs re-annotated pseudobulks — a Nicole ask) | reviewer reserve |
| S8 | IVW-vs-pseudobulk concordance | Validates the mock-up-era analyses | §B1 | methods sentence / reviewer reserve |
| S9 | Depth-dichotomy grouping (upper L1–3 vs lower L4–6, Xenium median depth, L2-3/L4 boundary 0.470) | (a) exogenous grouping — immune to the define-strata-from-own-composition circularity worry; (b) power-balanced arms: upper 57.7k cells / median 95 per donor vs lower 43.9k / 80.5, versus non-depleted's starved 13.6k / 23 per donor. Nesting is clean (upper = 5 depleted + Sst_23/11; lower = 5 non-depleted + Sst_9/12/13/19) so it coarsens, not crosscuts, the depletion axis | Rerun 27a/27b with two depth arms; expect the lower arm mildly negative (4 intermediates inside) — say so | sensitivity line / supp panel; NOT a replacement for the depletion strata |

Held in reserve (not in the main figure, per the v7 decision): the SCZ↔AD material —
divergence scatters, concordance classes (96% sign-concordance among jointly
|effect|>2 genes), A9 CPS strata GSEA + early-epoch version, exemplar boxplots
(scripts 21–26 outputs under `seaad_a9/`). One cautious discussion clause at most.

---

## D. Figure-finalization pass (cosmetic, after numbers stabilize)

- Panel c: decide the up-in-both gene labels (currently WDR19, CEP126; optionally add
  INSYN2B/RPS6KA5/GRIN2A from the shared-up list); keep cilium described neutrally.
- Panel d/e strips: final wording ("Graded: oxidative phosphorylation" is tight
  vertically — consider "Graded: OxPhos").
- Hand-chosen short pathway display names; final palette check against Fig 1–4 usage
  (#D55E00 doubles as depleted-stratum color and NES-positive — acceptable across
  panels, but confirm with Nicole's Fig 2 conventions).
- Legend rewrite from the draft in `figure5_section_DRAFT.md` (+ new panel-c cilium
  sentence, S2 power caveat sentence, layer-gutter provenance:
  `proposed_layer_boundaries.csv`).

### S11 — diagnosis × stratum interaction — DONE 2026-08-31, THE formal backstop

`transcriptomic/scripts/29_interaction_model.R` → `pseudobulk/interaction_{percohort,
meta,gsea}.csv`. Per cohort: joint-TMM donor×stratum pseudobulks, voom +
duplicateCorrelation(donor), ~ dx*stratum (ref = non_depleted) and ~ dx*score (trend);
metafor REML k≥5. RESULT: interaction GSEA (dep-vs-non coefficient) is dominated by
translation sets at extreme significance — SRP cotranslational NES −2.69 padj 4.4e-10,
elongation −2.69 padj 2.3e-9, cytosolic ribosome −2.55, initiation −2.39, Reactome
translation −2.02 padj 4e-6 — plus Hallmark OxPhos −2.0**; SYNAPTIC sets show NO
interaction (+1.2..+1.6 ns) and neither do SST (z = −0.84) or VGF (−0.55) → the
shared-vs-graded dissociation as tested quantities. Gene-level interaction FDR<0.05:
7 genes (dep-vs-non) / 5 (trend; top: CELF2 +5.1 padj .004, GRIK4 −4.9 padj .005,
CBLN4 −4.1 padj .088, COX5B −3.6); NDUFS8 z −2.7 / RPL36 −2.3 nominal. Within-donor
design ⇒ donor-level confounds (medication, PMI, age) cannot generate these
interactions. Results-draft sentence added. Consider quoting the two headline
interaction-GSEA padj values in the legend or results.

### S11+. Ordinal verdict & step localization (2026-08-31)

The ordinal linear-trend interaction (dx × score, 0/1/2) was run in script 29 (model B):
concordant but WEAKER than the categorical contrast (5 FDR genes; GSEA retains only
Hallmark OxPhos −2.1**). Reason, now demonstrated by GSEA on the intermediate-vs-non
coefficient: the two modules step at DIFFERENT boundaries — TRANSLATION's step is at
the intermediate/non-depleted boundary (int-vs-non interaction: elongation NES −2.78
padj 3.7e-11, ribosomal subunit −2.14 padj 5.7e-6 → dep ≈ int ≫ non), whereas OXPHOS's
step is at the depleted/intermediate boundary (int-vs-non ~null; dep-vs-non −2.0** →
dep ≪ int ≈ non). A single linear ordinal slope blurs two different step functions;
the categorical model contains the ordinal contrast as a nested special case and is
strictly more informative. REPORTING: categorical primary; one results/discussion
sentence on the step localization (translation extends through intermediate; OxPhos is
depleted-exclusive); linear-trend mentioned as a concordant sensitivity. OPTIONAL S12
(reviewer reserve): fully continuous dx × crumblr-β interaction at supertype
resolution — removes the FDR<0.20 binning decision entirely; costs supertype-pseudobulk
sparsity + β measurement error.

### DECISION 2026-08-31: no violin/boxplot panels in Figure 5 — full stop

Shreejoy: since no individual module gene is FDR-significant, gene-level violin
plots do not enter the figure; the presentation stays heatmap-centric (panels d/e
carry the gene-level content as z-value heatmaps, explicitly starless/illustrative).
The violin figures (figS_stratum_violin_exemplars, figS_interaction_gene_violins,
figS_gene_concordance_v3_boxplots) are retained as INTERNAL QC / reviewer-reserve
artifacts only. Figure 5 = a (definition + layers), b (burden), c (gene scatter,
family-colored), d (pathway NES blocks), e (exemplar-gene z heatmap).

### Figure v12 changes (2026-08-31)

a: strata named in-panel (colored, at top of each x-region); Spearman rho and legend
removed. b: Figure 2 two-tier significance (FDR<0.10 light / FDR<0.05 dark, Fig 2
colors + legend wording); counts 49/127/56/0 down and 7/18/5/0 up at FDR<0.10.
c: synaptic drawn small/translucent so the graded modules read; module centroids tried
and REMOVED (translation and OxPhos centroids nearly coincide, so they cluttered
without separating); labelled genes now drawn from panel e's exemplar lists via a
shared SEL/HL selection, ggrepel with segments. e: gene names rendered as a data-driven
text layer (ggtext unavailable) so the 8 c-labelled genes are BOLD; axis text suppressed.

UP-BLOCK REPLACEMENT: motile cilium -> K63-linked deubiquitination. The full up-set
inventory (padj<0.05, depleted stratum) is only ~3 distinct themes: (i) ciliary/axonemal
(motile cilium, cilium movement, cilium/flagellum motility, 9+2 cilium, AND the
gamete/reproduction sets - 40% leading-edge overlap with cilium, i.e. the same axonemal
genes, not independent options); (ii) K63-linked deubiquitination (NES 2.15, padj
0.0057, up in all three strata; LE = 16 bona fide DUBs incl. USP8, STAMBPL1, BRCC3,
ATXN3, OTUD1, VCP, PSMD14); (iii) TGF-beta/activin receptor kinase (NES 2.16, padj
0.015; LE = ACVR2B, TGFBR3, ACVR2A, ACVR1B, ACVR1C - only 5 genes, set size 12; also up
in all_sst, and echoes the BMP/TGF-beta theme in the project's earlier pathway work),
with anterior-posterior pattern specification sharing its ACVR genes; plus
GOCC_EXTERNAL_SIDE_OF_PLASMA_MEMBRANE (PECAM1, KLRD1, CD36 - looks non-neuronal,
avoid). CAVEATS ON THE CHOICE: (1) GOBP_REGULATION_OF_PROTEIN_UBIQUITINATION is DOWN
(-1.8, padj 0.0029), so the ubiquitin system is NOT uniformly up - the honest reading is
K63 deubiquitination up while ubiquitination regulation is down (reduced degradation
flux); do not label the block "ubiquitin up" in text without this nuance. (2) the second
row shown in d (GOMF_UBIQUITIN_CONJUGATING_ENZYME_BINDING, NES 1.72, padj 0.21, no star)
was added FOR LAYOUT - a one-row block cannot fit a vertical strip label. Revisit: either
keep (honest, unstarred) or drop it and shorten the strip label.

### CELL-COUNT CONFOUND: RESOLVED (2026-09-01) — definitive test passed

Cell-level Sst data arrived (`transcriptomic/data/stratum_sst_cells_export/`,
101,566 nuclei, 7 cohorts, reconciles exactly to the pseudobulks; the two donor
"FAIL" flags in `_sst_reconciliation.csv` are definitional -- 9 donors have zero
Sst nuclei and so cannot appear in a cell-level file, verified against the groups
CSVs). Scripts: `fig5/supp/subsample_cells.py` + `sensitivity_cellsubsample.R`.

**The test.** Subsample control donors' NUCLEI down to the case distribution
within each cohort and stratum (quantile matched, donors first restricted to
>= 10 nuclei so nothing drops out), rebuild the pseudobulks from cells, rerun the
unchanged pipeline. 5 independent draws.

**Matching worked.** The imbalance that motivated the whole exercise is gone:
depleted -30% (P = 1.6e-5) -> **-0.4% (P = 0.96)**; intermediate -18% -> +2%
(P = 0.80); non-depleted -4% -> +5% (P = 0.43).

**The result survives intact.**

| | unmatched | cell-matched, 5 draws |
|---|---|---|
| depleted, sets at FDR<0.10 | 130 | 101, 90, 117, 100, 98 |
| intermediate | 64 | 65, 39, 43, 51, 40 |
| non-depleted | 0 | 1, 0, 0, 0, 0 |
| translation, depleted | NES -2.6 | median -2.42, **20/20 tests significant** (down to 2e-11) |
| translation, intermediate | -2.2 | median -2.02, **20/20 significant** |
| translation, non-depleted | -0.7, FDR 1 | median -0.82, **0/20**, FDR 1 |
| oxphos, depleted | -1.9 to -2.0 | median -1.94, **10/10 significant** |
| oxphos, intermediate / non-depleted | n.s. | 0/10 and 0/10 |
| K63 deubiquitination, depleted | +2.15 | ~+2.05, significant in 5/5 |

The modest burden drop (130 -> ~100) is the expected power cost of discarding
~30% of control nuclei; the graded STRUCTURE is unchanged.

**Where this leaves the four tests.** Read-depth thinning and cell subsampling --
the two that remove the technical imbalance without removing the disease contrast
-- both preserve the signal. Covariate adjustment for log(nuclei) collapses it,
but that adjusts for a quantity that IS the disease effect (Fig 3), i.e. a
mediator. Donor matching also collapses it but discards 40% of donors and whole
cohorts. **The confound is excluded as an explanation**; cite the cell-matched
analysis, not the covariate model.

TEXT TO ADD (results or methods): one sentence noting that subsampling control
nuclei to the case distribution and rebuilding the pseudobulks leaves the graded
suppression unchanged (5 draws), so it is not a consequence of the depleted group
containing fewer cells in cases.

### Superseded first pass — kept for the reasoning

Script `fig5/supp/sensitivity_cellcount.R`; outputs `supp/sensitivity_cellcount*.csv`.

**The concern.** The groups are DEFINED by SCZ having fewer cells, so the
case-control nuclei imbalance is graded exactly like the result:
depleted -30% (P = 1.6e-5), intermediate -18% (P = 0.005), non-depleted -4% (n.s.).
Library size is worse: 36% lower in cases in the depleted group. The within-donor
interaction design does NOT remove this -- it cancels donor-level factors, but the
count difference is a donor-BY-group factor.

**The artifact is real.** Within CONTROL donors only (no disease involved),
log nuclei count predicts module score: translation beta +0.33 (P = 1.9e-14),
oxphos +0.31 (P = 2.6e-15), synaptic +0.22 (P = 4.9e-18). Naive projection of
that slope onto the -0.35 log-unit case-control gap would account for roughly
half of the observed depleted-group module effects.

**Three tests, which disagree:**
| Test | Depleted burden (FDR<0.10) | Verdict |
|---|---|---|
| Unadjusted | 130 | baseline |
| 1. Adjust for log(nuclei) | **1** | signal collapses |
| 2. Count-matched donors (110 pairs) | **8** (non-depleted goes 0 -> 39) | collapses, but unstable |
| 3. Depth-matched, controls thinned to the case library-size distribution | **109** (102/118 on other draws) | signal SURVIVES |

**Which to believe.** Test 1 adjusts for a quantity that IS the disease effect
(Fig 3), i.e. a mediator, so it asks a different and much harder question; it
cannot distinguish technical confounding from over-adjustment. Test 2 discards
~40% of donors, drops whole cohorts (6/5/4 of 7), needed k>=3 instead of k>=5,
and its non-depleted result (0 -> 39 sets) is not interpretable -- matching in
the depleted group forces pairing of atypically high-count cases with
atypically low-count controls. **Test 3 is the best-targeted**: it equalises the
technical quantity exactly, keeps every donor, and needs no confounder/mediator
assumption. Translation (ribosomal subunit NES -2.56, FDR 1.9e-12) and OxPhos
(-2.03, FDR 4.8e-4) are essentially unchanged, and the non-depleted group stays
empty. Stable across 3 draws.

**Still untested.** Thinning corrects sequencing depth, not cell number. In a
joint within-control model the two separate only partly (depleted group:
nuclei beta +0.217 P = 0.018; library beta +0.104 P = 0.18), so a cell-number
component remains uncorrected. THE DEFINITIVE TEST needs cell-level Sst matrices
(Shreejoy has offered them): subsample control cells down to the case
distribution, rebuild pseudobulks, rerun over ~10 draws.

**A second confound, checked and clean.** Supertype mixture within each pooled
group does not shift with diagnosis (no member supertype of the depleted pool
changes its within-pool fraction, all P > 0.07), so the pooled case pseudobulk
is not a different mixture of supertypes than the control one. The single
nominal hit is Sst_7 in the NON-depleted pool (P = 0.04, would not survive
correction), which cannot explain a depleted-group signal.

### Pipeline rerun VERIFIED end-to-end (2026-09-01)

Full rerun of the consolidated pipeline against a pre-refactor snapshot:
stratum_meta_de (45,107 rows), gsea_all_signatures (25,899), interaction_meta
(31,938), interaction_gsea (12,626) and module_scores_donor (3,444) are all
bit-identical (worst numeric diff 0.00e+00), and the rendered figure has the same
md5 as before the refactor. `09_verify.R` PASSES 0 failures.

Three text corrections the harness forced, all now applied to the draft:
- **Nuclei/median counts were computed over the wrong population.** The old
  numbers (44,950 / 43,000 / 13,616 nuclei; medians 69/77/21) summed over ALL
  donor-strata, while the donor counts beside them (395/411/342) applied the
  >=10-nucleus filter. Corrected to the population that actually enters the
  analysis: **44,744 / 42,814 / 13,157 nuclei, medians 85 / 90 / 31**.
- **"synaptic FDR >= 0.31" -> ">= 0.21".** 0.31 was the minimum over two synaptic
  gene sets; the BLOCKS table defines four, and GOBP_NEUROTRANSMITTER_SECRETION
  sits at 0.211. Still comfortably null, but the honest bound.
- **Xenium panel coverage "1 of 130" -> "2 of 130" translation genes** (oxphos
  stays 2 of 114). The pseudobulk signatures cover more panel genes than the IVW
  shortcut did (173 vs 144), so the overlap grew. SST/VGF z values unchanged.
- Also: "oxphos all FDR <= 0.0022" -> "<= 0.0023" (true max is 2.22e-3).

### Code consolidation (2026-09-01) — `transcriptomic/scripts/fig5/`

The Figure 5 code is now a single numbered pipeline with a shared helper module;
see `transcriptomic/scripts/fig5/README.md` for the DAG and runtimes.
01_strata -> 02_stratum_de -> 03_stratum_gsea -> 04_donor_pseudobulks ->
05_interaction -> 06_module_scores -> 07_xenium_stratum -> 08_figure5 -> 09_verify.
1,841 lines across 13 ad-hoc scripts became 1,154 lines across 10, and the
duplication that mattered is gone: module definitions were re-derived in NINE
scripts, strata membership hard-coded in THREE, MSigDB rebuilt in NINE. All now
come from `_common.R` (BLOCKS table -> modules(); strata_definition.csv ->
stratum_of(); cached msigdb_sets()).

VERIFIED: steps 02 and 03 rerun bit-identically against a pre-refactor snapshot
(maxAbsDiff = 0 on every numeric column of stratum_meta_de.csv, n = 45,107 rows,
and gsea_all_signatures.csv, n = 25,899 rows).

Three defects found and fixed during consolidation:
1. **Xenium provenance hole.** `xenium_stratum_concordance.csv` supplied the
   SST/VGF cross-platform z values and panel-coverage counts quoted in the
   results text, but had NO producer script (written ad hoc in an earlier
   session). Recipe recovered by testing candidates against the orphaned file --
   unweighted mean logFC across a stratum's supertypes, per-supertype z combined
   by Stouffer -- reproducing it to 0.0000 on both columns. Now
   `07_xenium_stratum.R`.
2. **S8 self-comparison.** The IVW-vs-pseudobulk sensitivity read the top-level
   CSVs, which had since been overwritten with pseudobulk output, so a rerun
   would have compared the pseudobulk result against itself and reported r = 1.00.
   The top-level copies are retired to `archive_exploratory/`; the genuine IVW
   results remain in `archive_ivw/`. RECOMMENDATION: drop this sensitivity from
   the paper entirely -- it compares our interim shortcut with the real analysis
   and is of no interest to a reader.
3. **Duplicate `27b` filenames.** Two different scripts were both named 27b; the
   swap-in variant was a superseded attempt (it would have written a signature
   label, `all_sst_pb`, that appears nowhere in the canonical output). Archived.

Dead code archived to `transcriptomic/scripts/archive/fig5_superseded/`: the
concordance drafts v1-v3 (24/25/26, the dropped SCZ-AD last row), the interaction
gene violins (31, violins excluded from the figure by decision), the violin half
of 28 (its pseudobulk-building half became 04), and the original exploration (17,
whose strata block became 01).

### SUPPLEMENT: FINAL — one figure, S10 (2026-09-01)

`fig5/supp/figS_strata_controls.R` -> `figS_strata_controls.(png|pdf)`, legend in
`manuscript/figure5_supp_legends_2026-09-01.md`. Three panels:
(a) case-control nuclei difference per group, observed vs cell-matched
(-29.5/-17.9/-4.0% -> -0.4/+1.8/+4.9%); (b) burden per group, unmatched bar +
5 matched draws; (c) interaction NES for the Fig 5d sets plus the two named in
the text.

Merged from two earlier drafts at Shreejoy's request; the top-15 interaction
ranking panel was dropped. NOTE what that costs: the word "concentrated" in the
interaction sentence is now supported only by the pre-selected Fig 5d sets, and
the disclosure that several top interaction sets are ribosome-heavy under
non-obvious names (selenoamino acid metabolism 69%, influenza infection 52%,
SLIT/ROBO 52%, GCN2 80%) is no longer anywhere in the paper. If a GSEA
supplementary table is included, a reader will meet those names there.
Predecessors archived in `scripts/archive/fig5_superseded/`.

### SUPPLEMENT SET, scoped to the current results text (2026-09-01, superseded)

Audit of the section as written. Only ONE supplement is explicitly cited, and one
more analysis is reported with numbers but no figure.

| | Called out in text? | Status |
|---|---|---|
| **S1 cell-count control** — the `(Fig. SXXX)` after the downsampling sentence | YES, the only explicit citation | **BUILT**: `supp/figS_cellcount_control.R` -> `figS_cellcount_control.(png\|pdf)`. Panels: (a) case-control nuclei difference per group, observed vs cell-matched (-29.5%/-17.9%/-4.0% -> -0.4%/+1.8%/+4.9%); (b) burden per group, unmatched bar + 5 matched draws; (c) the four blocks, unmatched vs 5 draws |
| **S2 interaction model** — "Formal models testing the interaction ... FDR = 4.4e-10 ... 0.053" | NO figure cited, though it is a whole analysis | NEEDS a citation and a focused figure. `supp/volcanoes.R` covers it but bundles per-stratum volcanoes the text never mentions |
| Xenium cross-platform | not mentioned any more | BUILT, now RESERVE |
| Module scores | not mentioned any more | BUILT, now RESERVE (but see below) |
| Baseline enrichment | not mentioned | BUILT, Figure 4 supplement if used |
| Jens replication, AD contrast | not mentioned | RESERVE |

Note the text's own threshold makes the oxphos interaction significant: at
FDR < 0.10 (the paper's GSEA criterion) the cited electron-transport-chain
interaction at FDR = 0.053 clears the bar, so "supported these findings" is if
anything understated. But the two oxphos sets PLOTTED in Fig 5d
(GOBP_OXIDATIVE_PHOSPHORYLATION 0.13, REACTOME_RESPIRATORY_ELECTRON_TRANSPORT
0.104) do NOT clear it -- the cited set is a third one not in the figure. A
reader checking panel d against that sentence will not find the cited set.

### Supplementary figure options (assessed 2026-09-01, superseded by the audit above)

Figure 5 asserts: (1) the strata are a real axis, (2) burden is graded,
(3) two components -- shared vs graded, (4) the dissociation is a formal
within-donor interaction, (5) it is not an artifact.

| # | Supplement | Supports | State | Recommendation |
|---|---|---|---|---|
| S1 | **Module scores** — forest of the diagnosis effect on each module score per stratum, all 7 cohorts shown, plus the paired within-donor interaction | 3, 4 | BUILT (`supp/module_scores_fig.R`) | **Include.** The only supplement that converts the aggregate claim into a per-donor, per-cohort measurement; without it the "no single gene is significant" objection stands unanswered |
| S2 | **Volcanoes** — 3 per-stratum + 3 interaction | 2, 4 | BUILT (`supp/volcanoes.R`) | **Include.** Discloses the gene level honestly (few hits, not tracking the gradient), which is what makes the "coordinated shifts" sentence credible rather than evasive |
| S3 | **Strata-definition sensitivity** — redefine strata at FDR<0.05 (drops Sst_20) and as depleted-vs-rest, show the graded pattern survives | 1, 2, 5 | NOT BUILT (~20 min compute per variant: rerun 02+03) | **Build.** The FDR<0.20 cut admitting Sst_20 is the likeliest reviewer target. NOTE: the leave-one-cohort-out half originally planned here is largely redundant with S1, whose per-cohort dots answer "is one dataset driving this?" more directly than seven leave-one-out metas would — build the strata variants, skip LOCO |
| S4 | **Xenium cross-platform** — stratum-level snRNA-seq vs Xenium z for the ~144 panel genes; panel-coverage bar | 5 | CSV built (step 07), figure NOT BUILT (~1 h) | **Build.** Supports the shared component AND discloses that the graded modules are untestable on a 300-gene panel — currently only a sentence |
| S5 | **Baseline enrichment** — ORA of the Fig 4g 580-gene signature + module baseline levels | "why vulnerable" (Etienne) | BUILT (`supp/baseline_enrichment*.R`) | **Include, but as a Figure 4 supplement** — it characterises neurotypical identity, not the disease state |
| S6 | Independent cohort (Jens) | 5 | BUILT, needs path updates | **Reserve.** 4/5 depleted supertypes map poorly there, so stratum assignment is porous; replicates the gradient and OxPhos but not the translation localisation. Requires disclosure |
| S7 | AD contrast (SEA-AD A9 CPS) | specificity | BUILT | **Reserve.** Strong answer if a reviewer asks whether this is generic neurodegeneration (the modules move in the SPARED strata along AD pseudo-progression) |
| S8 | IVW vs pseudobulk framework check | — | ARCHIVED | **Drop** (see defect 2) |

### Baseline-signature enrichment run (2026-09-01) — the "why vulnerable" supplement

Scripts `transcriptomic/scripts/34_baseline_marker_enrichment.R` + `34b_...fig.R`;
outputs `transcriptomic/results/sst_strata_gsea/baseline_marker_enrichment/`
(ora_580_baseline_signature.csv, gsea_baseline_log2fc.csv,
modules_baseline_log2fc.csv, figS_baseline_marker_enrichment.png).
Input = Fig 4g table `genetics/results/figures/r_panels/panel_volcano_vulnerable_vs_notdepleted.csv`
(580 genes at Wilcoxon FDR<0.05 & |log2FC|>0.25 verified; HCN1 +0.43 / CALB1 +0.55 anchors reproduce;
background = 11,473 genes detected in >=10% of nuclei in either group; 579/580 in background).
RESULT: higher-in-depleted = synaptic/connectivity identity (synapse organization
40/430 padj 8e-10; synaptic membrane; trans-synaptic signaling; cell-cell adhesion;
Ca2+ binding 30/310; GPCR; only positive GSEA set = pos. reg. of GABAergic synaptic
transmission padj 0.033). Higher-in-SPARED = ribosome/translation overwhelmingly
(elongation 30/89, 12.8x, padj 7.5e-22; cytosolic ribosome; SRP; initiation - same
~30 RP genes across sets). Module cross-check: translation module baseline median
log2FC -0.11 (7% positive, p 8.4e-28), OxPhos -0.03 (14% positive, p 2.4e-8),
synaptic ~0. INTERPRETATION: depleted types are synaptically over-invested and run
LOWER baseline translation/OxPhos than spared types; the disease-graded suppression
(Fig 5) lands on the programs they are lightest on ("thin margin"), and runs AGAINST
the detection-power gradient. Caveats: 3 neurotypical donors, cell-level Wilcoxon,
CP10K compositionality, identity-vs-vulnerability conflation.

### Results drafts v2 written (2026-09-01)

`manuscript/figure5_section_DRAFT.md` fully rewritten around the pseudobulk
numbers: LONG (~1,000 w, full refs) / PARED (~430 w, rationale-led) / SHORT
(~290 w, Doc-density paste candidate) + Discussion paragraph (goes between the
AD and L6b paragraphs; companion edit: Fig. 4 opener "Finally"->"Next") +
updated legend + methods stub. Every number machine-verified by
`transcriptomic/scripts/33_verify_fig5_text_numbers.R` and a 40-assertion
checker (0 failures). Key corrections vs the old draft: interaction GSEA has
only dep-vs-non + linear-trend contrasts (int-vs-non padj claims dropped);
module scores show all three modules modestly down even in non-depleted
(synaptic paired interaction p=1.3e-4), so the dissociation is phrased as
"exceeds the transcriptome-wide background only for translation/OxPhos";
stress-named Reactome sets are 55-83%% RP genes -> no ISR/eIF2 language for our
own results. Panel coverage updated to 1/130 and 2/114.

### Figure v14 (2026-09-01)

a: strata palette changed to purple / lilac / green, deliberately outside the
blue-orange (down/up) and magenta/teal (module) families used in b-e — this also
removes the collision where orange meant "depleted" in a but "up in SCZ" in b.
Panel a KEEPS the "(5)/(6)" group sizes because it is the definition panel; b, d, e
now use bare labels (All Sst / Depleted / Intermediate / Non-depleted).
b: count axis labelled "Significant gene sets (count)"; legend moved to the
top-right corner (the Non-depleted column is empty, so nothing is occluded).
c/e: EXEMPLAR RULE CHANGED for the two graded blocks — genes are now the most
DISCORDANT core members of each GO group (regex-restricted to cytosolic ribosomal
proteins + EIF3 for translation, and ETC complex subunits for OxPhos; require
depleted z < -1.2, then rank by non_depleted - depleted). This (i) puts them in the
upper-left quadrant of c, and (ii) fixes a real annotation problem: the previous
picks for "cytosolic translation" were mitoribosomal (MRPL33/MRPL55/MRPS25), and
the previous OxPhos picks were TCA/fatty-acid genes (IDH2, GOT1, ACADVL) rather
than ETC subunits. New picks: RPL10/RPL36/RPS21/RPL14/EIF3G/RPS25/EIF3I/RPL37A and
NDUFS8/CYCS/UQCRH/COX5B/NDUFA10/CYC1/NDUFA2/COX7B. Bolded in e / labelled in c:
SST, VGF, STAMBPL1, USP8, RPL36, EIF3G, NDUFS8, COX5B.
e: left gutter trimmed (expansion 2.7 -> 1.32 discrete units).

### S13. Module-score tests (2026-08-31) — the answer to "no FDR genes"

`transcriptomic/scripts/32_module_scores.R` → `pseudobulk/module_scores_donor.csv`.
Each module as ONE variable per donor-stratum (mean of within-cohort z-scored logCPM);
single tests, no gene-level multiplicity. (1) Score ~ dx per stratum, meta over 7
cohorts: OxPhos depleted −0.234 SD, p = 1.9e-4, NEGATIVE IN 7/7 COHORTS (translation
−0.216, p = 2.5e-4, 7/7; synaptic −0.227, p = 7e-8). All modules also nominally down
in intermediate and non-depleted with smaller magnitudes (graded). (2) Within-donor
paired interaction (depleted − non-depleted score ~ dx): OxPhos −0.105 p = .023 (6/7),
translation −0.091 p = .028, synaptic −0.124 p = 1.3e-4 (7/7). REMEDY for the
no-FDR-genes discomfort: quote GSEA + sign-bias binomial + module-score single test
together; same aggregate epistemology as Fig 4a's MAGMA enrichment.
NUANCE FOUND: at absolute score level the SYNAPTIC module's case-control gap is ALSO
graded (paired p = 1.3e-4) even though its genes show no competitive GSEA interaction
(rank-level null; background interaction median z = 0.01, so GSEA background is fair) —
i.e., synaptic interaction is a small uniform mean shift plus heterogeneous crossovers,
translation/OxPhos are uniformly sign-biased (63%/58% negative). WORDING ADJUSTMENT:
say the synaptic program is "downregulated across all three strata (largest in
depleted)" rather than "comparable effect size"; the graded blocks are distinguished
by per-gene coordination and competitive specificity, not by being the only modules
whose absolute magnitude grades.

### S11++. Interaction-gene violins (2026-08-31)

`pseudobulk/figS_interaction_gene_violins.png` (script 31): all 8 FDR interaction
genes are CROSSOVER interactions — diagnosis effect reverses sign across strata
(CELF2 +2.0/−0.1/−2.7; GRIK4 −3.5/−1.3/+1.7; FHOD3 −3.2/−1.0/+3.1), unlike the
modules' same-sign amplitude gradient. Crossovers are the easiest interactions to
detect at gene level, which is why these top the FDR list while module genes stay
sub-threshold. Motif: CELF2 and QKI (both neuronal splicing regulators) crossover
in the same orientation. Frame as gene-level complement, not the module story.

### D++. Violin exemplars (2026-08-31) — DONE, replaces the Jens placeholder

`transcriptomic/scripts/28_stratum_violin_exemplars.R` →
`pseudobulk/figS_stratum_violin_exemplars.png`: control-vs-SCZ violins by stratum
from the REAL 7-cohort per-donor pseudobulks (469 donors; values centered within
cohort × stratum; stratum meta z annotated). Shared row: SST, VGF (shift in all
strata); divergent: RPL36, NDUFS8 (shift in depleted only; NDUFS8 flat elsewhere).
Single-gene shifts are visually modest — the annotation carries the inference;
candidate Fig 5 third row (g) or companion supplement. Side product: per-donor
stratum pseudobulk matrices persisted under `pseudobulk/donor_stratum/` for
S2/S9/S10 and future panels.

### D+. Panel-e significance note (2026-08-31, pseudobulk meta)

Within-stratum gene-level FDR is sparse for the panel-e exemplars: only SST
(padj = 0.061) and RASGRF1 (0.079), both in the INTERMEDIATE stratum, pass
FDR < 0.10; in the depleted stratum the best is SST at 0.18 (RPL36 0.38,
NDUFS8 0.27). The pooled all-Sst analysis regains gene-level power: NMU (0.018),
ATP2B2 (0.040), TIMM13 (0.057), PICK1 (0.069), SLC4A8 (0.079), IDH2 (0.095).
Consequences: (i) panel e stays a z-value heatmap with NO significance stars;
(ii) legend sentence: "genes shown as illustrative leading-edge exemplars;
within-stratum gene-level FDR is limited by per-stratum power (Supp Table X) —
the stratum-level claims are carried by the gene-set statistics in (b, d)";
(iii) results text quotes z's, never within-stratum gene FDR (already the case).

## E. Text deliverables

1. Results section: `manuscript/figure5_section_DRAFT.md` SHORT version — fill
   brackets from the pseudobulk run; sharpen the opening to the independent-or-coupled
   framing; drop the panel-f sentence (row removed).
2. Discussion paragraph (drafted in the same file): keep; the cross-disorder sentence
   softens to "substantially shared response with module-level exceptions" or is cut.
3. Methods (~400–600 words): strata definition (crumblr FDR<0.20 on the 7-cohort
   composition + depth annotation), stratum pseudobulk DE (Nicole's framework, §A3),
   fgsea settings, module/leading-edge definitions, Xenium check, Jens cohort.
4. Verify [Table T2 — BDNF down across excitatory subclasses] programmatically before
   the discussion sentence stands.

## F. Decision gates (co-authors)

0. **Panel a vs Figure 3f — RESOLVED 2026-08-31: option A (keep the scatter).**
   Option C (compact key) was built and REJECTED by Shreejoy: same footprint, less
   information. Panel a therefore stays the depth-vs-crumblr-estimate scatter with
   stratum colors (STRAT_COLS), error bars, layer gutter, FDR<0.20 boundary and the
   rho annotation; the duplication with Fig 3f is accepted deliberately. Mitigations
   to apply at legend-writing time: state "supertypes and depths as in Fig. 3f" so the
   overlap reads as intentional; consider dropping the rho from ONE of the two panels
   so the statistic is reported once. STILL TO DO: add the layer gutter to Fig 3f
   (Shreejoy), and reconcile MEAN (Fig 3f legend) vs MEDIAN (5a) cortical depth — pick
   one for both panels. Original framing of the problem: Fig 3f is already "Mean
   cortical depth from Xenium ... versus snRNA-seq meta-analysis SCZ abundance change
   for the 16 Sst supertypes" — i.e. panel 5a is that exact plot plus a layer gutter
   and the stratum partition. Options: (A) keep as-is; (B) drop 5a, define strata in
   text/legend, cite Fig 3f; (C) RECOMMENDED — replace the scatter with a compact
   depth-ordered KEY (16 supertype labels positioned by depth, bracketed by stratum;
   no effect-size axis, no error bars, no rho), which keeps the partition and the
   laminar link but stops restating Fig 3f's result; (D) keep the scatter but demote
   it (smaller, no rho, no error bars) and caption "data as in Fig. 3f".
   ALSO TO RECONCILE regardless of option: Fig 3f legend says MEAN cortical depth,
   whereas 5a uses MEDIAN (`supertype_depth_platform_summary.csv: median_Xenium`).
   Inconsistent duplication is worse than duplication — pick one statistic for both.

1. Sst_20 in depleted (FDR<0.20) — S1 is the evidence either way.
2. Subclass reference column stays in b/d/e? (currently yes).
3. Figure 5 as main figure, last position (agreed here; Etienne/Laramie/Nicole to
   sign off) — RESULTS text already drafted lean.
4. Cilium up-in-both: highlight-only vs one results clause.
