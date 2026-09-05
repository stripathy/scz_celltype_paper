# Plan: responding to Etienne's review (Google Doc "Anonymous" comments, Aug 27–28, 2026)

**Source**: comments + tracked changes in the manuscript Google Doc, left as "Anonymous"
(Etienne, not signed in). Extracted 2026-08-29. Anchor letters ([ar], [as], …) refer to
the doc's comment anchors as of that date.

**Strategy in one paragraph.** Accept every rhetorical and clarity remedy — they are
cheap and they serve the paper's own thesis. Decline the new-analysis ask
(vulnerability-pathway / CMap analysis within the depleted supertypes) for *this* paper,
but make the scope decision visible in the text (a limitations addition pointing to
Tables T2/T2b) and in a comment reply (scope decision + selection-vs-regulation confound
+ Khodosevich follow-up as the designed vehicle). Reserve a GSEA of the already-validated
Sst subclass DE as revision ammunition, dry-run privately first. Do exactly one small new
analysis now: the Fig. 4b enrichment rerun without HCN1, which is a robustness check in
the paper's existing style, not a scope expansion.

---

## A. Discussion restructure (the core of the response)

### A1. Rewrite ¶1 to lead with the findings, not the adjudication ([ar], [as])

Etienne's novelty worry and his "altered state under-addressed" concern share one
remedy: the opening paragraph should carry the full findings package, with the
two-decade debate as frame rather than payload. This also serves the
"depletion-is-transformative" reading: the claim (loss) is old — Benes 1991 — but the
*evidence standard* is new, and ¶1 should say so.

Sketch (numbers as currently in the doc; re-verify at write time):

> Across seven harmonized snRNA-seq datasets (298 controls, 171 SCZ cases) and an
> independent Xenium spatial cohort (12 controls, 12 SCZ), we find that SCZ involves two
> distinguishable forms of Sst-interneuron pathology: a reduction of SST mRNA within Sst
> neurons that spans essentially all Sst subtypes, and a selective depletion of a small
> set of upper-layer Sst supertypes. The depleted population now has a specific
> identity: supertypes marked by CALB1 and high HCN1 expression, with pronounced
> HCN-dependent sag physiology; it includes Sst_25 (reduced ~27%), the likely
> transcriptomic correlate of the double-bouquet cell; it concentrates SCZ
> common-variant risk in proportion to depletion; and it overlaps the earliest neuronal
> populations lost in Alzheimer's disease. These findings bear on a two-decade dispute:
> the dominant view, based on in situ studies, holds that interneurons are
> transcriptionally altered without cell loss (Volk et al. 2000; Hashimoto et al. 2003;
> Dienel et al. 2022, 2023), whereas a competing line reports loss (Benes 1991; Beasley
> et al. 2002; Toker et al. 2018; Batiuk et al. 2022; Kiss et al. 2026). Our data
> indicate both accounts are correct — and explain why the dispute persisted: a ~27%
> depletion of individual supertypes is invisible at subclass resolution and to
> single-marker counting, and becomes detectable only with a granular taxonomy applied
> across hundreds of donors and two measurement platforms.

### A2. Compress ¶2 ("Distinguishing interneuron depletion…") ([at]: "redundant with intro")

He's right — granularity + pooling/power are already in the intro's final paragraph.
The one non-redundant idea (modest supertype depletion is invisible at subclass
resolution) moves into the rewritten ¶1 (above). Reduce ¶2 to a one-sentence bridge into
the detection-confound paragraph ("Because SCZ alters the very transcripts that define
cell identity…"), or delete it outright.

### A3. Limitations: add the scoping statement (answers [as] in the text)

Add as a fourth limitation (or extend the third — both work; the confound is related to
cross-sectionality but distinct):

> Fourth, this study identifies which cells are affected but not the molecular programs
> disrupted within them. We provide the complete subclass- and supertype-level
> differential-expression results as a resource (Tables T2, T2b), but we have
> deliberately not attempted pathway-level interpretation within the depleted
> supertypes: in cross-sectional post-mortem tissue, within-type expression differences
> conflate altered regulation with the selective survival of cell subpopulations.
> Distinguishing these — and identifying what renders these cells vulnerable — requires
> designs beyond this study and is the focus of ongoing work.

This sentence *is* the scope decision, made visible. It also pre-answers the reviewer
who will make Etienne's request.

### A4. Rework the final (therapeutics) paragraph ([ax] Fino nuance + tracked changes)

- Accept the tracked change: "have already advanced to clinical trials" → **"are
  advancing towards clinical trials"** (his domain; factual).
- Integrate the Fino point so the paragraph holds both logics — near-term compensation
  and longer-term protection — instead of the current "will not be enough," which
  undercuts the α5-PAM rationale he pioneered.

Sketch:

> …positive allosteric modulators of α5-subunit-containing GABA-A receptors, which
> potentiate the dendritic inhibition that Sst neurons supply (CITE Etienne's lab work —
> he is providing refs, see [aw]). Such compounds are advancing towards clinical trials,
> and our findings support their rationale: because Sst interneurons provide dense,
> largely unselective inhibition of pyramidal dendrites (Fino & Yuste 2011 — CONFIRM
> exact citation with Etienne), loss of a fraction of these cells should lower overall
> dendritic inhibitory tone rather than silence a discrete circuit — a deficit that
> enhancing the output of surviving Sst cells is well suited to offset. Compensation,
> however, treats the consequence of depletion rather than its cause. The longer-term
> goal is to protect the vulnerable cells themselves, and that population now has an
> identity: upper-layer Sst supertypes marked by CALB1 and high HCN1 expression.
> Understanding what makes these neurons vulnerable, and intervening before they are
> lost, is the natural next step.

### A5. Subplate/neurodevelopmental sentence ([au]: "refs for that?")

Add 1–2 citations for "altered migration or failed apoptosis of subplate neurons."
Candidates likely already in the bibliography's IWMN cluster (Akbarian 1993/1996;
Kubo 2020 review); otherwise [REFS TO FILL — ask co-authors]. Do not write from memory.

### A6. Accept the remaining Anonymous tracked changes in the Discussion

The 2:16 p.m. cluster (the "advancing towards" rewrite, A4 above) and the 2:02 p.m.
"Add: '-'" (hyphen; locate in doc — likely "bulk-tissue-based"). Review-then-accept.

---

## B. Clarity fixes elsewhere in the manuscript (all accepts)

| # | Location | Etienne's point | Change |
|---|----------|-----------------|--------|
| B1 | Abstract | 298/171 then 12/12 unclear which samples → which analysis | Attribute cohorts to analyses explicitly; reuse the sentence pattern Shreejoy already drafted in the Thomas DeLong thread |
| B2 | Abstract/intro | "lost or altered state remains unresolved" — say why it matters | Add clause: "…with direct implications for therapeutic strategy (compensating for dysfunction versus protecting or replacing cells)" — also seeds the final Discussion paragraph |
| B3 | First Results ¶ | Terminology note (Shreejoy already agreed in-thread) | Add 2–3 sentences: Allen nomenclature (Gabitto et al. 2024); class → subclass → supertype as successively finer levels of one hierarchy; "cell type" used generically for any level |
| B4 | Intro | Don't ignore Lewis-lab L3/4 pyramidal literature | One clause acknowledging L3 pyramidal dysfunction (refs TO FILL, e.g. Glausier & Lewis review — verify), then "here we focus on GABAergic cells" |
| B5 | Fig 4 results | "SCZ common-variant risk" — define | One clause at first use: gene-level GWAS enrichment (MAGMA convention), i.e., not PRS |
| B6 | Fig 4 results | "credible set" jargon | Parenthetical gloss: the set of variants most likely to contain the causal variant, from statistical fine-mapping |
| B7 | Results | "reference-based label transfer" unclear | Brief gloss at first use (each nucleus assigned to its best-matching reference type) |
| B8 | Fig 4 results | "pairs" — pairs of what? | Spell out "gene–cell-type pairs" at first use (he self-resolved but asked for clarification) |
| B9 | Fig 2 legend | FDR asterisk definitions unused; redundant with plusses/red italics | Reconcile legend and figure; keep ONE significance encoding. **Owner: Nicole** (bundle with the Fig 2 hand-off) |
| B10 | Global | "is" vs "was" for results | Decide once (recommend: past tense for what was done, present for what the data show), apply globally |
| B11 | Fig 2 header/text | "cell-intrinsic" word choice | Options: (a) keep but define at first use; (b) reword to "reduction of SST mRNA within Sst neurons" (recommended — it's what the section means and needs no definition). SHREEJOY DECISION |
| B12 | Methods/results | Fröhlich = OFC, others DLPFC — tested? | Point to the leave-one-dataset-out robustness (Fig. S6) in the comment reply; optionally add half-sentence in text citing S6 |
| B13 | Results (RNAscope, "Lastly…") | Keep as backup for reviewers; explain LCM-seq bias context | KEEP RNAscope in the paper; add one sentence noting the companion LCM-seq assessment biases against detecting cell loss, whereas the density counts are an independent assessment within the same study (Etienne can help word this) |

---

## C. Analyses: do / defer / decline

**C1. DO NOW — Fig. 4b enrichment without HCN1** (his 1:58 p.m. comment). Rerun the
Sst_3-specific gene-set / gene-driver analysis excluding HCN1. This is a one-gene
sensitivity check on an existing pipeline, in the same spirit as the S10 robustness
supplement. Either outcome is fine to report: signal survives → enrichment is not
single-gene-driven; signal collapses → confirms the text's existing claim that HCN1
carries the association (state it as such). Also clarify the Fig. 4b text pointer to
Methods (his "not clear how this analysis was performed").

**C2. DONE 2026-08-29 — stratified Sst GSEA (the better version of the dry run).**
`transcriptomic/scripts/17_sst_strata_gsea.R` → `transcriptomic/results/sst_strata_gsea/`.
Strata by SCZ crumblr (depleted = FDR<0.20: Sst_2/3/20/22/25; intermediate = negative ns;
non-depleted = β≥0), depth-validated (mean depth 0.31 / 0.50 / 0.72; ρ(β, depth)=0.74).
Headline: two-component structure. (1) SHARED pan-Sst state change — SST itself down in
all three strata (z −6.4 / −8.8 / −3.8) and the synaptic program at equal NES ≈ −1.5
everywhere (sig only where powered). (2) DEPLETION-GRADED signature — significant-pathway
burden 116 → 57 → 0 (padj<0.05); cytoplasmic translation/ribosome sets sign-flip
(NES ≈ −2.4 depleted vs +1.2 ns non-depleted) and OxPhos/mito (incl. PINK1, TOMM40)
grade monotonically. CALB1 itself is down only in the depleted stratum (ties to the
calbindin-IR literature). Nicole's March Unaffected/Vulnerable_AD/Vulnerable_SCZAD group
DEs match these strata (forensic z-correlations) and give the same picture (106 vs 5 sig
pathways). Result is clean and on-message → candidate supplement + results paragraph +
discussion sentence; if adopted, prefer proper stratum-level pseudobulk DE (Nicole's
group files are exactly that) over the IVW ranking, collapse the redundant translation
sets via leading-edge overlap, and keep the survivor-bias caveat sentence. C3's decline
still stands for the causal/CMap version of the ask.

Xenium checkability (2026-08-29, `fig5_xenium_stratum_concordance.png` +
`xenium_stratum_concordance.csv`): the pathway modules are NOT checkable — the 300-gene
panel holds 1/110 translation and 2/125 OxPhos leading-edge genes (quotable in the
limitation). What IS checkable checks out: stratum-level gene concordance on ~130 panel
genes, snRNA-seq vs Xenium (supertype DE aggregated to strata; 10.6k / 10.6k / 6.2k Sst
cells, 24 donors) gives r = 0.33 / 0.35 / 0.13 (all genes) and r = 0.65 / 0.79 / 0.21
among snRNA-seq signal genes (|z|>2; n = 26/30/12 — partly outlier-driven, Spearman is
more modest, so rerun properly before quoting). Single-gene replications: SST down in
ALL three strata on BOTH platforms (Xenium stratum z = −3.3 / −2.2 / −2.5); CALB1's
depleted-specific drop replicates (Xenium z = −3.1 in depleted vs +0.9 in intermediate);
VGF strongly down both platforms in intermediate. To make it legit for the paper:
(a) snRNA-seq side — per-cohort per-donor STRATUM pseudobulks → dream → metafor (same
pipeline one grain up; Nicole's March group files prove the machinery — NICOLE ASK,
bundle with the Fig 2 hand-off); (b) Xenium side — stratum-level pseudobulk edgeR rerun
via run_de.R (aggregate cells by stratum) instead of averaging supertype DE.

Literature hooks (`background/sibille_sst_literature.md`, 2026-08-29): the
translation/ribosome finding has direct Sibille-lab precedent — Lin & Sibille 2015
(EIF2-mediated translation suppressed in laser-captured SST neurons under chronic
stress), Tomoda et al. 2022 (SST neurons selectively enriched for PERK/eIF2α UPR arm),
Zhang et al. 2023 PNAS Nexus (ribosomal genes down in MDD cortex, stress-driven). OxPhos
was PREDICTED for SST cells by Lin & Sibille 2013 but never shown cell-resolved (Newton
2022 put it in pyramidal cells; Arbabi 2025 ATPase-genes-down in SST/SCZ is nearest) —
our graded OxPhos result is the first cell-resolved confirmation. Contrast to cite:
Seney 2015 concluded SST vulnerability is "general… independent of specific cell type" —
the graded/CALB1-linked result revises this. CAVEAT: Kiss 2026 reports SCZ Sst depletion
is age-dependent (present <~70y) — check cohort age-at-death before presenting strata
effects as stable.

Jens independent-cohort replication (2026-08-29, scripts 20a-c →
`transcriptomic/results/sst_strata_gsea/jens/`): Jens dataset (83 donors, 42 CTRL /
41 SCZ, ~30k Sst cells; brisc-mapped in scz_parse_taxonomy; NOT in the 7-cohort meta),
proper pseudobulk limma-voom + fgsea. VERDICT: (1) burden gradient REPLICATES —
172/102/14 down-sets across strata (conf≥0.5: 87/78/2; affected-vs-non: 112 vs 19,
non-depleted direction-balanced); (2) shared synaptic + SST/VGF REPLICATES (VGF t =
−3.3/−2.6/−1.6; SST down all strata; CALB1 depleted-only −2.35; GAD2 does not);
(3) OxPhos graded REPLICATES with blunted gradient (−2.5**/−2.6**/−1.9** GO:BP);
(4) translation does NOT localize to jens-depleted — it appears in jens-INTERMEDIATE
(−1.7**..−2.0**) and in the merged affected stratum (−1.6*/−1.8**), consistent with
the documented depleted-label unreliability in this cohort (4/5 depleted supertypes
CONFIDENTLY_WRONG, F1 0.33–0.54; runner-up audit: dep↔int porosity ~35%, dep↔non ~0%).
Cross-cohort NES r among meta-sig sets: intermediate 0.64, depleted 0.19 (the smearing
signature). Use as supplement-grade support: gradient + shared + OxPhos independent-
cohort confirmed; translation localization deferred to the 7-cohort pseudobulk rerun.

SEA-AD cross-disorder state check (2026-08-29, scripts 21a/21b →
`transcriptomic/results/sst_strata_gsea/seaad_a9/`): Gabitto 2024 (MTG) reports
vulnerable Sst supertypes "did not downregulate ETC and ribosomal genes" (which all
other neurons did), instead losing kinases/HECT E3 ligases/NGF/MME early in CPS. Our
recalculation on SEA-AD DFC (A9; our strata; pseudobulk ~ scale(CPS)+sex+age+PMI;
severely affected donors EXCLUDED per their DE convention, kept for composition):
CONFIRMS their claim and INVERTS our SCZ gradient — along CPS, translation/OxPhos
suppression is strongest in NON-depleted Sst (translation −2.0**, elongation −2.7**,
ribosomal −2.4**) and absent/positive in depleted (+0.5/+0.8/−0.6); synaptic program
goes UP along CPS (+1.8 to +2.0**, vs shared DOWN in SCZ); burden is up-dominated.
Conserved cross-disorder core: SST down all strata (−1.9/−2.4/−1.7) and VGF
(−2.3/−1.4/−2.3); MME down-trend in affected strata (echoes their claim). UPSHOT for
Fig 5/discussion: composition converges across disorders (Fig 4j) but STATE DIVERGES —
do not frame graded translation/OxPhos as a generic vulnerable-cell signature; frame as
SCZ-specific pathophysiology on a shared vulnerability axis. Candidate supplement panel
(SCZ vs AD strata heatmap, `key_blocks_scz_vs_a9.csv`). Caveats: linear CPS blends
their early/late epochs (early-epoch rerun = refinement); survivor-censoring could
also explain the AD flat-in-depleted pattern (suppressed cells die out of the data).

SEA-AD follow-ups (2026-08-29 overnight, scripts 21c/22/23): (1) EARLY-EPOCH sensitivity
(donors CPS ≤ 0.6, n = 31; our operationalization of their early epoch — exact spline in
their Suppl. Note): inversion holds and sharpens — a9early_depleted translation
+1.4/+1.0 (flat/up) vs non-depleted −1.8**/−2.8**; burden 10/68/91 down; SST/VGF down
everywhere; so the AD flat-in-depleted pattern is not a late-epoch artifact.
(2) SELECTION AUDIT of highlighted SCZ blocks: of 112 depleted-sig down sets, families
by ≥50% LE-overlap = translation 29 / OxPhos-mito 23 / synaptic 12 / other 48; top-20 by
padj = 18/20 in the three families (top-13 all translation-family); full annotated table
→ `suppl_table_strata_gsea_annotated.csv` (supp-table-ready). Selection rule to state:
"representative sets from the largest leading-edge families among depleted-significant
sets; full table in Supp Table X." (3) MARQUEE GENES (gabitto_marquee_genes_in_scz.csv):
MME tested and DOWN in SCZ Sst (subclass z = −3.15; intermediate −2.2) — cross-disorder
candidate; NGF flat in SCZ; NGFR/ATP5MPL not tested (expression-filtered); MAPK8 UP in
SCZ depleted (+2.6; opposite of their AD claim); RPL4 fits the sign-flip; none on Xenium
panel. (4) DIVERGENCE FIGURE `figS_scz_vs_ad_depleted_divergence.png` (pathway r = 0.21,
gene r = 0.39; SST/VGF isolated in the down-down quadrant) + Gabitto-style trajectories
`figS_a9_gene_trajectories.png` (SST/VGF decline in all strata; module inversion is a
module-level phenomenon, subtle per-gene). Figure-style borrowings from their Fig 5h/i:
family-colored two-axis effect scatters (adopted), LOESS trajectories (adopted for AD
supplement), epoch-split reporting (adopted as sensitivity).

Fig 5 mock-up v2 (2026-08-29, script 18 rebuilt): six panels — a strata definition,
b burden, c GENE-level depleted-vs-non-depleted scatter with module-family coloring
(Gabitto Fig 5h idiom; SST/VGF isolated off-diagonal), d pathway blocks, e exemplar
genes, f SCZ-vs-AD depleted-stratum NES divergence scatter (the "composition converges,
state diverges" close; r = 0.21). Panel f commits the main figure to the A9 CPS
reanalysis (methods + SEA-AD provenance; precedent exists via Fig 4j) — fallback if
co-authors prefer pure-SCZ: f moves to the companion supplement and the divergence
becomes one discussion sentence. NES pair scatters (old panel c) now supplement-only
(script 17 fig3). Not adopted from SEA-AD: CPS trajectories in the main figure (no
progression axis in SCZ; AD trajectories live in the supplement), epoch schematic.

Fig 5 iteration notes (2026-08-29 cont.): panel a now carries cortical-layer annotations
from `spatial/output/depth_proportions/proposed_layer_boundaries.csv` (xenium column) —
depleted stratum visibly L2/3. UP-IN-BOTH (panel c candidates): genes RPS6KA5, INSYN2B
(inhibitory postsynapse!), TMEM161B-DT, GHR, MEF2C-AS1, GRIN2A (~+2 all strata);
pathways = the cilium/Smoothened cluster (dep-sig, non-dep positive ns). Prospective
last-row gene-concordance panel (script 24, seaad_a9/figS_gene_concordance_DRAFTS.png +
gene_concordance_classes.csv): classes at stated thresholds — shared down n=124 (SST,
VGF, CBLN4, UNC13A, CIRBP, RASGRF2), SCZ-specific down n=61 (NTM, MEG3, HDAC4, RPL36,
RPS21), shared up n=272 (INSYN2B, GLCCI1, PTPRD, CELF2), SCZ-up/AD-down n=13 (TCF7L2,
ACAT1, SYTL2). Draft 1 = class heatmap, draft 2 = annotated scatter; ITERATING before
formal inclusion. Known wrinkle: NDUFS8 curated into SCZ-specific but falls to Other by
threshold (t_ad = −0.9) — resolve thresholds-vs-curation. Jens per-donor exemplar panel
(figS_gene_exemplars_jens_DRAFT.png; SST/VGF/RPL36/INSYN2B) stands in for raw-data
exemplars until Nicole's 7-cohort pseudobulks land.

Fig 5 v6 layout (2026-08-29, final mock-up state): row 1 a–c (strata definition with
layer gutter / burden / gene scatter); row 2 d–e ONLY (SCZ shared-vs-graded focus);
row 3 f–h = SCZ↔AD depleted-stratum comparison (f gene scatter with module centroids +
concordance annotation; g VGF boxplot pair; h RPL36 boxplot pair, SCZ side = Jens
placeholder, footered). FRAMING REVISION (Shreejoy's read, now quantified): DE is
largely SHARED between SCZ and AD in depleted cells — among genes with |effect| > 2 in
both, 96% sign-concordant (n = 327, r = 0.89; 99% at >2.5) — with the translation/
OxPhos modules as the conspicuous exceptions and SST/VGF the strongest shared elements.
Discussion sentence should read "substantially shared state, module-level exceptions,"
NOT the earlier blanket "state diverges."

Fig 5 v7 FINAL LAYOUT DECISION (2026-08-29): two rows, PURE SCZ — a (definition +
layer gutter), b (burden), c (gene scatter, now with the motile-cilium up-in-both
module in green, WDR19/CEP126 labeled), d–e (shared vs graded). SCZ↔AD row DROPPED
from the main figure (Shreejoy: DE similarity follows logically from Fig 4's shared
depletion; not worth making explicitly unless reviewers ask). The SCZ↔AD material
(scripts 22/24/25/26; divergence scatters, concordance classes, exemplar boxplots,
96%-concordance stat) is HELD IN RESERVE as supplement candidates / reviewer answers.
CELL-COUNT ANSWER (for panel-b caveats): non-depleted stratum has ~40–50% fewer cells
everywhere (Xenium controls 6.2k vs 10.6k; Jens 5.9k vs 11.2k; A9 16.3k vs 27.7k;
meta proxy: median IVW SE 0.0595 vs 0.046, fewer genes tested). This affects
significance counts in non-depleted but cannot explain the gradient: intermediate is
the BEST-powered stratum (SE 0.0425) yet has half the depleted burden, and the
translation sign-flip is directional, not attenuation. ROBUSTNESS TO-DO: cell-count-
matched downsample rerun (feasible on Jens/A9 now; meta once pseudobulks land).

Narrative-architecture note (2026-08-29): Shreejoy's worry — strata-DE as Fig 5 loops
back to DE and invites "why not abundance first?" Resolution: Fig 2 and Fig 5 ask
different questions (replicate the undisputed state change / test whether state and
loss are COUPLED); Fig 2 is the pipeline-credibility anchor and the abundance claim's
re-annotation defense presupposes it; Fig 5 must be FRAMED as the synthesis of the
paper's two axes, not as more DE. Litmus test for inclusion: if the Fig4→Fig5
transition sentence ("are the two debated phenomena independent or coupled?") reads
naturally after the pseudobulks land, keep as Fig 5; if it reads as assorted
stratified DE, demote to supplement. Do NOT reorder Figs 2/3.

Tomoda-model check (2026-08-29, `figS_sst_baseline_vs_depletion.png` + CSV, scripts
19a/19b): baseline SST expression per supertype vs crumblr depletion, three sources.
Directionally consistent with the preproSST-load model in all three (SEA-AD reference
ρ = −0.47, p = 0.068; SEA-AD DLPFC CP10K ρ = −0.26; Xenium controls ρ = −0.41,
p = 0.11) but never significant at n = 16, and depth-confounded in the SEA-AD sources
(partial | depth ≈ −0.11); only the Xenium-control baseline is depth-independent and
holds up (partial −0.48). Sst_10 is the standing counterexample (high SST, not
depleted). VERDICT: suggestive trend, back-pocket / Etienne-conversation material —
not a paper claim.

Module-gene check (fig6, `module_genes_xenium_concordance.csv`): the graded
translation/OxPhos signature gets NO meaningful Xenium test — its 3 panel genes are
mixed (MRPL27, NDUFB2 discordant; UQCRFS1 down in Xenium but in ALL strata, so not
stratum-specific) → snRNA-seq-only claim, either way. The SHARED synaptic program IS
partially checkable and largely concordant: 11/15 sn-significant module-gene×stratum
pairs sign-concordant, 8 also Xenium-nominal — anchored by VGF (down in all three
strata on BOTH platforms; second pan-Sst shared gene after SST, and an
Etienne-literature hook via the SST/VGF axis), RASGRF2, NTRK2, GAD2, CALB1
(depleted-specific). Discordant genes (GAD1, VAMP1, NEFL: down nuclear, flat/up
whole-cell) mirror the paper's PVALB nuclear-vs-cytoplasmic precedent.

**C3. DECLINE FOR THIS PAPER — vulnerability pathways within depleted supertypes +
CMap/L1000.** Rationale (for the reply and our own records): (a) cannot meet the paper's
evidentiary architecture — the 300-gene Xenium panel cannot replicate unbiased pathway
results even in principle, and within-supertype DE in rare populations is underpowered;
(b) conceptually premature in cross-sectional tissue — selection (survivor bias) vs
regulation cannot be separated when a subpopulation of the supertype is dying;
(c) CMap/L1000 is candidate generation, a different epistemic tier than the rest of the
paper. This is the centerpiece question for the Khodosevich collaboration, which can
generate data designed for it.

---

## D. Comment replies + asks back to Etienne

**Reply to [as]/[ar] (the scope question) — sketch:**

> Totally agree this is the key next question. We made a deliberate scope decision here:
> rigorous DE for all subclasses and supertypes with complete tables as a resource
> (T2/T2b), but no pathway-level interpretation in this paper — partly to keep the focus
> on the abundance result, and partly because within-supertype expression differences in
> cross-sectional tissue conflate altered regulation with selective survival of
> subpopulations, so we can't do this at the standard of the rest of the paper (discover
> in snRNA-seq, replicate in Xenium). In response we've (1) restructured the opening of
> the Discussion to foreground what we do establish about the depleted cells, and
> (2) added a limitation making the scope decision explicit. The vulnerability-pathway
> question (including therapeutic-signature analyses) is the centerpiece of follow-up
> work where we can generate data designed for it.

**On his novelty point ([ar])**: concede the remedy (¶1 rewrite, done), push back gently
on the premise — the citable record is a two-camp literature with the no-loss in-situ
view dominant; the claim is old but the evidence standard (supertype resolution,
marker-independent re-annotation, no compensatory relabeling, cross-platform, genetic
convergence) is new.

**Asks back to Etienne:**
1. References for the α5-PAM sentence ([aw] "CITE Etienne's lab work") — already
   requested Aug 19.
2. Exact Fino citation for the blanket-inhibition point (Fino & Yuste 2011 Neuron? Fino
   et al. 2013 review?).
3. Subplate-interpretation refs if he has favorites ([au]).
4. Funding statement ([ay]) — already requested.
5. His answer to Shreejoy's open [av] question (is it wise to spell out the
   double-bouquet counting experiment so explicitly?).

---

## E. Suggested order of operations

1. Review-and-accept the Anonymous tracked changes (typos + "advancing towards").
2. One writing pass for the Discussion: ¶1 rewrite (A1), ¶2 compression (A2),
   limitations addition (A3), final-paragraph rework (A4). These interlock — do them
   together so the "identity of the depleted cells" thread runs ¶1 → limitations →
   close.
3. B-item sweep (B1–B8, B10–B13) — one sitting; most are single clauses.
4. Hand B9 (Fig 2 legend encodings) to Nicole with the Fig 2 refresh bundle.
5. Run C1 (Fig 4b minus HCN1); drop into the robustness supplement if clean.
6. Post the comment replies (D), including the scope reply on [as]/[ar].
7. Scratch-run C2 (GSEA dry run); file results privately, not in the repo manuscript.

**Out of scope for this plan**: Marlen's and Thomas's comments (separate pass);
Shreejoy's own open to-dos in the doc (funding, data/code availability, provenance
notes).
