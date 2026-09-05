# Figure 5 integration audit — Google Doc read 2026-09-01

Full read of the live Doc (`1cO5ZStbp9b3cb6vb9s6H2YfbprYGikXBvFHQR7gROjo`,
288 paragraphs, 47 open / 81 resolved comment threads). Line numbers (L#) are
paragraph indices in today's pull and are only for locating text; quote the
anchor phrase when editing.

Every number below was checked against the pipeline today
(`transcriptomic/scripts/fig5/`): K63 FDR = 0.0066; ribosomal-subunit interaction
FDR = 3.40 × 10⁻⁶; electron-transport-chain FDR = 0.053; Fig 5a depth ρ = 0.74
(median Xenium depth); rebuilt Fig 4a ρ = 0.547, P = 0.031 (SEA-AD-only).

Principle: the Fig 5 section, legend and S10 are **already in the Doc** and are
numerically correct except for one rounding. What is missing is everything that
*surrounds* the figure — the abstract, the roadmap, the Discussion, the Methods,
the tables — plus a few collisions between Fig 5's vocabulary and Fig 4's.


> **Superseded in part (2026-09-01, later the same day).** The analysis is no longer a
> main figure: it is one paragraph closing the Figure 3 section plus supplement S8, and
> the controls supplement is dropped. Items A1, A2, A3, A5, A7, A8, A9 and A11 below
> refer to the main-figure configuration and are replaced by
> `manuscript/figure3_tail_paragraph_2026-09-01.md`. Sections B (with S8/S9/S10 now
> meaning strata / landscape / robustness), C, D and E still apply.

---

## A. Required edits to accommodate Figure 5

### A1. Abstract (L16) — add one sentence
Abstract is 145 words; Fig 5 is absent. Proposed sentence, placed **last**, after
"...intrinsically vulnerable population of upper-layer interneurons.":

> Consistent with this vulnerability, the depleted supertypes also showed the
> deepest transcriptional disruption, with suppressed protein-synthesis and
> oxidative-phosphorylation programs superimposed on the *SST* and synaptic
> deficits shared across Sst neurons.

(~32 words → ~177 total.) Alternative placement: directly after the "Second, we
find depletion..." sentence, which follows the lost-vs-altered logic of the
opening sentence more literally but breaks the First/Second/Notably cadence.

### A2. Introduction roadmap (L22) — extend the final sentence
Current: "Finally, we directly test whether the supertypes depleted in SCZ are
the same ones lost as Alzheimer's disease advances."

Proposed:
> We also test whether the supertypes depleted in SCZ are the same ones lost as
> Alzheimer's disease advances, and finally ask whether the transcriptional
> disease state of Sst neurons differs between the supertypes that are depleted
> and those that persist.

This is also the paper's most direct answer to the roadmap's *first* question
("lost, altered, or both"), which Etienne's comment #7 says is otherwise
unaddressed.

### A3. Results §4 opener (L64) — "Finally" → "Next"
"Finally, we asked whether the depletion of upper-layer Sst supertypes is
etiologically relevant..." → "**We next asked** whether...". Fig 5's section
(L72) keeps its "Finally".

### A4. Results §5 — one number (L74)
"K63-linked deubiquitination, NES = +2.15, **FDR = 0.006**" → **FDR = 0.007**
(pipeline: 0.0066; 0.006 is a truncation, not a rounding).

### A5. Results §5 — bridge the two "not/non-depleted" vocabularies (L72)
Fig 4 (L66, L70) and Methods (L150) use **not-depleted** for the *eleven*
supertypes outside the depleted five. Fig 5 (L72–78) uses **non-depleted** for a
*five*-supertype group. A reader will assume these are the same set. Add one
clause at L72 after "...non-depleted (n = 5)":

> ; the not-depleted supertypes of Fig. 4 are thus split here into intermediate
> and non-depleted groups

Also harmonize "not depleted" (L70, no hyphen) → "not-depleted". Checked: the
Fig 4e "not depleted" exemplars Sst_5 and Sst_1 both fall in Fig 5's
non-depleted group (non-depleted = Sst_1, Sst_4, Sst_5, Sst_7, Sst_10;
intermediate = Sst_9, Sst_11, Sst_12, Sst_13, Sst_19, Sst_23), so nothing
conflicts — it is only the naming.

### A6. Results §5 — Lin and Sibille 2015 (L72)
Bare citation, not Zotero-linked, absent from the reference list. Add to Zotero
(comment #4 attaches PMID 25600109) and confirm it is the paper reporting
suppressed translation under chronic stress rather than the companion review.

### A7. Results §2 — optional forward pointer (L50)
After "...left supertype-level DE underpowered for most supertypes (Figure S4,
full results in Supplementary Table T2b)" add ", motivating the grouped analysis
of Fig. 5". Optional; helps the reader see why Fig 5 pools supertypes.

### A8. Discussion — currently silent on Fig 5
The Discussion (L80–88) never mentions the graded translation/oxphos finding.
Three edits, the first two required:

**(a) ¶1 (L80), last sentence — add a clause.** "Namely, we find a broad
reduction of SST mRNA within Sst neurons that co-occurs with the depletion of
specific upper-layer Sst supertypes**, and the transcriptional disruption is
itself deepest within the supertypes being depleted**." This is the direct
answer to Etienne's comments #7 and #11.

**(b) New paragraph between the intrinsic-vulnerability ¶ (L85) and the L6b ¶
(L86).** A ~250-word version with literature exists in
`manuscript/figure5_section_DRAFT.md` (line 421). Trimmed ~170-word version
below, citation placeholders in brackets, no figure refs (Discussion
convention):

> Resolving the Sst disease state by depletion status separates two superimposed
> programs. The reduction of *SST*, *VGF* and synaptic transcripts is present in
> depleted and persisting supertypes alike and is most parsimoniously read as a
> circuit-level state imposed on all Sst neurons. The suppression of protein
> synthesis and oxidative phosphorylation, by contrast, marks the depleted
> supertypes specifically and connects to prior evidence that cortical Sst
> neurons reduce translation under chronic stress [Lin and Sibille 2015]. These
> are state differences, not causes: suppressed protein synthesis and energy
> metabolism in the depleted supertypes could precede and promote their loss,
> could reflect the response of surviving cells to an ongoing disease process, or
> could partly reflect which cells survive. Whichever holds, protein synthesis
> and oxidative ATP production are among the costliest activities of a neuron,
> and their coordinated withdrawal in the very cells declining in number is the
> signature expected of a population under sustained strain.

**(c) Limitations (L87), third point — add a clause.** "...cannot temporally
order transcriptional from compositional changes, **nor determine whether the
suppression of translation and oxidative phosphorylation in the depleted
supertypes precedes their loss or reflects the state of the cells that
remain**, nor separate the effects of disease from..."

**Consider (not required):**
- Loss-vs-erosion ¶ (L82). Fig 5 is compatible with both readings — under
  erosion, cells drifting *out* of a depleted label would leave the remaining
  labelled cells looking *less* altered, which is the opposite of what we see;
  under partial drift the graded DE is what erosion would predict. Net, it
  strengthens the paragraph's own conclusion ("either account entails
  substantial disruption of the same neurons"). One clause could say so: "...a
  disruption the depletion-graded analysis now shows directly".
- Closing ¶ (L88). Etienne's comment #8 argues that activating the surviving Sst
  cells should be therapeutic. Fig 5 supports this: the persisting supertypes
  carry only the shared *SST*/synaptic deficit, not the translation/oxphos
  collapse. Optional clause after "restoring their inhibitory output will not be
  enough": "although the persisting supertypes, which retain their
  protein-synthesis and energy-metabolism programs, remain plausible targets for
  it".

### A9. Methods — two insertions (text in `manuscript/figure5_methods_2026-09-01.md`)
- **Insert (A)** into "Cell type-specific differential gene expression", after
  the Meta-analysis paragraph (L130): depletion-group aggregation as a third
  aggregation level (~120 words).
- **New section (B)** "Gene-set enrichment and depletion-group interaction
  analyses", after "Patch-seq electrophysiology and morphology" (L153), before
  References (~560 words). The interaction model is now **equation (4)** — the
  Doc has (1) DE, (2) composition, (3) AD pseudo-progression.
- The draft's "Depletion groups" paragraph now **refers to the existing
  definition** in "Molecular characterisation of the depleted Sst supertypes"
  (L150) and splits its eleven not-depleted supertypes, rather than re-defining
  the groups from scratch (updated today).
- **Depth statistic mismatch.** Fig 3f (L62) plots *mean* Xenium depth against
  abundance change (ρ = 0.76, P < 0.001). Fig 5a (L78) plots *median* depth;
  the pipeline gives ρ = 0.74, P = 0.0016 for that. Same relationship, two
  summary statistics, two ρs. I removed the second ρ from the Methods draft and
  cite Fig. 3f instead. Decide whether Fig 5a should adopt Fig 3f's statistic
  (I could not find the Fig 3f input table in this repo — it is Nicole's panel).
- **Software versions.** The Doc has no software-version paragraph anywhere. The
  draft's Software line should either become a general Methods paragraph or be
  trimmed to the fgsea/msigdbr versions inline. MSigDB is now stated as release
  2025.1.Hs via msigdbr 25.1.1.
- **References to add:** Korotkevich et al. 2021 (fgsea), Liberzon et al. 2015
  (MSigDB), Lin and Sibille 2015. None is in the current list.

### A10. Supplementary tables (L267–276) and Data availability (L93)
The list holds T1, T2, T2b, T3. **T4 is cited in Methods (L147) but absent from
the list** — pre-existing gap, and it is also the table the SEA-AD-only rerun
replaces (see C). Fig 5 needs new tables; suggested:

| Table | Content | Source file | Rows |
|---|---|---|---|
| T5 | Depletion-group DE meta-analysis (gene × group: log₂FC, SE, P, FDR) | `stratum_meta_de.csv` | 45,107 |
| T6 | GSEA per depletion group and pooled Sst (NES, P, FDR, leading edge) | `gsea_all_signatures.csv` | 25,699 |
| T7 | Diagnosis × depletion-group interaction GSEA | `interaction_gsea.csv` | 12,526 |

(or one multi-sheet T5). Then fill "Supplementary Tables T1–T[n]" at L93.

### A11. S10 legend, panel c (L266) — stale phrase
"together with the two sets named in the text (diamonds)". The text now names
**ribosomal subunit** (already a Fig 5d row, so a regular point) and **electron
transport chain**; the figure's diamonds are **SRP-dependent cotranslational
protein targeting** (no longer named anywhere) and ETC. Either reword — "together
with two additional sets, SRP-dependent cotranslational protein targeting and
electron transport chain (diamonds)" — or rebuild S10c without the SRP row
(`transcriptomic/scripts/fig5/supp/figS_strata_controls.R`, `CITED`). Rewording
is the minimal fix. Also: title prefix "Fig. S10 |" vs "Figure S# |" elsewhere.

### A12. Table 1 legend (L45) — RNAscope residue
"The lower block lists the two independent datasets used for orthogonal
replication comprising the DLPFC Xenium spatial transcriptomics." → "The lower
block lists the independent Xenium spatial transcriptomics dataset used for
replication." Remove "sgACC, subgenual anterior cingulate cortex" from the
abbreviations.

---

## B. Supplementary-figure cross-references (off by one after the renumbering)

The legends were renumbered (S8 = enrichment landscape, S9 = genetics/AD
robustness, S10 = Fig 5 controls) but three in-text citations were not:

| Location | Reads | Should read |
|---|---|---|
| L64 "(Methods; Figure S9)" | S9 | **Figure S8** |
| L65 "robustness analysis ... in Figure S10" | S10 | **Figure S9** |
| L67 "MTG comparison in Figure S10" | S10 | **Figure S9** |
| L73, L74 "Fig. S10" | S10 | correct |

The archived RNAscope block (L277–287) still cites "Fig. S8", "Fig. S8c" and
"Figure SXXX"; recommend moving it to a separate holding Doc so it cannot be
mistaken for a live cross-reference. Likewise strip "Reviewer suggestions"
(L99–112) before export.

---

## C. Outstanding beyond Fig 5 — the SEA-AD-only taxonomy (Laramie, comment #5)

This is the largest unreconciled item. Laramie asked (2026-08-29) to drop the
combined SEA-AD + Siletti taxonomy. The repo landed the rerun on 2026-08-31
(uncommitted; `git status` shows the rebuilt Fig 4, new `*_seaad125.*` outputs
and retired `*_501.*` files). **The Doc has not been touched.** Shreejoy's own
"update me" comments (#1, #2) sit on the S8 and S9 images.

Verified today from the rebuilt panel data: **Fig 4a ρ = 0.55, P = 0.031**
(n = 16; Doc L65 still says ρ = 0.66, P = 0.0069). Top risk-carrying Sst are
Sst_2, Sst_20, Sst_23, Sst_3 — so "including Sst_2, Sst_3 and Sst_20" survives.

Cascade if adopted (numbers from the 08-31 memory note; re-verify the S8 counts
against `genetics/results/tables/scz_enrichment_seaad125_bigdeli_dlpfc.csv`
before pasting):

- Results L65: ρ/P → 0.55 / 0.031.
- Methods L145: drop "while still being ranked against a brain-wide background".
- Methods L146: delete the **Combined taxonomy** paragraph entirely.
- Methods L147: 16,668 × 501 → 16,544 × 125; 15,981 → 15,855 genes; "the 501
  types do not compete" → 125.
- Data sources L120 (Siletti atlas ¶) and Data availability L93 (Siletti link):
  delete.
- S8 legend (L259–260): single-panel, 125 supertypes, thresholds over 125 tests;
  image → `genetics/results/figures/supp_scz_enrichment_seaad125.png`. The
  submission folder currently holds `S08_scz_enrichment_501.*` — swap.
- S9 legend (L263a): the rebuilt figure shows all four reference × GWAS
  combinations; the MTG rows are non-significant without Siletti. Decide keep
  (honest) or drop MTG.
- Supplementary Table T4 → the seaad125 table.
- Check that the Fig 4 image in the Doc is the rebuilt one (panel a inset
  ρ = 0.55, P = 0.031; panel b 350 drivers, HCN1 at the 93.4th percentile).
- Commit the repo state once the decision is final.

Fig 5 is unaffected: its groups come from the compositional meta-analysis, not
from the genetics.

---

## D. Open comment threads that Fig 5 answers or touches

| # | Author, anchor | Fig 5 relevance / action |
|---|---|---|
| 7 | Etienne, Discussion ¶1: "you don't really address ... altered molecular state" | Answered by Fig 5. Reply after A8(a)/(b) are in. |
| 11 | Etienne: "both accounts are correct ... does this detract from novelty?" | Fig 5 adds the novelty: the state is graded by depletion. A8(a). |
| 8 | Etienne: "activating the surviving ones should have therapeutic efficacy" | Supported — persisting supertypes carry only the shared deficit. A8 closing-¶ option. |
| 12 | Etienne: redo Fig 4b without HCN1 | Still outstanding (Etienne response plan). |
| 3 | Shreejoy → Keon/Marlen: review §5 title | Pending their read. |
| 4 | Shreejoy: Lin & Sibille PMID | A6 — add to Zotero. |
| 5 | Laramie: remove combined taxonomy | Section C. |
| 22 | Anonymous, S5 legend: unused markers | Pre-existing; note S10 deliberately uses the extended scheme. |
| 1, 2 | Shreejoy: "update me" on S8/S9 images | Section C. |

The Etienne response plan (`manuscript/ETIENNE_REVIEW_PLAN.md`) recorded a
decision to *decline* the pathway-analysis ask. Fig 5 now supplies exactly that
analysis; the response should say so rather than decline.

---

## E. Other pre-existing items noticed on this read

- L48: "suggesting a subclass-wide rather than supertype-specific." — missing
  noun ("effect").
- L45 vs L115: "at or under 70 years" vs "younger than 70 years".
- S4 title (L247) is bare and lacks a terminal period.
- Acknowledgements (L90): LLM statement lists "Fable 5"; 5.1 was also used.
- Comments #35/#36 (data/code availability), #41 (GWAS provenance), #47
  (Patch-seq source), #21 (terminology note Shreejoy agreed to add to the first
  Results paragraph) — all still open.
- `spatial/code/analysis/plot_supertype_depth_casecontrol.R` still writes into
  the supplementary folder and would re-create the retired S10 orphan on re-run.

---

## Suggested order of operations

1. Paste the Methods insert and section (A9) — self-contained, no dependencies.
2. Apply the mechanical fixes together: A3, A4, A5, A12, and the three
   cross-references in B.
3. Abstract and roadmap (A1, A2).
4. Discussion (A8a, A8b, A8c), then reply to Etienne's #7/#11.
5. Tables (A10) and S10c legend (A11).
6. Decide Section C and, if adopted, do it as one pass with Laramie's thread
   open.
