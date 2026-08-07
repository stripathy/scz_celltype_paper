# Plan — coupling cell-type DE with compositional change at the person level

**Status:** design only (not yet run). The per-donor snRNA-seq matrices live on a
collaborator's machine and will be pulled later. What is local today: the snRNA
**DE meta** and the **crumblr compositional meta** (group-level summaries) — enough
for the Phase-1 cell-type-level check, *not* the person-level test.

---

## 0. The question

In schizophrenia (SCZ), are **transcriptional dysregulation** (cell-type DE) and
**cell-composition change** two faces of *one* disease axis, or *independent* axes
that vary from patient to patient?

Operationally (Shreejoy's framing): from the group analysis, define the *average
SCZ effect* along two modalities — a DE signature and a compositional signature.
Then **score each individual by how closely they resemble that average SCZ effect**
along (i) expression and (ii) composition, and ask whether the two per-person
scores **track each other across people, within cases** — or whether individuals
can have one (e.g. strong DE) without the other (little composition change), and
vice versa.

**Primary venue: the 7-cohort snRNA-seq meta dataset (~200 cases).** The
person-level question needs many people; Xenium's 24 donors are underpowered
per-person. Xenium instead serves as a **clean in-situ compositional validator**
(no dissociation bias).

---

## 1. Hypotheses / outcomes

The headline is the within-SCZ correlation of the two per-person scores:

| Within-SCZ corr(S_DE, S_comp) | Interpretation |
|---|---|
| **Positive** | A *unified* disease-severity axis — transcriptional and compositional dysregulation are coupled facets of one process. |
| **≈ Zero** | *Orthogonal* axes → patient **heterogeneity / subtypes** (some "transcriptional," some "compositional"). |
| **Negative** | A trade-off (unlikely; would be striking). |

Pre-state that the **independent (null) result is itself informative** — it would
argue SCZ molecular pathology is multi-dimensional, and yields a 2-D patient
phenotype for downstream subtyping. Ruzicka et al. 2024 (next section) already
showed the *transcriptional* axis alone splits SCZ into two subgroups; our novel
contribution is the **second, compositional axis** and whether it couples to the first.

---

## 1b. Precedent — Ruzicka et al. 2024 (Science): the DE-only version (TPS)

**"Single-cell multi-cohort dissection of the schizophrenia transcriptome"**
(McLean + MSSM, 140 individuals, 25 cell types) is the closest published precedent —
and **McLean + MSSM are two of our 7 meta cohorts**, so we can reuse / compare
against their result. They build the DE half of this plan as the **Transcriptional
Pathology Score (TPS)**.

**What they did**
- **TPS ≡ our `S_DE`.** Per individual × cell type, the consistency between the
  person's expression (relative to the cohort average) and the SCZ-vs-control DE
  direction, **aggregated by averaging across NEURONAL cell types** into one
  per-person score. (Validates our collapse-by-mean; refinement = neuronal types only.)
- **Out-of-sample via cohort split:** the SCZ direction was recomputed within McLean
  only, MSSM only, or the joint meta, and subgroup assignments were consistent across
  all three — a cross-cohort (independent-batch) alternative to within-cohort LOO.
- **Ranking TPS → four subgroups:** SZ (35), SZ_CON-like (13), CON_SZ-like (17),
  CON (36) — i.e. ~20% of clinically-SCZ individuals are *transcriptionally
  control-like*. Unstable assignments had **small-magnitude TPS** (low confidence,
  not irreproducible).
- **Un-collapsed companion = matrix decomposition.** An unsupervised decomposition
  over all cells/individuals recovered four **cellular states shared across multiple
  neuronal populations** (Ex_SZCS, In_SZCS positively correlated with TPS, R = 0.35 /
  0.56; Ex_SZTR R = −0.39). Ex_SZCS and In_SZCS were highly correlated across genes
  (R = 0.6) and across individuals (R = 0.7) — the excitatory and inhibitory disease
  states co-occur in the same people. Top shared genes: **DHFR** (one-carbon
  metabolism), **GRIN1** (NMDA), CNTNAP2, CHD5.
- **Genetics:** TPS correlated with PRS (r = 0.26, p = 9e-3), **but the subgroups were
  NOT separated by PRS** → the heterogeneity axis is largely not common-variant
  (they invoke one-carbon / folate, i.e. environmental).
- **No composition signal:** they explicitly report **no cell-type composition change**
  in their snRNA data — which is *why* they did the DE version only.

**Lessons folded into this plan**
1. `S_DE` is published-validated as TPS → keep the collapsed score as the headline,
   but **aggregate over neuronal cell types only** (§3); glia dilute the signal.
2. **Prefer cohort-split out-of-sample** (train DE on k−1 cohorts, score the k-th)
   over within-cohort LOO — our 7-cohort design supports it natively and it doubles
   as a reproducibility test (§4).
3. **Add a latent-state companion** (matrix decomposition / NMF over donor × cell
   expression) as the principled un-collapsed view — it recovers cross-cell-type
   shared states; the right answer to "don't just average across types" (§5).
4. **Composition is the orthogonal axis they could not access** — test whether our
   compositional axis aligns with, or is orthogonal to, the TPS transcriptional
   subgroups; that comparison *is* the novelty (§1, §5).
5. **Compare against / condition on PRS** (we have the GWAS set); prior = the
   heterogeneity axis is partly environmental, not common-variant (§5, §6).
6. **Down-weight near-origin donors** — small |TPS| = low confidence; likewise for
   donors near (0,0) in our (S_DE, S_comp) plane (§6).
7. **Reuse their TPS directly:** because McLean + MSSM are in our meta, correlate our
   `S_DE` against their published TPS as an external validity check.

---

## 2. Data

Per donor `d` (cohort `k`, diagnosis `dx_d`), on the **same cells** for both
modalities (cell-inclusion already harmonised across DE and crumblr):

| Item | What | Where | Local? |
|---|---|---|---|
| Expression | per-(donor × cell type) pseudobulk, log-normalised (logCPM) | snRNA cohort count matrices | ✗ (collaborator) |
| Composition | per-(donor × cell type) CLR proportions (crumblr input) | snRNA: collaborator; **Xenium: `spatial/output/crumblr/crumblr_input_subclass_corr.csv`** | partial |
| Covariates | cohort, diagnosis, age, sex, PMI, per-(donor×type) `n_cells` | cohort metadata | ✗ |
| **DE direction** β | snRNA meta + per-cohort DE | `data/DE_genes_all_cells_scz.csv`, `data/meta_results_cohorts_subclass.csv` | ✓ |
| **Comp direction** γ | crumblr 7-cohort meta (per supertype) | `…/nicole_scz_snrnaseq_betas/final_results_crumblr_7_cohorts.csv` (merged into `scz_cell_type_enrichment/results/tables/gwas_vs_casecontrol_composition.csv`) | ✓ |

**Granularity:** **subclass** as the primary unit (aggregate supertype composition
→ subclass); supertype as a sensitivity / confound check (§6).

---

## 3. Score definitions

All scores are computed **per cohort**, **centred on that cohort's controls**, and
**leave-one-donor-out (LOO)** so a donor is never scored against a signature their
own data helped define.

**Expression**
- Deviation: `δ[d,c,g] = logCPM[d,c,g] − mean_{control donors in cohort k}(logCPM[·,c,g])`
- LOO DE direction: `β⁻ᵈ[c,g]` = SCZ-vs-control logFC for cell type `c`, refit excluding donor `d`
- Per-type alignment: `s_DE[d,c] = proj(δ[d,c,·] onto β⁻ᵈ[c,·])` (restricted to the cell type's DE genes; |β|- or reliability-weighted). Projection = `⟨δ, β⟩ / ‖β‖`.
- **Per-donor DE score:** `S_DE[d] = reliability-weighted mean over c of z(s_DE[d,c])`
  — this is the Ruzicka **Transcriptional Pathology Score (TPS)**; **aggregate over
  neuronal cell types as primary** (glia dilute the signal), all-cell-types as a sensitivity

**Composition**
- Deviation: `η[d,c] = CLR[d,c] − mean_{control}(CLR[·,c])`
- LOO comp direction: `γ⁻ᵈ[c]` = crumblr SCZ CLR effect, refit excluding donor `d`
- **Per-donor comp score:** `S_comp[d] = proj(η[d,·] onto γ⁻ᵈ) = ⟨η, γ⁻ᵈ⟩ / ‖γ⁻ᵈ‖`

Both scores z-standardised within cohort before the scatter. (Cosine similarity is
an alternative to projection — see §5 robustness.)

> **Aggregate the per-cell-type *scores*, never average the gene *signatures* across
> cell types** — averaging signatures cancels genes that go up in one type and down
> in another, destroying the cell-type-specific signal that is the project's premise.

---

## 4. Three locked design principles (make-or-break)

1. **Out-of-sample scoring.** Avoids the circularity that would otherwise inflate
   every case's alignment and the coupling estimate. **Preferred: cohort-split**
   (train the DE / composition directions on k−1 cohorts, score the held-out cohort),
   following Ruzicka's cross-cohort TPS — it also doubles as a reproducibility test.
   Within-cohort leave-one-donor-out is the fallback when only one cohort is available
   (e.g. the Xenium prototype).
2. **Test *within* SCZ** (controls = negative control, expected ≈ 0). The pooled
   case+control correlation is trivially positive because both scores separate
   case from control ("both happen in SCZ"); that is not the question. Equivalent:
   partial out diagnosis.
3. **Cell-type-specific signatures, aggregated at the score level** (not the
   signature level) — see the box in §3.

---

## 5. Statistical analysis

- **Primary.** Per cohort `k`: `r_k = cor(S_DE, S_comp)` over **SCZ** donors.
  Meta-analyse `r_k` across the 7 cohorts (Fisher-z, DerSimonian–Laird random
  effects) → pooled coupling estimate + I² heterogeneity. (Re-uses the existing
  DE meta / forest machinery.)
- **Negative control.** Same correlation within **control** donors (expect ≈ 0;
  a non-zero value flags a shared technical confound).
- **Context.** Pooled partial correlation of (S_DE, S_comp) controlling for
  diagnosis + covariates (shows the gross alignment, distinct from the within-SCZ
  test).
- **Latent-state companion (the un-collapsed view).** Following Ruzicka, run a matrix
  decomposition / NMF over the donor × cell expression to recover cross-cell-type
  **cellular states** (their Ex_SZCS / In_SZCS / Ex_SZTR), then test whether
  composition couples to those data-driven states — the principled alternative to
  collapsing to one `S_DE`. RV coefficient / CCA between the donor×type
  expression-deviation matrix and the donor×type CLR matrix (within SCZ) gives the
  same "shared-structure" answer; a 2-factor model tests one-vs-two latent axes.
- **External check (leverages the cohort overlap).** Correlate our `S_DE` against
  Ruzicka's published TPS for the shared McLean + MSSM donors, and relate the
  (S_DE, S_comp) plane to their TPS subgroups (SZ / SZ_CON-like) — does composition
  track, or cross-cut, their transcriptional split?
- **Genetic-risk test.** Correlate both scores (and the coupling) with PRS; per
  Ruzicka, expect TPS↔PRS but the heterogeneity *not* fully PRS-explained.
- **Robustness.** Repeat across: (a) projection vs cosine vs correlation scores;
  (b) DE-gene set (all tested / FDR-thresholded / |β|-weighted); (c) subclass vs
  supertype; (d) with/without `n_cells` weighting; (e) neuronal-only vs all-cell-type
  `S_DE` aggregation.

---

## 6. Confounds & mitigations

| Confound | Risk | Mitigation |
|---|---|---|
| Circularity (self-scoring) | inflates scores + coupling | **LOO / cross-fit** |
| Diagnosis axis | pooled corr trivially positive | **within-SCZ** test; partial out dx |
| Power = cell number | `S_DE` noisy for rare types; couples to abundance | reliability weighting; standardise per-type; include `n_cells` |
| Within-subclass composition → pseudobulk DE | DE & composition entangled at supertype level | subclass primary; supertype sensitivity; condition DE on supertype mix |
| Cohort batch | both scores co-vary by cohort | within-cohort scoring + meta; centre on cohort controls |
| Shared demographics/technical (age, sex, PMI) | spurious coupling | residualise both scores |
| snRNA composition unreliability (dissociation bias) | weak/biased compositional axis | **Xenium in-situ composition as orthogonal validation** |
| Compositional closure (∑ proportions = 1) | spurious anti-correlations among types | CLR (crumblr) |
| Ambient RNA shifts with composition | biases pseudobulk | note; optional ambient correction |

---

## 7. Visualisations

- **Primary:** 2-D scatter `S_DE` vs `S_comp`; points = donors, coloured by
  diagnosis; within-SCZ fit + 95% CI; quadrant labels; control cloud for
  reference. (The scatter *shape* is the answer.)
- **Per-cohort small multiples** + a **forest plot** of `r_k` with the pooled meta
  (re-use the composite's forest style).
- **Per-cell-type contribution** bars — which cell types drive each score.
- **Phase-1:** cell-type-level scatter — per-type DE magnitude vs crumblr comp
  effect (n ≈ 23), sized by cell count.

---

## 8. Staged execution

**Phase 1 — now, no new data: cell-type-level coupling (snRNA meta).**
Per subclass: DE magnitude (effect-size based, not gene count) vs crumblr
compositional effect. Built from the local meta summaries; answers "do the
transcriptionally-dysregulated cell types also shift in abundance?" and frames
whether the person-level dig is warranted. Sibling of
`scz_cell_type_enrichment/scripts/13_gwas_vs_composition.py` (swap GWAS enrichment
→ DE magnitude). Harmonise supertype `comp_beta` → subclass.

**Phase 2 — needs the per-donor pull: person-level disease-score coupling (the core).**
Build per-donor `S_DE`, `S_comp` per cohort → within-SCZ correlation → meta across
cohorts, per §§3–5.

**Phase 3 — robustness + validation + (optional) subtyping.**
Sensitivity analyses (§5); confirm the compositional axis in Xenium; explore the
2-D `(S_DE, S_comp)` patient phenotype against clinical / genetic variables.

---

## 9. Deliverables

- `results/figures/de_composition_coupling.{png,pdf}` — 2-D scatter + per-cohort forest.
- `results/tables/person_disease_scores.csv` — `donor, cohort, dx, S_DE, S_comp, covariates`.
- `results/tables/coupling_meta.csv` — per-cohort `r` + pooled meta + I².
- Phase-1: `results/tables/celltype_de_vs_composition.csv` + scatter.
- A short findings note (coupled / independent / trade-off + which cell types drive it).

---

## 10. Decisions to lock before coding

- **Score type:** projection onto LOO direction (lean) vs cosine vs correlation.
- **DE gene set:** all tested vs FDR-thresholded vs |β|-weighted.
- **Unit:** subclass primary (supertype sensitivity).
- **Cross-cell-type aggregation weights:** reliability / `n_cells` vs equal.
- **Centring:** within-cohort controls (recommended).

---

## 11. Open questions for Shreejoy

1. Phase-1 (cell-type-level, doable now) first, or wait and go straight to the
   person-level once the per-donor data lands?
2. Is a clean per-donor **crumblr CLR input** (not just the group summary)
   retrievable from the collaborator's pipeline, alongside per-donor pseudobulk?
3. Subclass-only, or worth the supertype sensitivity pass from the start (it is
   also the cleanest test of the within-subclass-composition confound)?
