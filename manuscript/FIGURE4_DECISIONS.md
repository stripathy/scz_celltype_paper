# Figure 4 — decisions log and open items

Status as of 2026-08-07. Everything described here is built and rendered in the
session scratchpad (`magma_pgc3/fig4_bigdeli/`, PNG + PDF + SVG). **Nothing has
been written into the repo**: `genetics/results/figures/r_panels/` is unchanged,
so none of this is committed until we decide to land it.

---

## 1. Why the figure was rebuilt

An audit of the enrichment pipeline found two independent bugs, both in our code,
neither in the GWAS:

1. **The gene-property step was an in-house OLS reimplementation**, not MAGMA's own
   `--gene-covar` analysis. It omitted MAGMA's defaults — in particular the 5-SD
   truncation of specificity values — and disagreed with real MAGMA on **117 of 461**
   cell types in a like-for-like test.
2. **The Franken-503 specificity matrix was built differently from the published
   recipe**: CP10K-normalised means and symbol-space deduplication, where Duncan et al.
   use per-cell `ln(1+raw)` averaged within type, ENSEMBL identifiers, unique
   ENTREZ↔ENSEMBL mappings, and an MHC-excluded gene-loc join performed last.

Consequence of (1): a single gene, **AMY1B**, produced a spurious result. It carries the
highest Sst_25 specificity of all 503 types (0.4123) on a MAGMA statistic resting on
**one SNP**, and was present in Bigdeli but absent from PGC3. It alone drove Sst_25's
apparent collapse (β 4.47 → 31.72 when removed) and killed both SST sub-claims. MAGMA's
native truncation clips it automatically; the OLS version did not.

**Validation.** Regenerating the PGC3 gene-level analysis from scratch reproduces Duncan's
published output exactly (18,481/19,427 genes annotated; ZSTAT r = 0.99955). Running our
implementation on *their* taxonomy reproduces their published enrichment at **β r = 0.9915**,
with 17/461 cell types differing in FDR significance — residual attributable to
`org.Hs.eg.db` 3.21 vs their 3.12 changing which genes have unique ENTREZ↔ENSEMBL mappings.

---

## 2. Decisions made

### 2.1 Expression reference — SEA-AD **DLPFC**, not MTG

Specificity now comes from the SEA-AD A9/DLPFC release, using its **three neurotypical
reference donors** (H18.30.002, H19.30.001, H19.30.002; 90,579 nuclei).

*Rationale.* Every disease measurement in the paper is frontal — the seven snRNAseq
datasets (Fröhlich is OFC), the Xenium cohort, and the SEA-AD CPS analysis. Specificity
was the only temporal quantity, so panel a was correlating frontal-derived depletion
against temporal-derived specificity.

*Strength of the swap.* Three of the five MTG reference donors **are** these three
individuals, so this is close to a within-donor region swap rather than a change of
cohort. Every claim strengthens while using roughly **half** the Sst nuclei (6,829 vs
12,393) — signal, not power.

| | MTG | DLPFC |
|---|---|---|
| Depth gradient (PGC3) | ρ = −0.600, p = 0.0140 | **ρ = −0.759, p = 0.0007** |
| Convergence (PGC3) | ρ = +0.359, p = 0.173 | **ρ = +0.635, p = 0.0098** |
| GABAergic fraction of sig. types | 83% | 85% |

*Costs.* 125 SEA-AD supertypes instead of 137 (the 12 missing are all non-neuronal or
vascular, plus L5 ET_1 — none GABAergic). Panels g/h had to move from the donor-paired
test to cell-level, because 3 donors cannot support genome-wide donor-level FDR where 5 could.

### 2.2 GWAS — **Bigdeli 2026, European-ancestry**

Files are `bigdeli.*` (EUR, max effective N 114,827), **not** `multi.*` (AFR+EUR+EAS,
164,509).

*Rationale.* Bigdeli EUR edges PGC3 on both SST claims (convergence p = 0.0069 vs 0.0098;
depth p = 0.0001 vs 0.0007) and is the larger European meta-analysis. It also uses the
`g1000_eur` LD panel we have.

*Trans-ancestry was tested and rejected.* Run through the validated pipeline it is
**worse** on both claims and less specific — convergence ρ = 0.579, depth ρ = −0.721, and
78 significant types at only 72% GABAergic (vs 85% for PGC3). Doing it correctly would
need per-ancestry summary statistics and AFR/EAS LD panels, neither of which we have.

*Caveat to state.* "More people" does not mean more power here: Bigdeli EUR has 1.57×
PGC3's effective N but only 1.07× the polygenic signal — **0.68× signal per N** — because
the added cohorts are biobank/EHR-ascertained. HCN1's gene-level statistic actually
*drops* (−log₁₀P 11.85 → 10.86) even as the lead SNP becomes more significant.

*Supplement.* PGC3 EUR is retained as the robustness comparison (enrichment β correlates
at r = 0.97 between the two GWAS).

### 2.3 Fine-mapping — Bigdeli **SuSiE-R, EUR**

Credible set from Bigdeli Supplementary Table 13, LOCUS319: **4 variants, cumulative
PIP 0.985, lead rs10035564 PIP 0.556**. Against PGC3's 9 variants / cumulative 0.95 /
lead 0.525 — same lead SNP, tighter set.

*Why this analysis specifically.* Their tables report rs10035564 at PIP **0.001 to 0.841**
across six method/ancestry combinations. SuSiE-R EUR is ancestry-matched to the sumstats
used for enrichment and is the like-for-like swap for PGC3's EUR FINEMAP. The full range
belongs in the supplement. Quoting the trans-ancestry 0.841 while running EUR enrichment
would be selecting the best number after seeing all six.

*Incidental gain.* Bigdeli is natively hg38, so the hg19→hg38 liftover the PGC3 build
required is gone.

### 2.4 Taxonomy

500 types = 125 SEA-AD DLPFC supertypes + 375 non-redundant Siletti whole-brain clusters.
The SEA-AD↔Siletti matching was recomputed on the new normalisation and is **stable**:
still 95 reciprocal best hits, 88 identical to the previous set. Every downstream result
is invariant to the rebuild.

**A metric bug was found and must not be repeated in the methods.** The codebase's
"AUROC" is a row-wise rank statistic that equals **1.0 by construction for every best
hit** (verified: 137/137 at exactly 1.0). The documented "all mean AUROC ≥ 0.9" bar
therefore cannot fail and is not evidence of redundancy. The informative measure is the
underlying Spearman r: median 0.901, **minimum 0.336**. Eight merges fall below r = 0.8
and are **all microglial or vascular** — the worst being Micro-PVM_1_1-SEAAD ↔ Mono_3
(r = 0.336, microglia merged with monocytes). Un-merging those eight changes nothing
downstream. The metric is also a percentile-rank-of-correlation, **not** MetaNeighbor
neighbour-voting AUROC, so the Crow et al. 2018 citation overstates it.

### 2.5 Panel composition — nine panels

The HCN1-expression-vs-depletion panel (old **d**) was dropped to give the morphology
panel 2.4× more width (0.76 → 1.80 of its row). An HCN1 violin was built as a possible
replacement and **rejected** as duplicative of the sag panel, which is more comprehensive.

| Panel | Content | Key statistic |
|---|---|---|
| a | genetics ↔ depletion, 16 Sst | ρ = 0.659, p = 0.00685 |
| b | Sst_3 gene drivers | 345 drivers; HCN1 at 93.5th pct, P = 1.4 × 10⁻¹¹ |
| c | HCN1 locus zoom | 4-variant credible set, lead PIP 0.556 |
| d | HCN1 expression ↔ sag | ρ = 0.65, p = 0.00786 |
| e | five-cell morphology series | depleted / not-depleted bracketed |
| f | five-cell sag traces | 6.2 → 4.0 → 3.5 → 1.3 → 0.3 mV |
| g | depleted vs not-depleted volcano | 580 genes, FDR < 0.05, \|log₂FC\| > 0.25 |
| h | CALB1 violin | FDR < 0.001 |
| i | SCZ ↔ AD concordance | ρ = 0.87 |

### 2.6 Exemplar cells (panels e/f)

Depth-ordered, three depleted then two not:

| Supertype | Specimen | Depth | Sag | Group |
|---|---|---|---|---|
| Sst_25 | 819770858 | 406 µm | 0.564 | Depleted |
| Sst_22 | **907585117** | 517 µm | 0.363 | Depleted |
| Sst_3 | 1037461069 | 1109 µm | 0.357 | Depleted |
| Sst_5 | 758996755 | 1500 µm | 0.002 | Not depleted |
| Sst_1 | 797048104 | 1734 µm | 0.036 | Not depleted |

Sst_22 cell 907585117 was chosen over 894528079 because the latter sits at essentially the
same depth as Sst_25 (401 vs 406 µm), which broke the depth series. **Sst_2 has zero
reconstructions** and cannot serve as an e/f exemplar. NWBs absent from the upstream cache
were pulled from DANDI dandiset 000636.

### 2.7 Volcano threshold

Loosened from |log₂FC| > 0.5 to **> 0.25** so HCN1 (log₂FC = +0.434) is called
differentially expressed; it was failing on effect size alone, with FDR = 1.9 × 10⁻¹³⁸.
Cost: 183 → 580 called genes. Note **0.40 would also admit HCN1 at only 271 genes** if the
inflation is a concern.

Panels g/h use the **cell-level** Wilcoxon rather than the donor-paired test, because 3
DLPFC donors cannot support genome-wide donor-level FDR (only 1 gene survives, vs 234 in
MTG at n = 5). The donor-paired p-values remain nominally significant (CALB1 p = 0.0030,
HCN1 p = 0.0119) but do not survive correction across 36,601 genes at n = 3.

### 2.8 Cosmetic decisions

- **HCN1 recoloured** from blue (#1565C0) to CALB1's #693d07 in panels b, c and g — blue
  tied HCN1 to nothing else in the figure. HCN1 and CALB1 now read as a matched pair.
- **Panel c x-axis reversed** so HCN1's minus strand reads 5'→3' left-to-right with the
  3' end on the right. Implemented with `scale_x_reverse()` on both stacked subplots —
  `coord_cartesian` silently sorts its limits and will not flip the axis.
- **Transcription-direction chevrons** added to the gene model (7 marks, following the
  `strand` column, so they generalise).
- **PIP legend relabelled** "SuSiE-R PIP" (was "FINEMAP PIP", which named the wrong method).
- Panel e cells evenly spaced by centre-to-centre pitch; depleted/not-depleted brackets
  added beneath; panel f trace labels hand-placed with leaders.

---

## 3. Abandoned

**Mapping Siletti-461 enrichment onto Sst supertypes** for a taxonomy-independent
genetics↔depletion plot. The mapping is **16 → 10, not one-to-one** (MGE_242 is the best
match for Sst_11, Sst_20 *and* Sst_22; MGE_254 for Sst_5, Sst_4 and Sst_7). The all-16
version is pseudo-replication; the honest reciprocal-best-hit subset gives ρ = 0.60,
p = 0.073 at n = 10, and drops two of the most depleted supertypes. Not worth the space.

---

## 4. Still to decide

### 4.1 Text reconciliation — required, not optional

- **Panel letters shift**: old d dropped, then e→d, f→e, g→f, h→g, i→h, j→i.
- **"including Sst_2, Sst_25, and Sst_20"** is now wrong. Under Bigdeli/DLPFC the top Sst
  supertypes are Sst_23, Sst_3, Sst_2, Sst_20 — **Sst_25 is rank 6 of 16**.
- **Panel b is Sst_3**, not Sst_25. Either revert the panel or change the text.
- **"HCN1 mRNA expression was itself more highly expressed among the more depleted Sst
  supertypes (ρ = 0.79)"** has no panel behind it any more. Restore the panel, or cut the
  sentence and accept a weaker chain from HCN1 to depletion.
- **Numbers**: ρ = 0.66 / p = 0.0069 · HCN1 −log₁₀P = 10.86 · 4-variant credible set,
  cumulative PIP 0.985, lead 0.556 · 580 DE genes · CALB1 log₂FC = +0.546.
- **Check the Duncan citation**: the text says "specific cortical Sst subtypes carry the
  strongest genetic association of any brain cell type." In our reproduction of their
  461-cluster analysis the top hits are MGE/CGE-derived interneuron clusters. Confirm
  against their results text.

### 4.2 Which exemplar carries panel b

Sst_3 is currently in panel b while Sst_25 carries e/f. Under Bigdeli, Sst_2 (rank 1) and
Sst_3 (rank 2) both outrank Sst_25 (rank 6). Options: (i) keep Sst_25 throughout and stop
implying it is the most enriched; (ii) feature the top-enriched type in b and Sst_25 in
e–g, stating plainly that they are two different cells doing two different jobs. Either
works; what does not work is implying Sst_25 is both.

### 4.3 To declare in the legend

- **Three palette deviations**: Sst_22 → #9c5f14, Sst_3 → #c8862e, Sst_1 → #f5c77e.
  SEA-AD's native Sst tones cannot separate five cells (it assigns Sst_25/Sst_22 nearly
  identical dark browns and Sst_3/Sst_5 nearly identical light oranges).
- **Region mismatch**: patch-seq (panels d, e, f) is **temporal** MTG, while specificity is
  frontal. The DLPFC swap fixes the mismatch in panel a and creates a smaller one here.
- **Depleted/not-depleted grouping** in panel e comes from the Fig-3 compositional
  analysis, not from anything computed in that panel.
- Panels g/h are cell-level statistics over ~6,800 nuclei from **3 donors**.
- The credible set is **SuSiE-R EUR**; the enrichment is EUR MAGMA. Same GWAS, different
  analyses within it.

### 4.4 Code hygiene before landing

- `HCN1_COLOR` is defined **twice** — `fig4_style.R:26` and `scz_sst_hcn1_story.R:101`,
  where the second shadows the first. Collapse to one definition.
- The 125-type landscape hardcodes `bonf_y <- -log10(0.05 / 503)`; the taxonomy is now
  **500**. Read the denominator from the data.
- `export_for_R.py` is a script, not a module — importing it re-runs the full export over
  the live `r_panels/` directory. Do not import it.
- Panel h–j builders read `FIG4_PANELS` (hardcoded to the repo) rather than `DATA_DIR`,
  so they ignore any redirect applied to the main story script.

### 4.5 Supplementary figures

Two are built (`fig4_bigdeli/supp_landscape500.*`, `supp_siletti461.*`): the 500-type
combined landscape (238 FDR-significant, 80 Bonferroni at P < 1 × 10⁻⁴) and the standalone
Siletti-461 landscape (220 / 110). Both need decisions on: whether to keep the 91-cluster
"Splatter" block, how to handle **Misc_132** — the single strongest cluster in both
figures, sitting in an unresolved supercluster and matching only weakly to cortical types
(best r = 0.746 to L4 IT_2) — and the clipped supercluster labels along the bottom axis.

### 4.6 Irreducible

Exact reproduction of Duncan et al. requires `org.Hs.eg.db` from Bioconductor 3.12 (2020).
On the current 3.21, 1,323 gene-loc genes are dropped for multi-mapping where fewer were
in 2020, which is the source of the 282-gene difference and the 17/461 disagreements. Any
modern rerun of their pipeline — theirs or ours — lands a few hundred genes away from the
published numbers. Worth one sentence in the methods.
