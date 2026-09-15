# Known issues

**Only items that must be settled before the paper is submitted.** Everything
resolved is listed at the foot, one line each, so the ~20 component READMEs that
cite issue numbers still resolve. Numbers are stable and never reused.

Last worked: **2026-09-15**. Two items remain, and both need a person rather than
a script: one sentence of the Results is now wrong, and one governance decision.

---

## Must fix before submission

### 20. One Results sentence no longer matches the S7 figure

Commit `5133be3` (Nicole, 2026-09-15) aligned the sensitivity analysis's `Sex`
term with the primary analysis, closing issue 1, and re-ran it. Effect sizes
barely moved — max |Δβ| 0.034 across 109 supertypes — but two supertypes crossed
FDR 0.10, and one is Sst_25.

The five Sst supertypes depleted in the primary analysis, in the **current** S7:

| supertype | β | p | FDR (new) | FDR (old) |
|---|---|---|---|---|
| Sst_2 | −0.268 | 3e-05 | 0.0013 | 0.0013 |
| Sst_3 | −0.180 | 0.0020 | 0.0371 | 0.0307 |
| Sst_22 | −0.244 | 0.0028 | 0.0443 | 0.0307 |
| Sst_20 | −0.212 | 0.0035 | 0.0476 | 0.0459 |
| **Sst_25** | **−0.205** | **0.0120** | **0.1122** | **0.0899** |

The Results text says this analysis found "**all five Sst supertypes depleted at
FDR < 0.10**". Against the current file that is **four of five**. `L2/3 IT_7` also
flipped (0.0867 → 0.1122), so the figure's FDR < 0.10 count goes 10 → 8.

**The figure is already correct and self-consistent.** S7 was re-rendered from the
current file on 2026-09-15 into
`manuscript/figures/supplementary/S07_composition_no_sst_de_genes.png`, where
Sst_25 now draws in the italic-red FDR 0.10–0.20 tier rather than bold red. Only
the sentence is out of date.

Nothing about the conclusion changes: all five estimates stay negative, all five
are nominally significant, and all five remain under FDR 0.20, the threshold S7's
own legend uses for its trend tier.

**To do**, and it needs people, not a re-run:

1. Confirm with Nicole that the re-run is the version the paper reports — the
   current text was written against the previous numbers.
2. Reword. "All five depleted, four at FDR < 0.10 and Sst_25 at FDR = 0.11" or
   "all five at FDR < 0.20" both work and both stay true.
3. The committed S7 render is PNG only; the renderer has no `.pdf` output, so a
   vector version is still needed for submission.

### 6. The committed donor metadata is person-level

`Compositional_analysis/Files/7_cohorts_metadata_names.csv` holds donor ID, age,
sex, diagnosis and PMI for 469 donors. The repo is **private**, so nothing is
exposed today. This needs an explicit decision before it is shared with reviewers
or made public.

Assessed 2026-09-15, to make the decision an informed one rather than a worry:

- Donor IDs are the **source consortia's own de-identified identifiers**
  (`AMPAD_HBCC_*`, `AMPAD_MSSM_*`, `MB*`, `s*`, `CON*`, numeric for Multiome),
  not anything generated here, and they are the identifiers those consortia
  already publish against.
- Ages run **18 to 69** with none at or above 90, because the study's own
  `Age < 70` filter removes the range that HIPAA safe-harbour treats specially.
- The columns are age, sex, diagnosis and PMI — the standard demographic table
  that accompanies published post-mortem snRNA-seq.

So the residual question is not re-identification risk but **whether the data-use
agreements behind PsychAD/AMP-AD, Batiuk, Fröhlich and MultiomeBrain permit
redistributing donor-level demographics**. That is a governance call, and only
the PI can make it.

---

## Deferred — real, but not blocking submission

One line each; these keep their numbers because component READMEs cite them.

| # | Issue | Where |
|---|---|---|
| 8 | `Figure_1a_UMAP.r` is Python, not R. Rename to `.py`. | `Final_figures/` |
| 9 | `1_crumblr_meta_nonNeurons` has no file extension; it is R. | `Compositional_analysis/Non_neurons/` |
| 10 | Label-transfer reference is named `counts_hodge`/`meta_hodge` but is the **SEA-AD** taxonomy (Gabitto 2024), not Hodge 2019. Misleads anyone checking methods against code. | 6 files under `Label_transfer/` |
| 11 | 24 of 34 scripts hard-code `/scratch/nendresz/…` and open with `setwd()` on directories that do not exist here. Rendering S7 needs that line neutralised. | throughout `snrnaseq/` |
| 12 | Nicole's figure scripts write to cluster paths, not `manuscript/figures/` under an S-number. **S7 was fixed on 2026-09-15**; S1, S4 and S5 still have no rendered output here. | `Final_figures/` |
| 13 | Cell-type columns selected by hard-coded position (`colnames(meta)[1:24]`). Correct as written, silently wrong if a subclass set changes. | `snRNAseq_DE/*/2_DE.r` |
| 13b | Figure 3 panel variables are not named after their panel letters (`p3b` → panel c, etc.). Mapped in `Final_figures/README.md`. | `Final_figures/Figure_3.r` |
| 14 | A tracked symlink dangles in any fresh clone: `reserve/histology/coordinates`. | `reserve/histology/` |
| 17b | `2_Crumblr_analysis.r:78` hard-codes `coef = "DiagnosisSchizophrenia"`; fails on any matrix using the repo's `Control`/`SCZ` vocabulary. Derive it from the factor levels. | `Compositional_analysis/` |
| 22 | The Xenium **sensitivity** crumblr variants (`_corr`, `_hybrid`, `_margin_*`, `_pctl*`, `_no_high_umi`) are built without Br2039, 23 donors, while the four the paper uses (`_neuronal`, `_nonneuronal`, and the two pooled) have all 24. No paper figure reads a sensitivity variant, so nothing is wrong today; regenerate them before anyone does. | `spatial/output/crumblr/` |

---

## Closed

| # | Issue | Closed |
|---|---|---|
| 1 | `Sex` modelled as `(1\|Sex)` in the sensitivity analysis but `(Sex)` in the primary. All three crumblr analyses now use fixed-effect `Sex`. **Consequences in issue 20.** | 2026-09-15, `5133be3` |
| 2 | The subclass per-gene DE meta-analysis script was missing. Committed as `snRNAseq_DE/Subclass/3a_meta_per_gene.r`; the chain executes end to end. | 2026-09-13 |
| 3 | `Xenium_SCZ_R.rds` had no provenance and no producer, raising the risk that Figure 3 descended from the object with MGE and CGE transposed. **It does not.** Its committed metadata export holds 1,338,922 cells against the canonical 1,339,151 — a gap of exactly the 229 duplicate `obs_names` Seurat de-duplicates on conversion — with all 24 samples including Br2039 and subclass counts agreeing to within 0.02%; `xenium_crumblr_results_supertype_neuronal.csv` is byte-identical to the repo's own. Panels d, e and h use only metadata, so `Figure_3.r` now loads `Data/xenium_metadata.csv` when the `.rds` is absent and renders from a clone. | 2026-09-15 |
| 4 | `cluster_order_and_colors.csv` uncommitted and load-bearing (it defines which types are neuronal). Committed at `snrnaseq/`, byte-identical to source. | 2026-09-13 |
| 5 | A `log2_cells` covariate was computed but never entered the DE design, leaving it unclear whether DE was adjusted for cells per donor. **It was not, by design.** The dead line was removed in `c106186`, and re-running the Sst model without it reproduces the committed per-dataset table exactly for the three symbol-space datasets (McLean, MSSM 1, Multiome: same gene sets, max \|ΔlogFC\| 4.5e-14), while adding it drops agreement to r = 0.68–0.81. | 2026-09-15 |
| 15 | Figure 2 would not render from a clone and had lost its staleness guard. Inputs are committed snapshots with MANIFEST rows; the panel-i inset moved to snRNA-seq nuclei proportions. | 2026-09-14, `8364cc7` |
| 16 | Supplementary S8 read the git-ignored seam rather than a committed snapshot. Subclass falls back to the snapshot; the crumblr file was regenerated and tracked with a `PROVENANCE.md`. **Created issue 21.** | 2026-09-13 |
| 17 | Nicole's crumblr scripts assumed an NA-marked count matrix and would admit all-zero categories from a pivoted one. Fixed in all three; verified a no-op on her matrix. Residual hard-coded coefficient is 17b. | 2026-09-13 |
| 18 | Two committed docs contradicted each other on excluding Br2039, and a pre-canonical 23-donor Xenium snapshot sat in `figure_inputs/`. `all_samples_annotated_guide.md` now matches the code and the Methods — Br2039 is retained, and the cortical gate removes its white matter cell by cell — with the WM fraction stated per column (70.1% by `spatial_domain` against a 19.7% median, 41.8% by `layer` against 14.1%). The orphan snapshot, which nothing had read since 2026-09-14, was dropped from the refresh spec and deleted; the manifest rebuilt to 8 entries and Figure 2 still renders. Successor: issue 22. | 2026-09-15 |
| 19 | "The DE chain reproduces in substance but not bit-for-bit." Narrower than stated: it reproduces **exactly** wherever a dataset's gene universe is recoverable (the three symbol-space datasets, max \|ΔlogFC\| 4.5e-14). Datasets needing an identifier remap do not, because dropping unmappable genes moves the >= 80% expression filter and every library's TMM factor — a property of the rebuild, not the pipeline. | 2026-09-15 |
| 21 | `shared/verify_provenance.py` failed on a clean run: the issue-16 fix replaced the composition-betas seam, invalidating the checksum for the Figure 4 panel exported from it. Re-running the exporter reproduced the panel to 1e-15, confirming the regenerated seam. The exporter now sorts deterministically, since row order otherwise inherits the input's and churns the file for no reason. All four manifest sets pass. | 2026-09-15 |
| — | `Supertypes/3_meta_analysis.r` did not parse (stray `-`, missing paren), so the committed copy was not the copy that ran. | 2026-09-13 |
| — | *HCN1* credible-set variants were labelled on the wrong strand; `export_panels_abc.py` now strand-aware with a `dist_to_hcn1_kb` column. | 2026-09-13 |
| — | Figure 2 panels d, h redrawn for review comment G#123; panel j label collisions resolved. | 2026-09-14, `2623216` |
