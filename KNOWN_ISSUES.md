# Known issues

**Only items that must be settled before the paper is submitted.** Everything
resolved is listed at the foot, one line each, so the ~20 component READMEs that
cite issue numbers still resolve. Numbers are stable and never reused.

**No blocking items remain.** What is left below is real but does not stop
submission.

---

## Data that is deliberately not distributed

### 24. Donor-level clinical records

Per-cell metadata dumps must carry **only study variables**. Two classes of file
were removed from this repository and from its history because they carried the
donating brain banks' clinical records — cause and manner of death, psychiatric
and neuropathology narrative, medication and substance-use flags, and exact lab
dates — for identifiable donors:

- the 299-column per-cell metadata dump for the Batiuk cohort, which
  `Compositional_sensitivity_analysis/Code/1_Label_transfer_noSSTDE.r` produced
  by writing `Seu_sn@meta.data` whole. That writer now selects only the columns
  `2_Compile_metadata.r` consumes (`Donor, Age, Sex, Diagnosis, PMI,
  predicted.id`, plus Multiome's three spellings), and the committed file holds
  exactly those. The other six cohorts' objects never carried such columns.
- the two `reserve/histology/data/` inputs behind an analysis that is not in the
  paper. The code and the derived, de-identified analysis table remain; the raw
  inputs are available from the source lab on request.

**If you add a per-cell or per-donor export, select columns explicitly.** Never
write a metadata slot whole — what an upstream object carries is not under this
repository's control.

---

## Deferred — real, but not blocking submission

One line each; these keep their numbers because component READMEs cite them.

| # | Issue | Where |
|---|---|---|
| 8 | `manuscript/figures/main/README.md` names `Figure_1a_UMAP.r`; the file is `Final_figures/Code/Figure_1.py`. | `Final_figures/` |
| 10 | Label-transfer reference is named `counts_hodge`/`meta_hodge` but is the **SEA-AD** taxonomy (Gabitto 2024), not Hodge 2019. Misleads anyone checking methods against code. | 6 files under `Label_transfer/` |
| 11 | Of 41 scripts under `snrnaseq/`, 19 hard-code a `/scratch/` or `/project/` path and 31 open with `setwd()` on a directory that does not exist here — and they name **five different** project roots (`P1_SCZ_paper`, `scz_celltype_paper`, `P1_SCZ_DE_fresh`, `OFC_cohort`, `PsychAD`), so a reader cannot tell which tree a result came from. Rendering S7 needs the `setwd()` neutralised. | throughout `snrnaseq/` |
| 12 | Nicole's figure scripts write to cluster paths, not `manuscript/figures/` under an S-number. S7 is fixed; S1, S4 and S5 still have no rendered output here. | `Final_figures/` |
| 13 | Cell-type columns selected by hard-coded position (`colnames(meta)[1:24]`). Correct as written, silently wrong if a subclass set changes. | `snRNAseq_DE/*/2_DE.r` |
| 13b | Figure 3 panel variables are not named after their panel letters (`p3b` → panel c, etc.). Mapped in `Final_figures/README.md`. | `Final_figures/Figure_3.r` |
| 14 | A tracked symlink dangles in any fresh clone: `reserve/histology/coordinates`. | `reserve/histology/` |
| 17b | `coef = "DiagnosisSchizophrenia"` is hard-coded at `Compositional_analysis/Code/Neurons/1_Crumblr_analysis.r:84`, `Code/Non_neurons/1_Crumblr_analysis.r:53` and `Compositional_sensitivity_analysis/Code/3_Crumblr_analysis.r:52`; each fails on any matrix using the repo's `Control`/`SCZ` vocabulary. Derive it from the factor levels. | `snrnaseq/` |
| 22 | `spatial/output/crumblr/README.md` **is** tracked and documents variant files (`_corr`, `_hybrid`, `_margin_*`, `_pctl*`, `_no_high_umi`, `_all*`) that are **not** in the repo — only the four canonical inputs and their results are committed. A reviewer reading that table will look for files that are not there. Separately, on the authors' machines those variants are built without Br2039 (23 donors) while the committed four have all 24; regenerate them before anyone uses one. | `spatial/output/crumblr/README.md` |

---

## Closed

| # | Issue | Closed |
|---|---|---|
| 1 | `Sex` modelled as `(1\|Sex)` in the sensitivity analysis but `(Sex)` in the primary. All three crumblr analyses now use fixed-effect `Sex`. **Consequences in issue 20.** | 2026-09-15, `5133be3` |
| 2 | The subclass per-gene DE meta-analysis script was missing. Committed as `snRNAseq_DE/Subclass/3a_meta_per_gene.r`; the chain executes end to end. | 2026-09-13 |
| 3 | `Xenium_SCZ_R.rds` had no provenance and no producer, raising the risk that Figure 3 descended from the object with MGE and CGE transposed. **It does not.** Its committed metadata export holds 1,338,922 cells against the canonical 1,339,151 — a gap of exactly the 229 duplicate `obs_names` Seurat de-duplicates on conversion — with all 24 samples including Br2039 and subclass counts agreeing to within 0.02%; `xenium_crumblr_results_supertype_neuronal.csv` is byte-identical to the repo's own. Panels d and e use only metadata, so `Figure_3.r` now loads `Data/xenium_metadata.csv` when the `.rds` is absent and renders from a clone. | 2026-09-15 |
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
| 23 | S8 drew different numbers here than on a clone. `sst_strata/_common.R` preferred the git-ignored seam `shared/snrnaseq_de/DE_genes_all_cells_scz.csv` and fell back to the committed snapshot, on a comment's claim that the two were byte-identical. They were not: the seam left on this machine predates Nicole's 2026-09-13 rerun (222,028 rows, 345 Sst genes at FDR < 0.10) while the snapshot matches her current output (231,135 rows, 343 genes). Fixed by always reading the snapshot; S8 re-rendered (the change is confined to the `Sst_subclass` column of panel e) and `09_verify.R` passes. | 2026-09-15 |
| 6 | The committed donor metadata is person-level (469 donors: consortium ID, age, sex, diagnosis, PMI). Assessed as low re-identification risk — IDs are the source consortia's own, ages run 18–69 with none near the safe-harbour threshold, columns are the standard demographic table — leaving only whether the data-use agreements permit redistribution. **PI authorised release of a donor-level metadata sheet with this demographic content.** | 2026-09-15 |
| 20 | Nicole's `Sex`-alignment rerun moved Sst_25 to FDR 0.1122 in S7, contradicting a Results sentence that claimed all five depleted Sst supertypes at FDR < 0.10. Text updated to "all five Sst supertypes depleted at FDR < 0.20", verified true (max is Sst_25 at 0.1122), and the gene count corrected 345 → 343, also verified: the committed meta table gives exactly 343 unique Sst genes at padj < 0.10 (196 down, 147 up, *SST* included), and `1_Label_transfer_noSSTDE.r:8` derives the exclusion list from that same table, which was re-run in the same commit. S7 re-rendered from the current results to `manuscript/figures/supplementary/S07_composition_no_sst_de_genes.{png,pdf}`. | 2026-09-15 |
