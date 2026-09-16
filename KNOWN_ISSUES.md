# Known issues

**Only items that must be settled before the paper is submitted.** Everything
resolved is listed at the foot, one line each, so the ~20 component READMEs that
cite issue numbers still resolve. Numbers are stable and never reused.

**No blocking items remain.** What is left below is real but does not stop
submission.

---

## Deferred — real, but not blocking submission

One line each; these keep their numbers because component READMEs cite them.

| # | Issue | Where |
|---|---|---|
| 11 | Of 42 scripts under `snrnaseq/`, 31 call `setwd()` and 18 hard-code a `/scratch/` or `/project/` path, naming **8 different** top-level roots, so a reader cannot tell which tree a result came from. **Decided 2026-09-16: left as-is and documented** — rewriting the paths would make the archived code differ from the code that produced the numbers. The counts and the full root table are in `snrnaseq/README.md`. | throughout `snrnaseq/` |
| 12 | Nicole's figure scripts write to cluster paths, not `manuscript/figures/` under an S-number. S7 is fixed; S1, S4 and S5 still have no rendered output here. | `Final_figures/` |
| 25 | `Compositional_analysis/Files/7_cohorts_metadata_names.csv` spells Diagnosis three ways: `Control` (293), `Schizophrenia` (171) and **`control`** (5). All five lowercase rows are Multiome. Nothing is wrong today — every model is fitted one dataset at a time, so each fit sees a clean two-level factor, and `composition_sensitivity/00_prepare_counts.R` normalises case anyway. But any **pooled** fit over all 469 donors gets a three-level `Diagnosis`, and `model.matrix` then emits both `DiagnosisControl` and `DiagnosisSchizophrenia`, so a contrast picked by position silently becomes the wrong one. The issue-17b fix guards against this with `stopifnot(length(dx_coef) == 1)` rather than papering over it. Fixing the CSV itself would change a committed input, so it is left for a deliberate decision. | `snrnaseq/Compositional_analysis/Files/` |

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
| 8 | `manuscript/figures/main/README.md` named `Figure_1a_UMAP.r`; the row now names `Code/Figure_1.py` and `Code/Figure_1a_Barchart.r`. Every script the main README cites resolves. | 2026-09-16 |
| 10 | Label-transfer variables renamed `*_hodge` -> `*_seaad` across all 6 files (88 occurrences, all local variables — "hodge" never named a path), and each file gained a header naming the SEA-AD (Gabitto 2024) reference and why the old name was wrong. | 2026-09-16 |
| 13 | Positional cell-type selection replaced by `setdiff(colnames(meta), META_COLS)` at all 8 sites. Identical for a 24-subclass table; it also removes the latent bug where a 23-subclass dataset took `Donor` as a cell type, which the Multiome branch's `[1:24]` would have done. | 2026-09-16 |
| 13b | Figure 3 panel variables renamed to match their panel letters (`p3c`<->`p3b`, `p3h`->`p3e`). Verified the assembly draws the same plots in the same order; `Final_figures/README.md` updated. | 2026-09-16 |
| 14 | The dangling `reserve/histology/coordinates` symlink was removed. `COORD_DIR` in `code/config.py` is defined but never read; `reserve/histology/README.md` now records that and where the exports came from. | 2026-09-16 |
| 17b | The Diagnosis coefficient is now read off the fit rather than hard-coded — at **7** sites, not the 3 this issue listed (`snRNAseq_DE/{Subclass,Supertypes}/2_DE.r` have two each). Re-running both crumblr analyses reproduces the committed results to 4.3e-15, so no number moved. A `stopifnot(length(dx_coef) == 1)` guard was added: see issue 25. | 2026-09-16 |
| 22 | `spatial/output/crumblr/README.md` retitled and reworded so the variant table reads as a filename key for the pipeline machines, not an inventory of the directory, and points at `git ls-files` for what is actually present. | 2026-09-16 |
