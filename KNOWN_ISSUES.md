# Known issues

**Only items that must be settled before the paper is submitted.** Everything
the 2026-09-13 audit raised and has since been fixed was removed on 2026-09-15;
the two tables at the foot keep one line per removed item so the ~20 component
READMEs that cite issue numbers still resolve.

Numbers are stable and never reused. New items take the next free number.

---

## Must fix before submission

### 20. The S7 sensitivity rerun contradicts a stated result

**Found 2026-09-15, the day of the rerun. This is the one that changes text.**

Commit `5133be3` (Nicole, 2026-09-15) aligned the sensitivity analysis's `Sex`
term with the primary analysis — `(1|Sex)` → `(Sex)`, closing issue 1 — and
re-ran it. `Compositional_sensitivity_analysis/Files/meta_noSST_genes.csv` moved
with it. Effect sizes barely changed (max |Δβ| = 0.034 over 109 supertypes), but
two supertypes crossed FDR 0.10, and one of them is Sst_25.

The five Sst supertypes depleted in the primary analysis, in the **current** S7:

| supertype | β | p | FDR (new) | FDR (old) |
|---|---|---|---|---|
| Sst_2 | −0.268 | 3e-05 | 0.0013 | 0.0013 |
| Sst_3 | −0.180 | 0.0020 | 0.0371 | 0.0307 |
| Sst_22 | −0.244 | 0.0028 | 0.0443 | 0.0307 |
| Sst_20 | −0.212 | 0.0035 | 0.0476 | 0.0459 |
| **Sst_25** | **−0.205** | **0.0120** | **0.1122** | **0.0899** |

The Results text says, of this analysis, "**all five Sst supertypes depleted at
FDR < 0.10**". Against the current file that is **four of five**. `L2/3 IT_7`
also flipped (0.0867 → 0.1122), so the figure's FDR < 0.10 count goes 10 → 8.

Nothing about the conclusion changes: all five estimates stay negative, all five
are nominally significant, and all five remain FDR < 0.20 — which is the
threshold S7's own legend uses for its trend-level tier. Only the sentence is
now wrong.

**To do.** (a) Confirm with Nicole that the re-run is the version the paper
reports, since the previous S7 numbers are what the current text was written
against. (b) Then reword — "all five depleted, four at FDR < 0.10 and Sst_25 at
FDR = 0.11" or "all five at FDR < 0.20" both work. (c) Re-render S7, which has
no committed output (see issue 12), and re-check the legend's marker tiers
against the new file.

### 3. Which Xenium object Nicole's figures use — verified correct, provenance still unwritten

Flagged as the major outstanding item. **Checked on 2026-09-15: her analyses do
use the canonical object.** What remains is documentation, not a data risk.

`Final_figures/Code/Figure_3.r:166,360` reads `Data/Xenium_SCZ_R.rds`, a Seurat
conversion that no script in the repo produces and that is not committed. The
concern was that it might descend from the superseded object, which carries a
**transposed Lieber map** (cluster 9 ↔ 12, MGE and CGE swapped) and 373,514
rather than 356,313 neuronal cortical cells — an error that would silently
invert the S2a validation claim.

Evidence it does not, from the metadata export committed beside it
(`Data/xenium_metadata.csv`, git-lfs):

- **1,338,922 cells against the canonical object's 1,339,151.** The 229-cell gap
  is exactly the number of duplicate `obs_names` in the h5ad (1,339,151 total,
  1,338,922 unique), i.e. Seurat de-duplicating cell names on conversion, not a
  different cell population.
- **All 24 samples, Br2039 included.** No sample present in one and absent in the
  other.
- **Subclass counts agree to within 0.02%**, the residual tracking the same 229
  cells; obs schema matches `all_samples_annotated_clean.h5ad` column for column.
- `Data/xenium_crumblr_results_supertype_neuronal.csv`, which Figure 3 reads for
  its Xenium Sst_25 statistics, is **byte-identical** to the repo's own
  `spatial/output/crumblr/crumblr_results_supertype_neuronal.csv`
  (md5 `97c546a39016521db5ca58e92afd84ab`).

**To do.** Add a provenance line to `snrnaseq/Final_figures/README.md` recording
that `Xenium_SCZ_R.rds` is a Seurat conversion of
`SCZ_Xenium/output/all_samples_annotated_clean.h5ad` (md5
`d23458506f4f1d03a009b922d99cab2d`), that conversion drops 229 duplicate cell
names, and ideally commit the four-line conversion script. Correct the same
README's line 51, which currently warns the object may be the superseded one.
Decide whether Figure 3 panels d and e need to render from a clone; they cannot
today, because the `.rds` is not shipped.

### 21. The provenance verifier fails on a clean run

`shared/verify_provenance.py` reports `1 of 4 sets are stale`:

```
STALE  Figure 4 panel CSVs — panel_ad_concordance_sst.csv
       source: shared/snrnaseq_de/nicole_scz_snrnaseq_betas/final_results_crumblr_7_cohorts.csv
```

Self-inflicted by the issue-16 fix. The panel was exported on 2026-09-09 from the
original composition-betas file (`source_md5 c818add0…`); on 2026-09-13 that file
was replaced by the regenerated, git-tracked copy (`42df4fb2…`), so the recorded
checksum no longer matches.

The panel's numbers almost certainly do not move — `PROVENANCE.md` shows the
regenerated file reproducing three independent committed artefacts to ~1e-15 —
but that has not been demonstrated for this panel, and the repo currently ships a
verifier that fails.

**To do.** Re-run `genetics/scripts/figures/export_panel_ad_concordance.py`,
confirm the panel CSV changes only at floating-point noise, and commit the
refreshed MANIFEST row. Small, and it should be done last, after issue 20
settles whether any composition input moves again.

### 5. A cell-count covariate is computed but never enters the DE design

`snRNAseq_DE/Subclass/2_DE.r:85,202` and the supertype equivalent compute

```r
meta_subset$log2_cells <- scale(log2(meta_subset[[type]]))
```

but the design matrices (`2_DE.r:104,221`, `Supertypes/2_DE.r:112,242`) are
`~ scale(Age) + Sex + Diagnosis [+ scale(PMI)]`. `log2_cells` is not in them.

A reader will reasonably assume DE was adjusted for cells per donor. It was not.
A variant run evidently did adjust for it: `log2cells_Sst` appears in the subclass
exclusion list at `Supertypes/3_meta_analysis.r:15`.

**Not changed here** — deleting the line would be wrong if the adjustment was
intended, and adding it to the design would move every DE result, including
Figure 2. Nicole confirms which run the paper reports, then either the dead line
goes or the covariate comes back.

### 18. Two committed docs disagree about excluding Br2039

> The part of this issue that affected Figure 2 is gone: since 2026-09-14 the
> panel-i inset uses snRNA-seq nuclei proportions, not the Xenium snapshot. What
> is left is the contradiction, which a reviewer reading both files would see.

`spatial/methods_writeup.md:52` says Br2039 "is retained in all analyses (65%
white matter)". `spatial/all_samples_annotated_guide.md:312-322` says "**Exclude
Br2039 from cortical compositional analyses**" (41.8% WM) and gives a recipe that
does so. The code follows `methods_writeup`:
`spatial/code/modules/constants.py:31` is `EXCLUDE_SAMPLES = set()` and has been
since the monorepo was scaffolded. The two WM fractions differ because they count
different columns; by `spatial_domain` Br2039 is **70.1% WM** against a cohort
median of 19.7%.

**To do.** Pick one and make the other match, since the paper's Methods will
state the inclusion rule. Separately, `transcriptomic/data/figure_inputs/
crumblr_input_subclass_corr.csv` is still a pre-canonical 23-donor snapshot
(737,750 cells, no Br2039, source mtime 2026-03-12). Nothing in Figure 2 reads it
now; either refresh it from `spatial/output/crumblr/`, which is committed and
reproduces the canonical recipe exactly, or drop it from the manifest.

### 6. The committed donor metadata is person-level

`Compositional_analysis/7_cohorts_metadata_names.csv` holds donor ID, age, sex,
diagnosis and PMI for 469 donors. These are de-identified consortium IDs and the
repo is **private**, so nothing is exposed today. Needs an explicit decision
before the repo is shared with reviewers or made public.

---

## Deferred — real, but not blocking submission

One line each; these keep their numbers because component READMEs cite them.

| # | Issue | Where |
|---|---|---|
| 8 | `Figure_1a_UMAP.r` is Python, not R. Rename to `.py`. | `Final_figures/` |
| 9 | `1_crumblr_meta_nonNeurons` has no file extension; it is R. | `Compositional_analysis/Non_neurons/` |
| 10 | Label-transfer reference is named `counts_hodge`/`meta_hodge` but is the **SEA-AD** taxonomy (Gabitto 2024), not Hodge 2019. Misleads anyone checking methods against code. | 6 files under `Label_transfer/` |
| 11 | 24 of 34 scripts hard-code `/scratch/nendresz/…` and open with `setwd()` on directories that do not exist here. | throughout `snrnaseq/` |
| 12 | Nicole's figure scripts write to cluster paths, not `manuscript/figures/` under an S-number. Figs 1 and 3 and Supp S1, S4, S5, S7 have no rendered output in the repo. Blocks the S7 re-render in issue 20. | `Final_figures/` |
| 13 | Cell-type columns selected by hard-coded position (`colnames(meta)[1:24]`). Correct as written, silently wrong if a subclass set changes. | `snRNAseq_DE/*/2_DE.r` |
| 13b | Figure 3 panel variables are not named after their panel letters (`p3b` → panel c, etc.). Mapped in `Final_figures/README.md`. | `Final_figures/Figure_3.r` |
| 14 | A tracked symlink dangles in any fresh clone: `reserve/histology/coordinates`. | `reserve/histology/` |
| 17b | `2_Crumblr_analysis.r:78` hard-codes `coef = "DiagnosisSchizophrenia"`; fails on any matrix using the repo's `Control`/`SCZ` vocabulary. Derive it from the factor levels. | `Compositional_analysis/` |
| 19 | The DE chain reproduces in substance, not bit-for-bit: 12,490 of 12,492 Sst genes, r = 0.9968, *SST* at β −0.459/padj 0.048 vs committed −0.458/0.049. The residual is TMM drift from a different per-dataset gene universe. Exact reproduction needs Nicole's gene universes, not just her code. | `snRNAseq_DE/` |

---

## Closed

| # | Issue | Closed |
|---|---|---|
| 1 | `Sex` modelled as `(1\|Sex)` in the sensitivity analysis but `(Sex)` in the primary. Nicole aligned all three crumblr analyses on fixed-effect `Sex` and re-ran. **Consequences in issue 20.** | 2026-09-15, `5133be3` |
| 2 | The subclass per-gene DE meta-analysis script was missing. Committed as `snRNAseq_DE/Subclass/3a_meta_per_gene.r`; chain executes end to end. Exact-reproduction caveat is issue 19. | 2026-09-13 |
| 4 | `cluster_order_and_colors.csv` uncommitted and load-bearing (it defines which types are neuronal). Committed at `snrnaseq/`, byte-identical to source. | 2026-09-13 |
| 15 | Figure 2 would not render from a clone and had lost its staleness guard. All three cluster-only inputs are now committed snapshots with MANIFEST rows; the panel-i inset moved to snRNA-seq nuclei proportions. | 2026-09-14, `8364cc7` |
| 16 | Supplementary S8 read the git-ignored seam rather than a committed snapshot. Subclass falls back to the snapshot; the crumblr file was regenerated and tracked with a `PROVENANCE.md`. **Created issue 21.** | 2026-09-13 |
| 17 | Nicole's crumblr scripts assumed an NA-marked count matrix and would admit all-zero categories from a pivoted one. Fixed in all three; verified a no-op on her matrix. Residual hard-coded coefficient is 17b. | 2026-09-13 |
| — | `Supertypes/3_meta_analysis.r` did not parse (stray `-`, missing paren), so the committed copy was not the copy that ran. | 2026-09-13 |
| — | *HCN1* credible-set variants were labelled on the wrong strand; `export_panels_abc.py` now strand-aware with a `dist_to_hcn1_kb` column. | 2026-09-13 |
| — | Figure 2 panels d, h redrawn for review comment G#123; panel j label collisions resolved. | 2026-09-14, `2623216` |
