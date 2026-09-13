# Known issues

Open items found in a repo audit on **2026-09-13**, after the snRNA-seq pipeline
(`snrnaseq/`) was merged in. Kept here so they can be triaged in one place rather
than rediscovered by a reviewer.

Two categories:

- **Needs a decision** — something that could change a number, or a gap in a
  chain behind a figure. Resolve before submission.
- **Cosmetic / naming** — nothing depends on it being fixed; it just costs a
  reader time. Safe to batch.

Fixed items are listed at the bottom so the audit is auditable. Numbers are
stable and are cited from the component READMEs, so a resolved item leaves a gap
rather than causing a renumber.

---

## Needs a decision

### 1. `Sex` is modelled differently in the primary and the sensitivity composition analyses

The three crumblr analyses do not use the same formula:

| Analysis | Builds | `Sex` term | Line |
|---|---|---|---|
| `Compositional_analysis/2_Crumblr_analysis.r` | **Fig. 3a** (primary) | `+ (Sex)` — fixed effect | 66, 68 |
| `Compositional_analysis/Non_neurons/1_crumblr_meta_nonNeurons` | S5 | `+ (1\|Sex)` — random effect | 45 |
| `Compositional_sensitivity_analysis/3_Crumblr_analysis.r` | S7 | `+ (1\|Sex)` — random effect | 41 |

`(Sex)` in an R formula is just `Sex`, so the primary analysis fits Sex as a
fixed effect while the two supporting analyses fit it as a random intercept. The
bare parentheses look like a `(1|Sex)` that lost its `1|`.

Either convention is defensible for a two-level factor — a fixed effect is
arguably the better choice — so this is **not necessarily wrong**. What is a
problem is that S7 is presented as a sensitivity analysis of Fig. 3a while
differing from it in a second respect beyond the one it is meant to test. A
reviewer comparing the two scripts will see it.

**Not changed here**, because re-fitting would move every composition estimate in
Fig. 3 and everything downstream of it. Nicole's call: either align the three
formulas and re-run, or state the fixed-effect choice explicitly in the methods.

### 2. The subclass-level per-gene DE meta-analysis script is missing

`snRNAseq_DE/Subclass/3_meta_analysis.r:22` reads `meta_results_*.csv`, but
nothing in `Subclass/` writes them. The per-gene `rma()` pooling exists only in
`snRNAseq_DE/Supertypes/3_meta_analysis.r:42`, which explicitly *excludes*
subclasses (line 17, `setdiff(cell_types, subclasses)`).

So the step producing `DE_genes_all_cells_scz.csv` — the primary input to
**Figure 2** — is not in the repo. Every other chain here can be traced end to
end; this one cannot.

Most likely an earlier revision of the supertype script, run before the
subclass/supertype split. **Needs the actual script from Nicole.**

### 3. `Xenium_SCZ_R.rds` has no provenance, and no script creates it

`Final_figures/Figure_3.r:166,360` reads a Seurat conversion of the Xenium object
from `/scratch/nendresz/Xenium/Xenium_SCZ_R.rds`. Nothing in this repo produces
it, and it is not the canonical object.

This matters more than a usual missing input: `DATA_FLOW.md` records that a
superseded Xenium object carries a **transposed Lieber map** (cluster 9 ↔ 12,
i.e. MGE and CGE swapped) and a different cell count (373,514 vs 356,313
neuronal cortical cells). Anything derived from the wrong one is silently
inverted.

**Needs confirming** that `Xenium_SCZ_R.rds` derives from
`SCZ_Xenium/output/all_samples_annotated.h5ad` (md5
`763e6655fc55839f177d57aa98dc5453`, the 2026-04-01 state), and a one-line note
saying so — ideally the conversion script.

### 4. `cluster_order_and_colors.csv` is not committed, and it is load-bearing

Read by six scripts across `Compositional_analysis/`,
`Compositional_sensitivity_analysis/` and `Final_figures/`. It is not just a
palette: `2_Crumblr_analysis.r:44` uses its `class_label` column to **define
which cell types are neuronal**, and therefore which types enter the primary
composition analysis. `Non_neurons/1_crumblr_meta_nonNeurons:32` takes the
complement for S5.

It is a small CSV. Committing it under `snrnaseq/` would make the neuronal/
non-neuronal split inspectable instead of implicit.

### 5. A cell-count covariate is computed but never enters the DE design

`snRNAseq_DE/Subclass/2_DE.r:85,202` and the supertype equivalent compute

```r
meta_subset$log2_cells <- scale(log2(meta_subset[[type]]))
```

but the design matrices (`2_DE.r:104,221` and `Supertypes/2_DE.r:112,242`) are
`~ scale(Age) + Sex + Diagnosis [+ scale(PMI)]` — `log2_cells` is not in them.

A reader will reasonably assume DE was adjusted for cells per donor. It was not.
There is evidence a variant run *did* adjust for it: `log2cells_Sst` appears in
the subclass exclusion list at `Supertypes/3_meta_analysis.r:15`.

**Not changed here** — deleting the line would be wrong if the adjustment was
intended, and adding it to the design would change every DE result. Nicole
should confirm which run the paper reports, then either drop the dead line or
restore the covariate.

### 6. The committed donor metadata is person-level

`Compositional_analysis/7_cohorts_metadata_names.csv` holds donor ID, age, sex,
diagnosis and PMI for 469 donors. These are de-identified consortium IDs and the
repo is currently **private**, so nothing is exposed today. Worth an explicit
decision before the repo is made public for review.

---

## Cosmetic / naming

Batch these; nothing depends on them.

| # | Issue | Where |
|---|---|---|
| 8 | `Figure_1a_UMAP.r` is **Python**, not R (`from brisc import SingleCell`). Valid Python, wrong extension. Rename to `.py`. | `Final_figures/` |
| 9 | `1_crumblr_meta_nonNeurons` has **no file extension**; it is R. | `Compositional_analysis/Non_neurons/` |
| 10 | The label-transfer reference is loaded into variables named `counts_hodge` / `meta_hodge`, but it is the **SEA-AD** taxonomy (Gabitto 2024) — `refdata = ref$Supertype`, labels carry the `-SEAAD` suffix — not Hodge 2019, which is a different taxonomy of a different region. A reviewer cross-checking methods against code would read the wrong reference. Rename, or add a header comment. | 6 files under `Label_transfer/` and `Compositional_sensitivity_analysis/` |
| 11 | 24 of 34 scripts hard-code `/scratch/nendresz/…` or `/project/rrg-shreejoy/…`, and most open with `setwd("P1_SCZ_DE_fresh")` / `setwd("FINAL_FIGS")` / similar — directory names that do not exist in this repo. The absolute paths are expected (the analyses ran on the Alliance cluster); the `setwd()` calls are worse, because they make the scripts look like they belong to a tree that was not shipped. Dropping the `setwd()` lines would cost nothing. | throughout `snrnaseq/` |
| 12 | Her figure scripts write to `Figures/1a.png`, `Barchart_SCZ_meta_FINAL.png`, `Spatial3d.png` etc. on the cluster. Every other renderer in this repo writes its figure **into `manuscript/figures/` under its S-number**. Consequence: Figs 1 and 3 and Supp S1, S4, S5, S7 have no rendered output here, and no filename-to-panel mapping. | `Final_figures/` |
| 13 | Cell-type columns are selected by hard-coded position (`colnames(meta)[1:24]`, `[1:23]` for Multiome) rather than by name. Correct as written, but silently wrong if a cohort's subclass set changes. | `snRNAseq_DE/*/2_DE.r` |
| 13b | Figure 3's panel variables are not named after their panel letters — `p3b` becomes panel **c**, `p3c` becomes panel **b**, `p3h` becomes panel **e**. Only the `plot_grid` assembly at the foot of the script says so. Mapped in `Final_figures/README.md`. | `Final_figures/Figure_3.r` |
| 14 | One tracked symlink dangles in any fresh clone: `reserve/histology/coordinates` → `/Users/shreejoy/Github/sgACC_cell_depth_analysis/coordinates`. | `reserve/histology/` |

---

## Fixed on 2026-09-13

- **`snRNAseq_DE/Supertypes/3_meta_analysis.r` did not parse.** Two edit
  artifacts: a bare `-` on its own line after the meta-analysis loop, and a
  missing fourth closing parenthesis on the `sup_cohorts <- bind_rows(lapply(...`
  expression. Together they made the file un-runnable as committed, which means
  the committed copy is not the copy that produced the results — worth Nicole
  confirming the fixed version matches what she ran. Both fixed; the file now
  parses, and every other R script under `snrnaseq/` parses clean.

- **The *HCN1* credible-set variants were labelled on the wrong strand.**
  *HCN1* is transcribed on the minus strand (NM_021072.4), so its TSS is at the
  higher coordinate. `export_panels_abc.py` was labelling the three variants
  below the gene start `5'_upstream` when they are 3′ of the gene. Fixed
  strand-aware, and a `dist_to_hcn1_kb` column added; `panel_D_credible_set.csv`
  and `panel_D_meta.csv` regenerated and the manifest refreshed. No upstream
  source moved — every `source_md5` is unchanged.
