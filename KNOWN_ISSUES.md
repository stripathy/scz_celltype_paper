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

**Identified on 2026-09-13 (second audit pass).** It is
`Supertypes/3_meta_analysis.r` with the cell-type selection inverted: line 16
keeps `setdiff(cell_types, subclasses)`, the subclass run kept the intersection.
Everything else matches:

- lines 43-44 emit exactly the committed file's columns — `cell_type, genes,
  estimate, se, pval, ci.lb, ci.ub, k, tau2, I2` plus `padj`;
- line 53 writes `Files/meta_results_<ct>.csv`, the very files
  `Subclass/3_meta_analysis.r:20` globs;
- `k` in the committed output takes values {5,6,7}, matching `cohorts` (7) and
  the `nrow(df)>4` threshold; `tau2`/`I2` are non-zero, matching `method="REML"`.

Note the bootstrapping order this implies: `Supertypes/3_meta_analysis.r:10`
*derives* its cell-type list by globbing `meta_results_*.csv`, so the subclass
run must have come first and written them. That is also why line 13-15 carries a
hard-coded subclass exclusion list, `log2cells_Sst` included.

An ancestor of the same code is readable at
`/project/rrg-shreejoy/nendresz/Meta_DE/3_Meta_analysis_uns.r` — same 22
subclasses and per-gene `rma()`, but 6 datasets, `method="FE"`, and `.rds`
output with different column names, so it is not the run the paper reports.

**Closed on 2026-09-13.** `snrnaseq/snRNAseq_DE/Subclass/3a_meta_per_gene.r` is
committed, and the chain now executes end to end. The recipe is confirmed by
issue 19 — it recovers 12,490 of 12,492 committed Sst genes at r = 0.9968 — so
this is treated as the step, not as a stand-in for it.

**One caveat kept in the file's header:** re-running it does not reproduce the
published table to the last decimal, because exact recovery also needs Nicole's
**per-dataset gene universes** (issue 19 explains why). Worth asking her for
those if anyone needs to reconcile numbers exactly.

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

**Done on 2026-09-13**: copied from `/project/rrg-shreejoy/nendresz/` to
`snrnaseq/cluster_order_and_colors.csv` (139 supertypes, md5
`43649970c522eb232eb2ba8e0a126551`, byte-identical to the source). The eight
scripts still read the absolute path, which resolves on this cluster; repointing
them is a cosmetic follow-up (see item 11).

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

### 15. Figure 2 no longer renders from a clone, and lost its staleness guard

Commit `038490d` ("Update composite figure script", 2026-09-09) changed
`transcriptomic/scripts/09_composite_figure.R` in three ways that break the
repo's central reproducibility claim:

1. it added `setwd("scz_celltype_paper/transcriptomic")`, which only works if the
   cwd happens to be the repo's *parent*;
2. it replaced `source("scripts/_figure_inputs.R")` + `fig_input(...)` with raw
   `data/figure_inputs/...` paths, **dropping the MANIFEST staleness guard** —
   the guard whose own header records that its absence is "how Fig. 2 came to
   show 76% concordant while the manuscript said 72%";
3. it introduced three reads from `/scratch/nendresz/`, which no one else can
   read.

Items 1 and 2 are **fixed** (see below). Item 3 is not fixable here: the three
files are not in the repo and `/scratch/nendresz` is mode `drwx--x--x`.

Figure 2 therefore does **not** render from a clean clone, contrary to
`README.md` and `## Reproducing`. The three files are small and belong under
`transcriptomic/data/figure_inputs/` with MANIFEST entries:

| purpose | path |
|---|---|
| per-cohort donor n (forest labels) | `/scratch/nendresz/P1_Compositional_analysis/plotdata.csv` |
| Xenium Sst_25 donor n | `/scratch/nendresz/Xenium/xen_Sst_proportions.csv` |
| mean subclass proportion (panel i inset) | `/scratch/nendresz/FINAL_FIGS/Paper/df_mean_subclass_prop.csv` |

**Resolved 2026-09-14** (commit `8364cc7`). All three are now snapshots under
`transcriptomic/data/figure_inputs/` with MANIFEST entries: the first two are
copied from the tables Nicole committed under `snrnaseq/Final_figures/Data/`;
the third was **replaced, not copied** — the inset now uses the mean per-donor
share of nuclei per subclass from the Fig. 3a counts (469 donors) rather than
Xenium cell proportions, which disagreed with the nuclei by up to 25-fold for
Lamp5_Lhx6 and L5 ET (see the inset comment in `09_composite_figure.R`). The
two snapshots Nicole rewrote on 9 Sep are re-sourced from her in-repo files
(`snrnaseq/snRNAseq_DE/Files/`, the per-dataset table via git-lfs), so the guard
runs with nothing bypassed and a refresh can no longer revert her numbers.
Figure 2 renders from a clean clone again, into `manuscript/figures/main/`.

### 16. Supplementary S8 reads the git-ignored seam, not the committed snapshot

`transcriptomic/scripts/fig5/_common.R:40` sets
`subclass = "shared/snrnaseq_de/DE_genes_all_cells_scz.csv"`, and line 38 sets
`crumblr` to another path under the same git-ignored directory. But
`shared/snrnaseq_de/README.md` states that "downstream figure code does not read
this directory directly — it reads committed snapshots". S8 does read it
directly, so S8 does not render from a clone either.

For `subclass` the committed snapshot exists and is provably the same file:
`transcriptomic/data/figure_inputs/MANIFEST.tsv` records
`DE_genes_all_cells_scz.csv` as captured *from* that exact seam path, md5
`719da3d7519cbc5a3aec1b8c42d438da`. **Fixed** by falling back to the snapshot
when the seam file is absent.

`crumblr` (`nicole_scz_snrnaseq_betas/final_results_crumblr_7_cohorts.csv`) had
no committed snapshot anywhere in the repo and no copy on this cluster.

**Resolved on 2026-09-13.** The file was *regenerated* by re-running the repo's
own `2_Crumblr_analysis.r` and `3_meta_analysis.r` (paths changed, model
untouched) on a count matrix rebuilt from cleaned per-dataset h5ads, and is now
**tracked in git** at the seam path with a `PROVENANCE.md` beside it — the one
exception to `shared/snrnaseq_de/` being git-ignored, taken because the file is
17 KB and three components need it. It validates against three committed
artefacts independently: per-dataset betas to 1.8e-15, the FE meta to 8.1e-16,
and it reproduces `gwas_vs_casecontrol_composition.csv` to 1.0e-15 when fed to
`build_composition_table.py`. `genetics/`'s composition scripts, S6's
published-FE comparison and S8's `P$crumblr` all resolve from a clone now.

### 17. Nicole's crumblr scripts assume an NA-marked count matrix

`2_Crumblr_analysis.r:49` dropped structurally-empty cell types with
`colSums(!is.na(counts)) > 0`. Her committed matrix marks a supertype that was
never called in a cohort as `NA`, so this works. A matrix pivoted from per-cell
tables — as Shreejoy's clean h5ad exports produce — marks the same cells `0`,
where the test is a no-op and all-zero categories enter crumblr's CLR.

**Fixed** in all three crumblr scripts by also dropping all-zero columns.
Verified a no-op on her own matrix (it has no all-zero column in any dataset), so
no published number moves.

Still open, and **not** changed: `2_Crumblr_analysis.r:78` hard-codes
`coef = "DiagnosisSchizophrenia"`. Any matrix using the repo's own `Control`/`SCZ`
vocabulary (which `00_prepare_counts.R:37` and the clean h5ads both use) makes
that coefficient name not exist, and the script fails. Worth deriving the
coefficient name from the factor levels instead.

### 18. Figure 2's Xenium composition snapshot predates the object the paper uses, and drops Br2039

> **2026-09-14:** Figure 2 no longer reads this snapshot — the panel-i inset moved
> to snRNA-seq nuclei proportions (see #15). The discrepancy below still stands
> for the snapshot itself and for anything else that uses it, but it no longer
> affects any number in Figure 2.

`transcriptomic/data/figure_inputs/crumblr_input_subclass_corr.csv` (panel i
inset) holds **23 donors / 737,750 cells**. Every other Xenium-derived file in
the repo holds **24**, Br2039 included:

| file | donors | Br2039 |
|---|---|---|
| `spatial/output/crumblr/crumblr_input_{subclass,supertype}_{neuronal,nonneuronal}.csv` | 24 | yes |
| `transcriptomic/data/figure_inputs/pseudobulk_subclass.csv` | 24 | yes |
| `transcriptomic/results/tables/marker_norm_expr.csv` | 24 | yes |
| `transcriptomic/data/figure_inputs/crumblr_input_subclass_corr.csv` | **23** | **no** |

**No committed code drops Br2039.** `spatial/code/modules/constants.py:31` is
explicit and has been since the monorepo was scaffolded (`3c7b45a`, 2026-06-02):

```python
EXCLUDE_SAMPLES = set()  # No samples excluded; Br2039 (WM-heavy) included — improves snRNAseq concordance
```

The snapshot is simply **older than that decision**. Its MANIFEST row names the
source as `/Users/shreejoy/Github/SCZ_Xenium/output/crumblr/crumblr_input_subclass_corr.csv`
with `source_mtime` **2026-03-12** — three weeks before the canonical 2026-04-01
Xenium object, and from a laptop. Because that source is absent on any other
machine, the staleness guard skips its check by design and trusts the snapshot,
so nothing flags it.

Verified against the canonical object (via the clean export
`~/scratch/dataset_compliation/clean/Xenium.h5ad`), 2026-09-13:

- The canonical recipe is `corr_qc_pass & spatial_domain == "Cortical"`, all 24
  donors, split by class. It reproduces **exactly**: neuronal 356,313 and
  non-neuronal 385,790 cells, 419/419 rows identical to
  `crumblr_input_subclass_neuronal.csv`. Those are the numbers
  `spatial/output/crumblr/README.md` and `DATA_FLOW.md` report.
- Br2039 contributes **10,595** cells to that set (1,153 GABAergic /
  2,531 Glutamatergic / 6,911 non-neuronal).
- Canonical minus Br2039 is 731,508, still **6,242 short** of the snapshot's
  737,750 — so the snapshot is not the current object with one donor removed. It
  is a different, earlier annotation.

**Consequence:** Figure 2 panel i's inset is drawn from a pre-canonical Xenium
object with one donor missing, while every other Xenium number in the paper comes
from the 24-donor April object. **Needs a refreshed snapshot** taken from
`spatial/output/crumblr/` (which is committed, current, and reproduces exactly),
plus a MANIFEST row pointing at a path that exists outside the laptop.

Two committed docs also disagree and should be reconciled:
`spatial/methods_writeup.md:52` says "Br2039 is retained in all analyses
(65% white matter)", while `spatial/all_samples_annotated_guide.md:312-322` says
"**Exclude Br2039 from cortical compositional analyses**" (41.8% WM) and gives a
recipe that does so. The code follows `methods_writeup`. The WM fractions differ
because they count different columns; by `spatial_domain` Br2039 is **70.1% WM**
against a cohort median of 19.7%.

### 19. The DE chain reproduces in substance but not bit-for-bit

Related to issue 2. With the missing producer reconstructed from its recipe, the
subclass DE + per-gene meta was re-run for **Sst** from independently rebuilt
pseudobulks, following `1_Pseudobulk.r` and `2_DE.r` exactly (Age < 70,
>= 500 cells per donor, >= 1 count in >= 80% of samples, TMM, voom,
`~ scale(Age) + Sex + Diagnosis [+ scale(PMI)]`, then `rma(method = "REML")`).

Result: 12,490 of the committed 12,492 Sst genes recovered, effect sizes
**r = 0.9968** and 98.31% sign-concordant, 161 of 165 FDR < 0.05 genes shared,
and *SST* itself at beta -0.459 / padj 0.048 against the committed
-0.458 / 0.049. All 15 top committed hits reproduce.

The residual is a systematic beta shift (median -0.0025, IQR entirely negative)
that is **not** gene-symbol collapsing (3 of 12,490 genes). It is consistent with
**TMM normalisation drift** from a different per-dataset gene universe — the
rebuild maps to gene symbols and drops what will not map (Batiuk: 17,425 of
60,617), which moves the >=80% gene filter and hence every library's norm factor.

**Consequence for issue 2:** the recipe is confirmed correct, so recovering
Nicole's script matters less than recovering her **per-dataset gene universes**.
Exact reproduction of the published table needs the gene sets her objects
carried, not just her code. Worth asking her for those alongside the script.

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
| 13 | Cell-type columns are selected by hard-coded position (`colnames(meta)[1:24]`, `[1:23]` for Multiome) rather than by name. Correct as written, but silently wrong if a dataset's subclass set changes. | `snRNAseq_DE/*/2_DE.r` |
| 13b | Figure 3's panel variables are not named after their panel letters — `p3b` becomes panel **c**, `p3c` becomes panel **b**, `p3h` becomes panel **e**. Only the `plot_grid` assembly at the foot of the script says so. Mapped in `Final_figures/README.md`. | `Final_figures/Figure_3.r` |
| 14 | One tracked symlink dangles in any fresh clone: `reserve/histology/coordinates` → `/Users/shreejoy/Github/sgACC_cell_depth_analysis/coordinates`. | `reserve/histology/` |

---

## Fixed on 2026-09-14

- **Figure 2 renders from a clean clone again, with the staleness guard on**
  (#15). The three cluster-only inputs are committed snapshots; the two
  snapshots Nicole rewrote on 9 Sep are re-sourced from her in-repo pipeline
  files; the panel-i inset abundance axis is snRNA-seq nuclei, not Xenium cells;
  Lamp5_Lhx6 is classed Inhibitory. Canonical figure re-rendered; the duplicate
  under `results/Figures/` removed. Numbers that moved with her rerun: 343 (Sst)
  and 404 (Pvalb) DE genes at FDR < 0.10; panel j ρ = 0.74, 74% concordant.
- **Figure 2 panels d, h** redrawn for G#123 (counts worded, outlines named, 5 µm
  label legible and PDF-safe) — commits `da96f6d`, `479815a`.

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
