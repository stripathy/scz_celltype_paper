# Committed figure inputs — Figure 2

Real, self-contained snapshots of every CSV the Figure 2 scripts read, committed
here so the figure reproduces from a fresh clone with **no external files**. This
directory is a subdirectory of `data/`, so it escapes the single-level
`data/*.csv` git-ignore. Full provenance: `../../REPRODUCE.md` §2 and
`../../notes/figures_crossplatform_validation.md` §1.

## Every Figure 2 script reads from here

`09_composite_figure.R` and `12_marker_norm_expr.R` resolve their inputs through
`fig_input()` in `scripts/_figure_inputs.R`. No live paths to `SCZ_Xenium/output/`
or `shared/` remain in the figure code — that split is what previously let the
composite render stale numbers while the standalone panel scripts read fresh
ones.

## Refreshing after an upstream rerun

```bash
Rscript scripts/00_refresh_figure_inputs.R
```

That resyncs every file below and rewrites `MANIFEST.tsv`, which records each
source path, its checksum, and when the snapshot was taken. On every read,
`fig_input()` re-checksums the canonical source and **stops the script** if it no
longer matches the manifest, naming the file and the refresh command. When the
canonical source is absent (fresh clone, collaborator machine) the check is
skipped and the snapshot is used as-is.

After refreshing, re-run `12_marker_norm_expr.R` and then `09_composite_figure.R`.

## Contents

| File | What it is | Canonical source |
|---|---|---|
| `DE_genes_all_cells_scz.csv` | Meta-analytic snRNA-seq DE (7-cohort), full gene × subclass table | `shared/snrnaseq_de/` |
| `meta_results_cohorts_subclass_forest.csv` | Per-cohort snRNA-seq DE — **SST + PVALB rows only**; the full table is ~2.2M rows and stays external | `shared/snrnaseq_de/meta_results_cohorts_subclass.csv` |
| `de_results_subclass.csv` | Xenium pseudobulk DE, subclass level (panels b, f, j) | SCZ_Xenium `output/de/` |
| `crumblr_input_subclass_corr.csv` | Xenium per-donor subclass composition (panel i inset) | SCZ_Xenium `output/crumblr/` |
| `pseudobulk_subclass.csv` | Xenium pseudobulk counts (normalised marker expression, panels c, g) | SCZ_Xenium `output/de/` |
| `pseudobulk_subclass_samples.csv` | Xenium pseudobulk sample metadata (donor, diagnosis, sex, age, PMI) | SCZ_Xenium `output/de/` |
| `MANIFEST.tsv` | Checksums and provenance for the six files above — written by the refresh script, read by the guard | generated |

The subclass DE snapshot is the 2026-07-31 rerun: April-1 Xenium object,
`~ diagnosis + sex + age + PMI`, genes detected in ≥ 80% of retained donors.
