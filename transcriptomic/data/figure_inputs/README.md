# Committed figure inputs — Fig. 09 composite

Real, self-contained snapshots of the CSVs that `scripts/09_composite_figure.R`
reads, committed here so the composite reproduces from a fresh clone with **no
external files**. This directory is a subdirectory of `data/`, so it escapes the
single-level `data/*.csv` git-ignore. Full provenance: `../../REPRODUCE.md` §2 and
`../../notes/figures_crossplatform_validation.md` §1.

| File | What it is | Canonical source |
|---|---|---|
| `DE_genes_all_cells_scz.csv` | Meta-analytic snRNA-seq DE (7-cohort), full per-(cell type × gene) table | `shared/snrnaseq_de/` |
| `meta_results_cohorts_subclass_forest.csv` | Per-cohort snRNA-seq DE — **SST + PVALB rows only** (the forest subset; full per-cohort table is 313 MB and stays external) | `shared/snrnaseq_de/meta_results_cohorts_subclass.csv` |
| `de_results_subclass.csv` | Xenium spatial DE | SCZ_Xenium `output/de/` |
| `crumblr_input_subclass_corr.csv` | Xenium per-donor subclass composition (panel-i inset) | SCZ_Xenium `output/crumblr/` |

These are **snapshots** — to refresh after an upstream DE rerun, copy the new
source table in (re-subsetting the cohorts to `genes ∈ {SST, PVALB}`) and
re-render. See `notes/figures_crossplatform_validation.md` §7 (update checklist).
