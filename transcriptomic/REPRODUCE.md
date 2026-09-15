# Reproducing Figure 2

Checklist to regenerate `manuscript/figures/main/Fig2_cross_platform_de.{png,pdf}`
— the 7.1 × 6.625 in publication composite — on another machine.

## TL;DR

`scripts/09_composite_figure.R` runs from a clean clone: it `source()`s only two
tracked helpers beside it (`scripts/_figure_inputs.R`, `scripts/_exemplar_panels.R`)
and every CSV it reads is committed to the repo (the DE
inputs under `data/figure_inputs/`, the exemplar + CP1K tables under
`results/tables/`). You do NOT need the raw Xenium data, the Python pipeline, or
any external file — just R:

```bash
cd scz_celltype_paper/transcriptomic
Rscript scripts/09_composite_figure.R
# -> ../manuscript/figures/main/Fig2_cross_platform_de.{png,pdf}
# The script writes the submission figure in place, so there is no copy to drift.
```

## 1. The code (travels with git)

Clone the monorepo and enter the transcriptomic component:

```bash
git clone <remote-url> && cd scz_celltype_paper/transcriptomic
# (the composite-figure code lives in this transcriptomic/ tree; no separate tag)
```

Everything the figure needs is in git: the script, the **DE input CSVs**
(`data/figure_inputs/`, see §2), the panel d/h exemplar inputs
(`results/tables/exemplar_*.csv`, 18 files + `exemplar_cells_meta.csv`), the
per-donor CP1K table (`marker_norm_expr*.csv`), and all docs.

## 2. Data inputs — committed under `data/figure_inputs/`

The figure's four DE / composition inputs are committed as **real files** in
`data/figure_inputs/` (a subdirectory, so it escapes the single-level
`data/*.csv` ignore). Nothing has to be transferred separately — a fresh clone
has everything.

| File (in `data/figure_inputs/`) | Size | What it is |
|---|---|---|
| `DE_genes_all_cells_scz.csv` | 34 MB | Meta-analytic snRNA-seq DE (7-cohort), one row per (subclass × gene). **Full table.** Feeds the butterfly + inset (i), volcanoes (a,e), forest asterisks, scatter x-axis. |
| `meta_results_cohorts_subclass_forest.csv` | 20 KB | Per-cohort snRNA-seq DE — **SST + PVALB rows only** (141 of 2.16 M), the subset the forests (b,f) read. The full 313 MB per-cohort table stays external. |
| `de_results_subclass.csv` | 877 KB | Xenium spatial DE (snapshot of the SCZ_Xenium pipeline output). Feeds the forest Xenium rows + scatter y-axis. |
| `crumblr_input_subclass_corr.csv` | 21 KB | Xenium per-donor subclass composition. Feeds the panel-i DE-vs-proportion inset. |

These are **snapshots**. Canonical sources: `shared/snrnaseq_de/` (the two snRNA
tables) and the SCZ_Xenium repo (`output/de/`, `output/crumblr/`). To refresh
after an upstream DE rerun, regenerate them (the cohorts subset is just the
`genes ∈ {SST, PVALB}` rows) — see `notes/figures_crossplatform_validation.md`
§1 (column provenance) and §7 (update checklist).

## 3. R environment (match for an exact match)

R 4.5.1 with:

| package | version | | package | version |
|---|---|---|---|---|
| readr | 2.2.0 | | ggrepel | 0.9.6 |
| dplyr | 1.2.0 | | metafor | 4.8.0 |
| tidyr | 1.3.1 | | scales | 1.4.0 |
| ggplot2 | 4.0.2 | | tibble | 3.3.1 |
| cowplot | 1.2.0 | | | |

`ggplot2` and `ggrepel` matter most — they govern layout and label placement.

## 4. Run + verify

```bash
cd scz_celltype_paper/transcriptomic
Rscript scripts/09_composite_figure.R
# writes ../manuscript/figures/main/Fig2_cross_platform_de.{png,pdf}
```

The renderer writes the submission figure in place, so a successful run replaces
`manuscript/figures/main/Fig2_cross_platform_de.{png,pdf}`; `git diff --stat` on
that path is the check that it changed as expected. The figure legend — and a table verifying
every cited number against source — is in `notes/figure_composite_legend.md`.

## Exact vs functional reproduction

- **Functional** (identical panels, numbers, layout): the steps above suffice.
  All on-figure statistics (n = 166, r = 0.73, 72% concordant, the DE counts,
  the forest estimates) are computed from the committed `data/figure_inputs/`
  CSVs at run time, so a fresh clone reproduces them exactly.
- **Pixel-exact** also requires the same package versions (above) and the same
  device **font**. The script uses the default sans font; a different sans font
  changes text metrics, which can nudge the ggrepel labels. Note: the scatter
  (j) and volcano (a,e) labels both use ggrepel (fixed seed = 7), so they are
  the most font/version sensitive. For a guaranteed environment, capture a lockfile with
  `renv::init(); renv::snapshot()` or install the exact versions above.

## Optional — regenerate the panel d/h exemplar cells from scratch

NOT needed to reproduce the figure (the CSVs are committed). Do this only to
re-extract exemplar cells (e.g. for different genes — edit `PAIRS` in the
script). Requires Python (`numpy pandas scipy anndata matplotlib`) and these
Xenium exports from the SCZ_Xenium repo. Drawing sections are per gene — SST &
PVALB from `Br6432`/`Br5973`, FGFR3 from `Br5400`/`Br5973` — but the group-median
targets are computed over **all 24 donor h5ads**, so the full set is needed:

- `output/h5ad/*_annotated.h5ad`  (all 24 donors — for the pooled group medians)
- `output/deploy/boundaries/{Br6432,Br5400,Br5973}.json` + `*_nucleus.json`
- `output/deploy/transcripts/{Br6432,Br5400,Br5973}/{SST,FGFR3,PVALB}.json` + per-section `gene_index.json`

Then:

```bash
python scripts/10_xenium_exemplar_cells.py   # rewrites results/tables/exemplar_*.csv
Rscript scripts/09_composite_figure.R
```

See `notes/figures_crossplatform_validation.md` §6 for details.
