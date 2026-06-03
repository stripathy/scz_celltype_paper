# Reproducing the composite figure (Fig. 09)

Checklist to regenerate `results/figures/09_composite.{png,pdf}` — the
6.5 × 9.2 in publication composite — on another machine.

## TL;DR

`scripts/09_composite_figure.R` is self-contained (it does not `source()`
anything) and reads its **panel-K inputs from CSVs that are committed to this
repo**. So you do NOT need the raw Xenium data or the Python pipeline to
reproduce the figure — only three DE tables and R:

```bash
cd scz_pathway_enrichment
Rscript scripts/09_composite_figure.R     # -> results/09_composite.{png,pdf}
# compare to the committed reference: results/figures/09_composite.png
```

## 1. The code (travels with git)

Clone the repo and check out the pinned state:

```bash
git clone <remote-url> && cd scz_pathway_enrichment
git checkout composite-fig-v1        # the tagged commit this figure was built at
```

Everything code/text the figure needs is in git: the script, the panel-K
exemplar inputs (`results/tables/exemplar_*.csv`, 8 files), and all docs.

## 2. The three data files — NOT in git, transfer separately

These are git-ignored (large / external). They are the **only** things you must
hand over outside of git. Place each at the path shown:

| Place at this path | Size | What it is |
|---|---|---|
| `data/DE_genes_all_cells_scz.csv` | 34 MB | Meta-analytic snRNA-seq DE (7-cohort), one row per (subclass × gene). Feeds the butterfly (A), volcanoes (B,C), forest asterisks, and scatter x-axis. |
| `data/meta_results_cohorts_subclass.csv` | 313 MB | Per-cohort snRNA-seq DE (the 7 cohorts). Feeds the forest cohort rows + recomputed meta diamond. |
| `~/Github/SCZ_Xenium/output/de/de_results_subclass.csv` | 876 KB | Xenium spatial DE. Feeds the forest Xenium rows + scatter y-axis. Only this **one CSV** is needed (not the whole SCZ_Xenium repo); if you store it elsewhere, edit `INPUT_XENIUM` at the top of the script. |

Total to transfer ≈ 348 MB (the per-cohort file is the bulk). Column-level
provenance: `data/README.md` (first two) and
`notes/figures_crossplatform_validation.md` §1 (all three).

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
cd scz_pathway_enrichment
Rscript scripts/09_composite_figure.R     # writes results/09_composite.{png,pdf}
```

Compare `results/09_composite.png` against the committed reference
`results/figures/09_composite.png`. The figure legend — and a table verifying
every cited number against source — is in `notes/figure_composite_legend.md`.

## Exact vs functional reproduction

- **Functional** (identical panels, numbers, layout): the steps above suffice.
  All on-figure statistics (n = 166, r = 0.70, 74% concordant, the DE counts,
  the forest estimates) are computed from the data at run time, so they match as
  long as the three data files are identical.
- **Pixel-exact** also requires the same package versions (above) and the same
  device **font**. The script uses the default sans font; a different sans font
  changes text metrics, which can nudge the ggrepel labels. Note: the scatter
  (J) labels are placed deterministically (manual nudges) and are robust; the
  volcano (B,C) labels use ggrepel (fixed seed) and are the most font/version
  sensitive. For a guaranteed environment, capture a lockfile with
  `renv::init(); renv::snapshot()` or install the exact versions above.

## Optional — regenerate the panel-K exemplar cells from scratch

NOT needed to reproduce the figure (the CSVs are committed). Do this only to
re-extract exemplar cells (e.g. for different genes — edit `PAIRS` in the
script). Requires Python (`numpy pandas scipy anndata matplotlib`) and these
Xenium exports from the SCZ_Xenium repo:

- `output/h5ad/{Br6432,Br2039}_annotated.h5ad`
- `output/deploy/boundaries/{Br6432,Br2039}.json`
- `output/deploy/transcripts/{Br6432,Br2039}/{SST,RASGRF2}.json` + `gene_index.json`

Then:

```bash
python scripts/10_xenium_exemplar_cells.py   # rewrites results/tables/exemplar_*.csv
Rscript scripts/09_composite_figure.R
```

See `notes/figures_crossplatform_validation.md` §6 for details.
