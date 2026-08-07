# transcriptomic/ — Figure 2 and its supplements

Cell-type-specific differential expression in schizophrenia, across a 7-cohort
snRNAseq meta-analysis and the Xenium spatial dataset. This component builds
**Figure 2** and **Supplementary Fig. S7**, and nothing else.

Everything exploratory — the GSEA / gene-ontology / pathway analyses and the
narrative figures built on them — is parked in
[`archive/`](archive/README.md). None of it entered the manuscript.

## Figure 2

| Panels | What | Built by |
|---|---|---|
| a, e | Volcano, meta log₂FC vs −log₁₀ P, coloured by direction × FDR tier | `01_volcano_per_celltype.R` |
| b, f | Forest: 7 cohorts, pooled DerSimonian–Laird diamond, Xenium triangle | `07_forest_plots.R` |
| c, g | Per-donor CP1K expression, Control vs SCZ (12 vs 12) | `12_marker_norm_expr.R` |
| d, h | Representative Xenium cells at the pooled group-median transcript density | `10_xenium_exemplar_cells.py` |
| i | Up/down DE-gene counts per subclass, FDR < 0.10 with the FDR < 0.05 subset | `09_composite_figure.R` |
| i inset | DE-gene count vs mean per-donor Xenium proportion | `16_de_vs_proportion.R` |
| j | snRNAseq meta vs Xenium log₂FC, the 166 pairs at meta FDR < 0.10 | `08_meta_vs_xenium_scatter.R` |
| — | assembles all ten panels | `09_composite_figure.R` |

`09_composite_figure.R` rebuilds every panel internally, so it alone regenerates
the figure. The numbered scripts above are standalone, fully-labelled,
large-format renderers of the same panels — useful when a panel needs to be
inspected or redrawn on its own.

## Supplementary Fig. S7

Per-cell marker counts within the Xenium Sst and Pvalb populations, modelled with
negative-binomial mixed models under four normalisations. Cited in the Methods
for the four ratios it reports (SST 0.70× raw / 0.76× library-normalised;
PVALB 0.86× / 0.87×).

```
11_grain_density.py -> results/tables/percell_grain_density.csv
                    -> 13_supp_percell_metrics.R
                    -> results/tables/S_percell_stats.csv
                       results/figures/S_percell_metrics.{png,pdf}
```

## Layout

```
scripts/     00_refresh_figure_inputs.R   resync data/figure_inputs/ from canonical sources
             _figure_inputs.R             fig_input() — the staleness guard every figure reads through
             01, 07, 08, 09, 10, 12, 16   Figure 2 (see the table above)
             11, 13                       Supplementary Fig. S7
data/
  figure_inputs/                          committed input snapshots + MANIFEST.tsv
notes/
  figure_composite_legend.md              final legend wording + per-value provenance
  figures_crossplatform_validation.md     data provenance and update checklist
results/
  figures/, tables/                       curated snapshots
archive/                                  parked exploratory work — see its README
```

## Running it

```bash
cd transcriptomic
python3 scripts/10_xenium_exemplar_cells.py   # panel d/h exemplar tables — run BEFORE 09
Rscript  scripts/09_composite_figure.R        # the figure
```

The standalone panel renderers (`01`, `07`, `08`, `16`) and the S7 chain
(`11` then `13`) can be run in any order.

## Input snapshots and the staleness guard

Every figure script reads its inputs through `fig_input()`, which **refuses to
run** when a committed snapshot no longer matches the source it came from. That
guard exists because of a real incident: a figure once showed 76% cross-platform
sign concordance while the manuscript said 72%, because the snapshot and the
source had silently diverged.

After any upstream DE, crumblr or meta-analysis rerun:

```bash
Rscript scripts/00_refresh_figure_inputs.R    # resync + rewrite MANIFEST.tsv
python3 ../shared/verify_provenance.py        # confirm every set is consistent
```

`data/figure_inputs/DE_genes_all_cells_scz.csv` is Supplementary Table **T2**.

## Dependencies

R: `readr`, `dplyr`, `tidyr`, `purrr`, `stringr`, `ggplot2`, `cowplot`, `ggrepel`.
Python: `pandas`, `numpy`, `anndata`/`scanpy` (for the Xenium exemplar and
grain-density steps, which read the annotated h5ads).
