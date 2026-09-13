# transcriptomic/ — Figure 2 and Supplementary Fig. S8

Cell-type-specific differential expression in schizophrenia, across the
seven-dataset snRNA-seq meta-analysis and the Xenium spatial dataset. This
component builds **Figure 2** (cross-platform DE) and **Supplementary Fig. S8**
(the Sst depletion-strata analysis), and nothing else.

The exploratory GSEA / gene-ontology / pathway work that once lived here was
removed on 2026-09-04; none of it entered the manuscript. It is recoverable
from the git tag `pre-prune-2026-09-04`.

## Figure 2

Ten panels, all assembled by one script.

| Panels | What | Data from |
|---|---|---|
| a, e | Volcano, meta log₂FC vs −log₁₀ *P*, coloured by direction × FDR tier | `data/figure_inputs/DE_genes_all_cells_scz.csv` |
| b, f | Forest: seven datasets, pooled DerSimonian–Laird diamond, Xenium triangle | `meta_results_cohorts_subclass_forest.csv` + Xenium DE |
| c, g | Per-donor CP1K expression, Control vs SCZ (12 vs 12) | `12_marker_norm_expr.R` |
| d, h | Representative Xenium cells at the pooled group-median transcript density | `10_xenium_exemplar_cells.py` |
| i | Up/down DE-gene counts per subclass, FDR < 0.10 with the FDR < 0.05 subset | meta DE |
| i inset | DE-gene count vs mean per-donor Xenium proportion | Xenium crumblr input |
| j | snRNA-seq meta vs Xenium log₂FC, the pairs at meta FDR < 0.10 | meta DE ∩ Xenium DE |

```bash
cd transcriptomic
python3 scripts/10_xenium_exemplar_cells.py   # panels d, h   (needs the Xenium h5ads)
Rscript scripts/12_marker_norm_expr.R         # panels c, g   (needs Xenium pseudobulk)
Rscript scripts/09_composite_figure.R         # the figure
# -> manuscript/figures/main/Fig2_cross_platform_de.{png,pdf}
```

`09_composite_figure.R` builds every panel internally and reads only committed
CSVs, so **the figure regenerates from a clone with R alone** — see
[`REPRODUCE.md`](REPRODUCE.md). Steps 10 and 12 are only needed when their
inputs change; their outputs are committed under `results/tables/`.

`00_refresh_figure_inputs.R` resyncs `data/figure_inputs/` from the canonical
upstream sources and rewrites `MANIFEST.tsv`. `_figure_inputs.R` checks that
manifest on every read and halts the render if a source has moved on, so the
figure cannot quietly go stale.

## Supplementary Fig. S8 — Sst depletion strata

A nine-step pipeline in [`scripts/fig5/`](scripts/fig5/README.md) asking how the
SCZ transcriptional state of Sst interneurons differs between the supertypes
depleted in Figure 3 and those that persist. It was built as Figure 5 and moved
to the supplement on 2026-09-01, so it is still laid out as a main figure.

```bash
# from the repo root, in order; ~50 min total
Rscript transcriptomic/scripts/fig5/01_strata.R
# ... through ...
Rscript transcriptomic/scripts/fig5/08_figure5.R    # -> manuscript/figures/supplementary/S08_sst_strata
Rscript transcriptomic/scripts/fig5/09_verify.R     # must PASS
```

`09_verify.R` asserts every number quoted in the draft against its source table.
Run it after any change to the pipeline.

Its supplements and sensitivity analyses are not in the paper and live in
[`../reserve/sst_strata_supp/`](../reserve/sst_strata_supp/).

## Layout

```
data/figure_inputs/   committed Fig 2 inputs + MANIFEST.tsv (staleness guard)
data/stratum_*_export/  manifests for the cluster exports behind S8 (the
                        parquet/h5ad payloads are gitignored)
scripts/              Fig 2 chain (00, 09, 10, 12, _figure_inputs)
scripts/fig5/         the S8 pipeline (01-09 + _common.R)
results/tables/       committed Fig 2 panel inputs (exemplars, CP1K)
results/sst_strata_gsea/  committed S8 statistics
notes/                figure legend + cross-platform validation write-ups
```

## Cross-component seams

- **In:** the snRNA-seq meta-analysis and per-dataset tables, via
  `shared/snrnaseq_de/` (Endresz et al., in prep).
- **In:** Xenium DE, crumblr input and supertype depth from `spatial/output/`.
- **Out:** nothing. Both figures are terminal.
