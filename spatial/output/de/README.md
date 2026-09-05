# Xenium differential expression

## The canonical tables

Pseudobulk edgeR DE on the **2026-04-01** Xenium object (see `DATA_FLOW.md` for the
object's identity and cell-number accounting).

| File | Rows | Cell types | Donors | FDR < 0.10 | FDR < 0.05 |
|---|---|---|---|---|---|
| `de_results_subclass.csv` | 6,792 | 23 | 24 | 208 | 29 |
| `de_results_supertype.csv` | 33,562 | 127 | 24 | 450 | 188 |

**Model** (`code/analysis/run_de.R`): per cell type,
`DGEList → gene filter → calcNormFactors(TMM) → estimateDisp(robust) →
glmQLFit(robust) → glmQLFTest` on the SCZ coefficient, design
`~ diagnosis + sex + age + PMI` (age and PMI centred), keeping genes detected in
**≥ 80% of the donors retained for that cell type**.

Covariates and gene filter follow the snRNAseq meta-analysis, so the two platforms
are described by one sentence, and they match the crumblr composition model, which
already covaried PMI. Two snRNAseq conventions are deliberately **not** adopted,
because Xenium is a single 24-section cohort rather than seven pooled datasets:

- the ≥ 500 cells/donor/cell-type rule would drop 8 of 23 subclasses, **including
  L6b** (1/24 donors reach 500); `MIN_CELLS = 10` is used instead.
- BH across all gene × cell-type pairs is unattainable on a 300-gene panel
  (min P 3.6 × 10⁻⁵ × 6,640 pairs = 0.24), so **FDR is corrected within cell type**.

That last point matters when reading these tables next to the snRNAseq results: the
panel-wide FDR here is computed over far fewer genes and is **not** directly
comparable to the genome-wide snRNAseq FDR.

## Provenance, and what is superseded

The subclass table is also snapshotted at
`transcriptomic/data/figure_inputs/de_results_subclass.csv` for the Figure 2 scripts,
guarded by `MANIFEST.tsv` + `fig_input()`. Check both against their sources with:

```bash
python3 shared/verify_provenance.py
```

⚠️ **Xenium DE values quoted anywhere older than 2026-08-06 predate these tables.** They came
from a run whose object state no longer exists on disk and is not reproducible —
`de_results_subclass.archive_2026-06-04_preApril.csv` in the processing repo. Every
one of the seven named gene × cell-type results shifts slightly (e.g. SST in Sst:
−0.318, P 0.0522 then vs **−0.328, P 0.0524** now), sign-concordance across the 166
testable pairs moves 126/166 → **120/166**, and DE burden rises 94 → **208** hits at
FDR < 0.10. No conclusion reverses; the digits need updating. Reconcile against these tables,
not against the archived one. Figure 2 was re-rendered onto them on 2026-09-04.

Other archived DE tables in the processing repo (`archive_2026-02-24`,
`archive_2026-02-22`, `_hybrid`) are earlier 23-donor no-PMI runs. They are kept as
provenance records only.

## Regenerating

```bash
cd spatial
python3 code/analysis/build_de_input.py        # pseudobulk (both levels) + MANIFEST.tsv
Rscript  code/analysis/run_de.R subclass       # -> de_results_subclass_v2.csv
Rscript  code/analysis/run_de.R supertype      # -> de_results_supertype_v2.csv
```

`run_de.R` deliberately writes `_v2.csv` so the new fit can be diffed against the
committed table before replacing it. **Swap it in manually once checked** — that
review step is the guard, but note the swap itself is what lost the provenance of
the pre-April table, so record what you swapped and when.

The pseudobulk inputs these read are large (`pseudobulk_supertype.csv` is 21 MB) and
are not committed; regenerate them with the first command, which needs the annotated
h5ads.
