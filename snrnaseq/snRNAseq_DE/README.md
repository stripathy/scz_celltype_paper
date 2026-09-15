# snRNAseq_DE/ — differential expression in SCZ (Fig. 2 inputs, S4)

Stage 2b. Pseudobulk differential expression per cell type per dataset, pooled
across the seven datasets by per-gene meta-analysis. Run twice at two
granularities:

| Directory | Grouping | Feeds |
|---|---|---|
| `Subclass/` | 23–24 SEA-AD subclasses (`predicted.id` with the supertype suffix stripped) | **Figure 2**, Supplementary Fig. **S8**, and the Sst DE gene list used by `../Compositional_sensitivity_analysis/` |
| `Supertypes/` | the full SEA-AD supertypes | Supplementary Fig. **S4** |

Both follow the same three steps and differ only in the grouping variable.

```
1_Pseudobulk.r    -> per-dataset pseudobulk matrices + per-donor metadata
2_DE.r            -> DE_results_{cohort}_{celltype}.rds     (limma-voom)
3_meta_analysis.r -> meta_results_{celltype}.csv -> pooled tables
```

In `Subclass/` the pooling half of step 3 is a separate file,
`3a_meta_per_gene.r`, which runs before `3_meta_analysis.r`. It is a
**reconstruction**, not the original — see step 3 below.

## 1 — Pseudobulk

`Seurat::AggregateExpression(group.by = c(<celltype>, "Donor"))` on each dataset,
saved as one matrix per dataset. Alongside it, a per-donor metadata CSV carrying
the `Donor × celltype` cell counts plus `Age, Sex, Diagnosis, PMI`.

Subclass labels are derived from `predicted.id` by stripping the trailing
supertype index (`gsub("_[0-9].*$", "", …)`), with `Lamp5_Lhx6` collapsed to
`Lamp5Lhx6` so the underscore rule does not split it.

## 2 — DE per dataset per cell type

For each cell type in each dataset:

- **donors:** keep those with **≥ 500 total cells**,
- **genes:** keep those with ≥ 1 count in ≥ **80%** of samples,
- **normalisation:** TMM (`edgeR::calcNormFactors`),
- **model:** `voom` → `lmFit` → `eBayes` on
  `~ scale(Age) + Sex + Diagnosis + scale(PMI)` (Multiome drops PMI — it has
  none),
- **contrast:** `coef = "DiagnosisSchizophrenia"`, BH-adjusted within cell type.

Cell types where `Diagnosis` has fewer than two levels after filtering are
skipped.

> **DE is not adjusted for cells per donor**, by design. A `log2_cells` variable
> was computed here but never entered the design matrix; it was removed in
> `c106186`. Confirmed empirically on 2026-09-15: re-running the Sst model
> without it reproduces the committed per-dataset table **exactly** (McLean,
> MSSM 1 and Multiome; same gene sets, max |ΔlogFC| 4.5e-14), while adding it
> drops agreement to r = 0.68-0.81. Closed as issue 5 in
> [`../../KNOWN_ISSUES.md`](../../KNOWN_ISSUES.md).

> Cell-type columns are selected by hard-coded position (`colnames(meta)[1:24]`,
> and `[1:23]` for Multiome, which lacks one subclass) rather than by name.
> Issue 13.

## 3 — Meta-analysis across datasets

Per cell type, per gene: `metafor::rma(yi = logFC, sei = SE, method = "REML")`,
requiring the gene to be testable in **more than 4** of the 7 datasets. Standard
errors reconstructed as `|logFC / t|`. FDR (BH) applied across genes within cell
type.

`Supertypes/3_meta_analysis.r` runs this loop and writes, in order:

| Output | What |
|---|---|
| `meta_results_{celltype}.csv` | pooled result, one file per cell type |
| `DE_results_cohorts_supertype.csv` | all per-dataset rows, long |
| `DE_meta_results_supertype.csv` | all pooled rows, long |
| `DE_all_cohorts_meta_supertype.csv` | the two above stacked, datasets renamed to paper labels, `"Meta-analysis"` as a pseudo-dataset — this is what S4 reads |
| `DE_all_cohorts_meta_supertype_by_subclass_cohort.xlsx` | the same split one sheet per subclass × dataset (supplementary table) |

`Subclass/3_meta_analysis.r` *consumes* `meta_results_*.csv` and assembles:

| Output | What |
|---|---|
| `sig_genes_all_cells_scz.csv` | union of genes at FDR < 0.05 in any subclass |
| `sig_genes_all_cells_scz_w_estimates.RDS` | the same, with effect sizes, per subclass |
| `DE_genes_all_cells_scz.csv` | **all genes × all subclasses** — the primary DE export |

Its last line writes `DE_genes_all_cells_scz.csv` directly into
`transcriptomic/data/figure_inputs/`, which is the live seam into **Figure 2**.

`Subclass/3a_meta_per_gene.r` runs the same loop for the subclasses and writes
their `meta_results_{celltype}.csv`. ⚠️ **It is a reconstruction committed on
2026-09-13, not the script that produced the published numbers** — the original
was never in the repo and has not been recovered. It is
`Supertypes/3_meta_analysis.r` with the cell-type selection inverted; the file's
own header lists the evidence. It is now committed as `Subclass/3a_meta_per_gene.r`.

Re-running the per-dataset DE reproduces the committed table **exactly** wherever
the dataset's gene universe is recoverable: McLean, MSSM 1 and Multiome, which are
already in gene-symbol space, match gene for gene at max |ΔlogFC| 4.5e-14
(2026-09-15). Datasets needing an identifier remap do not, because dropping
unmappable genes moves the >= 80% expression filter and hence every library's TMM
factor — that is the whole of the earlier r = 0.9968 gap, and it is a property of
the rebuild, not of the pipeline. Issues 2 and 19 in
[`../../KNOWN_ISSUES.md`](../../KNOWN_ISSUES.md).
