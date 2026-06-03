# Data

`DE_genes_all_cells_scz.csv` is the symlinked input — a meta-analytic DE table from a 7-cohort SCZ snRNA-seq integration, one row per (SEA-AD subclass × gene). Columns:

- `cell_type` — SEA-AD subclass (23 categories: 9 excitatory, 8 inhibitory, 6 glial)
- `genes` — HGNC symbol
- `estimate` — meta-analytic SCZ-vs-control effect (log fold change-like)
- `se` — standard error of estimate
- `pval`, `padj` — nominal and BH-adjusted p
- `ci.lb`, `ci.ub` — 95% CI on estimate
- `k` — number of contributing studies (max 7)
- `tau2`, `I2` — between-study variance and heterogeneity index

The file (~35 MB, 222K rows) is symlinked from `/Users/shreejoy/Downloads/DE_genes_all_cells_scz.csv` and is excluded from git via `.gitignore`. Replace the symlink if the source location changes.

External data referenced by scripts (not redistributed here):
- `~/Github/scz_cell_type_enrichment/data/gwas/scz_gwas_gene_set_no_mhc.csv` — 484 PGC3 SCZ GWAS genes (Trubetskoy 2022), MHC excluded.
- MSigDB collections via the `msigdbr` R package (cached locally on first use).
