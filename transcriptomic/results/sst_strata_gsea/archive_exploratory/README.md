# Superseded outputs (not used by Figure 5)

- `fig1..fig6*.png`, `gsea_stratum_comparison.csv`, `module_genes_xenium_concordance.csv`
  Exploratory figures from the original single-script exploration (archived as
  `transcriptomic/scripts/archive/fig5_superseded/17_sst_strata_gsea.R`).
- `figS_sst_strata_mockup.(png|pdf)` — the figure's filename before it became a
  main figure; current output is `figure5_sst_strata.(png|pdf)`.
- `*_TOPLEVEL_COPY.csv` — copies of the pseudobulk GSEA/signature tables that
  once sat at the top level so the figure script could read them unchanged. The
  canonical versions live in `pseudobulk/`. Keeping a third copy is what let an
  earlier sensitivity analysis silently compare the pseudobulk result against
  itself; do not restore them.

`../archive_ivw/` is different: it holds the genuine inverse-variance-weighted
results and is the correct comparator for the framework-robustness supplement.
