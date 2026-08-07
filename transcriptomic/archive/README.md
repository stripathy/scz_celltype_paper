# transcriptomic/archive/ — parked exploratory work

**Archived 2026-08-06. Nothing here is deleted, and nothing here is in the paper.**

`transcriptomic/` now holds only what is needed to build **Figure 2** and its
supplements. Everything that was exploratory moved here, keeping its original
relative layout so the scripts still resolve their own inputs — run any of them
with `transcriptomic/archive/` as the working directory and the paths work as
they did.

## Why these were parked

The gene-set enrichment, gene-ontology and pathway analyses never entered the
manuscript. Checked against the authoritative Google Doc, **"GSEA", "gene set
enrichment", "gene ontology" and "leading edge" return zero hits** — the only
"pathway" match is "causal pathway" in an unrelated sentence. Two narrative
figures built on top of them (the OxPhos and cholesterol stories) likewise have
no presence in the paper.

This is a scope decision, not a judgement about the analyses. They ran, they
produced real output, and that output is preserved here in full.

## What is here

```
scripts/     02_gsea_pipeline.R              GSEA across 23 cell types -> gsea_cache.rds
             03_pathway_summary_heatmap.R    pathway summary heatmap (reads the ORA snapshot below)
             04_gwas_overlap_and_leading_edge.R
             05_multi_testing_correction.R   Stouffer / Fisher combination across cell types
             06_explainer_leading_edge.R     leading-edge explainer panel
             15_coupling_prototype_xenium.py PROTOTYPE, self-described; DE-vs-composition coupling
exploratory/ archive_ora_*.R                 early over-representation analyses
             story_oxphos_inhibitory.R       5-panel OxPhos narrative figure
             story_cholesterol_glia_exc.R    5-panel cholesterol narrative figure
             mechanism_oxphos_tests.R
             de_butterfly_bar.py             superseded by 09_composite_figure.R's build_butterfly()
R/           shared.R                        dead: nothing sourced it; 09 re-declares its palettes inline
notes/       findings.md, literature_context.md, literature_per_celltype.md,
             per_celltype_top_hits.json, plan_de_composition_coupling.md,
             _agent_outputs/, _agent_prompts/
results/     tables/  gsea_*, early_ora_*, story_*_gsea, 00_de_*, coupling_prototype_*,
                      08_meta_vs_xenium_pairs*, dotfig_candidates.csv
             figures/ gsea_*, early_ora_*, story_*, 06_explainer_*, 07_forest_panel*,
                      00_de_butterfly*, coupling_prototype_xenium.png
             exploratory/archive_ora_sst/    the ORA snapshot script 03 reads
```

Three scripts here (`story_*`, `mechanism_oxphos_tests.R`) hard-code an
out-of-repo absolute path, `~/Github/scz_cell_type_enrichment/data/gwas/…`, so
they were already not runnable from a clean clone.

## Two things worth knowing

**`08_meta_vs_xenium_pairs.csv` is archived because it is stale, not exploratory.**
It is a June 2 render giving 74.1% sign concordance and r = 0.70, against the
canonical 72.3% / r = 0.73. The live version is written to `results/` by
`scripts/08_meta_vs_xenium_scatter.R` on each run. Keeping the stale copy in the
live tree was a genuine trap — the numbers are close enough that a substitution
would not be noticed.

**`dotfig_candidates.csv` is untracked and has no producer**, so archiving rather
than deleting it is the only way to keep it at all.

## Restoring something

```bash
cd transcriptomic
git mv archive/scripts/02_gsea_pipeline.R scripts/
```

and restore its inputs alongside. Everything here is in git history from this
commit forward, so `git log --follow` on any path will find it.
