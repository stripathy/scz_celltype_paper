# SCZ Pathway Enrichment

Pathway-level interpretation of differential-expression results from a 7-cohort meta-analysis of schizophrenia (SCZ) postmortem snRNA-seq, across 23 SEA-AD subclasses.

The input is a single meta-analytic DE table (`data/DE_genes_all_cells_scz.csv`) with one row per (cell type × gene) and meta-analytic effect, SE, p-value, and BH-adjusted p-value. The pipeline tests pathway-level enrichment against multiple gene-set databases, then layers in genetic (PGC3 SCZ GWAS) and cell-type-specific interpretations.

## Headline findings

1. **Synaptic / vesicle machinery is down-regulated across most cell types** — pan-cell-type signal, robust to all multi-testing corrections (Stouffer Z = -15 for trans-synaptic signaling regulation).
2. **OxPhos / mitochondrial collapse is selective to inhibitory neurons** — Sst (NES = -2.14, padj = 3e-5), Vip, Lamp5, Pvalb. Gene-level dissection points to **mitochondrial assembly + import + quality control failure** (TOMM40 down, AFG3L2 down, UQCC2 down, PINK1 down) rather than transcriptional biogenesis suppression (PGC-1α axis only mildly affected).
3. **Cholesterol biosynthesis collapses in excitatory neurons + macroglia** — coordinated suppression of the entire mevalonate pathway (HMGCR, HMGCS1, FDPS, ... → SREBF2) in L5/L6 IT excitatory + Astro + Oligo; OPCs go the opposite direction (NES > 0, potential compensation).
4. **GWAS convergence**: PGC3 SCZ GWAS genes are broadly down-shifted in DE z-scores, strongest in Oligo (***), L2_3 IT (**), OPC (**), Vip (**), Astro (*). **SREBF2** and **GOT2** are GWAS hits that double as leading-edge drivers of the cholesterol and OxPhos signatures respectively.
5. **Hypothesis dismissals worth noting**: glycolytic compensation (ρ = +0.44 not -1), activity-driven OxPhos (ρ(IEG, OxPhos) = -0.15 n.s.), and pure PGC-1α suppression all fail as explanations of the OxPhos signal.

See [`notes/findings.md`](notes/findings.md) for the full take-stock with caveats and open questions.

## Directory layout

```
.
├── data/                                 # Input DE table (symlinked, not committed)
│   └── DE_genes_all_cells_scz.csv      → /Users/shreejoy/Downloads/...
├── R/shared.R                            # Palettes, cell-class groupings, helpers
├── scripts/                              # Formalized, numbered pipeline
│   ├── 01_volcano_per_celltype.R         # Per-cell-type volcano (Sst exemplar)
│   ├── 02_gsea_pipeline.R                # Core GSEA across 23 cell types
│   ├── 03_pathway_summary_heatmap.R      # Curated theme heatmap + ORA-vs-GSEA
│   ├── 04_gwas_overlap_and_leading_edge.R # PGC3 GSEA + leading-edge × GWAS
│   ├── 05_multi_testing_correction.R     # Cross-cell-type Fisher / Stouffer / global BH
│   ├── 06_explainer_leading_edge.R       # GSEA running-ES walkthrough (Sst OXPHOS)
│   ├── 07_forest_plots.R                 # Per-gene cohort + meta + Xenium forest plots
│   ├── 08_meta_vs_xenium_scatter.R       # snRNA-seq vs Xenium concordance scatter
│   ├── 09_composite_figure.R             # Publication composite (butterfly+volcano+forest+scatter+exemplars)
│   ├── 10_xenium_exemplar_cells.py       # Xenium exemplar-cell extraction (panel l inputs)
│   ├── 11_grain_density.py               # Per-cell grain density, 24 donors → exemplar selection + suppl.
│   ├── 12_marker_norm_expr.R             # Per-donor CP1K + edgeR p (panel k; SST/FGFR3 + PVALB suppl.)
│   ├── 13_supp_percell_metrics.R         # Suppl.: per-cell SST/PVALB across normalisations
│   └── 14_supp_pvalb.R                   # Suppl.: PVALB forest + CP1K + exemplar cells
├── exploratory/                          # Deep-dives and archived early work
│   ├── story_oxphos_inhibitory.R         # 5-panel OxPhos figure
│   ├── story_cholesterol_glia_exc.R      # 5-panel cholesterol figure
│   ├── mechanism_oxphos_tests.R          # PGC-1α / Warburg / IEG hypothesis tests
│   ├── archive_ora_sst.R                 # First-pass ORA (superseded by GSEA)
│   ├── archive_ora_sst_figure.R
│   ├── archive_ora_sst_nominal.R
│   ├── archive_ora_all_celltypes.R       # Cross-CT ORA (superseded)
│   └── de_butterfly_bar.py               # Initial Python bar chart
├── results/
│   ├── gsea_cache.rds                    # GSEA results cache (regeneratable)
│   ├── figures/                          # Curated PNG/PDF snapshots from this session
│   ├── tables/                           # CSVs (full GSEA results, leading-edge, etc.)
│   └── exploratory/                      # Per-script subdirs for re-run outputs
└── notes/
    ├── findings.md                       # Full annotated findings
    ├── literature_context.md             # OxPhos / cholesterol novelty vs literature
    ├── literature_per_celltype.md        # Per-cell-type DE replication vs literature
    ├── figures_crossplatform_validation.md  # Docs + data provenance for figure scripts 07–14
    └── figure_composite_legend.md        # Nature Neuroscience legend + per-value provenance (fig 09)
```

> **Figures & cross-platform validation (scripts 07–14)**: the standalone forest
> plots and concordance scatter, the publication composite (script 09), and the
> Xenium exemplar-cell extraction (script 10) are documented with full data
> provenance and an update checklist in
> [`notes/figures_crossplatform_validation.md`](notes/figures_crossplatform_validation.md);
> the composite's figure legend (final wording + per-value provenance) is in
> [`notes/figure_composite_legend.md`](notes/figure_composite_legend.md).

## How to run

All scripts assume `pwd` is the repo root. Order:

```bash
cd ~/Github/scz_celltype_paper/transcriptomic
Rscript scripts/02_gsea_pipeline.R          # builds results/gsea_cache.rds (~10 min)
Rscript scripts/03_pathway_summary_heatmap.R # uses cache
Rscript scripts/04_gwas_overlap_and_leading_edge.R
Rscript scripts/05_multi_testing_correction.R
Rscript scripts/06_explainer_leading_edge.R  # GSEA explainer
Rscript scripts/01_volcano_per_celltype.R    # Sst volcano (edit CELL_TYPE constant to change)
```

Cross-platform validation figures + the publication composite (see
`notes/figures_crossplatform_validation.md`):

```bash
Rscript scripts/07_forest_plots.R            # standalone forest composite
Rscript scripts/08_meta_vs_xenium_scatter.R  # standalone concordance scatter
python  scripts/10_xenium_exemplar_cells.py  # panel-l exemplar tables (run BEFORE 09)
Rscript scripts/09_composite_figure.R        # publication composite (7.1 x 6.7 in)
```

Exploratory deep-dives (reuse the cached GSEA):

```bash
Rscript exploratory/story_oxphos_inhibitory.R
Rscript exploratory/story_cholesterol_glia_exc.R
Rscript exploratory/mechanism_oxphos_tests.R
```

## Dependencies

R: `readr`, `dplyr`, `tidyr`, `purrr`, `stringr`, `ggplot2`, `cowplot`, `ggrepel`, `fgsea`, `msigdbr`, `gprofiler2` (for archived ORA scripts only).

External:
- `~/Github/scz_cell_type_enrichment/data/gwas/scz_gwas_gene_set_no_mhc.csv` — PGC3 SCZ GWAS gene set
- MSigDB via `msigdbr` (downloads on first use, cached locally)

## Provenance

Generated during a Claude Code session on 2026-05-27 from the meta-analysis output. The pipeline ordering reflects how the work actually evolved: initial DE inspection → Sst volcano → ORA (revealed false negatives) → GSEA (rescued the signal) → cross-cell-type comparison → multi-testing → GWAS layer → mechanistic deep-dives. Each script's docstring records what question it was built to answer.
