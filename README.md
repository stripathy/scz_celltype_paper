# scz_celltype_paper

Code and figures for a study of cell-type-specific changes in schizophrenia,
combining a seven-cohort snRNA-seq meta-analysis, Xenium spatial
transcriptomics, common-variant genetics and patch-seq physiology.

**Headline:** upper-layer somatostatin (Sst) interneuron subtypes are
selectively depleted in schizophrenia, the same subtypes carry the most SCZ
common-variant risk, and they are distinguished by *HCN1* and a
sag-dominated intrinsic physiology.

The manuscript itself is not in this repo. This is the analysis code, the
committed figure inputs, and the rendered figures.

---

## Finding the code behind a figure

| Figure | Rendered file | Built by | Component |
|---|---|---|---|
| **Fig 1a** | not rendered here | `Final_figures/Figure_1a_UMAP.r`, `Figure_1a_Barchart.r` | [`snrnaseq/`](snrnaseq/README.md) |
| **Fig 2** | `manuscript/figures/main/Fig2_cross_platform_de` | `scripts/09_composite_figure.R` | [`transcriptomic/`](transcriptomic/README.md) |
| **Fig 3** | not rendered here | `Final_figures/Figure_3.r` | [`snrnaseq/`](snrnaseq/README.md) |
| **Fig 4** | `manuscript/figures/main/Fig4_scz_genetics_hcn1_sst` | `scripts/figures/scz_sst_hcn1_story.R` | [`genetics/`](genetics/README.md) |
| **S1** | not rendered here | `Final_figures/Supplemental/Markerplot.py` | [`snrnaseq/`](snrnaseq/README.md) |
| **S2** | `manuscript/figures/supplementary/S02_xenium_celltype_annotation` | `plot_markers_resolvability_combined.R` | [`spatial/`](spatial/README.md) |
| **S3** | `.../S03_xenium_merfish_concordance` | `plot_xenium_merfish_composite.R` | [`spatial/`](spatial/README.md) |
| **S4** | not rendered here | `Final_figures/Supplemental/Supertype_DE.r` | [`snrnaseq/`](snrnaseq/README.md) |
| **S5** | not rendered here | `Final_figures/Supplemental/NonNeuron_supplement.r` | [`snrnaseq/`](snrnaseq/README.md) |
| **S6** | `.../S06_composition_pooling` | `02_plot_pooling_heatmap.R` | [`snrnaseq/composition_sensitivity/`](snrnaseq/composition_sensitivity/README.md) |
| **S7** | not rendered here | `Final_figures/Supplemental/Sensitivity_analysis_Barchart.r` | [`snrnaseq/`](snrnaseq/README.md) |
| **S8** | `.../S08_sst_strata` | `scripts/fig5/08_figure5.R` | [`transcriptomic/`](transcriptomic/scripts/fig5/README.md) |
| **S9** | `.../S09_scz_enrichment_seaad125` | `plot_supp_enrichment_seaad125.R` | [`genetics/`](genetics/README.md) |
| **S10** | `.../S10_genetics_ad_robustness` | `plot_fig4_robustness.R` | [`genetics/`](genetics/README.md) |
| **T6** | `genetics/results/tables/supp_T6_patchseq_sst_annotations.csv` | `build_supp_table_patchseq_labels.py` | [`genetics/`](genetics/README.md) |

Supplementary Methods **SM1** describes the Xenium pipeline in
[`spatial/`](spatial/README.md). Figures 1b–f are assembled by hand from the
annotated Xenium objects that `spatial/` produces.

**Two conventions, not one.** The renderers in `transcriptomic/`, `spatial/`,
`genetics/` and `snrnaseq/composition_sensitivity/` write into
`manuscript/figures/` in place and run from committed inputs, so those figures
regenerate from a clean clone without the raw data. The renderers in
`snrnaseq/Final_figures/` ran on the Alliance cluster against the full
per-cohort objects and wrote to a working directory there — their figures are
**not** in this repo and their chains are archived for reading, not
re-execution. Each component's README gives the detail.

Open items from the 2026-09-13 audit are collected in
[`KNOWN_ISSUES.md`](KNOWN_ISSUES.md).

---

## Components

| Directory | What it does |
|---|---|
| [`transcriptomic/`](transcriptomic/README.md) | Cross-platform differential expression (Fig 2) and the Sst depletion-strata analysis (S8) |
| [`spatial/`](spatial/README.md) | The Xenium pipeline: cell typing, cortical depth, laminar segmentation, plus Xenium DE and composition (SM1, S2, S3; inputs to Figs 1–3) |
| [`genetics/`](genetics/README.md) | MAGMA cell-type enrichment of SCZ common-variant risk, fine-mapping, patch-seq physiology (Fig 4, S9, S10, T6) |
| [`crossdisorder/`](crossdisorder/README.md) | Whether the SCZ-depleted Sst supertypes also decline in Alzheimer's (the AD axis of Fig 4i) |
| [`snrnaseq/`](snrnaseq/README.md) | The seven-cohort snRNA-seq pipeline: label transfer, composition and DE meta-analysis (Fig 1a, Fig 3, S1, S4, S5, S7), plus a composition robustness check (S6) |
| [`shared/`](shared/README.md) | Cross-component interfaces: the snRNA-seq DE seam, and the figure-input staleness guard |
| [`manuscript/figures/`](manuscript/figures/main/README.md) | The rendered figures, main and supplementary |
| [`reserve/`](reserve/README.md) | Analyses built and deliberately left out of the paper |

The snRNA-seq DE and composition meta-analysis (Endresz et al., in prep) is the
root of the dependency graph — everything else is downstream of it. It was
merged into [`snrnaseq/`](snrnaseq/README.md) on 2026-09-09; the interface
through which the other components consume it is documented in
[`shared/snrnaseq_de/`](shared/snrnaseq_de/README.md).

See [`DATA_FLOW.md`](DATA_FLOW.md) for who produces what and who consumes it.

---

## Data

**Code is tracked; large data is not.** Raw Xenium images, sequencing data, GWAS
summary statistics and reference atlases are git-ignored and wired in via
symlinks or documented download URLs. What *is* committed is the small set of
figure inputs each renderer needs, so no figure depends on data you have to
fetch first. Each component's `data/README.md` lists what to obtain and from
where.

## Reproducing

```bash
git clone <remote-url> && cd scz_celltype_paper
```

Then follow the component README for the figure you want. Figures 2 and 4 and
supplementary figures S2, S3, S6, S8, S9 and S10 all render from committed
inputs. Re-running the *analyses* behind them needs the external data.

Figure 1a, Figure 3 and supplementary figures S1, S4, S5 and S7 are built by
[`snrnaseq/Final_figures/`](snrnaseq/Final_figures/README.md) and do **not**
render from a clone — those scripts ran on the Alliance cluster against the
full per-cohort objects, with cluster-absolute paths. They are here to be read.
The one committed input from that pipeline is the per-donor cell-count matrix
behind the composition analysis,
[`snrnaseq/Compositional_analysis/7_cohorts_metadata_names.csv`](snrnaseq/Compositional_analysis/7_cohorts_metadata_names.csv).

## Repository history

The tree was pruned on 2026-09-04 down to the code behind the paper, so that
finding a figure's source does not mean walking past a much larger body of work
that did not make it in. The state before that pass is the git tag
**`pre-prune-2026-09-04`**; nothing was lost, and each pruning commit says what
it removed and why.
