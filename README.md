# scz_celltype_paper

Code and figures for a study of cell-type-specific changes in schizophrenia,
combining a seven-dataset snRNA-seq meta-analysis, Xenium spatial
transcriptomics, common-variant genetics and patch-seq physiology.

**Headline:** upper-layer somatostatin (Sst) interneuron subtypes are
selectively depleted in schizophrenia, the same subtypes carry the most SCZ
common-variant risk, and they are distinguished by *HCN1* and a
sag-dominated intrinsic physiology.

The manuscript itself is not in this repo. This is the analysis code, the
committed figure inputs, and the rendered figures.

---

## Finding the code behind a figure

| Figure | Built by | Component |
|---|---|---|
| **Fig 1** | [`Final_figures/Code/Figure_1.py`](snrnaseq/Final_figures/Code/Figure_1.py), [`Figure_1a_Barchart.r`](snrnaseq/Final_figures/Code/Figure_1a_Barchart.r) | [`snrnaseq/`](snrnaseq/README.md) |
| **Fig 2** | [`scripts/09_composite_figure.R`](transcriptomic/scripts/09_composite_figure.R) | [`transcriptomic/`](transcriptomic/README.md) |
| **Fig 3** | [`Final_figures/Code/Figure_3.r`](snrnaseq/Final_figures/Code/Figure_3.r) | [`snrnaseq/`](snrnaseq/README.md) |
| **Fig 4** | [`scripts/figures/scz_sst_hcn1_story.R`](genetics/scripts/figures/scz_sst_hcn1_story.R) | [`genetics/`](genetics/README.md) |
| **S1** | [`Supplemental/Code/Markerplot.py`](snrnaseq/Final_figures/Supplemental/Code/Markerplot.py) | [`snrnaseq/`](snrnaseq/README.md) |
| **S2** | [`code/analysis/plot_markers_resolvability_combined.R`](spatial/code/analysis/plot_markers_resolvability_combined.R) | [`spatial/`](spatial/README.md) |
| **S3** | [`code/analysis/plot_xenium_merfish_composite.R`](spatial/code/analysis/plot_xenium_merfish_composite.R) | [`spatial/`](spatial/README.md) |
| **S4** | [`Supplemental/Code/Supertype_DE.r`](snrnaseq/Final_figures/Supplemental/Code/Supertype_DE.r) | [`snrnaseq/`](snrnaseq/README.md) |
| **S5** | [`Supplemental/Code/NonNeuron_supplement.r`](snrnaseq/Final_figures/Supplemental/Code/NonNeuron_supplement.r) | [`snrnaseq/`](snrnaseq/README.md) |
| **S6** | [`code/02_plot_pooling_heatmap.R`](snrnaseq/composition_sensitivity/code/02_plot_pooling_heatmap.R) | [`snrnaseq/composition_sensitivity/`](snrnaseq/composition_sensitivity/README.md) |
| **S7** | [`Supplemental/Code/Sensitivity_analysis_Barchart.r`](snrnaseq/Final_figures/Supplemental/Code/Sensitivity_analysis_Barchart.r) | [`snrnaseq/`](snrnaseq/README.md) |
| **S8** | [`scripts/fig5/08_figure5.R`](transcriptomic/scripts/fig5/08_figure5.R) | [`transcriptomic/`](transcriptomic/scripts/fig5/README.md) |
| **S9** | [`scripts/figures/plot_supp_enrichment_seaad125.R`](genetics/scripts/figures/plot_supp_enrichment_seaad125.R) | [`genetics/`](genetics/README.md) |
| **S10** | [`scripts/figures/plot_fig4_robustness.R`](genetics/scripts/figures/plot_fig4_robustness.R) | [`genetics/`](genetics/README.md) |
| **T6** | [`scripts/figures/build_supp_table_patchseq_labels.py`](genetics/scripts/figures/build_supp_table_patchseq_labels.py) | [`genetics/`](genetics/README.md) |

Rendered figures, where they are in the repo, are under
[`manuscript/figures/`](manuscript/figures/main/README.md); that folder's README
says which are present.

Supplementary Methods **SM1** describes the Xenium pipeline in
[`spatial/`](spatial/README.md). Figures 1b–f are assembled by hand from the
annotated Xenium objects that `spatial/` produces.

**Two conventions, not one.** The renderers in `transcriptomic/`, `spatial/`,
`genetics/` and `snrnaseq/composition_sensitivity/` write into
`manuscript/figures/` in place and run from committed inputs, so those figures
regenerate from a clean clone without the raw data. Most renderers in
`snrnaseq/Final_figures/` ran on the Alliance cluster against the full
per-dataset objects and wrote to a working directory there; those are here to be
read rather than re-executed. Each component's README gives the detail.

Items that must be settled before submission are collected in
[`KNOWN_ISSUES.md`](KNOWN_ISSUES.md), with resolved ones listed at its foot.

---

## Components

| Directory | What it does |
|---|---|
| [`transcriptomic/`](transcriptomic/README.md) | Cross-platform differential expression (Fig 2) and the Sst depletion-strata analysis (S8) |
| [`spatial/`](spatial/README.md) | The Xenium pipeline: cell typing, cortical depth, laminar segmentation, plus Xenium DE and composition (SM1, S2, S3; inputs to Figs 1–3) |
| [`genetics/`](genetics/README.md) | MAGMA cell-type enrichment of SCZ common-variant risk, fine-mapping, patch-seq physiology (Fig 4, S9, S10, T6) |
| [`crossdisorder/`](crossdisorder/README.md) | Whether the SCZ-depleted Sst supertypes also decline in Alzheimer's (the AD axis of Fig 4i) |
| [`snrnaseq/`](snrnaseq/README.md) | The seven-dataset snRNA-seq pipeline: label transfer, composition and DE meta-analysis (Fig 1a, Fig 3, S1, S4, S5, S7), plus a composition robustness check (S6) |
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
fetch first (Figure 2 excepted — see below). Each component's `data/README.md`
lists what to obtain and from where.

**Independently reproduced.** On 2026-09-13 the composition chain behind Figure
3a and S6 was rebuilt from a separate set of cleaned per-dataset h5ads and
compared against the committed artefacts. The Figure-3a input matrix
(`snrnaseq/Compositional_analysis/Files/7_cohorts_metadata_names.csv`) came back
**identical in all 61,908 donor x supertype counts** (2,399,784 cells), with Age,
Sex and PMI matching exactly; re-fitting the crumblr models reproduced the
per-dataset, leave-one-dataset-out and subclass estimates to **1e-15**. The
Xenium composition set reproduces exactly too (356,313 neuronal and 385,790
non-neuronal cortical cells, all 24 donors). Two seams that the comparison
exposed are recorded as issues 17 and 18 in
[`KNOWN_ISSUES.md`](KNOWN_ISSUES.md).

## Reproducing

```bash
git clone <remote-url> && cd scz_celltype_paper
```

Then follow the component README for the figure you want. Figure 4 and
supplementary figures S2, S3, S6, S8, S9 and S10 render from committed inputs;
verified by running all of them on 2026-09-13. Re-running the *analyses* behind
them needs the external data.

Figure 2 renders from a clone as well, since 2026-09-14: the three inputs that
had lived under `/scratch/nendresz/` are committed snapshots with `MANIFEST.tsv`
rows (issue 15, closed).

Two prerequisites, both learned the hard way:

- **R packages.** Beyond the usual tidyverse/ggplot2 stack the renderers need
  `ragg` and `ggsignif`, and Figure 4 needs `svglite`. Missing any of them fails
  only at the final `ggsave`, after the whole figure has been computed.
- **Working directory.** Most renderers locate their own root from `--file=` and
  run from anywhere. The two in `spatial/` do not: run them **from `spatial/`**
  (`cd spatial && Rscript code/analysis/<script>.R`), or their relative
  `source()` and input paths miss.

Figures 1 and 3 and supplementary figures S1, S4 and S5 are built by
[`snrnaseq/Final_figures/`](snrnaseq/Final_figures/README.md). Most do **not**
render from a clone — those scripts ran on the Alliance cluster against the full
per-dataset objects, with cluster-absolute paths, and are here to be read. Two
exceptions: **S7** renders anywhere (both its inputs are committed), and
**Figure 3 panels d, e and h** do since 2026-09-15, because they need only the
committed Xenium metadata rather than the Seurat object.
The one committed input from that pipeline is the per-donor cell-count matrix
behind the composition analysis,
[`snrnaseq/Compositional_analysis/Files/7_cohorts_metadata_names.csv`](snrnaseq/Compositional_analysis/Files/7_cohorts_metadata_names.csv).

## Repository history

The tree was pruned on 2026-09-04 down to the code behind the paper, so that
finding a figure's source does not mean walking past a much larger body of work
that did not make it in. The state before that pass is the git tag
**`pre-prune-2026-09-04`**; nothing was lost, and each pruning commit says what
it removed and why.
