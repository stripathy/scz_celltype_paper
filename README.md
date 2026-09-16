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
| **S8** | [`scripts/sst_strata/08_figure.R`](transcriptomic/scripts/sst_strata/08_figure.R) | [`transcriptomic/`](transcriptomic/scripts/sst_strata/README.md) |
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
root of the dependency graph — everything else is downstream of it. The
interface through which the other components consume it is documented in
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

### Published data (Zenodo)

The data this study generated is deposited at
**[Zenodo record 22801230](https://zenodo.org/records/22801230)**
(DOI [10.5281/zenodo.22801230](https://doi.org/10.5281/zenodo.22801230),
CC-BY-4.0, published 2026-09-16): harmonized per-cell cell-type annotations for
all eight datasets, and the cell-type-specific differential expression results at
subclass and supertype resolution.

| File | Size | md5 | Contents |
|---|---|---|---|
| `cell_metadata.parquet` | 53.5 MB | `4d51070c637070ec4ab2075d175afbc7` | 3,738,935 cells x 14 columns; 493 donors, 8 datasets (7 snRNA-seq + Xenium); `class` / `subclass` / `supertype`, `diagnosis`, `sex`, `age`, `n_counts`, `n_genes`, `qc_pass`, `spatial_domain` |
| `DE_subclass.parquet` | 73.7 MB | `5da0d816f713191645dd9bf82ece8dd1` | 2,471,813 rows: `gene, logFC, SE, PValue, FDR, cell_type, cohort`; 48,252 genes x 24 subclasses x 9 cohorts (the 7 datasets, `Meta-analysis`, and `Xenium`) |
| `DE_supertype.parquet` | 223.5 MB | `4d6029d131eb51cb3846e5a9cc1c4796` | the same at supertype resolution, for the 7 snRNA-seq datasets and the meta-analysis |

**This is the citable route to the large tables that are otherwise in Git LFS.**
Verified 2026-09-16: the `Meta-analysis` rows of `DE_subclass.parquet` are the
same data as the committed
`transcriptomic/data/figure_inputs/DE_genes_all_cells_scz.csv` — 231,135 rows,
16,361 genes, 24 cell types, agreeing to 1e-16 — and its per-dataset rows are the
public form of `snrnaseq/snRNAseq_DE/Files/DE_genes_all_cohorts_subclass.csv`
(the 301 MB LFS object). `cell_metadata.parquet` is the harmonized annotation set
behind the composition analyses.

Two schema notes for anyone joining Zenodo to repo files: Zenodo uses the
`Control` / `SCZ` diagnosis vocabulary and the dataset names `MSSM1` / `MSSM2` /
`Frohlich`, whereas the in-repo `snrnaseq/` tables use `Control` /
`Schizophrenia` and `MtSinai` / `MSSM` / `OFC`. The paper's labels are MSSM 1,
MSSM 2 and Fröhlich.

Supplementary Tables T1–T4 are published with the article, not here.

## Reproducing

```bash
git clone <remote-url> && cd scz_celltype_paper
git lfs install && git lfs pull        # REQUIRED — see below
```

**This repository uses Git LFS.** Twenty-three files are stored there, including
the large per-cohort DE tables under `snrnaseq/snRNAseq_DE/Files/` and
`Supertypes/Files/`. Without `git lfs pull` those paths hold a 134-byte pointer
stub instead of the data, and **nothing errors**: `read.csv()` on a pointer
returns a two-row frame of pointer text, so a script can run to completion on
garbage. The one guard that does catch it is the figure-input manifest, which
refuses to draw Figure 2 and reports the source as changed — the symptom that
led here on 2026-09-16. If a chain behaves strangely, check for stubs first:

```bash
head -c 40 <path> | grep -q git-lfs && echo "pointer stub, run: git lfs pull"
```

If you would rather not fetch LFS at all, the two largest DE tables are also on
Zenodo in parquet form — see **Published data** above.

Then follow the component README for the figure you want. Figure 2, Figure 4
and supplementary figures S2, S3, S6, S8, S9 and S10 render from committed
inputs. Re-running the *analyses* behind them needs the external data.

Two prerequisites:

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
exceptions: **S7** renders anywhere (both its inputs are committed), and so do
**Figure 3 panels d and e**, which need only the committed Xenium metadata
rather than the Seurat object.
The one committed input from that pipeline is the per-donor cell-count matrix
behind the composition analysis,
[`snrnaseq/Compositional_analysis/Files/7_cohorts_metadata_names.csv`](snrnaseq/Compositional_analysis/Files/7_cohorts_metadata_names.csv).
