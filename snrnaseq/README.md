# snrnaseq/ — the seven-dataset snRNA-seq meta-analysis

The **root of this paper's dependency graph** (Endresz et al., in prep). Every
other component consumes its exports: cell-type composition betas and
meta-analytic differential expression across seven SCZ snRNA-seq datasets,
469 donors, labelled against the SEA-AD supertype taxonomy.

It builds **Figure 1a**, **Figure 3**, and Supplementary Figs. **S1, S4, S5 and
S7**; it supplies the inputs behind **Figure 2**, **Figure 4** and
Supplementary Figs. **S6 and S8** in the other components.

> The subdirectory READMEs were reconstructed from the code on 2026-09-13, not
> written by their author — they describe what the scripts do, and Nicole should
> correct anything that misreads intent. Open items are collected in
> [`../KNOWN_ISSUES.md`](../KNOWN_ISSUES.md).

## The datasets

Seven datasets, assembled from five sources. Two naming points that otherwise
cause confusion: Ruzicka's *MtSinai* is the paper's **MSSM 1**, while PsychAD's
*MSSM* is the paper's **MSSM 2** (they share donors — see the composition
sensitivity work); and `OFC` throughout the code is the **Fröhlich** dataset.

| Paper label | In the code | Source | Object built by |
|---|---|---|---|
| HBCC | `HBCC` | PsychAD | `Label_transfer/PsychAD/` |
| MSSM 2 | `MSSM` | PsychAD | `Label_transfer/PsychAD/` |
| MSSM 1 | `MtSinai`, `Ruz_MtSinai` | Ruzicka (PsychENCODE) | `Label_transfer/Ruzicka/` |
| McLean | `McLean`, `Ruz_McLean` | Ruzicka (PsychENCODE) | `Label_transfer/Ruzicka/` |
| Fröhlich | `OFC` | Fröhlich | `Label_transfer/Frohlich/` |
| Batiuk | `Batiuk`, `Bat` | Batiuk | `Label_transfer/Batiuk/` |
| Multiome | `Multi` | PsychENCODE2 / brainSCOPE | `Label_transfer/Multiome/` |

**"Dataset", not "cohort".** The paper calls these seven units *datasets*,
following Kiss et al.: the Mount Sinai brain bank contributes two of them, so a
dataset is the right unit and "cohort" is reserved for a single brain bank's
donor group (the Xenium/LIBD cohort, the RNAscope cohort). The code predates
that convention and calls them cohorts throughout — a `Cohort` metadata column,
`cohorts_use`, `7_cohorts_metadata_names.csv`. These READMEs follow the paper in
prose and keep `cohort` only when naming an actual column, variable or file.

Donor inclusion is applied at load time, not in a separate filtering step:
age < 70 everywhere, plus age > 20 for HBCC, and Multiome restricted to
`Disorder ∈ {control, Schizophrenia}`.

## Build order

Each stage is a numbered chain inside its own directory; see that directory's
README for the per-script detail.

| Stage | Directory | Produces |
|---|---|---|
| 1 | [`Label_transfer/`](Label_transfer/README.md) | SEA-AD supertype labels on all seven datasets (`predicted.id`), as per-dataset Seurat objects |
| 2a | [`Compositional_analysis/`](Compositional_analysis/README.md) | per-donor cell counts → crumblr per dataset → fixed-effect meta-analysis (**Fig. 3a**, **S5**) |
| 2b | [`snRNAseq_DE/`](snRNAseq_DE/README.md) | pseudobulk → limma-voom DE per dataset → per-gene meta-analysis, at subclass and supertype level (**S4**; inputs to **Fig. 2**) |
| 3 | [`Compositional_sensitivity_analysis/`](Compositional_sensitivity_analysis/README.md) | stage 1–2a repeated with Sst DE genes withheld from the reference (**S7**) |
| 4 | [`Final_figures/`](Final_figures/README.md) | **Fig. 1a**, **Fig. 3**, **S1**, **S4**, **S5**, **S7** |
| — | [`composition_sensitivity/`](composition_sensitivity/README.md) | a *downstream* robustness check on stage 2a, run from this repo (**S6**) |

`composition_sensitivity/` is not part of the upstream pipeline — it is a
separate check asking whether the Fig. 3a result depends on how the seven
per-dataset estimates are pooled, and it is the one chain here that runs from
committed inputs:

```bash
Rscript snrnaseq/composition_sensitivity/code/02_plot_pooling_heatmap.R
# -> manuscript/figures/supplementary/S06_composition_pooling.{png,pdf}
```

## The reference

All label transfer uses one neutral reference — `raw_counts_ref.rds` plus
`Neurotypical_ref_metadata.rds` — and transfers its `Supertype` column via
Seurat `FindTransferAnchors` / `TransferData`. This is the **SEA-AD**
neurotypical taxonomy (Gabitto et al. 2024); the transferred labels carry the
`-SEAAD` suffix and are the same supertype names used by `genetics/` and
`spatial/`.

The reference is loaded into variables named `counts_hodge` / `meta_hodge`
throughout. That naming is a leftover and is **misleading** — it is not the
Hodge et al. 2019 MTG taxonomy. Logged as issue 10 in
[`../KNOWN_ISSUES.md`](../KNOWN_ISSUES.md).

## Running this code

**These scripts do not run from a clone.** They ran on the Alliance cluster
against the full per-dataset objects: paths are cluster-absolute
(`/scratch/nendresz/…`, `/project/rrg-shreejoy/…`) and most scripts open with a
`setwd()` naming a working directory that is not part of this repo. They are
archived here so the analysis can be read and checked, not re-executed.

Environment is selected per script by a `conda activate` comment in the header
(`crumblr_env`, `de_env`, `r_env_meta_analysis`, `py_anndata_env`, `brisc`).

The one committed input is
[`Compositional_analysis/Files/7_cohorts_metadata_names.csv`](Compositional_analysis/Files/7_cohorts_metadata_names.csv)
— the per-donor cell-count matrix with age, sex, diagnosis and PMI for all 469
donors, which is what stage 2a actually fits.

## What this exports downstream

Consumed through [`../shared/snrnaseq_de/`](../shared/snrnaseq_de/README.md):

| Export | Produced by | Consumed by |
|---|---|---|
| `DE_genes_all_cells_scz.csv` — subclass-level meta DE | `snRNAseq_DE/Subclass/4_Compile_results.r` | `transcriptomic/` (Fig. 2, S8) |
| `meta_results_cohorts_subclass.csv` — per-dataset DE | `snRNAseq_DE/Subclass/` | `transcriptomic/` (Fig. 2 forests) |
| composition betas (crumblr, 7-dataset meta) | `Compositional_analysis/Code/Neurons/2_meta_analysis.r` | `genetics/` (Fig. 4a, 4i), `transcriptomic/` (S8 strata), `composition_sensitivity/` (S6) |

`snRNAseq_DE/Subclass/4_Compile_results.r` writes `DE_genes_all_cells_scz.csv`.
Its last line targets `/scratch/nendresz/scz_celltype_paper/transcriptomic/data/
figure_inputs/`, a cluster-absolute path, so the seam is live only on the machine
that ran it; elsewhere the committed snapshot is refreshed by
`transcriptomic/scripts/00_refresh_figure_inputs.R`.

See [`../DATA_FLOW.md`](../DATA_FLOW.md) for the whole graph.
