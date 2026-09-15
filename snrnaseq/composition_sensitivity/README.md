# snRNA-seq composition: pooling-strategy and single-dataset sensitivity

Does the cell-type abundance result of Fig. 3a depend on how the seven per-dataset crumblr
estimates are pooled, or on any one dataset? This is **Supplementary Fig. S6**. For the
109 neuronal supertypes this folder
compares the fixed-effect meta-analysis used in the paper against a random-effects
meta-analysis and two mega-analyses that fit all 469 donors in a single model (dataset as a
random intercept, without and with a dataset-specific random SCZ slope), reports each
dataset's own estimate alongside the Xenium estimate, and runs a leave-one-dataset-out
fixed-effect meta-analysis. Output: the supplementary figure
(`manuscript/figures/supplementary/S06_composition_pooling.{png,pdf}`, written in place by
`02`; legend in `results/supp_composition_pooling_legend.md`) and the tables behind it.

A dataset-as-fixed-covariate mega-analysis is also fit and kept in the results CSV, but is
not plotted: it is numerically indistinguishable from the random-intercept version
(β Pearson r = 0.9996, identical FDR sets), so only one of the two is shown.

## Run it

Scripts locate the repo from their own path, so they can be run from any working directory.

```
Rscript code/00_prepare_counts.R  [path/to/7_cohorts_metadata_names.csv]  # verify (add --write to rebuild)
Rscript code/01_fit_pooling_models.R      # per-dataset + pooled models, LODO   (~2 min)
Rscript code/01b_fit_subclass_models.R    # the same models on subclass-level counts
Rscript code/02_plot_pooling_heatmap.R    # the figure
```

R 4.5.1 with crumblr 1.0.0, variancePartition 1.38.1 (dream/lmer backend), metafor 4.8.0,
ggplot2 4.0.2, cowplot 1.2.0, ggrepel 0.9.6, readr, dplyr, tidyr, jsonlite.

Steps 1–3 need nothing outside this folder except the two tracked Xenium tables
(`spatial/output/crumblr/crumblr_{results,input}_supertype_neuronal.csv`), which supply the
Xenium row of panel a. Step 0 needs Nicole's raw export, which is **not** tracked here.

## Inputs

| File | Tracked | Notes |
|---|---|---|
| `data/neuron_counts_469donors.csv` | yes | Analysis table: one row per donor, `dataset, donor, dx, age, sex, PMI` + 109 supertype counts. Built by `00` from the raw export. |
| `data/neuronal_supertypes_109.csv` | yes | The 109 neuronal supertypes, in Fig. 3a order. Defines what `00` selects and what `02` plots. |
| `data/seaad_supertype_colors.json` | yes | Supertype palette of Figs 1b/3a/4. Same file as `genetics/data/seaad_supertype_colors.json`. |
| `7_cohorts_metadata_names.csv` | no | Nicole's raw per-donor counts + metadata (donor rows, cell-type columns, both neuronal and non-neuronal). Pass its path to `00`; the default is `~/Downloads/`. |
| `shared/snrnaseq_de/nicole_scz_snrnaseq_betas/final_results_crumblr_7_cohorts.csv` | yes | The published FE meta-analysis. Optional: when present, `01` cross-checks the rebuilt FE estimates against it and uses it for the "FE meta-analysis (paper)" row; when absent, `01` uses its own FE meta-analysis for that row and says so. |

**Running it on your own copy of the data.** Point `00` at your export and add `--write`, then
re-run `01`, `01b` and `02`. Without `--write`, `00` only reports whether it reproduces the
tracked table, which is the quick way to diff a new export against the one used here.

`00` harmonizes what differed between datasets in the raw export: `Sex` and `Diagnosis` case
(`Schizophrenia` → `SCZ`) and dataset naming (`MSSM` → MSSM 2, `MtSinai` → MSSM 1, `OFC` →
Fröhlich). PMI is in hours everywhere; the Multiome dataset has none and one MSSM 2 donor is
coded 0 h, which `01` treats as missing (PMI is dropped from the Multiome model and
mean-imputed in the mega models). Supertypes with no counts in a dataset are dropped from
that dataset's model, so the meta-analyses pool 5–7 datasets per supertype.

## Outputs

| File | Contents |
|---|---|
| `results/composition_per_dataset_estimates.csv` | β, SE, P per supertype × dataset |
| `results/composition_pooling_sensitivity.csv` | long table: supertype × method (FE, RE, three mega models, published FE), with FDR, k, I², Cochran's Q |
| `results/composition_lodo.csv` | FE meta-analysis with each dataset omitted in turn |
| `results/composition_subclass_estimates.csv`, `composition_subclass_lodo.csv` | the same at subclass level (18 neuronal subclasses), for the composite Sst statement in the text |
| `manuscript/figures/supplementary/S06_composition_pooling.{png,pdf}` | the figure (written in place; legend and verified numbers in `results/supp_composition_pooling_legend.md`) |

The rebuilt fixed-effect meta-analysis reproduces Nicole's published
`final_results_crumblr_7_cohorts.csv` exactly (Spearman ρ = 1.00, identical β and P).
