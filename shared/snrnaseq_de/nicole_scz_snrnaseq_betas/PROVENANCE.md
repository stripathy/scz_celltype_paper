# final_results_crumblr_7_cohorts.csv — provenance

Fixed-effect meta-analysis of cell-type composition across the seven snRNA-seq
datasets: one row per neuronal supertype, `estimate` being the SCZ-vs-control
CLR abundance change. Read by `genetics/` (Fig. 4a, 4i), `transcriptomic/`
(the S8 strata definition) and `snrnaseq/composition_sensitivity/` (S6).

Unlike the rest of `shared/snrnaseq_de/`, this file is **tracked in git**. The
original was never staged down from the cluster, which left every consumer above
unable to resolve the seam from a clone. It is 17 KB, so it is committed rather
than symlinked.

## How this copy was produced

Generated 2026-09-13 by running the repo's own scripts, unmodified except for
their input and output paths:

    snrnaseq/Compositional_analysis/2_Crumblr_analysis.r   per-dataset crumblr/dream
    snrnaseq/Compositional_analysis/3_meta_analysis.r      fixed-effect rma pooling

The input was a donor x supertype count matrix rebuilt independently from cleaned
per-dataset h5ads, not from `7_cohorts_metadata_names.csv`.

## Checks

Every comparison that can be made against a committed artefact passes:

| check | result |
|---|---|
| the rebuilt input vs `7_cohorts_metadata_names.csv` | 61,908 / 61,908 counts identical |
| per-dataset betas vs `composition_per_dataset_estimates.csv` | 752/752 rows, max abs diff 1.8e-15 |
| this file vs the `FE meta` rows of `composition_pooling_sensitivity.csv` | 109/109, max abs diff 8.1e-16 |
| fed to `genetics/scripts/figures/build_composition_table.py`, reproduces `gwas_vs_casecontrol_composition.csv` | 109/109 supertypes, max abs diff 1.0e-15 |

Quick sanity value: **Sst_25 estimate = -0.258993, p = 1.719e-03**.

Full log: `/scratch/shreejoy/celltype_repro/out/REPRODUCTION_LOG.md`.
