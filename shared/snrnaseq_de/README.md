# shared/snrnaseq_de/ — upstream snRNA-seq DE interface

The exports of the upstream snRNA-seq meta-analysis pipeline (Endresz et al.,
in prep), consumed by the downstream components. Data here is git-ignored — this
directory is the designated home for large out-of-git inputs — until the
upstream pipeline itself lives in [`../../snrnaseq/`](../../snrnaseq/).

Entries are either **real files held here** or **symlinks** to a canonical
location in another repo. A file is held here whenever the alternative would be
depending on a volatile path; it is symlinked when the source is another git
repo (stable) or too large to duplicate.

| file | what | storage | consumed by |
|---|---|---|---|
| `DE_genes_all_cells_scz.csv` | meta-analytic DE (23 subclasses × genes; estimate, se, padj), 35.6 MB | **real file** (moved out of `~/Downloads` on 2026-07-31 — a Downloads cleanup would have broken every DE figure) | `transcriptomic/` (butterfly, volcano, scatter, forest stars) |
| `meta_results_cohorts_subclass.csv` | per-cohort DE (7 cohorts), 328 MB | symlink → `scz_pathway_enrichment/data/` (another git repo; too large to duplicate) | `transcriptomic/` (forest rows + meta diamond) |
| `nicole_scz_snrnaseq_betas/` | cell-type COMPOSITION betas (crumblr, 7-cohort meta) | symlink → `SCZ_Xenium/data/nicole_scz_snrnaseq_betas/` | `genetics/` (Fig 4a, 4i), `transcriptomic/` (the S8 strata definition), `snrnaseq/composition_sensitivity/` (S6) |

Downstream figure code does not read this directory directly. It reads committed
snapshots taken from here, guarded by a checksum manifest — see
`../../transcriptomic/data/figure_inputs/README.md`.

Two export specs also live here — `EXPORT_SPEC_stratum_pseudobulks.md` and
`EXPORT_SPEC_sst_celllevel.md`. They are the instructions used to produce the
per-cohort pseudobulk and Sst cell-level exports behind Supplementary Fig. S8,
run on the cluster that holds the full per-dataset objects. The exports
themselves are too large to track; their manifests are committed under
`transcriptomic/data/stratum_*_export/`.

**Owner:** Endresz et al. (in prep). When that pipeline is added under
`../../snrnaseq/`, repoint these symlinks to its outputs (one place) instead of
the scattered locations above — that is the single seam between it and this
paper.
