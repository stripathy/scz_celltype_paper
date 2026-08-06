# shared/snrnaseq_de/ — upstream snRNA-seq DE interface

The exports of the student's snRNA-seq meta-analysis pipeline (Endresz et al.,
in prep), consumed by the downstream subdirs. Data here is git-ignored — this
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
| `nicole_scz_snrnaseq_betas/` | cell-type COMPOSITION betas (crumblr) | symlink → `SCZ_Xenium/data/nicole_scz_snrnaseq_betas/` | `genetics/` (GWAS × composition), `spatial/` (concordance) |

Downstream figure code does not read this directory directly. It reads committed
snapshots taken from here, guarded by a checksum manifest — see
`../../transcriptomic/data/figure_inputs/README.md`.

**Owner:** Endresz et al. (in prep). When her pipeline is added under
`../../snrnaseq/`, repoint these symlinks to its outputs (one place) instead of
the scattered locations above — that is the single seam between her work and
this paper.
