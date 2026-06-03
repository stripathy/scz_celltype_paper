# shared/snrnaseq_de/ — upstream snRNA-seq DE interface

The exports of the student's snRNA-seq meta-analysis pipeline (Endresz et al.,
in prep), consumed by the downstream subdirs. Data is git-ignored and symlinked
to the current canonical locations until the upstream pipeline itself lives in
[`../../snrnaseq/`](../../snrnaseq/).

| file | what | current source (symlink) | consumed by |
|---|---|---|---|
| `DE_genes_all_cells_scz.csv` | meta-analytic DE (23 subclasses × genes; estimate, se, padj) | `~/Downloads/` export | `transcriptomic/` (butterfly, volcano, scatter, forest stars) |
| `meta_results_cohorts_subclass.csv` | per-cohort DE (7 cohorts) | `scz_pathway_enrichment/data/` | `transcriptomic/` (forest rows + meta diamond) |
| `nicole_scz_snrnaseq_betas/` | cell-type COMPOSITION betas (crumblr) | `SCZ_Xenium/data/nicole_scz_snrnaseq_betas/` *(link when genetics/+spatial/ wired)* | `genetics/` (GWAS × composition), `spatial/` (concordance) |

**Owner:** Endresz et al. (in prep). When her pipeline is added under
`../../snrnaseq/`, repoint these symlinks to its outputs (one place) instead of
the scattered locations above — that is the single seam between her work and
this paper.
