# shared/snrnaseq_de/ — upstream snRNA-seq DE interface

The exports of the upstream snRNA-seq meta-analysis pipeline (Endresz et al.,
in prep), consumed by the downstream components. Data here is git-ignored — this
directory is the designated home for large out-of-git inputs.

The pipeline's **code** lives in [`../../snrnaseq/`](../../snrnaseq/README.md)
as of 2026-09-09. Its **outputs** still arrive here as staged copies and
symlinks rather than from `snrnaseq/` directly, because those analyses ran on
the Alliance cluster and the results were staged down by hand. So this directory
remains the seam; what changed is that each export can now be traced to the
script that wrote it.

Entries are either **real files held here** or **symlinks** to a canonical
location in another repo. A file is held here whenever the alternative would be
depending on a volatile path; it is symlinked when the source is another git
repo (stable) or too large to duplicate.

| file | what | written by | storage | consumed by |
|---|---|---|---|---|
| `DE_genes_all_cells_scz.csv` | meta-analytic DE (23 subclasses × genes; estimate, se, padj), 35.6 MB | `snrnaseq/snRNAseq_DE/Subclass/3_meta_analysis.r` (last line) | **real file** (moved out of `~/Downloads` on 2026-07-31 — a Downloads cleanup would have broken every DE figure) | `transcriptomic/` (butterfly, volcano, scatter, forest stars) |
| `meta_results_cohorts_subclass.csv` | per-dataset DE (7 datasets), 328 MB | `snrnaseq/snRNAseq_DE/Subclass/2_DE.r`, collated | symlink → `scz_pathway_enrichment/data/` (another git repo; too large to duplicate) | `transcriptomic/` (forest rows + meta diamond) |
| `nicole_scz_snrnaseq_betas/` | cell-type COMPOSITION betas (crumblr, 7-dataset meta) | `snrnaseq/Compositional_analysis/3_meta_analysis.r` | **regenerated 2026-09-13** (the original was never staged down); see `PROVENANCE.md` beside it | `genetics/` (Fig 4a, 4i), `transcriptomic/` (the S8 strata definition), `snrnaseq/composition_sensitivity/` (S6) |

⚠️ The script that writes the **subclass-level** per-gene meta-analysis feeding
`DE_genes_all_cells_scz.csv` was never in the repo. A **reconstruction** is now
committed as `snrnaseq/snRNAseq_DE/Subclass/3a_meta_per_gene.r`, so the chain
executes end to end — but it is not the original, and it is not bit-exact
against the published table. Nicole should still confirm the copy she ran, and
supply her per-dataset gene universes, which is what exact recovery needs. See
issues 2 and 19 in [`../../KNOWN_ISSUES.md`](../../KNOWN_ISSUES.md).

Downstream figure code should not read this directory directly. It reads
committed snapshots taken from here, guarded by a checksum manifest — see
`../../transcriptomic/data/figure_inputs/README.md`.

Supplementary S8 is the exception. `transcriptomic/scripts/fig5/_common.R` still
points `P$crumblr` at `nicole_scz_snrnaseq_betas/final_results_crumblr_7_cohorts.csv`
here, and that file has no committed snapshot anywhere in the repo, so S8 cannot
fall back to one. `P$subclass` had the same problem and was fixed on 2026-09-13:
it now prefers this directory and falls back to
`transcriptomic/data/figure_inputs/DE_genes_all_cells_scz.csv`, which the
manifest records as a byte-identical capture of the seam file (md5
`719da3d7519cbc5a3aec1b8c42d438da`). See issue 16 in
[`../../KNOWN_ISSUES.md`](../../KNOWN_ISSUES.md).

Two export specs also live here — `EXPORT_SPEC_stratum_pseudobulks.md` and
`EXPORT_SPEC_sst_celllevel.md`. They are the instructions used to produce the
per-dataset pseudobulk and Sst cell-level exports behind Supplementary Fig. S8,
run on the cluster that holds the full per-dataset objects. The exports
themselves are too large to track; their manifests are committed under
`transcriptomic/data/stratum_*_export/`.

**Owner:** Endresz et al. (in prep), whose code is now in
[`../../snrnaseq/`](../../snrnaseq/README.md). The symlinks above still resolve
to the scattered staging locations rather than to one place under `snrnaseq/`,
because that pipeline's outputs live on the cluster. Repointing them is only
worth doing if those outputs are ever brought down locally; until then this
directory stays the single seam between her pipeline and this paper.
