# shared/snrnaseq_de/ — upstream snRNA-seq DE interface

The exports of the upstream snRNA-seq meta-analysis pipeline (Endresz et al.,
in prep), consumed by the downstream components. Data here is git-ignored — this
directory is the designated home for large out-of-git inputs. The single
exception is `nicole_scz_snrnaseq_betas/final_results_crumblr_7_cohorts.csv`,
which is small, load-bearing for three components, and had no copy on any machine
but the cluster; it is tracked so a clone resolves the seam.

The pipeline's **code** lives in [`../../snrnaseq/`](../../snrnaseq/README.md).
Its **outputs** arrive here as staged copies and symlinks rather than from
`snrnaseq/` directly, because those analyses run on the Alliance cluster and the
results are staged down by hand. This directory is the seam, and each export
here can be traced to the
script that wrote it.

Entries are either **real files held here** or **symlinks** to a canonical
location in another repo. A file is held here whenever the alternative would be
depending on a volatile path; it is symlinked when the source is another git
repo (stable) or too large to duplicate.

| file | what | written by | storage | consumed by |
|---|---|---|---|---|
| `DE_genes_all_cells_scz.csv` | meta-analytic DE (23 subclasses × genes; estimate, se, padj), 35.6 MB | `snrnaseq/snRNAseq_DE/Subclass/3_meta_analysis.r` (last line) | **real file** (moved out of `~/Downloads` on 2026-07-31 — a Downloads cleanup would have broken every DE figure) | `transcriptomic/` (butterfly, volcano, scatter, forest stars) |
| `meta_results_cohorts_subclass.csv` | per-dataset DE (7 datasets), 328 MB | `snrnaseq/snRNAseq_DE/Subclass/2_DE.r`, collated | symlink → `scz_pathway_enrichment/data/` (another git repo; too large to duplicate) | `transcriptomic/` (forest rows + meta diamond) |
| `nicole_scz_snrnaseq_betas/` | cell-type COMPOSITION betas (crumblr, 7-dataset meta) | `snrnaseq/Compositional_analysis/3_meta_analysis.r` | **tracked in git** — the one exception in this directory; regenerated 2026-09-13 because the original was never staged down, see `PROVENANCE.md` beside it | `genetics/` (Fig 4a, 4i), `transcriptomic/` (the S8 strata definition), `snrnaseq/composition_sensitivity/` (S6) |

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

Supplementary S8 reads two entries here directly. `P$crumblr` points at
`nicole_scz_snrnaseq_betas/final_results_crumblr_7_cohorts.csv`, which is
tracked, so it resolves from a clone. For subclass DE, S8 always reads the
committed snapshot `transcriptomic/data/figure_inputs/DE_genes_all_cells_scz.csv`
(md5 `f7ff12756d1ffcc143e2f723ea4fd249`), **not** this directory.

That distinction is load-bearing, because a superseded subclass DE table is
still in circulation on some machines. It has 222,028 rows and yields **345**
Sst genes at FDR < 0.10; the canonical table has 231,135 rows and yields
**343**. Reading the seam in preference to the snapshot draws one set of numbers
locally and another on a clone. **343 is the current number.**

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
