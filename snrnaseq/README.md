# snrnaseq/ — RESERVED slot for the upstream snRNA-seq pipeline

Placeholder. This subdir will hold the **snRNA-seq DE + cell-type composition
meta-analysis** pipeline (7 SCZ cohorts; Endresz et al., in prep) — the
student's repo, which is the **root of this paper's dependency graph** and is
not yet accessible.

- Everything else in this monorepo (`transcriptomic/`, `spatial/`, `genetics/`,
  `histology/`) is **downstream**: it consumes this pipeline's exports.
- Those exports (the interface) are documented + symlinked in
  [`../shared/snrnaseq_de/`](../shared/snrnaseq_de/).
- **Merge direction is undecided and the layout supports either:** drop her repo
  here, OR merge this repo into hers (this whole tree becomes a subdir there).
  Keep paths repo-relative and naming neutral so both work.

When access is granted: `git archive HEAD` her repo into here (code only), point
its outputs at `../shared/snrnaseq_de/`, and delete this placeholder.
