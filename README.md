# scz_celltype_paper — convergent cell-type-specific changes in human schizophrenia

Consolidated code + paper layer for one integrative SCZ paper combining four
evidence streams, each a former standalone repo namespaced as a subdir:

| subdir | source repo | evidence stream |
|---|---|---|
| `snrnaseq/` | **student's upstream pipeline (TBD)** | snRNA-seq DE + composition meta-analysis (Endresz et al., in prep) |
| `genetics/` | scz_cell_type_enrichment | SCZ GWAS → cell-type enrichment; owns `franken_taxonomy/` + GWAS set |
| `spatial/` | SCZ_Xenium | Xenium spatial composition (Kwon et al. 2026); owns Xenium DE |
| `transcriptomic/` | scz_pathway_enrichment | DE pathway/GSEA + cross-platform composite figure (smoke-tested ✓) |
| `histology/` | sgACC_cell_depth_analysis | RNAscope FISH SST density (standalone) |

Shared narrative: **SST / inhibitory (and deep-layer excitatory) dysregulation
in SCZ**, seen convergently across genetics, spatial composition, transcriptomic
DE, and histology.

See **[`DATA_FLOW.md`](DATA_FLOW.md)** for the full dependency graph (who
produces what, who consumes it) and the realized symlink seams.

## Status: SCAFFOLD (4 components in)
All four downstream components are merged as **tracked code only** (`git archive`
— no data, no git history; ~250 MB). The composite figure regenerates inside the
repo (smoke-tested). **Data is git-ignored and currently symlinked back to the
original `~/Github/` repos**, so the monorepo is code-complete but not yet
self-contained for data (see DATA_FLOW.md → "Current state"). Pending: the
`snrnaseq/` upstream pipeline, the convergence figure (`figures/`), and moving
data in-repo.

## The snRNA-seq dependency (root of the graph)
The snRNA-seq DE + composition meta-analysis is produced by a separate
**upstream pipeline owned by the student** (Endresz et al., in prep), not yet
accessible. Everything here is **downstream** of it. Its outputs are the
interface in [`shared/snrnaseq_de/`](shared/snrnaseq_de/). When the pipeline is
available it slots into [`snrnaseq/`](snrnaseq/) — or, equivalently, this repo
merges into hers. **The subdir layout is direction-agnostic on purpose.**

## Data & conventions
Code only is tracked; data is git-ignored and wired via symlinks. Large shared
reference assets stay in the EXTERNAL `~/Github/shared_data/` hub (it serves
non-SCZ projects too) and are symlinked, not absorbed — following the
conventions in `~/Github/DATA_LAYOUT.md`.

## Originals
The source repos under `~/Github/` are untouched and serve as backups; this
monorepo was built by copying their current state (no history needed).
