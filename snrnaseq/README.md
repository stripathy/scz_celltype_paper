# snrnaseq/ — reserved slot for the upstream pipeline, plus one robustness check

## Reserved slot

This directory will hold the **snRNA-seq DE and cell-type composition
meta-analysis** pipeline (7 SCZ cohorts, 469 donors; Endresz et al., in prep) —
the root of this paper's dependency graph, and not yet accessible.

Everything else in the repo is downstream of it and consumes its exports, which
are documented and symlinked in
[`../shared/snrnaseq_de/`](../shared/snrnaseq_de/README.md).

Merge direction is undecided and the layout supports either: drop that repo in
here, or make this whole tree a subdirectory of it. Paths are repo-relative and
naming is neutral so both work. When access is granted, archive the code in
here, point its outputs at `../shared/snrnaseq_de/`, and delete this section.

## What is here now

[`composition_sensitivity/`](composition_sensitivity/README.md) — a downstream
robustness check on the composition meta-analysis, asking whether the Fig. 3a
abundance result depends on how the seven per-dataset estimates are pooled or on
any one dataset. It builds **Supplementary Fig. S6** and runs from its own
README; it is not part of the reserved slot above.

```bash
Rscript snrnaseq/composition_sensitivity/code/02_plot_pooling_heatmap.R
# -> manuscript/figures/supplementary/S06_composition_pooling.{png,pdf}
```
