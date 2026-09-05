# reserve/ — built, checked, and deliberately not in the paper

Analyses that were carried far enough to trust and then left out. They are kept
because each answers a question a reviewer is likely to ask, and because the
decision to cut them was editorial rather than a finding that they were wrong.

**Nothing in this directory produces a figure, table or number in the
manuscript.** If you are looking for the code behind a paper figure, you want
`genetics/`, `spatial/`, `transcriptomic/`, `snrnaseq/` or `crossdisorder/` —
see the map in the root `README.md`.

| Directory | What it is | Why it is here |
|---|---|---|
| [`histology/`](histology/) | RNAscope re-analysis of SST⁺ interneuron density in sgACC (Arbabi, Newton et al. 2025), layer-stratified with a VIP laminar quality filter | Cut 2026-09-01 on Etienne's advice. The effect is in the predicted direction but does not reach significance (SCZ β = −0.60 per frame, *P* = 0.071; L2/3 −0.76 vs L5/6 −0.43), and it conflates cell number with per-cell SST mRNA. Revisit if reviewers ask for protein- or RNA-level corroboration. |
| [`sst_strata_supp/`](sst_strata_supp/) | Supplements and sensitivity analyses for the Sst depletion-strata result (Supplementary Fig. S8) | The S8 section was shortened when the analysis moved from main figure to supplement. Holds the controls figure that was dropped, the independent-cohort replication, the Alzheimer's contrast, and the cell-count confound tests. |
| [`gse158516/`](gse158516/) | Reiner et al. reprocessed as an independent eighth SCZ cohort | Not in the meta-analysis and not cited. Reprocessed 2026-08-14 to ask whether an eighth dataset agrees; it does. **Its cell-type labels are provisional** and must not be used without re-annotation — see its README. |
| [`percell_normalisation/`](percell_normalisation/) | SST and PVALB effects under four per-cell normalisation schemes | Was proposed as Supplementary Fig. S7; that slot now holds the re-annotation control instead. The robustness claim it supports is not currently made in the paper. |

Two rendered figures that are likewise out of the current version live beside
the submission set, in
[`../manuscript/figures/not_in_current_version/`](../manuscript/figures/not_in_current_version/),
so the supplementary folder stays exactly the submission set.

## Running any of it

The strata supplements are run **from the repo root** and source the live
pipeline's `_common.R`, so they inherit the same palettes, module definitions
and stratum membership as Supplementary Fig. S8:

```bash
Rscript reserve/sst_strata_supp/volcanoes.R
```

They require the core S8 pipeline (`transcriptomic/scripts/fig5/01`–`07`) to
have run first. The scripts numbered `19`–`23` predate the pipeline
consolidation, still use the old paths, and need updating before reuse; the
rest were rewired. `sst_strata_supp/README.md` records what each one showed and
why it was held back.

`histology/` and `gse158516/` are self-contained and run from their own
directories, per their own READMEs.

## Provenance

Everything here was tracked in the main tree until 2026-09-04, when the
codebase was pruned so a reader can find the code behind a paper figure without
walking past a much larger body of work that did not make it. The state before
that pass is the git tag `pre-prune-2026-09-04`.
