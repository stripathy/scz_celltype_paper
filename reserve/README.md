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

## One analysis is held back from this repository

A polygenic-score analysis (`reserve/genotype/`) exists on the authors' machines
and is **deliberately not published here**. Its donor table carries, per donor,
genetic ancestry principal components, a superpopulation label, polygenic risk
scores and AD neuropathology scores. That is individual-level genotype-derived
data — a different category from the demographic donor sheet cleared in
[`../KNOWN_ISSUES.md`](../KNOWN_ISSUES.md) issue 6 — and the upstream genotypes
are governed by data use agreements that restrict redistribution of derived
individual-level data.

It is listed in `.gitignore` so it cannot be committed by accident. Un-ignoring
it needs a data-use-agreement check and explicit PI sign-off on that specific
content, not the demographic clearance. The manuscript makes no polygenic-score
claim, so nothing in the paper depends on it.

## The data is not on this machine

The large inputs behind three of these are not held locally, to save disk: the GSE158516 count matrices, the Jens cohort object, and the
SEA-AD 2026-06-22 DFC nuclei file that the AD contrast reads. All three are
public downloads — see the READMEs below and `crossdisorder/data/SEAAD_INPUTS.md`
for the URLs. The committed results, manifests and figures are unaffected; only
re-running these analyses needs a re-download.

## Running any of it

The strata supplements are run **from the repo root** and source the live
pipeline's `_common.R`, so they inherit the same palettes, module definitions
and stratum membership as Supplementary Fig. S8:

```bash
Rscript reserve/sst_strata_supp/volcanoes.R
```

They require the core S8 pipeline (`transcriptomic/scripts/sst_strata/01`–`07`) to
have run first. The scripts numbered `19`–`23` predate the pipeline
consolidation, still use the old paths, and need updating before reuse; the
rest were rewired. `sst_strata_supp/README.md` records what each one showed and
why it was held back.

`histology/` and `gse158516/` are self-contained and run from their own
directories, per their own READMEs.

## Why these are here and not in the main tree

Keeping them under `reserve/` means a reader looking for the code behind a paper
figure does not have to walk past analyses that did not make it in. Nothing
downstream of the paper reads this directory.
