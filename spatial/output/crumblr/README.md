# Xenium compositional (crumblr) results

## The canonical set — these eight files, nothing else

Four results, one per **{neuronal, non-neuronal} × {subclass, supertype}**, plus the
four inputs they were fitted on. These are the numbers the paper reports.

| Level | Population | Input | Results | Cells | Types in → out |
|---|---|---|---|---|---|
| subclass | neuronal | `crumblr_input_subclass_neuronal.csv` | `crumblr_results_subclass_neuronal.csv` | 356,313 | 18 → 17 |
| supertype | neuronal | `crumblr_input_supertype_neuronal.csv` | `crumblr_results_supertype_neuronal.csv` | 356,313 | 108 → 106 |
| subclass | non-neuronal | `crumblr_input_subclass_nonneuronal.csv` | `crumblr_results_subclass_nonneuronal.csv` | 385,790 | 6 → 6 |
| supertype | non-neuronal | `crumblr_input_supertype_nonneuronal.csv` | `crumblr_results_supertype_nonneuronal.csv` | 385,790 | 29 → 27 |

356,313 + 385,790 = 742,103 — the full analysis set of the canonical
**2026-04-01** Xenium object (see `DATA_FLOW.md`). All 24 donors, 12 Control /
12 SCZ. "Types in → out" differs because crumblr drops types it cannot fit.

Model (`code/analysis/run_crumblr.R`): `~ diagnosis + sex + scale(age_num) + scale(pmi)`
— the same covariates as the snRNAseq composition model and as the Xenium DE.

Built 2026-05-03 from the April object. Verified: the inputs rebuild from that
object with **zero differing rows**, and re-running crumblr reproduces the
committed results to **max |diff| = 1e-15**.

## Sensitivity variants — not in this repository, and not the paper's numbers

The pipeline also produces 65 alternative QC/label-confidence cuts. **None of
them are committed here**; the table below is a key to their filenames on the
machines that ran the pipeline, not an inventory of this directory. `git ls-files
spatial/output/crumblr/` lists everything that is actually present: this README,
the four canonical inputs and their four results.

| Suffix | What it varies |
|---|---|
| *(none)* | all cells pooled, neuronal + non-neuronal together |
| `_corr` | correlation-classifier labels instead of the deployed hybrid |
| `_hybrid` | hybrid QC pass |
| `_margin_{strict,moderate,permissive}` | label-confidence margin threshold |
| `_pctl{01,05,10}` | per-donor transcript-count percentile filters |
| `_no_high_umi` | drops high-UMI (likely doublet) cells |
| `_all*` | subclass and supertype concatenated — redundant with the four above |
| `crumblr_depth_*` | depth-stratified model, a separate question |

> **Only the four canonical files are in this repository** — `crumblr_input_{subclass,supertype}_{neuronal,nonneuronal}.csv`
> and their `crumblr_results_*` counterparts. Every variant suffix in the table
> above (`_corr`, `_hybrid`, `_margin_*`, `_pctl*`, `_no_high_umi`, `_all*`) and
> the `crumblr_depth_*` set exist only on the machines that ran the pipeline;
> they are git-ignored. Rebuild them with `build_crumblr_input.py` and
> `run_crumblr.R` if you need one. Note also that as built on those machines the
> variants exclude Br2039 (23 donors) while the four committed files include all
> 24 — see issue 22 in [`../../../KNOWN_ISSUES.md`](../../../KNOWN_ISSUES.md).

⚠️ **The variants do not all agree, and the disagreement is not obvious from the
filenames.** L6b at subclass level is the clearest trap:

| File | L6b logFC | P | FDR |
|---|---|---|---|
| **`_neuronal` (canonical)** | **+0.330** | **0.176** | **0.470** |
| *(unsuffixed)* | +0.401 | 0.073 | 0.330 |
| `_pctl10` | +0.317 | 0.0017 | 0.042 |

Reading `_pctl10` instead of `_neuronal` turns a null result into a significant
one. The L6b abundance increase is supported at **supertype** resolution
(L6b_2 FDR 0.028, L6b_4 0.028, L6b_5 0.041), not at subclass level.

## Regenerating

```bash
cd spatial
python3 code/analysis/build_crumblr_input.py   # writes inputs + MANIFEST.tsv
Rscript  code/analysis/run_crumblr.R           # writes results
python3 ../shared/verify_provenance.py         # confirm inputs match the object
```

`snrnaseq_vs_xenium_comparison*.csv` are downstream joins of these results
against the snRNAseq composition betas, not crumblr outputs; regenerate them
from the canonical results rather than trusting the copies here.
