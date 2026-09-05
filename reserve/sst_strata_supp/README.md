# Sst depletion-strata supplements (reserve)

Supplements and sensitivity analyses for **Supplementary Fig. S8**, none of
which are in the paper. Run from the repo root; each sources
`transcriptomic/scripts/fig5/_common.R`, so palettes, module definitions and
stratum membership match S8 automatically.

Moved here from `transcriptomic/scripts/fig5/supp/` on 2026-09-04. The headings
below are the state as of 2026-09-01, when the S8 section was shortened: the
"proposed" set was built for a version of the section that is no longer in the
paper, so in current terms every script on this page is reserve.

## Was proposed for the paper

| Script | Figure | Supports |
|---|---|---|
| `module_scores_fig.R` | `figS_module_scores.png` | Each module as one measurement per donor-stratum, all 7 cohorts shown, plus the paired within-donor interaction. Answers "no individual gene is significant" |
| `volcanoes.R` | `figS_volcanoes.png` | Gene level for the three strata and the three interaction contrasts. Discloses that significant genes are few and do not track the gene-set gradient |
| `xenium_crossplatform.R` | `figS_xenium_crossplatform.png` | Shared component replicates on an independent platform/cohort; graded modules are absent from the 300-gene panel |
| `baseline_enrichment.R` → `baseline_enrichment_fig.R` | `baseline_marker_enrichment/` | What distinguishes the depleted supertypes in NEUROTYPICAL tissue. Arguably a Figure 4 supplement |

Requires the core S8 pipeline (`transcriptomic/scripts/fig5/01`–`07`) to have run first.

## Reserve — built, not proposed for the initial submission

| Script(s) | Why held back |
|---|---|
| `20a/20b/20c` Jens independent cohort | 4 of the 5 depleted supertypes map poorly in that dataset, so stratum assignment is porous. Replicates the gradient and OxPhos but not the translation localisation. Needs explicit disclosure; paths still point at the pre-consolidation layout |
| `21a/21b/21c`, `22`, `23` SEA-AD A9 / AD contrast | The AD row was cut from the main figure by decision. Strong rebuttal material if a reviewer asks whether this is generic neurodegeneration — along AD pseudo-progression the same modules decline in the SPARED strata, not the vulnerable ones |
| `19a/19b` baseline SST expression vs depletion | Tests the preproSST-load hypothesis (higher baseline SST in depleted types). Positive trend only, not significant at n = 16 supertypes (rho 0.26–0.47, p >= 0.07). Reportable if asked; too thin to stand alone |

Reserve scripts have NOT been rewired to `_common.R` and still use the
pre-consolidation paths. Update them before reuse.

## Not built

- **Strata-definition sensitivity** (S3): redefine strata at FDR<0.05 (drops
  Sst_20) and as depleted-vs-rest. ~20 min compute per variant (rerun 02 + 03).
  The likeliest reviewer target, since the FDR<0.20 cut is what admits Sst_20.
- Leave-one-cohort-out was considered and dropped: `module_scores_fig.R` already
  shows all seven per-cohort estimates, which answers the same question more
  directly than seven leave-one-out meta-analyses would.
