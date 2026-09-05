# Built, but not in the current version of the paper

Figures that exist and are regenerable but are cited nowhere in the manuscript.
Kept out of `../supplementary/` so that folder is exactly the submission set.
The depth-by-diagnosis renderer writes here directly, so re-running it cannot
drop a stray file into the submission folder.

| File | Why it is here | Renderer |
|---|---|---|
| `S09_scz_enrichment_501` | SCZ enrichment landscape over the 501-type combined SEA-AD + Siletti taxonomy. Replaced 2026-09-02 by `../supplementary/S09_scz_enrichment_seaad125` after L. Duncan asked that Fig. 4a and its supplement use the SEA-AD taxonomy alone. | Renderer and MAGMA run removed 2026-09-04; recover from tag `pre-prune-2026-09-04` (`plot_supp_enrichment_501.R` + `T_a9rbh_bigdeli.gsa.out`). The rendered figure below is the record. |
| `S10_supertype_depth_by_diagnosis` | Per-supertype cortical depth, control vs SCZ (null control: composition is not a depth artifact). Has no legend and no citation in the Doc; its old S10 slot went to the genetics robustness figure. | `spatial/code/analysis/plot_supertype_depth_casecontrol.R` (run from `spatial/`) |

RNAscope (`histology/`) is also out of this version, cut on Etienne's advice
2026-09-01, to be revisited if reviewers ask. Its figure was never moved into the
submission folder.
