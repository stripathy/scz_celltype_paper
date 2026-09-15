# Built, but not in the current version of the paper

Figures that exist and are regenerable but are cited nowhere in the manuscript.
Kept out of `../supplementary/` so that folder is exactly the submission set.
The depth-by-diagnosis renderer writes here directly, so re-running it cannot
drop a stray file into the submission folder.

| File | Why it is here | Renderer |
|---|---|---|
| `S09_scz_enrichment_501` | SCZ enrichment landscape over the combined SEA-AD + Siletti taxonomy (501 types). Superseded by `../supplementary/S09_scz_enrichment_seaad125`, which puts Fig. 4a and its supplement on the SEA-AD taxonomy alone. | No renderer in the tree; the rendered figure here is the record. |
| `S10_supertype_depth_by_diagnosis` | Per-supertype cortical depth, control vs SCZ (null control: composition is not a depth artifact). Has no legend and no citation in the Doc; its old S10 slot went to the genetics robustness figure. | `spatial/code/analysis/plot_supertype_depth_casecontrol.R` (run from `spatial/`) |

RNAscope (`histology/`) is also out of this version, to be revisited if
reviewers ask. Its figure was never moved into the submission folder.
