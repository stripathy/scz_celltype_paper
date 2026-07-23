# Figure 4 (SCZ × HCN1 × Sst) — port provenance

Ported 2026-07-22 from the original `scz_cell_type_enrichment` repo into this
monorepo. Two figures, one R script, one Python extractor.

## Outputs
| File (in `genetics/results/figures/`) | What |
|---|---|
| `scz_sst_hcn1_multipanel_R_v3.pdf` / `.png` | Figure 4 **6-panel** variant — A depletion↔genetics, B HCN1 locus, C HCN1 expr↔genetics (ρ=0.69), D sag↔HCN1, E morphologies, F sag traces |
| `scz_sst_hcn1_multipanel_withdriver.pdf` / `.png` | Figure 4 **7-panel** variant — A full 137-type enrichment landscape, B depletion↔genetics, **C Sst_25 gene-driver scatter (HCN1)**, D HCN1 locus, E sag↔HCN1, F morphologies, G sag traces |
| `supp_scz_enrichment_landscape_R.pdf` / `.png` | Supplementary — SCZ enrichment landscape, neurons-only |

The two variants share one script (`scz_sst_hcn1_story.R`) and one `Rscript` run
produces all three. The 7-panel variant reuses the 6-panel's builders plus
`build_gene_driver_plot()` (panel C) and a full-taxonomy landscape
(`build_enrichment_landscape(neurons_only = FALSE, dense_labels = TRUE)`).
Note the 7-panel variant reuses the v3 depletion panel (**−β / depletion**,
ρ=+0.56) and the v3 sag panel (points sized by **−log10 P_MAGMA**); the original
figure it reproduces used raw β (ρ=−0.56) and n-cells sizing — cosmetic only.

Panel map of the main figure (final lettering after the enrichment-landscape
panel was demoted to the supplement and the Sst_25 gene-driver panel moved to a
*separate* supplement): **A** cell depletion ↔ GWAS enrichment (ρ=0.56); **B**
HCN1 locus zoom (FINEMAP PIP); **C** HCN1 expression ↔ GWAS enrichment (ρ=0.69);
**D** sag ratio ↔ HCN1 expression (ρ=0.65); **E** two exemplar Sst morphologies
(Sst_25·L2, Sst_5·L4); **F** patch-seq sag traces.

## How to regenerate
```bash
# 1. raw sources -> flat panel CSVs (needs the external inputs below to exist)
python3 genetics/scripts/figures/export_for_R.py
# 2. panel CSVs -> rendered figures (R 4.5.1 + tidyverse/cowplot/ggrepel/scales)
Rscript genetics/scripts/figures/scz_sst_hcn1_story.R
```
Step 2 alone reproduces the figures from the **committed** `r_panels/` CSVs and
needs nothing external. Step 1 additionally needs the external raw sources.

## Inputs — committed here vs. referenced in place
The R script reads only panels A/B/D/E/F/G; the gene-driver panel C is generated
by the extractor but is **not** used by either shipped figure.

**Committed into this repo:**
- `genetics/results/figures/r_panels/*.csv` — the 19 extracted panel CSVs (3.1 MB); step 2 uses only these
- `genetics/results/tables/{rbh_combined_enrichment,gwas_vs_casecontrol_composition,sst_supertype_ephys_summary}.csv` (panels A/B/E)
- `genetics/data/seaad_supertype_colors.json`
- `genetics/data/fine_mapping/pgc3_finemap_credible_sets.csv` (panel D, 1.2 MB)
- `genetics/data/gwas/ncbiRefSeq_hg38.txt.gz` (panel D gene track, 7 MB)

**Referenced in place (LARGE / cross-repo — deliberately NOT vendored, ~325 MB):**
| Source | Size | Panel | Location (see extractor config block) |
|---|---|---|---|
| PGC3 SCZ sumstats `.tsv.gz` | 229 MB | D | `ENRICH_REPO/data/gwas/` |
| `seaad_supertype_mean_expression.csv` | 43 MB | E | `ENRICH_REPO/results/intermediates/` |
| patch-seq NWBs (819770858, 758996755) | 53 MB | F/G | `PATCHSEQ_REPO/data/nwb_cache/` |
| SWC morphologies (`*_upright.swc`) | small | F | `PATCHSEQ_REPO/data/morphology/swc/` |
| `rbh_combined_specificity.csv` + MAGMA out + gene.loc | 289 MB | C (unused) | `ENRICH_REPO/results/intermediates/`, `MAGMA_REPO/` |

`ENRICH_REPO = ~/Github/scz_cell_type_enrichment`,
`PATCHSEQ_REPO = ~/Github/human_int_patch_seq` (paths set at the top of
`export_for_R.py`). To make the repo fully self-contained, vendor those files in
and repoint the config block.

## Verification (2026-07-22)
- Extractor ran end-to-end from raw in 7.6 s; all 7 panels emitted.
- Regenerated `r_panels/` CSVs are **byte-identical** (19/19) to the versions
  shipped from the original repo.
- R re-render is a pixel match to the original `scz_sst_hcn1_multipanel_R_v3`.

## Downstream TODO (not part of the port)
The manuscript's current Figure 4 Results + legend describe the *old* matplotlib
panels (depth-gradient r=−0.50, genetics↔depletion r=0.645, a main-figure Sst_25
driver panel). This ported figure supersedes that layout: panels/ρ-values differ
(A ρ=0.56, C ρ=0.69, D ρ=0.65) and the driver panel is demoted. Reconcile the
Fig 4 text and legend to this A–F layout; this likely resolves the
"Fig 4c Sst_25 driver pending regeneration" TODO.
