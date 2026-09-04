# Supplementary figures — submission set

The single place the supplementary figures live. Each generator writes its figure
**directly into this folder** under its **S-number**, so the whole supplement can
be handed to a journal from one folder rather than assembled from six module
directories at deadline.

## Conventions

- **Naming:** `S<nn>_<short_slug>.{png,pdf}` — zero-padded so the folder sorts in
  figure order.
- **Formats:** `.png` (400 dpi, for preview and for pasting into the draft) and
  `.pdf` (vector, editable text, for submission).
- **These are originals, not copies.** Each renderer's output path points here, so
  re-running a script updates the submission figure in place. There is no copy
  step and no second copy to drift. (Until 2026-08-06 these were copies of files
  in `spatial/supplemental_figures/`; that directory is retired and gitignored.)
- **Tracked in git**, unlike the module output directories, so the submission set
  has version history and is available to collaborators on clone.
- **The S-number lives in the renderer**, as a `FIGSTEM` constant next to the
  output directory. Renumbering means editing that one line per script.

## Contents

Numbering follows the manuscript Doc (`1cO5ZSt…`) as of 2026-09-01. **S1, S4, S5
and S7 are Nicole's** (from her snRNA-seq DE and composition analyses); they are
not in this repo and are not in this folder yet.

| S# | File | Shows | Built by | Run from |
|----|------|-------|----------|----------|
| S1 | — | Marker dot plots, integrated + SEA-AD reference | Nicole | — |
| S2 | `S02_xenium_celltype_annotation` | Agreement with Kwon/Lieber annotations (a); subclass and Sst-supertype marker genes (b, d); panel-vs-transcriptome classification F1 (c, e) | `plot_markers_resolvability_combined.R` | `spatial/` |
| S3 | `S03_xenium_merfish_concordance` | Xenium vs MERFISH proportions and depth (a–c); per-supertype depth distributions (d, e) | `plot_xenium_merfish_composite.R` | `spatial/` |
| S4 | — | Supertype-level differential expression | Nicole | — |
| S5 | — | Non-neuronal compositional changes | Nicole | — |
| S6 | `S06_composition_pooling` | Abundance results across pooling strategies and leave-one-dataset-out | `02_plot_pooling_heatmap.R` | repo root |
| S7 | — | Composition robust to re-annotation excluding DE genes | Nicole | — |
| S8 | `S08_sst_strata` | Sst depletion groups: definition (a), gene-set burden (b), gene-level z (c), NES and exemplar-gene heatmaps (d, e). Built as a main figure (was Fig. 5 until 2026-09-01) so it can be promoted | `fig5/08_figure5.R` | repo root |
| S9 | `S09_scz_enrichment_seaad125` | SCZ genetic-risk enrichment across the 125 supertypes of the SEA-AD DLPFC taxonomy (single panel; replaced the 501-type combined-taxonomy version 2026-09-02) | `plot_supp_enrichment_seaad125.R` | repo root |
| S10 | `S10_genetics_ad_robustness` | Genetic-risk and AD comparisons across GWAS version and reference region | `plot_fig4_robustness.R` | repo root |

Regenerate (each writes into this folder in place):

```bash
# from spatial/
Rscript code/analysis/plot_markers_resolvability_combined.R      # S2
Rscript code/analysis/plot_xenium_merfish_composite.R            # S3
# from the repo root
Rscript snrnaseq/composition_sensitivity/code/02_plot_pooling_heatmap.R   # S6
Rscript transcriptomic/scripts/fig5/08_figure5.R                         # S8
Rscript genetics/scripts/figures/plot_supp_enrichment_seaad125.R         # S9
Rscript genetics/scripts/figures/plot_fig4_robustness.R                  # S10
```

Not in this version: the RNAscope figure (cut on Etienne's advice 2026-09-01);
the 501-type enrichment landscape and `S10_supertype_depth_by_diagnosis`, both
cited nowhere in the Doc and moved to `../not_in_current_version/` (their renderers
now write there, so re-running them cannot drop a stray file into this folder);
and the Sst depletion-group controls figure
(nuclei imbalance, cell-matched draws, interaction NES), dropped 2026-09-01 when the
analysis moved to the supplement -- its renderer `fig5/supp/figS_strata_controls.R`
is kept and now writes `reserve_strata_controls.*` to `transcriptomic/results/`.

The small figure-input CSVs are committed (force-added past the `spatial/output/`
ignore rule), so the three renderers run from a clean clone with no external data.
The two ~20 MB per-cell tables behind the S3 violins and the retired depth-by-diagnosis figure are **not** committed —
regenerate them with the `build_*` steps below, which need the Xenium h5ads.

**All three were regenerated from the reinstated 2026-04-01 Xenium dataset** (see
`SCZ_Xenium/archive/README.md`). S2's resolvability panels (c, e) are
dataset-independent in practice — they read only the 300-gene panel list from the
Xenium object and are otherwise leave-one-donor-out CV on the snRNAseq reference.

**Numbering follows order of first citation in SM1**: S2 is cited in the
cell-type classification and resolvability sections, S3 in the validation section. S2 supersedes three earlier figures, which are now its panels: the former
S4 (panel marker genes), S5 (panel resolvability) and S6 (Lieber/Kwon agreement).
Their standalone generators still exist and still work —
`plot_xenium_marker_figure.R`, `plot_panel_resolvability_supplement.R` and
`plot_lieber_celltype_diagonal.py` — but are no longer part of the submission set.
S3 likewise merges the former Xenium-vs-MERFISH scatter figure with the standalone
per-supertype depth figure; its supertype median-depth scatter was dropped, since
panels d and e make the same comparison per supertype with the full distribution.
`plot_celltyping_validation_supplement.R` and `plot_supertype_depth_platform.R`
still build those two as standalones but are no longer part of the submission set.
**Three S-numbers are freed by the two merges; renumber S4–S14 before submission.**

S2 and S3 are both rendered on a 7.1 in canvas.

To regenerate (from `spatial/`):

```bash
# S2  (the lieber data step stages changes; pass --apply to write h5ad columns)
python3 code/analysis/integrate_lieber_annotations.py
python3 code/analysis/plot_sst_supertype_depth_dotplot.py
python3 code/analysis/plot_subclass_marker_dotplot.py
Rscript code/analysis/plot_markers_resolvability_combined.R

# S3
python3 code/analysis/build_celltyping_validation_data.py
python3 code/analysis/build_supertype_depth_platform_data.py
Rscript code/analysis/plot_xenium_merfish_composite.R

# depth by diagnosis -- NOT in this version; writes to ../not_in_current_version/
python3 code/analysis/build_supertype_depth_casecontrol_data.py
Rscript code/analysis/plot_supertype_depth_casecontrol.R
```

## Not in this folder

S1, S4, S5 and S7 are Nicole's figures and live outside the repo (see
`manuscript/SUPPLEMENT_PLAN.md`).

S8 is rendered by `transcriptomic/scripts/fig5/08_figure5.R` straight into this
folder; see `transcriptomic/scripts/fig5/README.md` for the pipeline behind it. To
promote it back to a main figure, change `FIGSTEM` and the output directory at the
bottom of that script.
