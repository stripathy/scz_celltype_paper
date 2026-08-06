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

All generators live in `spatial/code/analysis/` and are run from `spatial/`.

| S# | File | Shows | Built by | Data inputs |
|----|------|-------|----------|-------------|
| S2 | `S02_xenium_celltype_annotation` | Five panels: agreement with the Kwon/Lieber authors' own annotations (a); subclass and Sst-supertype marker genes (b, d); classification F1 from the 300-gene panel vs the full transcriptome (c, e) | `integrate_lieber_annotations.py`, `plot_subclass_marker_dotplot.py`, `plot_sst_supertype_depth_dotplot.py`, `resolvability_report.py` (data) → `plot_markers_resolvability_combined.R` | `output/depth_validation/lieber_layers/` + `output/celltyping_supplement/data/` |
| S3 | `S03_xenium_merfish_concordance` | Five panels: subclass and neuronal-supertype proportions and subclass median depth, Xenium vs MERFISH (a–c); per-supertype depth distributions, glutamatergic (d) and GABAergic (e) | `build_celltyping_validation_data.py` + `build_supertype_depth_platform_data.py` (data) → `plot_xenium_merfish_composite.R` | `output/celltyping_supplement/data/` + `output/depth_platform/` |
| S10 | `S10_supertype_depth_by_diagnosis` | Per-supertype cortical depth, control vs SCZ, with mixed-effects model | `build_supertype_depth_casecontrol_data.py` (data) → `plot_supertype_depth_casecontrol.R` | `output/depth_casecontrol/` |

The small figure-input CSVs are committed (force-added past the `spatial/output/`
ignore rule), so the three renderers run from a clean clone with no external data.
The two ~20 MB per-cell tables behind the S3 and S10 violins are **not** committed —
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

S2 and S3 are both rendered on a 7.1 in canvas (S10 is not, and is due a pass).

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

# S10
python3 code/analysis/build_supertype_depth_casecontrol_data.py
Rscript code/analysis/plot_supertype_depth_casecontrol.R
```

## Still to add

S1, S7–S9, S11–S14. Most already exist in module directories (see
`manuscript/SUPPLEMENT_PLAN.md` for the full inventory and readiness state);
S8 (supertype-level DE power) has not been built yet.
