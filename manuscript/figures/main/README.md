# Main figures — submission set

The two main figures **rendered into this folder**. Each renderer writes
directly here, so re-running a script updates the submission figure in place;
there is no copy step and no second copy to drift.

Figures 1 and 3 have renderers in the repo as of 2026-09-09, but they write
elsewhere — see the note below the table.

| Fig | File | Shows | Built by | Run from |
|----|------|-------|----------|----------|
| 1a | — | Cohort UMAPs and cell/donor counts | `Figure_1a_UMAP.r`, `Figure_1a_Barchart.r` | `snrnaseq/Final_figures/` |
| 1b–f | — | Xenium sections, cell-type annotation, laminar segmentation | Assembled by hand; `spatial/` supplies the annotated h5ads it draws on. | — |
| 2 | `Fig2_cross_platform_de` | Cross-platform differential expression: volcanoes, seven-cohort forests for SST and PVALB, per-donor CP1K, exemplar Xenium cells, DE burden by subclass, snRNA-seq vs Xenium concordance | `09_composite_figure.R` | `transcriptomic/` |
| 3 | — | Compositional depletion of upper-layer Sst supertypes: abundance bar chart (a), Sst_25 forest (b) and per-donor proportions (c), spatial rendering (d), abundance vs depth (e), snRNA-seq vs Xenium concordance (f) | `Figure_3.r` | `snrnaseq/Final_figures/` |
| 4 | `Fig4_scz_genetics_hcn1_sst` | SCZ common-variant enrichment vs depletion, Sst_2 gene drivers, the *HCN1* locus, *HCN1* expression vs patch-seq sag, exemplar morphologies and traces, marker volcano, *CALB1*, and the Alzheimer's comparison | `scz_sst_hcn1_story.R` | `genetics/` |

```bash
# Figure 2
cd transcriptomic && Rscript scripts/09_composite_figure.R

# Figure 4
cd genetics && Rscript scripts/figures/scz_sst_hcn1_story.R
```

Both run from committed inputs, so they regenerate from a clean clone without
the raw data. Figure 2's inputs are under `transcriptomic/data/figure_inputs/`
and `transcriptomic/results/tables/`; Figure 4's are the panel CSVs under
`genetics/results/figures/r_panels/`. Both chains carry a `MANIFEST.tsv` that
records the upstream file behind every input, and both renderers halt rather
than draw stale numbers if a source has moved on.

**Formats:** `.png` (400 dpi, preview and pasting into the draft) and `.pdf`
(vector, for submission). Figure 4's renderer also writes an `.svg` for
hand-editing; it is git-ignored because the patch-seq traces make it ~25 MB.

**Figures 1 and 3 are not in this folder, even though their code is in the
repo.** `snrnaseq/Final_figures/` follows the other convention: those scripts
ran on the Alliance cluster against the full per-cohort objects and wrote to a
working directory there, under names that do not carry figure numbers
(`1a.png`, `Figure3_composite.png`). They will not run from a clone. See
[`snrnaseq/Final_figures/README.md`](../../../snrnaseq/Final_figures/README.md)
for the panel map, and `KNOWN_ISSUES.md` issue 12 for the convention gap. What
Figure 1's Xenium panels need from `spatial/` is documented in
`spatial/README.md`.

Supplementary figures are in [`../supplementary/`](../supplementary/); figures
built but not in the current version are in
[`../not_in_current_version/`](../not_in_current_version/).
