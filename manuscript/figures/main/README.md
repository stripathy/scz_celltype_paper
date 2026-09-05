# Main figures — submission set

The two main figures this repo builds. Each renderer writes **directly into this
folder**, so re-running a script updates the submission figure in place; there
is no copy step and no second copy to drift.

| Fig | File | Shows | Built by | Run from |
|----|------|-------|----------|----------|
| 1 | — | Xenium sections, cell-type annotation, laminar segmentation | Assembled outside this repo (Nicole). `spatial/` supplies the annotated h5ads it draws on. | — |
| 2 | `Fig2_cross_platform_de` | Cross-platform differential expression: volcanoes, seven-cohort forests for SST and PVALB, per-donor CP1K, exemplar Xenium cells, DE burden by subclass, snRNA-seq vs Xenium concordance | `09_composite_figure.R` | `transcriptomic/` |
| 3 | — | Compositional depletion of upper-layer Sst supertypes | Assembled outside this repo (Nicole). `spatial/output/crumblr/` supplies the Xenium replication. | — |
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

Figures 1 and 3 are not in this repo. What they need from here is documented in
`spatial/README.md`.

Supplementary figures are in [`../supplementary/`](../supplementary/); figures
built but not in the current version are in
[`../not_in_current_version/`](../not_in_current_version/).
