# Compositional_sensitivity_analysis/ — is the depletion an artefact of labelling? (S7)

Builds **Supplementary Fig. S7**.

The concern this answers: cells are assigned to supertypes by transferring
labels from a reference, and the same genes that are differentially expressed in
schizophrenia also help drive that assignment. If SCZ changes Sst expression,
Sst cells in cases might simply be *labelled* differently — producing an
apparent loss of a supertype with no cells lost at all.

The test: withhold the Sst DE genes from the reference, re-transfer labels,
and re-run the whole composition analysis. If the depletion survives, it is not
a labelling artefact.

## Chain

```
1_Label_transfer_noSSTDE.r -> per-dataset counts, re-labelled  ({dataset}_*SST_DEgenes.csv)
2_Compile_metadata.r       -> 7_cohorts_metadata_NoSST_DEgenes_names.csv
3_Crumblr_analysis.r       -> crumblr_results_noSSTgenes.csv  (per dataset)
4_Meta_analysis.r          -> meta_noSST_genes.csv            (pooled)
```

Steps 2–4 mirror `../Compositional_analysis/` exactly; only step 1 differs.

**The step-1 outputs are not tracked.** The 13 per-dataset
`{dataset}_*SST_DEgenes.csv` files are per-cell dumps running to 2.8 GB
together; they are git-ignored, and regenerating them means re-running step 1
against the full Seurat objects. What *is* committed is the compiled table step
2 produces, `Files/7_cohorts_metadata_NoSST_DEgenes_names.csv` (0.18 MB), which
is what `3_Crumblr_analysis.r` reads — so steps 3–4 run from a clone.

Step 1 writes only the columns step 2 consumes (`Donor, Age, Sex, Diagnosis,
PMI, predicted.id`, plus Multiome's three spellings) rather than the whole
`meta.data` slot. Keep it that way: what an upstream object carries is not
under this repository's control, and one cohort's arrived with the donating
brain banks' full clinical record attached.

### 1 — Re-transfer with the Sst DE genes removed

Reads the subclass meta-DE table from `../snRNAseq_DE/Subclass/`, takes every
gene with **`cell_type == "Sst"` and `padj < 0.1`**, and excludes that gene set
from the feature space before anchor-finding. Otherwise the transfer is the same
as [`../Label_transfer/`](../Label_transfer/README.md) — same reference, same
`TransferData(refdata = ref$Supertype)`.

It processes eleven objects rather than seven, because the PsychAD datasets are
held as multiple files: `MSSM2` arrives as four parts and `HBCC` as two, and
they are recombined into their seven dataset labels in step 2.

### 3 — Crumblr

Same neuronal restriction (via `cluster_order_and_colors.csv`) and the same
per-dataset `dream` fit, including the same fixed-effect `Sex` term.

> `Sex` is modelled exactly as in the primary analysis, so S7 differs from
> Fig. 3a only in the thing it is meant to test. In the current fit Sst_25 sits
> at FDR 0.1122, above 0.10; four of the five depleted Sst supertypes clear
> FDR < 0.10 here. The Results text reflects this.

### 4 — Meta-analysis

`metafor::rma(method = "FE")` across datasets, ≥ 2 datasets required, BH FDR.
Same as the primary chain.

Rendered by `../Final_figures/Supplemental/Sensitivity_analysis_Barchart.r`,
which draws the same bar chart as Fig. 3a so the two can be read side by side.
