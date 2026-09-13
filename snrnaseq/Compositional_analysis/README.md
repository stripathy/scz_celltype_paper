# Compositional_analysis/ — cell-type abundance in SCZ (Fig. 3a, S5)

Stage 2a. Asks whether each cell type makes up a different fraction of the
neuronal population in schizophrenia, per dataset, then pools the seven datasets.
This is the analysis behind **Figure 3a** and the abundance axis of **Figure 4a
and 4i** in `genetics/`.

Compositional data, so the model is [crumblr](https://github.com/GabrielHoffman/crumblr)
(counts → CLR with precision weights) fitted by `dreamlet::dream`, one fit per
dataset, then a fixed-effect meta-analysis across datasets.

## Chain

```
1_bind_metadata.r    -> 7_cohorts_metadata_names.csv        [committed here]
2_Crumblr_analysis.r -> crumblr_results_final_one_doc.csv   (per dataset)
3_meta_analysis.r    -> crumblr_results_final_meta.csv      (pooled)
```

### 1 — `1_bind_metadata.r`

Loads the seven labelled objects from `../Label_transfer/`, applies the donor
filters (age < 70; age > 20 as well for HBCC; Multiome restricted to
`Disorder ∈ {control, Schizophrenia}`), collapses cell-level metadata to one row
per donor, cross-tabulates `Donor × predicted.id` into per-donor cell counts,
and standardises the supertype names (`.SEAAD` → `-SEAAD`, `L2.3` → `L2/3`,
`Sst.Chodl` → `Sst Chodl`, remaining `.` → space) so they match the names used
by `genetics/` and `spatial/`.

Output: **`7_cohorts_metadata_names.csv`** — one row per donor (469), columns
`Cohort, Donor, Age, Sex, Diagnosis, PMI` followed by one count column per
supertype. **This file is committed**, so the model in step 2 can be read
against its actual input.

### 2 — `Code/Neurons/1_Crumblr_analysis.r`

Per dataset: restrict the count columns to **neuronal** supertypes, fit

```r
crumblr(counts) |> dream(~ scale(Age) + scale(PMI) + Diagnosis + (Sex), meta) |> eBayes()
```

and take `topTable(coef = "DiagnosisSchizophrenia")`. Multiome has no PMI, so it
drops that term.

"Neuronal" is defined by an **external** file,
`cluster_order_and_colors.csv`, via `class_label ∈ {Neuronal: GABAergic,
Neuronal: Glutamatergic}`. That file is not committed — see issue 4 in
[`../../KNOWN_ISSUES.md`](../../KNOWN_ISSUES.md).

> The `(Sex)` term is a plain fixed effect. The two supporting analyses
> (`Non_neurons/`, `../Compositional_sensitivity_analysis/`) use `(1|Sex)`
> instead. Logged as issue 1 in [`../../KNOWN_ISSUES.md`](../../KNOWN_ISSUES.md).

### 3 — `3_meta_analysis.r`

Per cell type, pools the seven dataset estimates with
`metafor::rma(yi = logFC, sei = SE, method = "FE")`, requiring ≥ 2 datasets.
Standard errors are reconstructed as `|logFC / t|` from the crumblr output.
Reports `estimate, se, zval, pval, I2, ci.lb, ci.ub, k`, with Benjamini-Hochberg
FDR across cell types.

## `Non_neurons/` — Supplementary Fig. S5

`1_crumblr_meta_nonNeurons` runs the same fit-then-pool on the **complement**
set (everything not `Neuronal: GABAergic` / `Neuronal: Glutamatergic`), reading
the same committed counts file. Both steps are in the one script.

Two differences from the primary analysis, both worth knowing: `Sex` enters as
`(1|Sex)`, and the file has **no extension** (it is R). Issues 1 and 9 in
[`../../KNOWN_ISSUES.md`](../../KNOWN_ISSUES.md).

Rendered by `../Final_figures/Supplemental/NonNeuron_supplement.r`.
