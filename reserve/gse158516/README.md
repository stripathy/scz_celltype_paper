# GSE158516 — an 8th SCZ cohort, reprocessed and parked

> **STATUS: parked, not in the paper. The count matrices were removed from
> this machine on 2026-09-04; re-download from GEO to re-run.** This cohort is **not** in the meta-analysis and is
> not cited in the manuscript. It was reprocessed (2026-08-14) to answer one question —
> *would an independent 8th dataset agree with our findings?* — and the answer was yes.
> It is kept as a reviewer-response asset and as a ready-made input if we later decide to
> include it.
>
> **The cell-type labels here are PROVISIONAL and must not be used.** They come from a
> lightweight brisc kNN label transfer, not the Seurat reference-based pipeline the
> paper's seven datasets went through. To use this cohort for anything real, re-annotate
> `GSE158516_counts_qc.h5ad` with the canonical pipeline. The provisional labels are
> deliberately stored *outside* the h5ad so they cannot be picked up by accident.

Reprocessing of **Reiner et al.**, *Single-nuclei transcriptomics of schizophrenia
prefrontal cortex primarily implicates neuronal subtypes* (bioRxiv 2020.07.29.227355;
Reiner, Crist, Stein, Weller, Doyle, Arauco-Shapiro & Berrettini, UPenn), deposited as
[GSE158516](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158516). GEO lists no
linked publication, so it appears still to be a preprint.

**Not one of the seven** (Batiuk, HBCC, MSSM 1, MSSM 2, McLean, Multiome, Fröhlich), so
it is a genuinely independent cohort.

Every number below is extracted programmatically from `output/`; `code/13_verify_readme.py`
re-derives 67 of them and passes.

---

## Read this before using the cohort

Four properties constrain what this dataset can be used for. None are obvious from the
GEO record.

1. **Only 26 of the 32 GEO samples are usable.** CT_03, CT_08, SZ_02, SZ_12, SZ_15 and
   SZ_16 are absent from the authors' supplement and have no age or PMI, so they cannot
   enter a covariate-adjusted model. Two of them (CT_03, SZ_15) also fail QC badly (~60%
   of nuclei removed), which is presumably why the authors dropped them. All 32 are built
   and shipped; `obs['in_authors_analysis']` flags the 26.
2. **Every donor is male** (26/26). Sex cannot be a covariate, and this cohort cannot
   contribute to any sex-stratified analysis. This is the sharpest limitation.
3. **Controls are 12.8 years older than cases** (54.4 vs 41.7, Welch *p* = 0.064); PMI runs
   the other way (22.0 vs 29.1 h, *p* = 0.12). Applying the paper's **≤ 70 years** rule
   drops 5 donors, *all controls* (ages 71, 73, 76, 80, 87), which balances age exactly but
   leaves only **9 controls vs 12 cases**.

   | donor set | n (CT/SCZ) | meeting ≤ 70 | CT age | SCZ age | Welch *p* |
   |---|---|---|---|---|---|
   | all | 14 / 12 | CT 9/14 = 64.3%, SCZ 12/12 = 100% | 54.4 | 41.7 | 0.064 |
   | ≤ 70 | 9 / 12 | — | 41.7 | 41.7 | 1.000 |

4. **Cell calling was ours to do.** GEO ships CellRanger **raw** matrices (33,538 features
   × 6,794,880 barcodes per sample), not filtered ones. Keep the rule established below
   rather than deriving a new one.

---

## The deliverable

`~/Github/shared_data/GSE158516/`

| file | what |
|---|---|
| **`GSE158516_counts_qc.h5ad`** | **the artifact.** 294,974 nuclei × 33,538 genes, 8.0 GB. Raw counts, donor metadata, QC metrics, `passed_QC`, `in_authors_analysis`. **No cell-type labels** — this is a label-transfer input. |
| `GSE158516_provisional_labels.parquet` | provisional labels + confidences + UMAP, per cell. **Not for analysis.** |
| `GSE158516_provisional_harmony.npy` | provisional Harmony embedding (cell order matches the h5ad) |
| `raw/` | the 9.9 GB GEO download; byte-exact reproducible via `download.sh` + `filelist.txt` |
| `barcode_totals/` | per-barcode UMI totals (pass 1 of cell calling); basis of the knee plots |
| `per_sample_qc/` | per-sample QC'd h5ads, the input to the concatenation |
| `pseudobulk_subclass/`, `de/` | provisional pseudobulk and DE outputs |

`per_sample/` and `reference/` were deleted after packaging — both regenerate in ~10 min
(`01c_build_h5ad.py`, `02b_subsample_reference.py`).

---

## How it was QC'd

1. **Cell calling** (`01a`, `01b`, `01c`). Barcode-rank curves cached and three rules
   compared (`output/knee_plots.png`, `cell_calling_comparison.csv`). CellRanger-v2 ordmag
   recovers only 0.56× the authors' nuclei. **Top-N by UMI, with N = the authors' reported
   per-sample nuclei count, reproduces their cell set almost exactly**: our recovered
   median UMI/nucleus matches their published Supplementary Table 2 value to < 1% in 23/26
   samples (max 3.2%, median ratio 1.0001), and the total, 361,681 nuclei, equals the sum
   of their reported counts. **Keep this rule.**
2. **QC** (`03`, brisc). `max_mito_fraction` 5%, `min_genes` 500, `nonzero_MALAT1`, cxds
   doublet removal, per sample. 28.5% of nuclei removed on average → **259,187 nuclei
   across the 26** (the preprint reports ~275,000), 294,974 across all 32.
   **QC attrition does not track diagnosis**: 28.2% control vs 28.9% SCZ (*p* = 0.80);
   mito failure *p* = 0.71.
3. **Provisional annotation** (`02b`, `04`, brisc) — *to be redone*. Joint HVG → normalize
   → PCA → Harmony against the SEA-AD MTG reference, then kNN label transfer. Subclass
   confidence > 0.8 for 93.9% of nuclei; supertype median 0.75. De novo Leiden clusters are
   97.9% pure by transferred subclass, markers form a clean diagonal
   (`output/fig_markers.png`), and Harmony mixes samples and diagnoses
   (`output/fig_annotation.png`).

**One trap worth recording:** when subsampling a reference for kNN transfer, sample
**proportionally**. Capping at *N* cells per supertype flattens the reference's class
priors, which kNN voting leaks into the calls — a 400-cell cap inflated GABAergic from
32.7% to 60.6% of reference neurons and pushed the query to 37.3% instead of 32.3%.
Subclass calls were 95.6% stable between the two schemes; supertype only 67.7%.

---

## What we learned (not in the paper)

Run on the provisional labels, so directional rather than definitive — but the pipeline
itself is validated: against Reiner et al.'s own published DE (their Supplementary Table
4, a MAST per-nucleus model on their own 20-cluster taxonomy) our reprocessing is
**95.6% sign-concordant** (4,647/4,863 gene × cluster pairs), and **29/29** on their SST+
interneuron genes (Spearman ρ = 0.81). So disagreements would be biology, not handling.

| paper's claim | GSE158516 |
|---|---|
| *SST* down in Sst cells (meta log₂FC = −0.458, FDR = 0.049) | −0.175, *p* = 0.61 — same direction, underpowered |
| *PVALB* **not** DE in Pvalb cells (meta −0.056) | −0.021, *p* = 0.94 — reproduced |
| DE effects generalize | *r* = 0.376, 68.6% sign-concordant over 6,864 meta-significant pairs |
| named Sst/Pvalb/Chandelier genes | 15/16 directionally concordant |
| upper-layer Sst depleted, L6b increased | **8/8** headline supertypes concordant |
| composition generalizes (Xenium: neuronal subclass *r* = 0.70) | neuronal subclass *r* = **0.702**; supertype *r* = 0.367 |

Adding it to the composition meta-analysis (`14`, `15`; **not done in the paper**) would
give it ~7.5% of the inverse-variance weight, move β by a median of 0.012, and take
FDR<0.10 from 6 to 8 supertypes — Sst_20 from FDR 0.19 to 0.065, Sst_3 the other way
(0.046 → 0.065). The two supertypes that newly appear (L5/6 NP_4, Sst_11) **fail to
replicate in Xenium** and are poorly resolved on the 300-gene panel (F1 at the 31st and
28th percentile), so nothing new here survives orthogonal checking.
See `output/status_changes.csv`.

---

## Findings independent of this cohort — worth acting on

Three things surfaced from *paper* data while doing this. None depend on GSE158516.

1. **The published compositional pool is fixed-effect, not the REML the Methods
   describe.** FE reproduces the published rows to 8.9 × 10⁻¹⁶; REML/DL/ML/HE/SJ/EB/PM are
   all off by ~2 × 10⁻². Corroborating: the composition table has no `tau2`/`I2`/`k`
   columns while the DE table does. Under REML two of the five vulnerable Sst types weaken
   materially (BH within the 16 SST types: Sst_3 0.010 → 0.060, Sst_20 0.062 → 0.21).
   Heterogeneity is low (I² ≤ 41%), so FE is defensible — but text and analysis must be
   made to agree. Evidence: `code/16_validate_pooling_method.R`,
   `output/pooling_method_validation.csv`. **Needs a decision.**
2. **Sst depletion vs cortical depth is directly quantifiable**: Spearman ρ = **+0.72**
   (*p* = 0.0016) across the 16 Sst supertypes, on the *existing seven datasets*. The
   manuscript asserts the depletion is upper-layer-concentrated but demonstrates it only
   transitively (composition↔genetics, genetics↔depth). One sentence would make the
   central claim direct.
3. **Cross-platform replication tracks panel resolvability**: well-resolved supertypes
   agree between snRNA-seq and Xenium at *r* = 0.60 (n = 53) versus *r* = 0.41 for
   poorly-resolved ones (n = 52) — a quantitative version of the caveat the Limitations
   state in words. The per-type F1 difference between concordant and discordant types is
   itself not significant (MWU *p* = 0.33), so treat it as graded, not a threshold.

---

## If we decide to include this cohort

1. Re-annotate `GSE158516_counts_qc.h5ad` with the canonical Seurat label transfer;
   discard `*_provisional_*`.
2. Decide the donor set: 26 (age covaried) or 21 (≤ 70, age balanced, 9 controls).
3. Re-run crumblr with `~ diagnosis + age + PMI` — **sex must be dropped**, it is constant.
4. Pool with the other seven. Inverse-variance weighting is associative, so a fixed-effect
   pool can be combined with the published 7-dataset summary directly; a random-effects
   pool would need the per-cohort estimates (available for the 16 SST supertypes only, in
   `~/Downloads/SCZ_SST_7_cohort_estimates (3).csv`).

---

## Files

| file | what |
|---|---|
| `code/01a`–`01c` | barcode totals → cell-call comparison → per-sample h5ad |
| `code/02b` | SEA-AD reference subsample (proportional / capped) |
| `code/03`–`05b` | QC, provisional label transfer, pseudobulk + DE |
| `code/06` | crumblr composition (optional counts prefix) |
| `code/07`, `09` | comparison vs our findings; reproduction of Reiner's own DE |
| `code/08`, `11`, `15`, `17` | QC/annotation, results, meta-8 and status-change figures |
| `code/10`, `18` | concatenate, then package as a label-transfer input |
| `code/12` | ≤ 70 y sensitivity analysis |
| `code/13` | re-derives every number in this README (67 checks) |
| `code/14`, `16` | 8-dataset composition update; pooling-method validation |
| `code/brisc_setup.py` | version-aware brisc shim (see `~/Github/brisc_test/BRISC_BUGS.md`) |

## brisc

Run on brisc 0.1.0 and re-verified on 0.1.3; the DE table is bit-identical between them
(max |Δ logFC| = 0 over 377,011 rows). Two integer-handling defects found in
`Pseudobulk._create_design_matrix` are fixed and submitted upstream as
[briscverse/brisc#2](https://github.com/briscverse/brisc/pull/2). `brisc_setup.py` applies
the 0.1.0 workarounds only when the running version needs them.
