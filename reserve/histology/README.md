# sgACC SST interneuron density (RNAscope re-analysis)

Re-analysis of RNAscope cell-count data from human subgenual anterior cingulate cortex (sgACC; Arbabi, Newton et al., *Molecular Psychiatry* 2025) asking whether the density of SST-positive interneurons is lower in schizophrenia, and whether any reduction is larger in upper (L2/3) than deep (L5/6) layers. Feeds one supplementary figure and one supplementary table of the SCZ cell-type paper.

## Background

The original study sampled sgACC sections from 76 post-mortem subjects (Control, MDD, Bipolar, Schizophrenia tetrads) across 20 stereological counting frames per section, sites 1–10 in L2/3 and 11–20 in L5/6, and reported a non-significant trend toward lower SST density in schizophrenia when averaging over layers. Our snRNA-seq meta-analysis predicts that Sst depletion is concentrated in upper layers, so we re-analyzed the counts layer-stratified.

## Pipeline (run from this directory, in order)

| Step | Script | Output |
|---|---|---|
| 1. Build the analysis dataset | `code/build_analysis_dataset.py --write` | `data/sst_analysis_data.csv` |
| 2. Fit the mixed models | `code/fit_sst_mixed_models.py` | `results/sst_mixed_models.csv`, `results/sst_subject_density.csv` |
| 3. Extract micrograph panels | `code/extract_micrograph_panels.py` | `results/microscopy/panels/*.png` |
| 4. Render the figure | `Rscript code/plot_figS_rnascope.R` | `results/figS_rnascope_sst.{png,pdf}` |

Step 1 parses the raw count sheets (`code/config.py`, `code/data_loading.py`: subject-ID corrections, exclusions, diagnosis and demographics merge), assigns binary layers, and applies the **VIP laminar quality filter**: VIP interneurons are concentrated in superficial layers, so for each subject we test VIP(L2/3) > VIP(L5/6) (one-sided t-test); subjects with P ≥ 0.05 have layer annotations that do not reflect laminar biology and are excluded (`vip_pass`), retaining 54 of 68 subjects (15 Control, 11 SCZ, 15 Bipolar, 13 MDD). Without `--write` the script only verifies that it reproduces the tracked CSV and that per-subject means match the original study's summary counts (r = 1.000).

Step 2 fits, on VIP-filtered per-frame counts,

```
SST ~ diagnosis + age + sex + PMI + (1 | subject)
```

once on all frames and once per layer (statsmodels MixedLM, REML; Control as reference). **These are the numbers quoted in the manuscript**: SCZ β = −0.597 per frame (−5.4 cells/mm²), P = 0.071 overall; −0.762, P = 0.126 in L2/3; −0.434, P = 0.180 in L5/6; Bipolar −0.379, P = 0.214; MDD +0.075, P = 0.813. Frame area = 0.110889 mm².

Steps 3–4 build the figure in the main-figure conventions (7-pt cowplot theme, Control blue / SCZ orange as in Fig. 2): (a) per-subject SST⁺ density by layer, (b) representative lipofuscin-suppressed composites with the manually annotated cells, (c) diagnosis coefficients by layer. Legend draft and provenance: `results/figS_rnascope_sst_legend.md`.

Supplementary-table material: `results/rnascope_subjects_final.csv` (the 54 analysed subjects with demographics and per-subject densities) and `results/rnascope_demographics_summary.csv`.

## Data

| File | Description |
|---|---|
| `data/Cell_counts_NU.csv` | Per-site 488/568-nm cell counts (20 sites × 2 sections per subject). R section: 488 = SST, 568 = VIP. |
| `data/full cell counts(Excel).xlsx` | Subject-to-diagnosis mapping and the original study's per-subject summary counts. |
| `data/pTable with correct med info.csv` | Demographics (age, sex, PMI). |
| `data/sst_analysis_data.csv` | Consolidated per-frame dataset built by step 1 (68 subjects, `vip_pass` column). |
| `results/microscopy/` | Representative-image processing (lipofuscin suppression, curated annotation SVG, marker coordinates); see its README. |

## Requirements

Python 3: pandas, numpy, scipy, statsmodels, openpyxl, Pillow. R: readr, dplyr, tidyr, ggplot2, cowplot, ggsignif, png.

## Citation

> Arbabi K, Newton DF, Oh H, Davie MC, Lewis DA, Wainberg M, Tripathy SJ, Sibille E. Transcriptomic pathology of neocortical microcircuit cell types across psychiatric disorders. *Molecular Psychiatry*. 2025;30:1057-1068. doi:10.1038/s41380-024-02707-1
