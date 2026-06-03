# sgACC SST Interneuron Density Analysis

Analysis of RNAscope cell density data from human subgenual anterior cingulate cortex (sgACC), testing whether SST interneuron density is reduced in schizophrenia using layer-stratified mixed-effects models.

## Background

Arbabi, Newton et al. (*Molecular Psychiatry*, 2025) performed RNAscope fluorescent in situ hybridization on sgACC sections from 76 postmortem subjects organized as 19 matched tetrads (Control, MDD, Bipolar, Schizophrenia). Each section was sampled across 20 stereological counting frames assigned to upper (L2/3, sites 1–10) or deep (L5/6, sites 11–20) cortical layers. The original study reported a trend toward reduced SST density in schizophrenia that did not reach significance when averaging across all sites.

Our snRNAseq meta-analysis predicts that SST reductions in schizophrenia should be concentrated in upper cortical layers. If so, averaging across layers would dilute a real layer-specific effect. This repository re-analyzes the original data to test that prediction.

## Approach

### VIP-based section quality filter

VIP interneurons are concentrated in superficial cortical layers, providing an internal control for laminar structure. For each subject, we test whether VIP counts are significantly higher in L2/3 than L5/6 (one-sided t-test). Subjects where this expected laminar pattern is absent (p ≥ 0.05) are excluded from the primary analysis, as their layer annotations may not reflect reliable laminar biology. This retains 54 of 68 subjects.

### Statistical models

We fit mixed-effects models with subject as a random intercept, controlling for age, sex, and postmortem interval:

```
SST ~ diagnosis + age + sex + PMI + (1 | subject)
```

Models are fit on all 20 ROIs per subject (not collapsed to subject-level means), giving ~1,080 observations for the aggregate analysis and ~540 for each layer-stratified analysis.

## Data

| File | Description |
|------|-------------|
| `data/Cell_counts_NU.csv` | Per-site cell counts (20 sites × 2 sections × 76 subjects) |
| `data/full cell counts(Excel).xlsx` | Subject-to-diagnosis mapping |
| `data/pTable with correct med info.csv` | Subject demographics (age, sex, PMI) |
| `data/sst_analysis_data.csv` | Consolidated analysis dataset with `vip_pass` filter column |

### Channel mapping

| Section | 488nm | 568nm |
|---------|-------|-------|
| R | SST | VIP |
| L | PYR (SLC17A7) | PV (PVALB) |

## Results

### VIP-filtered subjects (N = 54)

| Analysis | SCHIZ β | p |
|----------|---------|---|
| Aggregate | −0.60 | 0.071 |
| L2/3 stratified | −0.76 | 0.126 |
| L5/6 stratified | −0.43 | 0.180 |

Bipolar disorder shows intermediate effects in the same direction (β = −0.38, p = 0.21). MDD shows no SST reduction (β = +0.08, p = 0.81).

## Repository structure

```
code/
  config.py             # Paths, exclusions, channel mappings
  data_loading.py       # Data parsing and merging
  depth.py              # Coordinate extraction (used for supplementary analyses)
  models.py             # Mixed-effects model fitting
  plotting.py           # Visualization functions
  run_analysis.py       # Main analysis pipeline
  subject_tracking.py   # Subject inclusion/exclusion tracking
data/
  sst_analysis_data.csv     # Consolidated dataset (68 subjects, vip_pass column)
  Cell_counts_NU.csv        # Raw per-site cell counts
results/
  plot_main_results.Rmd     # Main figure generation (R Markdown)
  fig_main_result.R         # Standalone main figure script
  fig_main_result.png       # Main 4-panel figure
coordinates/                # Stereology Excel files (site XY positions)
archive/                    # Prior exploratory and continuous depth analyses
```

## Requirements

**Python 3**: numpy, pandas, scipy, statsmodels, matplotlib, openpyxl

**R**: tidyverse, cowplot, ggsignif, patchwork

## Author

Shreejoy J. Tripathy — [Tripathy Lab](https://triplab.org), Krembil Centre for Neuroinformatics, CAMH & University of Toronto.

## Citation

> Arbabi K, Newton DF, Oh H, Davie MC, Lewis DA, Wainberg M, Tripathy SJ, Sibille E. Transcriptomic pathology of neocortical microcircuit cell types across psychiatric disorders. *Molecular Psychiatry*. 2025;30:1057-1068. doi:10.1038/s41380-024-02707-1
