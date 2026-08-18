# Data Directory

This directory contains the raw input data for the sgACC SST-density re-analysis. All data originate from the RNAscope FISH densitometry experiment described in [Arbabi, Newton et al. (2025)](https://doi.org/10.1038/s41380-024-02707-1).

## Primary files used in the analysis

| File | Description |
|------|-------------|
| `Cell_counts_NU.csv` | Per-site cell counts for all subjects. Each row is one subject-section-site combination (e.g., `1234 - R - 7`), with columns for 488nm and 568nm channel counts and field-of-view area. 20 sites per section, two sections (L and R) per subject. |
| `full cell counts(Excel).xlsx` | Authoritative subject-to-diagnosis mapping and Dwight's aggregated cell count summaries (sheet: `full cell counts`). Diagnosis groups: Control, MDD, Bipolar, SCHIZ. |
| `pTable with correct med info.csv` | Subject demographics and clinical covariates: age, sex, PMI, and medication history. Multiple rows per subject (one per cell type); deduplicated during loading. |

## Derived file

| File | Description |
|------|-------------|
| `sst_analysis_data.csv` | Consolidated per-frame dataset built by `code/build_analysis_dataset.py` from the three files above: subject, section, site, layer (L2/3 = sites 1-10, L5/6 = 11-20), SST and VIP counts, diagnosis, age, sex, PMI, and the VIP laminar-QC filter (`vip_ttest_p`, `vip_pass`). 68 subjects, 54 pass the filter. |

## Notes

- Subject IDs in `Cell_counts_NU.csv` contain known typos that are corrected during loading (683 &rarr; 863, 1159 &rarr; 1157, 1449 &rarr; 1444, 1381 &rarr; 1391).
- Channel mapping depends on section: R-section 488nm = SST, 568nm = VIP; L-section 488nm = PYR, 568nm = PV.
- See `code/config.py` for the full list of subject exclusions and ID corrections applied during analysis.
