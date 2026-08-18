# Legend draft for results/figS_rnascope_sst.{png,pdf}

**Supplementary Fig. S7 | SST-positive cell density in subgenual anterior cingulate cortex by RNAscope.**
Re-analysis of RNAscope counts from Arbabi et al. (2025), restricted to donors passing the VIP laminar quality filter (15 control, 11 SCZ; the 13 MDD and 15 bipolar donors enter the model in c; Methods). (**a**) SST-positive cell density per donor (mean over counting frames), control versus SCZ, for all frames and for L2/3 and L5/6 frames separately; P, diagnosis effect from a mixed-effects model of per-frame counts (SST ~ diagnosis + age + sex + PMI + (1 | donor)) fitted to all frames or to each layer. (**b**) Representative lipofuscin-suppressed composites from one control and one SCZ donor in L2/3 (top) and L5/6 (bottom); channel colors and cell markers as keyed on the first image (manually annotated SST+ and VIP+ cells); scale bar, 100 µm. (**c**) Diagnosis coefficients from the same models (difference from control in cells/mm², 95 % CI) for MDD, bipolar disorder and SCZ, for all frames (circles), L2/3 (up triangles) and L5/6 (down triangles); P values are for the all-frames estimates. Lipofuscin suppression is for display only; all counts were made on the original images (Methods).

## Provenance

- Statistics: `code/fit_sst_mixed_models.py` → `results/sst_mixed_models.csv`, `results/sst_subject_density.csv` (statsmodels MixedLM, REML; reproduces β = −0.597 / P = 0.071 overall, −0.762 / 0.126 in L2/3, −0.434 / 0.180 in L5/6).
- Micrographs: `code/extract_micrograph_panels.py` pulls the four cleaned composites out of `results/microscopy/RNAscope_fig_representative_inkscape_SJT.svg`; markers from `results/microscopy/marker_coordinates.csv`.
- Render: `code/plot_figS_rnascope.R` (theme_cowplot 7 pt, Control `#0072B2` / SCZ `#D55E00` as in Fig. 2, 7.1 × 3.95 in, 400 dpi; a over c at left, 2 × 2 micrographs (b) at right with the channel/marker key drawn on the first image).
