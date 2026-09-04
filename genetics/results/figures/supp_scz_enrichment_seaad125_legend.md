# Legend draft for manuscript/figures/supplementary/S09_scz_enrichment_seaad125.{png,pdf}

**Supplementary Fig. S9 | SCZ genetic-risk enrichment across the 125 supertypes of the SEA-AD DLPFC taxonomy.**
The analysis of Fig. 4a extended to every cell type tested: MAGMA gene-property enrichment (one-sided) of SCZ common-variant association (Bigdeli et al. 2026) in the genes specific to each supertype, with expression specificity computed across the 125 SEA-AD DLPFC supertypes (Methods). Each point is one supertype, and all 125 tested supertypes are shown, so the figure and the multiple-testing universe are the same. Supertypes are grouped by subclass (bars beneath the axis, order as in Fig. 3a: GABAergic, then glutamatergic, then non-neuronal) and ordered within subclass by decreasing significance; fill colors are the SEA-AD supertype palette of Fig. 1b. Purple dot-dashed line, Bonferroni *P* = 0.05 over 125 tests (−log₁₀ *P* = 3.40); gray dashed line, FDR = 0.05 (−log₁₀ *P* = 1.74). Supertypes passing Bonferroni are labeled: 18 of 125, all GABAergic interneurons. Of the 51 supertypes at FDR < 0.05, 43 (84%) are GABAergic, eight are glutamatergic and none are non-neuronal. The five Sst supertypes depleted in SCZ (Fig. 3a) are outlined in black; Sst_2, Sst_20 and Sst_3 pass Bonferroni, whereas Sst_25 (*P* = 0.044) and Sst_22 (*P* = 0.055) do not reach FDR < 0.05.

## Provenance

- Table: `scripts/figures/export_supp_enrichment_seaad125.py` → `results/tables/scz_enrichment_seaad125_bigdeli_dlpfc.csv` (from `results/intermediates/T_a9only_bigdeli.gsa.out`, the run behind Fig. 4a; BH FDR and Bonferroni over the 125 tests; subclass/class/palette annotations).
- Figure: `scripts/figures/plot_supp_enrichment_seaad125.R` (7.1 × 3.2 in, cowplot 7 pt; single panel, SEA-AD supertype palette, 24-subclass strip beneath the axis).
- Counts: 18 of 125 supertypes pass Bonferroni — Lamp5 ×5, Pvalb ×5, Sncg ×4, Sst ×4, i.e. every one GABAergic; 51 pass FDR < 0.05 (43 GABAergic, 8 glutamatergic, 0 non-neuronal). Strongest signals: Lamp5_4 (−log₁₀ *P* = 5.59), Sncg_8 (5.58), Sst_2 (5.46).
- Two glutamatergic L6b supertypes reach FDR < 0.05 — L6b_4 (β = 2.96, *P* = 0.0053, FDR = 0.017) and L6b_3 (β = 2.06, *P* = 0.015, FDR = 0.039). L6b_4 is also one of the two supertypes *increased* in SCZ (β = +0.26, FDR = 0.022; Fig. 3a) and is FDR-significant under PGC3 and under the retired combined taxonomy as well; L6b_3 is Bigdeli-specific (PGC3 FDR = 0.28). L6b_1, the strongest compositional increase (β = +0.48, FDR = 2×10⁻⁴), narrowly misses here (*P* = 0.052, FDR = 0.10) but passes under PGC3 (FDR = 0.039). Not currently mentioned in the manuscript text.

## Supersedes

`supp_scz_enrichment_501_legend.md` and the two-panel 501-type figure it described. The
combined SEA-AD + Siletti taxonomy was retired on L. Duncan's advice (2026-08-29); the
Siletti panel is therefore gone and the significance lines are corrected over the 125
supertypes actually tested rather than over 501. The retired figure is kept at
`manuscript/figures/not_in_current_version/S09_scz_enrichment_501.{png,pdf}`.
