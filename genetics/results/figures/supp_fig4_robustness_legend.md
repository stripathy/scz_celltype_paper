# Legend draft for manuscript/figures/supplementary/S10_genetics_ad_robustness.{png,pdf}

> Numbering note: this is **S10** in the manuscript Doc and in
> `manuscript/figures/supplementary/README.md`; **S9** is the enrichment landscape
> (`S09_scz_enrichment_seaad125`). Renumber if the supplement is renumbered before
> submission — the README flags that S4–S14 are due for renumbering.

**Supplementary Fig. S10 | Robustness of the genetic-risk and Alzheimer's disease comparisons to the choice of GWAS and reference region.**
Each point is one of the 16 Sst supertypes; thick outlines mark the five depleted in SCZ (FDR < 0.20 in Fig. 3a); dashed lines, linear fits with 95 % confidence bands; ρ, Spearman correlation. (**a**) SCZ GWAS enrichment per supertype (−log₁₀ *P*) versus depletion in SCZ (−β from Fig. 3a), as in Fig. 4a, with the enrichment recomputed substituting each input in turn: the SCZ GWAS (PGC3, Trubetskoy et al. 2022, in place of Bigdeli et al. 2026) and the expression reference (SEA-AD MTG in place of SEA-AD DLPFC, each region contributing its own supertype set — 125 supertypes in DLPFC, 137 in MTG). Specificity calculation, MAGMA settings and compositional data are held fixed. (**b**) Depletion in SCZ (−β) versus decline along the SEA-AD Alzheimer's disease pseudo-progression score (−β per SD of CPS), as in Fig. 4i, computed in DLPFC (main figure) and recomputed in MTG from the same donors with the same model.

## Not in the legend

Kept out to match the brevity of the rest of the supplement, but useful for the Results text
or a rebuttal:

- Row a, per facet: DLPFC + Bigdeli ρ = 0.55, *P* = 0.031; DLPFC + PGC3 ρ = 0.66,
  *P* = 0.0072; MTG + Bigdeli ρ = 0.42, *P* = 0.11 (a trend of the same sign, in the region
  mismatched to the frontal disease data). Values as drawn in the insets.
- Row b: DLPFC ρ = 0.87, *P* = 1.3 × 10⁻⁵ (inset renders *P* < 10⁻⁴); MTG ρ = 0.79,
  *P* = 4.9 × 10⁻⁴. Fifteen of 16 supertypes agree in direction in both regions, Sst_12 the
  exception; the DLPFC and MTG CPS slopes themselves correlate at ρ = 0.94.
- The panel plots uncorrected −log₁₀ *P*, so the differing DLPFC/MTG correction universes
  (125 vs 137) do not affect what is drawn.

## Provenance

- Tables: `scripts/figures/export_fig4_robustness.py` → `results/tables/fig4_robustness_sst16.csv` (panel a; from the four `results/intermediates/T_{a9only,mtgonly}_{bigdeli,pgc3}.gsa.out` runs joined to `r_panels/panel_B_genetics_vs_depletion.csv`) and `results/tables/fig4_robustness_ad_mtg.csv` (panel b; from `r_panels/panel_ad_concordance_sst.csv` and `crossdisorder/results/crumblr_results_mtg_supertype_neurons.csv`).
- Specificity and enrichment: `scripts/figures/build_spec_seaad_only.py` builds both region-own specificity matrices and runs all four MAGMA gene-property jobs (one-sided, `direction=greater`) off the existing `.genes.raw` files.
- Figure: `scripts/figures/plot_fig4_robustness.R` (10.6 × 8.0 in, cowplot 12 pt; three facets in row a, two in row b, shared geometry and y-axis).
- **The fourth cell of the 2 × 2 is computed but not drawn.** Swapping both inputs at once (MTG + PGC3) gives ρ = 0.12, *P* = 0.65 — available as a column of `fig4_robustness_sst16.csv`. Row a is deliberately a one-at-a-time design, so the both-swapped cell is omitted from the figure; quote it if a referee asks for the full factorial. Under the retired combined taxonomy that same cell was also non-significant (ρ = 0.35, *P* = 0.18), so it is not an artifact of the taxonomy change.
- **Quote the figure's own *P* values.** The insets use R `cor.test` (exact Spearman); scipy's t-approximation of the same data differs in the third decimal (main-figure facet: 0.031 vs 0.028). The earlier "stale 0.0069" scare in the genetics section was this difference, not a stale number.

## Changed from the previous version

- Row a keeps its three facets, but they are now framed as a one-at-a-time sensitivity
  design rather than as "each combination", which the three facets never were. The
  both-swapped cell was computed during the rebuild (see Provenance) and left out by
  decision, not by omission — the earlier version had simply never run it.
- The MTG facet previously reused the retired 501-type combined SEA-AD + Siletti taxonomy
  and its RBH removals, which isolated the reference region. Now that the combined taxonomy
  is retired (L. Duncan's advice, 2026-08-29), each region is computed on its own SEA-AD
  supertype set, so the MTG facet varies reference region and supertype set together. The 16
  Sst supertypes plotted are the same 16 in every facet.
- MTG attenuation is not new to the SEA-AD-only taxonomy — the temporal-lobe reference was
  always the fragile axis, and the whole-brain context partially masked it.
