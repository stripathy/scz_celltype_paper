# Legend draft for results/figures/supp_scz_enrichment_501.{png,pdf}

**Supplementary Fig. SX | SCZ genetic-risk enrichment across the 501 cell types of the combined taxonomy.**
The analysis of Fig. 4a extended to every type: MAGMA gene-property enrichment of SCZ common-variant association (Bigdeli et al. 2026) in the genes specific to each type, with expression specificity recomputed from the SEA-AD DLPFC reference for the SEA-AD supertypes and from the Siletti et al. 2023 whole-brain atlas for the remaining clusters (Methods). Purple dot-dashed line, Bonferroni P = 0.05 over 501 tests; gray dashed line, FDR = 0.05. Each point is one cell type. (**a**) The 125 SEA-AD supertypes retained in the combined taxonomy, grouped by subclass (order as in Fig. 3a) and colored by supertype as in Fig. 1b; supertypes passing Bonferroni are labeled; the five Sst supertypes depleted in SCZ (Fig. 3a) are outlined in black. (**b**) The 376 Siletti clusters without a SEA-AD counterpart, grouped by Siletti supercluster (alternating gray; superclusters ordered by their strongest member); clusters with −log₁₀ P ≥ 8 are labeled.

## Provenance

- Table: `scripts/figures/export_supp_enrichment_501.py` → `results/tables/scz_enrichment_501_bigdeli_dlpfc.csv` (from `results/intermediates/T_a9rbh_bigdeli.gsa.out`, the run behind Fig. 4a; BH FDR and Bonferroni over the 501 tests; source/group/palette annotations).
- Figure: `scripts/figures/plot_supp_enrichment_501.R` (7.1 × 6.3 in, cowplot 7 pt; points in both rows, SEA-AD palette in a, alternating grays by supercluster in b).
- Counts: 81 of 501 types pass Bonferroni (18 SEA-AD, 63 Siletti); 244 pass FDR < 0.05. Top SEA-AD supertype is Sst_23 (−log₁₀ P = 7.2, rank 31 overall); the strongest signals overall are Siletti clusters (Misc_132 20.8, ULIT_133 15.1, then Eccentric-MSN, MGE, LAMP5-LHX6/Chandelier and deep-layer IT clusters).
