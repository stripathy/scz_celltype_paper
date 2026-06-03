# Composite figure legend (script 09)

Nature Neuroscience style (terse; description only, no interpretation).
Figure file: `results/figures/09_composite.{png,pdf}` (6.5 × 9.2 in).
Wording below is Shreejoy's edited version; the only change applied on top is
the `•` symbol definition (precision fix — see provenance notes).

---

**Fig. X | Cell-type-resolved differential expression in schizophrenia (SCZ)
with cross-platform replication.**

(A) Differentially expressed (DE) gene counts per cell type subclass (23
subclasses) from a random-effects meta-analysis of seven frontal cortical
snRNA-seq cohorts. (B,C) Volcano plots for Sst (B) and L2/3 IT (C): SCZ log₂
fold change versus −log₁₀ FDR; points colored as in A; selected genes labeled.
(D–I) Forest plots for individual mRNA DE genes: SST mRNA (Sst cells, D),
PVALB mRNA (Pvalb cells, E), BDNF (L2/3 IT, F), FKBP5 (OPC, G), CX3CR1
(Micro-PVM, H) and RASGRF2 (Pvalb, I). Individual cohorts (grey squares; up to
seven), pooled meta-analytic estimate (black diamond) and independent Xenium
spatial replication (green triangle); whiskers, 95% CI; \*FDR < 0.10, \*\*FDR <
0.05, \*\*\*FDR < 0.01; •, nominal P < 0.05 (FDR ≥ 0.10); n.s., pooled estimate
not significant. (J) snRNA-seq meta-analytic log₂ versus Xenium log₂ fold
change for all gene × cell-type pairs with meta FDR < 0.10 (n = 166); line,
linear fit (shaded 95% CI); dashed line, identity. Pearson r = 0.70; 74%
sign-concordant. (K) Exemplar Xenium cells, SST mRNA in Sst cells (top) and
RASGRF2 in Pvalb cells (bottom), one control (Br6432) and one SCZ (Br2039)
cell each; grey outline, segmented cell boundary; red dots, individual
marker-gene transcript molecules within the boundary (count, top left). Scale
bar, 5 µm.

---

## Provenance of cited values (verified programmatically, do not edit by hand)

| Claim | Value | Source |
|---|---|---|
| Subclasses (panel A) | 23 | `data/DE_genes_all_cells_scz.csv` (unique `cell_type`) |
| snRNA-seq cohorts | 7 (Bat, HBCC, Mclean, MSSM, MtSinai, Multi, OFC) | `data/meta_results_cohorts_subclass.csv` (unique `cohort`) |
| Meta model | random-effects DerSimonian–Laird (`rma(method="DL")`) — not named in current legend text | `scripts/09_composite_figure.R` |
| Forest cohorts available | SST 7, PVALB 7, BDNF 7, FKBP5 6, CX3CR1 5, RASGRF2 7 | per-cohort table (non-NA logFC, t≠0) |
| Sig. symbols | \*\*\* FDR<0.01, \*\* FDR<0.05, \* FDR<0.10 | `ast()` in script 09 |
| `•` symbol | nominal P<0.05 AND FDR≥0.10 (drawn only when NOT starred; mutually exclusive with stars) | `is_dot()` in script 09 |
| Panel J n | 166 gene × cell-type pairs (meta FDR<0.10 ∩ Xenium) | inner join, script 09 |
| Panel J Pearson r | 0.7049 → 0.70 | `cor(meta_est, xen_logFC)` |
| Panel J concordance | 74.10% → 74% | mean(sign(meta)==sign(xen)) |
| Exemplar molecule counts | SST 47 / 10; RASGRF2 24 / 5 | `results/tables/exemplar_cells_meta.csv` (`n_dots_in_poly`) |
| Xenium samples | Br6432 = Control, Br2039 = SCZ | `exemplar_cells_meta.csv` (`sample`,`dx`); SCZ_Xenium `config.SAMPLE_TO_DX` |
| Scale bar | 5 µm | `build_exemplar()` `sb <- 5` |
| **"frontal cortical"** | **UNVERIFIED — no brain-region/tissue field in the meta or cohort tables** | external (source studies); OFC = orbitofrontal is consistent, confirm others (e.g. DLPFC/PFC) |

Regenerate the verification by loading the three CSVs above and recomputing.
The figure renders these same values, so any edit to the underlying data
requires re-running `scripts/09_composite_figure.R` and re-checking this table.
