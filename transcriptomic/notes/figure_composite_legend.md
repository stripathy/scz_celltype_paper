# Composite figure legend (script 09)

Nature Neuroscience style (terse; description only, no interpretation).
Figure file: `results/09_composite.{png,pdf}` (7.1 × 6.625 in, 400 dpi; 10 panels, a–j).

Layout: one interneuron marker per row across four views (volcano → forest →
per-donor expression → exemplar cells), then a transcriptome-wide overview.
Row 1 (a–d) SST in Sst cells; Row 2 (e–h) PVALB in Pvalb cells; Row 3 (i)
DE-gene burden per cell type with a DE-vs-proportion inset, (j) cross-platform
concordance.

---

**Fig. X | Cross-platform schizophrenia (SCZ) differential expression for the
two canonical interneuron markers, in the context of transcriptome-wide DE
burden.**

(**a–d**) **SST in Sst cells.** (a) Volcano of SCZ log₂ fold change versus
−log₁₀ FDR for genes tested in Sst cells; points coloured by direction × FDR
tier (up orange, down blue; dark FDR < 0.05, light FDR < 0.10; NS grey); SST,
NAT16, SMAD1 and AFG3L2 labelled; data-driven axes. (b) Forest plot for SST:
seven frontal cortical snRNA-seq cohorts (grey squares), pooled random-effects
(DerSimonian–Laird) meta-analytic estimate (black diamond; log₂FC −0.46,
FDR = 0.049) and independent Xenium spatial estimate (green triangle; log₂FC
−0.32, P = 0.052); whiskers, 95% CI. (c) Per-donor library-normalised expression
(counts per 1,000 transcripts, CP1K) of SST in Sst cells, Control versus SCZ
(12 vs 12 donors); p, Xenium edgeR quasi-likelihood test (P = 0.052; the DE
defining the Xenium estimate in b). (d) Representative Xenium Sst cells — the
cell at the pooled group-median transcript density per diagnosis (Control
Br6432; SCZ Br5973), matched for size and typical outline; grey outline, cell
boundary; dashed, nucleus; red dots, SST transcript molecules (count, top-left);
scale bar, 5 µm. (**e–h**) **PVALB in Pvalb cells**, same four views: (e) volcano
(PVALB, SMAD1, ANXA2, SCN3A, NAT16 labelled); (f) forest (meta log₂FC −0.06,
FDR = 0.86, n.s.; Xenium −0.22, P = 0.044); (g) CP1K (edgeR P = 0.044); (h)
exemplar Pvalb cells (Control Br6432; SCZ Br5973). (**i**) Up- (orange) and
down-regulated (blue) DE-gene counts per cell-type subclass (23 subclasses) from
the meta-analysis, FDR < 0.10 (light) with the FDR < 0.05 subset overlaid (dark);
counts annotated. *Inset:* number of DE genes (FDR < 0.10) versus mean per-donor
cell-type proportion (Xenium, log₁₀ axis) across the 22 subclasses with ≥ 1 DE
gene; points coloured by cell class; Astro, L5 IT, Vip and L6b labelled; Spearman
ρ = 0.83. (**j**) snRNA-seq meta-analytic versus Xenium log₂ fold change for all
gene × cell-type pairs with meta FDR < 0.10 (n = 166); points coloured by cell
class (Excitatory, Inhibitory, Glia) and sized by meta FDR (larger, FDR < 0.05;
smaller, 0.05–0.10); dashed line, identity. Pearson r = 0.73; 76% sign-concordant.
Significance throughout: \*FDR < 0.10, \*\*FDR < 0.05, \*\*\*FDR < 0.01; •, nominal
P < 0.05 (FDR ≥ 0.10); n.s., pooled estimate not significant.

**Supplementary.** *Per-cell SST/PVALB reduction across normalisations* (raw
dots, grain density, library-normalised, library size; negative-binomial mixed
model) is in **Supplementary Fig. S_percell**. The standalone PVALB supplement
(`scripts/14_supp_pvalb.R`, *S_pvalb*) is now superseded by main-figure row 2
(e–h); the fully-labelled DE-vs-proportion scatter (all subclasses) is
`scripts/16_de_vs_proportion.R`.

---

## Provenance of cited values (verified programmatically, do not edit by hand)

| Claim | Value | Source |
|---|---|---|
| Panels / canvas | 10 (a–j); 7.1 × 6.625 in; 400 dpi | `scripts/09_composite_figure.R` `ggsave` |
| Subclasses (panel i) | 23 | `data/DE_genes_all_cells_scz.csv` (unique `cell_type`) |
| snRNA-seq cohorts | 7 (Bat, HBCC, Mclean, MSSM, MtSinai, Multi, OFC) | `data/meta_results_cohorts_subclass.csv` (unique `cohort`) |
| Meta model | random-effects DerSimonian–Laird (`rma(method="DL")`) | script 09 |
| Volcano a (Sst) genes | SST, NAT16, SMAD1, AFG3L2 | `build_volcano("Sst", …)`, script 09 |
| Volcano e (Pvalb) genes | PVALB, SMAD1, ANXA2, SCN3A, NAT16 | `build_volcano("Pvalb", …)`, script 09 |
| Forest b — SST/Sst | meta log₂FC −0.458, FDR 0.0487 (\*\*); Xenium log₂FC −0.318, P 0.0521 (FDR 0.304); 7 cohorts | `DE_genes_all_cells_scz.csv` (`padj`); `de_results_subclass.csv` (`logFC`,`PValue`,`FDR`) |
| Forest f — PVALB/Pvalb | meta log₂FC −0.056, FDR 0.8556 (n.s.); Xenium log₂FC −0.219, P 0.0438 (FDR 0.392); 7 cohorts | same |
| CP1K c — SST/Sst (edgeR) | logFC −0.318, P = 0.0521 | `results/tables/marker_norm_expr_stats.csv` |
| CP1K g — PVALB/Pvalb (edgeR) | logFC −0.219, P = 0.0438 | `results/tables/marker_norm_expr_stats.csv` |
| Exemplar d — SST dots | Control Br6432 35 (grain density 18.7); SCZ Br5973 28 (13.2) | `results/tables/exemplar_cells_meta.csv` (`n_dots_in_poly`,`disp_grain_density`) |
| Exemplar h — PVALB dots | Control Br6432 10 (grain density 5.2); SCZ Br5973 9 (4.2) | `results/tables/exemplar_cells_meta.csv` |
| Exemplar eccentricity | SST 0.58 / 0.51; PVALB 0.60 / 0.62 (Ctrl / SCZ) | `results/tables/exemplar_cells_meta.csv` (`eccentricity`) |
| Scale bar | 5 µm | `build_exemplar()` `sb <- 5` |
| Inset i — relationship | # DE genes (FDR<0.10) vs Xenium mean per-donor proportion; subclasses with ≥1 DE | `build_de_prop_inset()`; `../spatial/output/crumblr/crumblr_input_subclass_corr.csv` |
| Inset i — n / ρ | n = 22 subclasses; Spearman ρ = 0.8272 → 0.83 | as above (`cor(prop, n_de, method="spearman")`) |
| Inset i — labelled | Astro, L5 IT, Vip, L6b | `build_de_prop_inset()`, script 09 |
| Panel j — n | 166 gene × cell-type pairs (meta FDR<0.10 ∩ Xenium) | inner join, script 09 |
| Panel j — Pearson r | 0.7310 → 0.73 | `cor(meta_est, xen_logFC)` |
| Panel j — concordance | 75.90% → 76% | `mean(sign(meta_est)==sign(xen_logFC))` |
| Panel j — size split | FDR<0.05 n = 109; 0.05–0.10 n = 57 | `meta_padj` |
| Panel j — labelled pairs | SST/Sst, BDNF/L2_3 IT, FKBP5/OPC, CX3CR1/Micro-PVM, SMAD1/Pvalb, SERPING1/Astro, FGFR3/Astro | `lab_pairs`, script 09 |
| Panel j — fit line | none (identity dashed line only; geom_smooth removed) | `build_scatter()`, script 09 |
| Sig. symbols | \*\*\* FDR<0.01, \*\* FDR<0.05, \* FDR<0.10 | `ast()` in script 09 |
| `•` symbol | nominal P<0.05 AND FDR≥0.10 (drawn only when NOT starred) | `is_dot()` in script 09 |
| **Suppl. S_percell** (NB mixed model, grain density `m_area`) | SST 0.72× (p = 0.002), PVALB 0.87× (p = 0.086) | `results/tables/S_percell_stats.csv`; `scripts/11` + `scripts/13` |
| **"frontal cortical"** | **UNVERIFIED — no brain-region/tissue field in the meta/cohort tables** | external (source studies); OFC = orbitofrontal is consistent |

Regenerate by loading the CSVs above and recomputing (e.g. the verification block
in the session that produced this file). Re-run `scripts/09` (composite) after any
data change, then re-check this table.
