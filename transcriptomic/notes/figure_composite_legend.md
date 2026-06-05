# Composite figure legend (script 09)

Nature Neuroscience style (terse; description only, no interpretation).
Figure file: `results/figures/09_composite.{png,pdf}` (7.1 × 6.7 in; 12 panels, a–l).

---

**Fig. X | Cell-type-resolved differential expression in schizophrenia (SCZ)
with cross-platform replication.**

(a) Differentially expressed (DE) gene counts per cell type subclass (23
subclasses) from a random-effects meta-analysis of seven frontal cortical
snRNA-seq cohorts; up-regulated (orange) and down-regulated (blue) at FDR < 0.10
(light) with the FDR < 0.05 subset overlaid (dark). (b–e) Volcano plots for Sst
(b), L2/3 IT (c), Astro (d) and Micro-PVM (e): SCZ log₂ fold change versus
−log₁₀ FDR; points colored by the up/down × FDR tiers in a (NS, grey); selected
genes labeled; each panel on independent, data-driven axes. The glial panels
show a reactive shift — astrocytes (d) up-regulate SERPING1 and CHI3L1 while
down-regulating the identity genes FGFR3 and NOTCH1; microglia (e) up-regulate
complement (C1QA, C1QB) and down-regulate homeostatic markers (CX3CR1, P2RY12)
and SORL1. (f–i) Forest plots for individual mRNA DE genes: SST (Sst cells, f),
BDNF (L2/3 IT, g), FGFR3 (Astro, h) and FKBP5 (OPC, i). Individual cohorts (grey
squares; up to seven), pooled meta-analytic estimate (black diamond) and
independent Xenium spatial replication (green triangle); whiskers, 95% CI;
\*FDR < 0.10, \*\*FDR < 0.05, \*\*\*FDR < 0.01; •, nominal P < 0.05 (FDR ≥ 0.10);
n.s., pooled estimate not significant. (j) snRNA-seq meta-analytic log₂ versus
Xenium log₂ fold change for all gene × cell-type pairs with meta FDR < 0.10
(n = 166); points colored by cell class (Excitatory, Inhibitory, Glia; as in a)
and sized by meta FDR (larger, FDR < 0.05; smaller, 0.05–0.10); line, linear fit
(shaded 95% CI); dashed line, identity. Pearson r = 0.73; 76% sign-concordant.
(k) Library-size-normalised expression (counts per 1,000 transcripts) of SST in
Sst cells (top) and FGFR3 in astrocytes (bottom), Control versus SCZ; boxplots,
per-donor pseudobulk values (12 control, 12 SCZ); p, Xenium edgeR
quasi-likelihood test (the same DE that defines the Xenium replication in f–j).
(l) Exemplar Xenium cells, one Control and one SCZ cell each — the representative
cell at the pooled group-median grain density for its diagnosis, matched for cell
size and typical outline (median eccentricity); grey outline, segmented cell
boundary; dashed outline, nucleus; red dots, individual marker-gene transcript
molecules within the cell boundary (count, top left); scale bar, 5 µm. Sst cells
are drawn from Br6432 (Control) / Br5973 (SCZ); astrocytes from Br5400 (Control)
/ Br5973 (SCZ). The panel-k titles (SST mRNA in Sst cells, FGFR3 mRNA in
Astrocytes) also name the aligned exemplar rows in l.

**Supplementary figures.** *PVALB mRNA in Pvalb cells* across the same three views
(forest, per-donor CP1K, exemplar cells) is in **Supplementary Fig. S_pvalb**
(meta not significant, but Xenium-replicated: edgeR p = 0.044; exemplars 10 vs 9
molecules). *Per-cell SST/PVALB reduction across normalisations* (raw dots, grain
density, library-normalised, library size; negative-binomial mixed model) is in
**Supplementary Fig. S_percell**.

---

## Provenance of cited values (verified programmatically, do not edit by hand)

| Claim | Value | Source |
|---|---|---|
| Subclasses (panel a) | 23 | `data/DE_genes_all_cells_scz.csv` (unique `cell_type`) |
| snRNA-seq cohorts | 7 (Bat, HBCC, Mclean, MSSM, MtSinai, Multi, OFC) | `data/meta_results_cohorts_subclass.csv` (unique `cohort`) |
| Meta model | random-effects DerSimonian–Laird (`rma(method="DL")`) | `scripts/09_composite_figure.R` |
| Volcano highlight genes (b–e) | b Sst: SST/NAT16/SMAD1/AFG3L2; c L2/3 IT: BDNF/SMAD1/VWA5B2/ADAMTS9-AS2/ST6GAL2; d Astro: SERPING1/CHI3L1/FGFR3/NOTCH1; e Micro-PVM: C1QA/C1QB/CX3CR1/P2RY12/SORL1 | `build_volcano()` calls, script 09 |
| Forests (f–i) cohorts available | SST 7, BDNF 7, FGFR3 7, FKBP5 6 | per-cohort table (non-NA logFC, t≠0) |
| Forest meta FDR (f–i) | SST 0.049 (\*\*), BDNF 0.034 (\*\*), FGFR3 1.5e-4 (\*\*\*), FKBP5 0.091 (\*) | `DE_genes_all_cells_scz.csv` `padj`; `ast()` |
| Forest Xenium (f–i) | SST −0.32 (p=0.052), BDNF −0.72 (p=0.001), FGFR3 −0.15 (p=0.075), FKBP5 +0.79 (p=0.003) | `de_results_subclass.csv` (`logFC`,`PValue`) |
| Sig. symbols | \*\*\* FDR<0.01, \*\* FDR<0.05, \* FDR<0.10 | `ast()` in script 09 |
| `•` symbol | nominal P<0.05 AND FDR≥0.10 (drawn only when NOT starred) | `is_dot()` in script 09 |
| Panel j n | 166 gene × cell-type pairs (meta FDR<0.10 ∩ Xenium) | inner join, script 09 |
| Panel j Pearson r | 0.7310 → 0.73 | `cor(meta_est, xen_logFC)` |
| Panel j concordance | 75.90% → 76% | mean(sign(meta)==sign(xen)) |
| Panel j labelled pairs | SST/Sst, BDNF/L2_3 IT, FKBP5/OPC, CX3CR1/Micro-PVM, SMAD1/Pvalb, SERPING1/Astro, FGFR3/Astro | `lab_pairs`, script 09 |
| Panel k CP1K (edgeR DE) | SST/Sst logFC −0.32, p = 0.052; FGFR3/Astro logFC −0.15, p = 0.075 | `marker_norm_expr_stats.csv`; per-donor CP1K (counts/1,000 tx) via `scripts/12` |
| Panel l exemplar dots (in cell) | SST 35 / 28; FGFR3 14 / 12 (Ctrl / SCZ) | `exemplar_cells_meta.csv` (`n_dots_in_poly`) |
| Panel l exemplar grain density | SST 18.7 → 13.2; FGFR3 7.1 → 6.5 grains/100 µm² (Ctrl → SCZ) ≈ pooled group-median | `exemplar_cells_meta.csv` (`disp_grain_density`) |
| Panel l exemplar eccentricity | SST 0.58 / 0.51; FGFR3 0.57 / 0.52 (≈ median ~0.57) | `exemplar_cells_meta.csv` (`eccentricity`, `ecc_target`) |
| Panel l Xenium sections | SST Br6432/Br5973; FGFR3 Br5400/Br5973 (Ctrl/SCZ) | `exemplar_cells_meta.csv` (`sample`,`dx`); SCZ_Xenium `config.SAMPLE_TO_DX` |
| Scale bar | 5 µm | `build_exemplar()` `sb <- 5` |
| **Suppl. S_pvalb** (PVALB/Pvalb) | meta padj 0.86 (n.s.), Xenium −0.22 (p=0.044); CP1K p=0.044; exemplars 10 / 9 dots, gd 5.2 → 4.2; Br6432/Br5973 | `scripts/14_supp_pvalb.R`; `marker_norm_expr_stats.csv`; `exemplar_cells_meta.csv` |
| **Suppl. S_percell** (NB mixed model, grain density) | SST 0.72× (p=0.002), PVALB 0.88× (p=0.086) | `results/tables/S_percell_stats.csv` (metric `m_area`); `scripts/11` + `scripts/13` |
| **"frontal cortical"** | **UNVERIFIED — no brain-region/tissue field in the meta/cohort tables** | external (source studies); OFC = orbitofrontal is consistent |

Regenerate by loading the CSVs above and recomputing. Re-run `scripts/09`
(composite) and `scripts/14` (PVALB supplement) after any data change, then
re-check this table.
