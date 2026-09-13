# genetics/ — SCZ common-variant risk across cell types

Produces **Figure 4**, **Supplementary Figs. S9 and S10**, and **Supplementary
Table T6**. Everything else that once lived here (the combined SEA-AD + Siletti
"Franken" taxonomy, the Siletti whole-brain enrichment, conditional/forward
selection, gene-driver scatters, the interactive web app) was retired and
removed on 2026-09-04; it is recoverable from the git tag
`pre-prune-2026-09-04`.

## The analysis in one paragraph

SCZ common-variant association (**Bigdeli et al. 2026**, European ancestry) is
tested for enrichment in the genes specific to each cell type, by MAGMA
gene-property regression. Specificity is computed across the **125 supertypes
of the SEA-AD neurotypical DLPFC (A9) taxonomy** — the region matched to the
paper's frontal-cortex datasets. Enrichment is then compared with the
compositional depletion measured in Figure 3, and the convergence is followed
down to one gene (*HCN1*), one locus, and the intrinsic physiology of the
patch-seq cells that express it.

PGC3 (Trubetskoy et al. 2022) and the SEA-AD **MTG** taxonomy (137 supertypes)
appear only as robustness checks in S10. The combined SEA-AD + Siletti taxonomy
was dropped on L. Duncan's advice (2026-08-29).

## What builds what

Run from the repo root. Steps 1–2 need the large inputs in `data/`
(see `data/README.md`); step 3 onward runs from committed files.

```bash
# 1. Expression reference -> per-supertype means      (needs the 3 A9 h5ads)
python3 genetics/scripts/figures/build_dlpfc_specificity.py    # DLPFC 125
python3 genetics/scripts/figures/seaad_supertype_log1p.py      # MTG 137, S10 only

# 2. Specificity matrices + the four MAGMA runs       (needs MAGMA + sumstats)
python3 genetics/scripts/figures/build_spec_seaad_only.py
#    -> results/intermediates/T_{a9only,mtgonly}_{bigdeli,pgc3}.gsa.out  [tracked]

# 3. Figure 4 panel data -> results/figures/r_panels/  [tracked]
python3 genetics/scripts/figures/build_composition_table.py     # crumblr betas
python3 genetics/scripts/figures/build_sst_ephys_summary.py     # patch-seq sag
python3 genetics/scripts/figures/export_panels_abc.py           # a, b, c
python3 genetics/scripts/figures/export_panel_d_genetrack.py    # c gene track
python3 genetics/scripts/figures/export_panels_dehi.py          # d, g, h
python3 genetics/scripts/figures/export_panels_ef.py            # e, f
python3 genetics/scripts/figures/export_panel_ad_concordance.py # i
python3 genetics/scripts/figures/r_panels_provenance.py         # MANIFEST.tsv

# 4. Render
Rscript genetics/scripts/figures/scz_sst_hcn1_story.R    # Figure 4  (9 panels)

# 5. Supplements and the patch-seq table
python3 genetics/scripts/figures/export_supp_enrichment_seaad125.py
Rscript genetics/scripts/figures/plot_supp_enrichment_seaad125.R   # S9
python3 genetics/scripts/figures/export_fig4_robustness.py
Rscript genetics/scripts/figures/plot_fig4_robustness.R            # S10
python3 genetics/scripts/figures/build_supp_table_patchseq_labels.py  # T6
```

`export_panel_d_genetrack.py` must run **after** `export_panels_abc.py`, which
writes the plot window it reads. Everything else in step 3 is independent.

## Figure 4, panel by panel

The renderer builds panels under their historical letters and
`fig4_assemble_nod.R` relabels them a–i (an earlier HCN1-expression-vs-depletion
panel was dropped, so the old e–j became d–i).

| Panel | Shows | Data | Built by |
|---|---|---|---|
| a | SCZ enrichment vs compositional depletion, 16 Sst supertypes | `panel_B_genetics_vs_depletion.csv` | `export_panels_abc.py` |
| b | Per-gene drivers of Sst_2 enrichment; *HCN1* highlighted | `panel_C_*.csv` | `export_panels_abc.py` |
| c | *HCN1* locus zoom, credible set, gene track | `panel_D_*.csv` | `export_panels_abc.py` + `export_panel_d_genetrack.py` |
| d | *HCN1* expression vs patch-seq voltage sag | `panel_E_hcn1_vs_sag.csv` | `export_panels_dehi.py` |
| e | Five exemplar Sst reconstructions, ordered by soma depth | `panel_F_*.csv` | `export_panels_ef.py` |
| f | Voltage responses of the same five cells | `panel_G_*.csv` | `export_panels_ef.py` |
| g | Depleted vs not-depleted Sst marker volcano | `panel_volcano_*.csv` | `export_panels_dehi.py` |
| h | *CALB1* by depletion group | `panel_violin_*.csv` | `export_panels_dehi.py` |
| i | SCZ depletion vs SEA-AD Alzheimer's CPS slope | `panel_ad_concordance_sst.csv` | `export_panel_ad_concordance.py` |

Layout is `fig4_assemble_nod.R`; shared fonts, colours and `theme_panel()` are
in `fig4_style.R`, which also styles panels g–i in `fig4_new_panels.R`.

**Why the panel CSVs are committed.** The renderer draws from them, not from
the analyses that produced them, so the figure rebuilds from a clone without
the 329 MB of summary statistics. That also means it can go quietly stale, so
`r_panels/MANIFEST.tsv` records the upstream file behind every panel and
`shared/figure_inputs.R` halts the render when one has moved on.

Output: `results/figures/scz_sst_hcn1_figure4.{png,pdf}`. The renderer also
writes an SVG for hand-editing; it is git-ignored because the patch-seq traces
make it ~25 MB.

## Supplements and tables

| Item | Output | Chain |
|---|---|---|
| Supp Fig S9 | `manuscript/figures/supplementary/S09_scz_enrichment_seaad125.{png,pdf}` | `export_supp_enrichment_seaad125.py` -> `results/tables/scz_enrichment_seaad125_bigdeli_dlpfc.csv` -> `plot_supp_enrichment_seaad125.R` |
| Supp Fig S10 | `manuscript/figures/supplementary/S10_genetics_ad_robustness.{png,pdf}` | `export_fig4_robustness.py` -> `results/tables/fig4_robustness_{sst16,ad_mtg}.csv` -> `plot_fig4_robustness.R` |
| Supp Table T6 | `results/tables/supp_T6_patchseq_sst_annotations.csv` | `build_supp_table_patchseq_labels.py` |

The S9 table doubles as the paper's supplementary table of SCZ enrichment: one
row per SEA-AD DLPFC supertype, so the figure and the multiple-testing universe
are the same 125 types.

Legend drafts sit beside the figures they describe, in `results/figures/`.

## Cross-component seams

- **In:** the 7-dataset crumblr composition meta-analysis, via
  `shared/snrnaseq_de/nicole_scz_snrnaseq_betas/` (panels a and i).
- **In:** the Alzheimer's pseudo-progression slopes from `crossdisorder/`
  (panel i, and S10b).
- **Out:** `results/figures/r_panels/panel_volcano_vulnerable_vs_notdepleted.csv`
  and `results/intermediates/*_mean_expression*.csv` are read by the reserve
  analyses in `reserve/sst_strata_supp/`.

## Citations

- Bigdeli et al. 2026 — SCZ GWAS, European ancestry (primary)
- Trubetskoy et al. 2022, *Nature* — PGC3 SCZ GWAS (robustness)
- Gabitto et al. 2024, *Nature Neuroscience* — SEA-AD taxonomy and reference
- de Leeuw et al. 2015, *PLoS Comput Biol* — MAGMA
