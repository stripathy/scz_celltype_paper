# Figure generation & cross-platform validation (scripts 07–10)

Documentation + data provenance for the figures that present the snRNA-seq SCZ
DE meta-analysis and validate it against an independent Xenium spatial dataset.
Scripts 07 and 08 are standalone components; **script 09 is the publication
composite** that combines them (forests + concordance scatter) with a butterfly
summary, per-cell-type volcanoes, and Xenium exemplar cells; script 10 produces
the exemplar-cell inputs that panel K of the composite reads.

| Script | Figure(s) / output | What it shows |
|---|---|---|
| `scripts/07_forest_plots.R` | `results/figures/07_forest_composite.{png,pdf}` | Per-gene forest plots: each snRNA-seq cohort + pooled meta + Xenium replication |
| `scripts/08_meta_vs_xenium_scatter.R` | `results/figures/08_meta_vs_xenium_scatter*.{png,pdf}`, `08b_*_uniform*.{png,pdf}` | Concordance scatter of meta logFC vs Xenium logFC across all testable gene×cell-type pairs |
| `scripts/10_xenium_exemplar_cells.py` | `results/tables/exemplar_*.csv` | Per-cell boundary + marker-molecule coordinates for the panel-K exemplar cells (run **before** 09) |
| `scripts/09_composite_figure.R` | `results/figures/09_composite.{png,pdf}` | **Publication composite (6.5 × 9.2 in)**: butterfly (A), volcanoes (B,C), forests (D–I), concordance scatter (J), Xenium exemplar cells (K) |

The validation figures (07, 08, and panels D–K of 09) all ask the same question
from different angles: **does the snRNA-seq discovery DE replicate on an
independent spatial-transcriptomics platform?** The composite figure legend
(final wording + per-value provenance) lives in
[`notes/figure_composite_legend.md`](figure_composite_legend.md).

---

## 1. Data sources (UPDATE THESE IF THE DATA CHANGES)

Three input files feed these figures. If any DE analysis is re-run upstream,
update the corresponding file (or the path constant at the top of each script)
and re-render.

### A. Per-cohort snRNA-seq DE — `data/meta_results_cohorts_subclass.csv`
- **What**: limma differential expression, run separately in each of 7 SCZ
  snRNA-seq cohorts, for each SEA-AD subclass. One row per (cohort × cell_type × gene).
- **Cohorts (7)**: `Bat`, `HBCC`, `MSSM`, `Mclean`, `MtSinai`, `Multi`, `OFC`.
- **Columns used**: `genes`, `logFC`, `t`, `P.Value`, `adj.P.Val`, `cell_type`, `cohort`.
  (Also present, unused here: `AveExpr`, `B`, a leading row-name column.)
- **Size**: ~310 MB, ~2.16M rows.
- **Used by**: script 07 (the individual cohort rows + the recomputed meta diamond position).
- **Where it lives**: symlinked into `data/` (large; git-ignored). Source: the
  per-cohort outputs that the meta-analysis pools.

### B. Meta-analytic snRNA-seq DE — `data/DE_genes_all_cells_scz.csv`
- **What**: the random-effects meta-analysis (metafor) across the 7 cohorts,
  one row per (cell_type × gene). This is the project's canonical DE table.
- **Columns used**: `genes`, `cell_type`, `estimate`, `pval`, `padj`.
  (Also present: `se`, `ci.lb`, `ci.ub`, `k`, `tau2`, `I2`.)
- **Size**: ~35 MB, 222,028 rows, 16,368 unique genes, 23 cell types.
- **Used by**:
  - script 07 — ONLY for the canonical meta FDR (`padj`) that sets the diamond's asterisks.
  - script 08 — defines the testable set (filtered to `padj < PADJ_THR`) and the x-axis.
- **Where it lives**: symlinked into `data/` (git-ignored). See `data/README.md`.

### C. Xenium spatial DE — `~/Github/SCZ_Xenium/output/de/de_results_subclass.csv`
- **What**: pseudobulk DE on Xenium spatial transcriptomics, one row per
  (cell_type × gene). Produced by `SCZ_Xenium/code/.../run_de_edgepython.py`:
  raw counts summed per (sample × subclass) → edgeR quasi-likelihood GLM
  F-test for SCZ vs control, adjusting for sex + age.
- **Dataset**: Kwon et al. 2026 (GSE307404), 24 DLPFC sections (12 SCZ, 12 control),
  **300-gene panel**. Cell types annotated by label transfer from SEA-AD.
- **Columns used**: `gene`, `logFC`, `F`, `PValue`, `FDR`, `celltype`.
  (Also present: `logCPM`, `n_scz`, `n_ctrl`, `n_samples`, `n_genes_tested`, `cell_class`.)
- **Used by**: both scripts (the Xenium replication row / y-axis).
- **IMPORTANT — only 300 genes**: most snRNA-seq DE genes are NOT on the panel.
  A gene/cell-type pair can only be validated in Xenium if its gene is one of the 300.
- **Subclass-name harmonization**: Xenium uses `Astrocyte`, `L2/3 IT`, `L5/6 NP`,
  `Microglia-PVM`, `Oligodendrocyte`, `Endothelial`; the SEA-AD convention used
  everywhere else is `Astro`, `L2_3 IT`, `L5_6 NP`, `Micro-PVM`, `Oligo`, `Endo`.
  Scripts 07–09 remap via a `ct_map` lookup. **If the Xenium pipeline changes its
  subclass labels, update `ct_map` in all three scripts.**

> If the Xenium repo moves, edit `INPUT_XENIUM` at the top of scripts 07, 08
> **and 09**. If a newer Xenium DE table changes column names, update the
> `transmute()` calls and `ct_map`.

### D. Panel-K exemplar tables (derived) — `results/tables/exemplar_*.csv`
- **What**: per-cell boundary polygons + marker-gene transcript-molecule
  coordinates (µm, recentred), plus `exemplar_cells_meta.csv`. Produced by
  `scripts/10_xenium_exemplar_cells.py` from the Xenium **h5ad + boundary +
  transcript** exports (NOT the DE table). See §6 for full provenance.
- **Used by**: script 09 panel K only. **Run script 10 before script 09.**

---

## 2. Shared methods

> **Script 09 reuses everything in this section.** Its forests (D–I), scatter
> (J), volcanoes (B,C) and butterfly (A) read the same three DE inputs and use
> the same SE derivation, meta-diamond computation, significance markers,
> `ct_map`, and class colours described below. The constants are re-declared at
> the top of `09_composite_figure.R` (`INPUT_*`, `EXC`/`INH`/`GLI`, `CLASS_COL`).

### Standard error derivation
Neither input stores a usable SE directly, so both are derived from the test statistic:
- **limma cohorts**: `SE = logFC / t`.
- **Xenium edgeR QL F-test** (1-df numerator): `t = sign(logFC) * sqrt(F)`,
  then `SE = logFC / t`. This is a Wald-equivalent approximation; it is mildly
  biased for very large/small F but fine for the moderate effects shown.

### The meta diamond (script 07)
The pooled estimate shown as the black diamond is **recomputed** in-script via
`metafor::rma(yi = cohort logFC, sei = cohort SE, method = "DL")` (DerSimonian–
Laird random effects) from the 7 cohorts. Its **position** (estimate + 95% CI)
comes from this recomputation; its **asterisks** come from the canonical meta
table's `padj` (input B), so the significance shown matches the rest of the
project rather than a re-derived p-value. (The two agree closely.)

### Xenium is held OUT of the meta
Xenium is **never** pooled into the meta-analysis. It is shown only as
independent replication (the green triangle in 07; the y-axis in 08). This
preserves the discovery/replication separation — see
`notes/findings.md` discussion. (An 8-cohort sensitivity pool was explored
earlier but is not the primary analysis.)

### Significance markers (script 07)
Applied per row:
| Marker | Meaning |
|---|---|
| `***` | FDR < 0.01 |
| `**`  | FDR < 0.05 |
| `*`   | FDR < 0.10 |
| `•` (dot) | uncorrected p < 0.05 in an individual dataset (cohort or Xenium) not reaching FDR < 0.10 |
| `n.s.` | shown only on the meta diamond row when it is not FDR-significant |

- Cohort rows: asterisks from per-cohort `adj.P.Val`, dot from `P.Value`.
- Meta diamond: asterisks from canonical meta `padj`; never a dot (it is a synthesis).
- Xenium row: asterisks from Xenium `FDR`, dot from Xenium `PValue`.
- Consequence to expect: individual datasets usually show **dots** (nominally
  significant, under-powered per study); the pooled meta shows **asterisks**.
  Note Xenium's FDR is across only the 300-gene panel.

### Cell-class colour scheme (script 08)
`CLASS_COL`: Excitatory `#117733` (green), Inhibitory `#882255` (magenta),
Glia `#DDCC77` (yellow). Class assignment is hard-coded in the `EXC`/`INH`/`GLI`
vectors at the top of script 08 — update these if the subclass set changes.

### Inclusion threshold (script 08)
`PADJ_THR` (default 0.10) filters input B to the meta-DE genes. Passing a CLI
arg overrides it and adds an `_fdrNN` suffix to the output filenames so multiple
cutoffs coexist (e.g. `Rscript scripts/08_meta_vs_xenium_scatter.R 0.05` →
`08_meta_vs_xenium_scatter_fdr05.png`).

---

## 3. Figure 07 — forest plots

**How to read**: one panel per (gene, cell type). Top→bottom: the 7 snRNA-seq
cohorts (grey squares, sorted by effect), the pooled `snRNA-seq meta` diamond
(black), a dotted separator, then the `Xenium` replication (green triangle).
Bars are 95% CIs. Panel letters top-left, gene name (italic) as the subtitle,
markers per the table above.

**Composite panels** (curated in the `panels` tribble near the bottom of the
script): A. SST/Sst, B. PVALB/Pvalb, C. BDNF/L2_3 IT, D. FKBP5/OPC,
E. CX3CR1/Micro-PVM, F. SMAD1/Pvalb. Edit that tribble to change the panel set.

**Run**:
```bash
cd ~/Github/scz_pathway_enrichment
Rscript scripts/07_forest_plots.R                  # 6-panel composite
Rscript scripts/07_forest_plots.R GAD1 Sst         # any single (gene, cell type)
```

---

## 4. Figure 08 — concordance scatter

**How to read**: each point is a (gene, cell type) pair where the gene is
meta-DE (padj < PADJ_THR) AND on the Xenium panel AND tested in that subclass.
x = snRNA-seq meta logFC, y = Xenium logFC. Dashed line = y=x (perfect
agreement); solid line = OLS fit. Colour = cell class. In the primary figure,
**dot size = snRNA-seq meta FDR bin** (large < 0.05, small 0.05–0.10); the
`08b_*_uniform` variant makes all dots one size. A curated set of genes is
labelled with ggrepel (full-data repel + directional nudges so labels avoid all
points). The view is cropped to ±~0.96 via `coord_cartesian` so the bulk of the
data fills the panel — points outside the view remain in the fit and stats.

**On-plot stats**: Pearson `r` and `% concordant` only. The quadrant breakdown
(both-up / both-down / discordant) and binomial p are printed to the console
for the caption, not drawn.

**Labelled genes** (edit `FOREST_PAIRS` + `EXTRA_PAIRS` near the top of the
plotting section): the 5 significant forest-panel pairs + SERPING1/Astro,
ABCG2/Endo, VGF/Lamp5. The script warns if a requested pair is absent from the
testable set (e.g. PVALB/Pvalb, whose meta padj > 0.1, cannot appear).

**Run**:
```bash
Rscript scripts/08_meta_vs_xenium_scatter.R        # padj < 0.10 (primary)
Rscript scripts/08_meta_vs_xenium_scatter.R 0.05   # stricter sensitivity
```

**Current numbers** (padj < 0.10): 166 pairs, 60 both-up / 63 both-down /
43 discordant, 74% directionally concordant, Pearson r = 0.70, slope = 0.77,
binomial p = 2e-10. At padj < 0.05: 109 pairs, 76% concordant, r = 0.77.

---

## 5. Figure 09 — publication composite

The single multi-panel figure for the paper
(`results/figures/09_composite.{png,pdf}`, **6.5 × 9.2 in, 400 dpi**). Width is
**locked at 6.5 in** (`FIG_W`); height (`FIG_H = 9.2`) is free and was grown to
give the enlarged text room. It reuses the inputs (§1) and methods (§2) of
07/08 and adds the butterfly, volcanoes, and exemplar cells.

| Panel | Builder | Content |
|---|---|---|
| A | `build_butterfly()` | Up/down DE-gene counts per cell type; FDR<0.10 (light) with the FDR<0.05 subset overlaid (dark); 4-entry legend at `c(0.70,0.16)`. |
| B,C | `build_volcano()` | Volcanoes for Sst (B) and L2/3 IT (C); points coloured by the same up/down × FDR tiers as A (NS = grey); dashed line at FDR=0.10; selected genes via ggrepel. |
| D–I | `build_forest()` | Six (gene,cell) pairs (the `FOREST` tribble); cohort squares + DL meta diamond + Xenium triangle; markers per §2. |
| J | `build_scatter()` | meta logFC vs Xenium logFC over meta FDR<0.10 ∩ Xenium; colour = class, size = meta FDR tier; OLS fit + identity line; `r`/`% concordant` text. |
| K | `build_exemplar()` | Xenium exemplar cells (inputs from script 10): boundary + marker molecules. |

**Layout** (cowplot): row 1 = butterfly + volcano column (B over C), `rel_widths
c(1.85, 1)`; row 2 = 6 forests (3×2); row 3 = scatter J + exemplar matrix K,
`rel_widths c(1, 1)`; `rel_heights c(1.15, 1.05, 1.0)`. Panel letters 12 pt bold,
top-left.

**Conventions locked during the cosmetic pass** (change these together, not
piecemeal):
- **Text size**: master `BASE = 9` (nearly all text is sized relative to it).
  This is about the ceiling for the locked 6.5-in width before the panel-J
  labels collide; enlarging further would require shortening the J labels.
- **Forests**: the `SCZ log₂ FC` x-axis title is drawn on the **bottom row only**
  (`show_xlab = c(F,F,F,T,T,T)` in the `Map` call) so each column is labelled
  once. Panel set + order = the `FOREST` tribble (SST/Sst, PVALB/Pvalb,
  BDNF/L2_3 IT, FKBP5/OPC, CX3CR1/Micro-PVM, RASGRF2/Pvalb).
- **Volcano highlights**: edit the `build_volcano()` calls — B = SST, NAT16,
  SMAD1, AFG3L2; C = BDNF, SMAD1, VWA5B2, ADAMTS9-AS2, ST6GAL2.
- **Scatter J labels are placed deterministically**: `force = 1`, `force_pull =
  0.1`, plus explicit per-gene `nx`/`ny` nudges (each = `desired_label_pos −
  point`). The legend is moved **out of the plot to the bottom**
  (`legend.box="vertical"`, 2 rows) so the freed top-left holds the wide
  top-right labels (SERPING1, FKBP5); the in-plot `r` / `% concordant` text sits
  in the lower-right. Labelled pairs = the `lab_pairs` tribble (SST/Sst,
  BDNF/L2_3 IT, FKBP5/OPC, CX3CR1/Micro-PVM, SMAD1/Pvalb, SERPING1/Astro,
  RASGRF2/Pvalb). To move a label edit its `nx`/`ny`; to add/remove one edit
  `lab_pairs` **and** its nudge rows.
- **Butterfly legend** sits at `c(0.70, 0.16)` (shifted right of the bar-tip
  count labels), keys `unit(9,"pt")`.
- **Panel K is a labelled matrix**: `Control` / `SCZ` column headers + rotated
  row labels (`SST in Sst`, `RASGRF2 in Pvalb`). All four cells share one
  coordinate limit (`ex_lim` = max boundary extent × 1.15) so the zoom is
  identical and the **5 µm scale bar is the same physical length in every
  panel**; the `5 µm` text is on the SCZ/RASGRF2 cell only.

**Run** (script 10 must have produced the exemplar tables first):
```bash
cd ~/Github/scz_pathway_enrichment
python scripts/10_xenium_exemplar_cells.py   # panel-K inputs (if not present)
Rscript scripts/09_composite_figure.R         # -> results/09_composite.{png,pdf}
cp results/09_composite.png results/09_composite.pdf results/figures/  # snapshot
```

---

## 6. Script 10 — Xenium exemplar cells (panel K inputs)

Extracts, for each (gene, subclass) pair, one **Control** and one **SCZ**
exemplar cell, writing the cell boundary + the marker-gene transcript molecules
inside it, for panel K.

- **Pairs** (`PAIRS`): `SST`/Sst and `RASGRF2`/Pvalb. Chosen because both are
  (1) strongly meta-DE-down in the snRNA-seq meta, (2) expressed highly enough
  per cell to render as countable molecule dots, and (3) depth-matchable between
  Control and SCZ cells. Rejected: BDNF, SMAD1 (~1 molecule/cell — no visible
  dot range) and CX3CR1, TF (total-count / segmentation depth-confounded).
- **Samples**: `Br6432` (Control), `Br2039` (SCZ) — the two Xenium sections with
  transcript-level molecule export; diagnosis from SCZ_Xenium `config.SAMPLE_TO_DX`.
- **Raw inputs** (SCZ_Xenium repo):
  - `output/h5ad/<sample>_annotated.h5ad` — raw counts (`.X`), `subclass_label`,
    `qc_pass`, `total_counts`.
  - `output/deploy/boundaries/<sample>.json` — per-cell polygons (25 verts;
    decode `micron = quant*scale + offset`); polygon index = h5ad obs index.
  - `output/deploy/transcripts/<sample>/<GENE>.json` (+ `gene_index.json`) —
    per-gene molecule coordinates.
- **Cell selection** (`pick_cell`): restrict to the subclass, `qc_pass`, marker
  count > 0; keep typically-sized cells (total counts in the 20–90th percentile);
  prefer cells whose marker count is near the group target (`group_target` =
  subclass median, or the expressing-cell median when the group median is 0);
  among those pick the **roundest convex** cell (highest circularity `4πA/P²`,
  solidity `A/hull ≥ 0.93`, area ≥ 40 µm²) to avoid segmentation artefacts from
  neighbouring cells. Molecules are clipped to the polygon by point-in-polygon
  (`matplotlib.path.Path`).
- **Outputs**: `exemplar_<gene>_<dx>_{boundary,dots}.csv` (recentred µm coords)
  and `exemplar_cells_meta.csv` (sample, `marker_count`, `n_dots_in_poly`,
  `total_counts`, `circularity`, `solidity`, `area_um2`).
- **Current exemplars** (from `exemplar_cells_meta.csv`): SST/Sst 47 (Ctrl) vs
  10 (SCZ); RASGRF2/Pvalb 24 vs 5 molecules; all circularity 0.97–0.98, solidity
  ≈ 1.0. Panel K shows `n_dots_in_poly` (molecules inside the boundary).
- **To change the exemplar genes**: edit `PAIRS` (and, for cross-referencing,
  `lab_pairs`/`FOREST` in script 09), then re-run 10 → 09.

**Run**:
```bash
python scripts/10_xenium_exemplar_cells.py
```

---

## 7. Update checklist

When upstream DE results change, to refresh the figures:

1. Replace/repoint the relevant input file(s):
   - new meta-analysis → `data/DE_genes_all_cells_scz.csv` (check cols `estimate`, `pval`, `padj`)
   - new per-cohort DE → `data/meta_results_cohorts_subclass.csv` (check cols `logFC`, `t`, `P.Value`, `adj.P.Val`, `cohort`)
   - new Xenium DE → `INPUT_XENIUM` path (check cols `gene`, `logFC`, `F`, `PValue`, `FDR`, `celltype`)
2. If subclass labels changed in any input, update `ct_map` (scripts 07–09) and
   the `EXC`/`INH`/`GLI` vectors (scripts 08 and 09).
3. If the cohort set changed, no code change needed (cohorts are read from the
   `cohort` column), but check the forest layouts still fit.
4. If the Xenium h5ad / boundary / transcript exports change (or you want
   different exemplar genes), edit `PAIRS` / sample constants in script 10 and
   **re-run script 10 first** — it regenerates `results/tables/exemplar_*.csv`
   that panel K reads.
5. Re-run the figures and copy outputs into `results/figures/`:
   - `Rscript scripts/07_forest_plots.R`
   - `Rscript scripts/08_meta_vs_xenium_scatter.R`
   - `Rscript scripts/09_composite_figure.R`   (after step 4)
6. Re-verify every cited number against source data (per the project's
   numerical-precision rule) and refresh the provenance table in
   `notes/figure_composite_legend.md` before using any value in text.

---

## 8. Caveats

- **300-gene Xenium panel** limits replication to ~a third of meta-DE genes;
  key novel candidates (NAT16, AFG3L2, TOMM40, SREBF2, …) are not on the panel
  and cannot be cross-validated here.
- **Xenium FDR is panel-wide** (only 300 genes), so it is not directly
  comparable to the genome-wide snRNA-seq FDR. Under the marker scheme this
  makes most Xenium rows show as dots even when nominally strong.
- **Cell-type-marker circularity**: for genes that define a Xenium cell type
  (e.g. PVALB in Pvalb), the Xenium effect can be mildly biased by the
  label-transfer step. Concordance for such genes is supportive but not fully
  independent.
- **SE approximations** (logFC/t, logFC/√F) are Wald-equivalents; exact pooling
  would back out SE from p-value + df.
- **Panel-K exemplars are illustrative, not quantitative**: each is a single
  hand-vetted cell per group (Control/SCZ), depth-matched and chosen for round
  morphology to minimise segmentation artefacts. The molecule counts shown
  depend on the specific cell; the quantitative SST/RASGRF2 down-regulation is
  carried by the DE / forest / scatter panels, not by these two cells.
- **Only two Xenium sections** (`Br6432`, `Br2039`) have transcript-level
  molecule export, so exemplars are necessarily drawn from those.
