# Figure generation & cross-platform validation (scripts 07–14)

Documentation + data provenance for the figures that present the snRNA-seq SCZ
DE meta-analysis and validate it against an independent Xenium spatial dataset.
Scripts 07 and 08 are standalone components; **script 09 is the publication
composite**, laid out one interneuron marker per row — volcano → forest → CP1K
boxplot → exemplar cells — for SST (a–d) and PVALB (e–h), then a transcriptome-wide
DE-burden butterfly with a DE-vs-proportion inset (i) and the concordance scatter
(j); script 10 produces the exemplar-cell inputs that panels d/h read.

| Script | Figure(s) / output | What it shows |
|---|---|---|
| `scripts/07_forest_plots.R` | `results/figures/07_forest_composite.{png,pdf}` | Per-gene forest plots: each snRNA-seq cohort + pooled meta + Xenium replication |
| `scripts/08_meta_vs_xenium_scatter.R` | `results/figures/08_meta_vs_xenium_scatter*.{png,pdf}`, `08b_*_uniform*.{png,pdf}` | Concordance scatter of meta logFC vs Xenium logFC across all testable gene×cell-type pairs |
| `scripts/10_xenium_exemplar_cells.py` | `results/tables/exemplar_*.csv` | Per-cell boundary + marker-molecule coordinates for the panel d/h exemplar cells (run **before** 09) |
| `scripts/09_composite_figure.R` | `results/figures/09_composite.{png,pdf}` | **Publication composite (7.1 × 6.625 in, 10 panels a–j)**, one marker per row: **Row 1 SST in Sst** — volcano (a), SST/Sst forest (b), CP1K boxplot (c), exemplar cells (d); **Row 2 PVALB in Pvalb** — volcano (e), PVALB/Pvalb forest (f), CP1K boxplot (g), exemplar cells (h); **Row 3** — DE-gene-count butterfly per subclass with a DE-genes-vs-proportion inset (i), concordance scatter (j) |
| `scripts/11_grain_density.py` | `results/tables/percell_grain_density.csv` | Per-cell grain-density input (canonical cells, 24 donors) → panel d/h exemplar selection + Supplementary |
| `scripts/12_marker_norm_expr.R` | `results/tables/marker_norm_expr.csv` (+ `_stats`) | Per-donor CP1K (counts/1,000 transcripts) + edgeR p for SST (**panel c**) and PVALB (**panel g**); FGFR3 also computed (former composite gene, now spare/supplement) |
| `scripts/13_supp_percell_metrics.R` | `results/figures/S_percell_metrics.{png,pdf}` | **Supplementary**: per-cell SST/PVALB across normalisations (raw, grains/cell-area, lib-norm, library size) |
| `scripts/14_supp_pvalb.R` | `results/figures/S_pvalb.{png,pdf}` | **Supplementary (superseded)**: PVALB mRNA in Pvalb cells — forest + CP1K + exemplar cells; now redundant with composite Row 2 (e–h) |
| `scripts/16_de_vs_proportion.R` | `results/figures/de_vs_proportion_subclass.{png,pdf}` | Standalone fully-labelled DE-genes (FDR<0.10) vs cell-type-proportion scatter (subclass; Xenium proportion proxy) — minimal version is the panel-i inset |

The validation figures (07, 08, and the cross-platform panels of 09 — forests
b/f, CP1K c/g, scatter j) all ask the same question
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

### C. Xenium spatial DE — `../spatial/output/de/de_results_subclass.csv`
- **What**: pseudobulk DE on Xenium spatial transcriptomics, one row per
  (cell_type × gene). Read via the monorepo-internal path `../spatial/output/de/`
  (a git-ignored symlink to `~/Github/SCZ_Xenium/output/de/`). Produced by
  `spatial/code/analysis/{build_de_input.py,run_de.R}` (pseudobulk + edgeR QL):
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

### D. Panel d/h exemplar tables (derived) — `results/tables/exemplar_*.csv`
- **What**: per-cell boundary polygons + marker-gene transcript-molecule
  coordinates (µm, recentred), plus `exemplar_cells_meta.csv`. Produced by
  `scripts/10_xenium_exemplar_cells.py` from the Xenium **h5ad + boundary +
  transcript** exports (NOT the DE table). See §6 for full provenance.
- **Used by**: script 09 panels d/h only. **Run script 10 before script 09.**

---

## 2. Shared methods

> **Script 09 reuses everything in this section.** Its forests (b, f), scatter
> (j), volcanoes (a, e) and butterfly + inset (i) read the same three DE inputs and use
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
cd ~/Github/scz_celltype_paper/transcriptomic
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

**Current numbers** (padj < 0.10): 166 pairs, 61 both-up / 65 both-down /
40 discordant, 76% directionally concordant, Pearson r = 0.73, slope = 0.79,
binomial p = 7e-12. At padj < 0.05: 109 pairs, 78% concordant, r = 0.78.

---

## 5. Figure 09 — publication composite

The single multi-panel figure for the paper
(`results/figures/09_composite.{png,pdf}`, **7.1 × 6.625 in, 400 dpi**). Width
`FIG_W = 7.1` in (NN double-column) and height `FIG_H = 6.625` in. It reuses the
inputs (§1) and methods (§2) of 07/08 and adds the butterfly + inset and exemplar
cells. Layout = **one interneuron marker per row** (volcano → forest → CP1K boxplot
→ exemplar cells) for SST (a–d) and PVALB (e–h), then a transcriptome-wide row (i, j).

| Panel | Builder | Content |
|---|---|---|
| a / e | `build_volcano()` | Volcano for Sst (a) and Pvalb (e); **y = −log₁₀ P** (raw meta p), points coloured by up/down × FDR tiers (NS = grey); dashed line at the FDR<0.10 p-boundary (`max(pval[padj<0.10])` — BH-monotone, so every coloured point sits above it); selected genes via ggrepel (subtitle `<cell> cells`; the cell-type marker SST/PVALB is **bold**-italic, others italic). **Tight, per-panel data-driven x/y limits** (not symmetric, not shared); labels kept inside via ggrepel `xlim`/`ylim` + `coord_cartesian`. a Sst = SST, NAT16, SMAD1, AFG3L2, STAC, KCTD4, SLC9A9, DRD3; e Pvalb = PVALB, SMAD1, ANXA2, SCN3A, NAT16, VGF, CIRBP, TCAF2, FGF10. |
| b / f | `build_forest()` | Forest for SST/Sst (b) and PVALB/Pvalb (f) (the `FOREST` tribble); 7 cohort squares + DL meta diamond + Xenium triangle; markers per §2; `SCZ log₂ FC` x-axis on both. |
| c / g | `build_normexpr()` | Per-donor library-normalised expression (CP1K, counts/1,000 tx) boxplots, Control vs SCZ; SST-in-Sst (c), PVALB-in-Pvalb (g); y-axis `<gene> expr (CP1K)`; **edgeR** p via ggsignif (same DE as b/f, j). Input from `scripts/12`. |
| d / h | `build_exemplar()` | Xenium exemplar cells (inputs from script 10): cell boundary + dashed nucleus + red marker molecules (count top-left); SST/Sst (d), PVALB/Pvalb (h), Control \| SCZ per row; shared `ex_lim` zoom + 5 µm scale bar. |
| i | `build_butterfly()` + `build_de_prop_inset()` | Up/down DE-gene counts per subclass; FDR<0.10 (light) + FDR<0.05 (dark); 4-entry legend at `c(0.70,0.16)`. **Inset** in the empty bottom-left (`draw_plot`, x=0.135 y=0.085, 0.36×0.378): # DE genes (FDR<0.10) vs mean per-donor Xenium proportion (log10 %, ≥1-DE subclasses); points by class; Astro/L5 IT/Vip/L6b labelled; Spearman ρ. |
| j | `build_scatter()` | meta logFC vs Xenium logFC over meta FDR<0.10 ∩ Xenium (n=166); colour = class, size = meta FDR tier (legends removed; in figure legend); **identity line only — no fit line**; 10 ggrepel-labelled pairs; `r` / `% concordant` text bottom-right. |

**Layout** (cowplot): rows 1–2 each = [volcano | forest | CP1K boxplot]
(`align="h", axis="tb"`, `rel_widths RW3 = c(1.05, 0.95, 0.70)`) then the exemplar
Control\|SCZ pair as a 4th column (`rel_widths c(sum(RW3), 1.25)`); row 3 = butterfly
i + scatter j (50/50). `rel_heights c(2, 2, 2.625)` are **inches** (rows 1–2 = 2.0 in
each, row 3 = 2.625 in), summing to `FIG_H = 6.625`. Panel letters 8 pt bold
lowercase (a–j), top-left.

**Conventions locked (Nature / Nature Neuroscience print spec)** — change these
together, not piecemeal:
- **Size**: 7.1 × 6.625 in (NN double-column width). Rows 1–2 = 2.0 in each, row 3 = 2.625 in (`rel_heights` are inches).
- **Text size**: master `BASE = 7`; all figure text falls in **5–7 pt** (the
  Nature/NN requirement). `geom_text` multipliers are tuned so nothing exceeds
  7 pt or drops below 5 pt; panel letters are 8 pt bold lowercase. Font =
  Helvetica (PDF device default).
- **One marker per row**: row 1 = SST in Sst cells (a–d), row 2 = PVALB in Pvalb
  cells (e–h); the volcano/forest/CP1K triplet is x-aligned (`align="h", axis="tb"`).
- **Forests**: 2 panels (b SST/Sst, f PVALB/Pvalb); `SCZ log₂ FC` x-axis on both
  (`show_xlab = TRUE`). Panel set = the `FOREST` tribble. (Earlier drafts used a
  single four-forest row — SST/Sst, BDNF/L2_3 IT, FGFR3/Astro, FKBP5/OPC — before
  the reorg to one marker per row; BDNF/FGFR3/FKBP5 now appear only as scatter labels.)
- **Volcano highlights**: edit the `build_volcano()` calls — a Sst = SST, NAT16,
  SMAD1, AFG3L2; e Pvalb = PVALB, SMAD1, ANXA2, SCN3A, NAT16. **Tight, per-panel
  data-driven limits** (not symmetric, not shared); to widen for labels, raise the
  `xpad` / `yhi` multipliers in `build_volcano()`. (The earlier glial volcanoes —
  Astro SERPING1/CHI3L1/FGFR3/NOTCH1 and Micro-PVM C1QA/C1QB/CX3CR1/P2RY12/SORL1 —
  are no longer in the main figure after the reorg.)
- **Scatter j labels use ggrepel** (`force = 5`, fixed `seed = 7`; ATP2B4/Sst is
  `nudge`d up-left so it does not collide with CALB1/L6 IT); cell-class + meta-FDR
  legends **removed from the panel** (in the figure legend); `r` / `% concordant`
  lower-right. **No fit line** (dashed identity line only). Labelled pairs (10) =
  the `lab_pairs` tribble: SST/Sst, BDNF/L2_3 IT, FKBP5/OPC, CX3CR1/Micro-PVM,
  SMAD1/Pvalb, SERPING1/Astro, FGFR3/Astro, VGF/Chandelier, CALB1/L6 IT, ATP2B4/Sst.
  To add/remove one, edit `lab_pairs`.
- **Panel-i inset**: `build_de_prop_inset()` (reads `INPUT_CRUMBLR`) draws # DE
  genes (FDR<0.10) vs mean per-donor Xenium proportion across ≥1-DE subclasses;
  minimal styling (semi-transparent fit, 1%/10% log stops, ρ, no point labels except
  Astro/L5 IT/Vip/L6b). Standalone fully-labelled version = `scripts/16`.
- **Butterfly legend** sits at `c(0.70, 0.16)` (right of the bar-tip count labels),
  keys `unit(9,"pt")`.
- **Exemplar pair (d/h)**: `Control` / `SCZ` column headers on **row 1 only** (a
  blank spacer of equal height keeps the row-2 cells the same size). Each cell draws
  the solid boundary, the **dashed nucleus**, and red marker dots (count top-left).
  All four cells share one coordinate limit (`ex_lim` = max boundary extent × 1.05)
  so the zoom and the **5 µm scale bar** are identical; the `5 µm` text is on the
  h (PVALB/SCZ) cell only.

**Run** (script 10 must have produced the exemplar tables first):
```bash
cd ~/Github/scz_celltype_paper/transcriptomic
python scripts/10_xenium_exemplar_cells.py   # panel d/h exemplar inputs (if not present)
Rscript scripts/09_composite_figure.R         # -> results/09_composite.{png,pdf}
cp results/09_composite.png results/09_composite.pdf results/figures/  # snapshot
```

---

## 6. Script 10 — Xenium exemplar cells (panel d/h inputs)

Extracts, for each (gene, subclass) pair, one **Control** and one **SCZ**
exemplar cell, writing the cell boundary + the marker-gene transcript molecules
inside it, for panels d/h.

- **Pairs** (`PAIRS`): `(gene, subclass, {dx: section})` triples — `SST`/Sst (panels
  c/d) and `PVALB`/Pvalb (panels g/h) for the composite, plus `FGFR3`/Astrocyte (a
  former composite gene, kept as a spare / supplement candidate). SST and PVALB =
  down-regulated interneuron markers; FGFR3 = down-regulated astrocyte-identity gene.
  All render as countable molecule dots. (Caveat: PVALB partly *defines* the Pvalb
  type, so its Xenium effect is mildly circular with label transfer — see §8.)
- **Samples** (per pair, in `PAIRS`): the donor whose per-section canonical-cell
  grain-density median is closest to that diagnosis's pooled group median. SST =
  `Br6432` (Control) / `Br5973` (SCZ); FGFR3 = `Br5400` / `Br5973` (`Br6432` is
  atypically LOW in FGFR3); PVALB = `Br6432` / `Br5973`. Cell + nucleus boundaries
  are exported for all 24 donors; transcript molecules for six (Br2039, Br5400,
  Br5746, Br5973, Br6432, Br8667), so drawing sections are chosen among those —
  `Br5400` was re-exported to serve as the representative FGFR3 control. Diagnosis
  from SCZ_Xenium `config.SAMPLE_TO_DX`.
- **Raw inputs** (SCZ_Xenium repo):
  - `output/h5ad/<sample>_annotated.h5ad` — raw counts (`.X`), `subclass_label`,
    `qc_pass`, `total_counts`, `predicted_norm_depth`, `layer`.
  - `output/deploy/boundaries/<sample>.json` + `<sample>_nucleus.json` — per-cell
    and per-nucleus polygons (25 verts; decode `micron = quant*scale + offset`;
    same index = h5ad obs index).
  - `output/deploy/transcripts/<sample>/<GENE>.json` (+ `gene_index.json`) —
    per-gene molecule coordinates.
- **Canonical cells** (`canonical_mask`): `corr_subclass` + cortical +
  (`qc_pass` & `corr_qc_pass`) — identical to the crumblr / edgeR-DE cell set.
- **Group target** (`pooled_grain_median`): median GRAIN DENSITY (marker dots per
  cell area, grains/100 µm²) over **all 24 donors'** canonical cells, per
  diagnosis — the Dienel-style target the exemplar aims at (SST 18.7→12.7,
  FGFR3 7.3→6.3, PVALB 5.2→4.3 Ctrl→SCZ). Cell area from the deploy boundary polygons.
- **Shape anchors** (`shape_anchors`): median cell area (size-match) and median
  eccentricity (typical outline) over the pair's two drawing sections' canonical cells.
- **Cell selection** (`pick_cell`): among canonical, **size-matched** (area within
  ±20% of the median) and convex (solidity `A/hull ≥ 0.93`) cells, pick the cell
  whose grain density is closest to the group target AND whose **eccentricity is
  closest to the median outline** (a typical-shaped cell, not the roundest); a
  final stage matches the *displayed* in-polygon grain density to the target so
  the drawn dots are representative and the Control>SCZ direction holds. Depth/
  layer is **reported, not used**. Molecules clipped by point-in-polygon.
  Grain-density input is produced by `scripts/11_grain_density.py`.
- **Outputs**: `exemplar_<gene>_<dx>_{boundary,nucleus,dots}.csv` (recentred µm
  coords) and `exemplar_cells_meta.csv` (sample, `cell_index`, `grain_density_target`,
  `ecc_target`, `marker_count`, `n_dots_in_poly`, `grain_density`,
  `disp_grain_density`, `eccentricity`, `total_counts`, `norm_depth`, `layer`,
  `circularity`, `solidity`, `area_um2`).
- **Current exemplars** (from `exemplar_cells_meta.csv`): dots shown
  (`n_dots_in_poly`) SST/Sst 35 (Ctrl) vs 28 (SCZ) — **panel d** — and PVALB/Pvalb
  10 vs 9 — **panel h**; FGFR3/Astrocyte 14 vs 12 (spare). Grain density 18.7→13.2
  (SST), 5.2→4.2 (PVALB), 7.1→6.5 (FGFR3) grains/100 µm²; eccentricity 0.51–0.62
  (≈ median ~0.57); size-matched (area within ±20%).
- **To change the exemplar genes**: edit `PAIRS` (and, for cross-referencing,
  `FOREST` / `lab_pairs` in script 09 and the `ex_pair` block), re-run 10 → 09.

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
   that panels d/h read.
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
- **Panel-l exemplars are illustrative, not quantitative**: each is the single
  representative cell per group at the pooled group-median grain density
  (size-matched, typical eccentricity). At the median the per-cell contrast is
  modest — SST 35 vs 28, FGFR3 14 vs 12 dots (PVALB 10 vs 9 in the supplement) —
  because the median downregulation is itself modest; the quantitative
  down-regulation is carried by the DE / forest / scatter panels, not by these
  cells. The
  molecule count shown (`n_dots_in_poly`) runs slightly above the cell's assigned
  `.X` count because the boundary captures some unassigned/neighbour molecules.
- **PVALB defines the Pvalb type**, so its Xenium signal is mildly circular with
  label transfer (see the marker-circularity caveat above). A per-cell
  negative-binomial mixed model (donor random intercept) shows the PVALB
  reduction is modest and borderline — raw dots/cell 0.86× (p=0.065), and it does
  NOT strengthen after library-size normalisation (0.87×, p=0.060, n.s.) —
  consistent with the modest edgeR effect (logFC −0.22, p=0.044). The SST
  reduction, by contrast, is robust (raw 0.70×, p=0.002; library-normalised
  0.76×, p=0.011). Values: `results/tables/S_percell_stats.csv` (`scripts/13`).
- **Five Xenium sections** have transcript-level molecule export (Control:
  `Br6432`, `Br8667`; SCZ: `Br2039`, `Br5746`, `Br5973`), so exemplars are drawn
  from those; `Br6432`/`Br5973` are used as the most group-representative pair.
