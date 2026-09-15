# Sst depletion strata — the Supplementary Fig. S8 pipeline

Self-contained pipeline for the paper's final main figure: how the SCZ
transcriptional state of Sst interneurons differs between the supertypes that are
depleted (Fig. 3) and those that persist.

Run every script from the **repo root**. Each sources `_common.R`.

```bash
Rscript transcriptomic/scripts/sst_strata/01_strata.R            #  <1 min
Rscript transcriptomic/scripts/sst_strata/02_stratum_de.R        # ~10 min  (8 cores)
Rscript transcriptomic/scripts/sst_strata/03_stratum_gsea.R      #  ~8 min
Rscript transcriptomic/scripts/sst_strata/04_donor_pseudobulks.R #  ~3 min
Rscript transcriptomic/scripts/sst_strata/05_interaction.R       # ~25 min
Rscript transcriptomic/scripts/sst_strata/06_module_scores.R     #  ~4 min
Rscript transcriptomic/scripts/sst_strata/07_xenium_stratum.R    #  <1 min
Rscript transcriptomic/scripts/sst_strata/08_figure.R           #  <1 min
Rscript transcriptomic/scripts/sst_strata/09_verify.R            #  <1 min  ← must PASS
```

## Dependency graph

```
                     external inputs
   crumblr 7-cohort meta ──┐
   Xenium supertype depth ─┴─► 01_strata ──► strata_definition.csv
                                                 │ (single source of truth for
                                                 │  stratum membership)
   stratum_pseudobulks_export/ ──┬──────────────►┤
     (7 cohorts, parquet+groups) │               │
                                 ├─► 02_stratum_de ──► stratum_percohort_de.csv
                                 │                     stratum_meta_de.csv ──┐
                                 │                                          │
                                 │                     03_stratum_gsea ◄─────┤
                                 │                       ├─► gsea_all_signatures.csv
                                 │                       └─► stratum_gene_signatures.csv
                                 │                                          │
                                 └─► 04_donor_pseudobulks ──► donor_stratum/ │
                                          │                                 │
                                          ├─► 05_interaction ──► interaction_{percohort,meta,gsea}.csv
                                          └─► 06_module_scores ─► module_scores_donor.csv
                                                                  module_score_tests.csv
   spatial supertype DE ──────────────────► 07_xenium_stratum ──► xenium_stratum_concordance.csv

   08_figure  ◄── strata_definition + gsea_all_signatures + stratum_gene_signatures + subclass DE
       └─► manuscript/figures/supplementary/S08_sst_strata.(png|pdf)

   09_verify   ◄── everything above; asserts every number quoted in the draft
```

`02` and `03` collapse each stratum to **one pseudobulk per donor**; `04` keeps
**one pseudobulk per donor per stratum**, which is what makes the within-donor
interaction test in `05` possible. Both aggregations come from the same export.

## Outputs

All under `transcriptomic/results/sst_strata_gsea/`, **except the figure**, which
is written directly into its submission location so there is no copy to drift:

- `manuscript/figures/supplementary/S08_sst_strata.{png,pdf}` -- laid out at
  main-figure size, so it can be promoted; the S-number is the `FIGSTEM`
  constant at the bottom of `08_figure.R`
- the controls figure (nuclei imbalance, cell-matched draws, interaction NES) is
  **not in the current manuscript**; `reserve/sst_strata_supp/figS_strata_controls.R`
  still renders it as `reserve_strata_controls.{png,pdf}` under results/ if a
  reviewer asks. The supplements and sensitivity analyses built on this pipeline
  live in `reserve/sst_strata_supp/`; they source this directory's `_common.R`
  and run from the repo root.

| File | Written by | Used for |
|---|---|---|
| `strata_definition.csv` | 01 | stratum membership everywhere; panel a |
| `pseudobulk/stratum_meta_de.csv` | 02 | gene-level meta statistics (canonical) |
| `pseudobulk/stratum_percohort_de.csv` | 02 | leave-one-cohort-out sensitivity |
| `pseudobulk/gsea_all_signatures.csv` | 03 | panels b, d; module definitions |
| `pseudobulk/stratum_gene_signatures.csv` | 03 | panels c, e |
| `pseudobulk/donor_stratum/*` | 04 | inputs to 05 and 06 |
| `pseudobulk/interaction_{meta,gsea}.csv` | 05 | the "graded" claim |
| `pseudobulk/module_score_tests.csv` | 06 | per-donor module statistics |
| `xenium_stratum_concordance.csv` | 07 | cross-platform check; panel coverage |
| *(figure output moved — see below)* | 08 | the figure |

## Figure sizing (matches Figure 4)

Figure 4 is built at **8.0 in wide and placed in a 7.1 in column**, so its text
prints at its nominal size x 7.1/8.0. S8 uses Figure 4's constants exactly,
via `FIG_SCALE <- 8.0/7.1` in `_common.R`. Note Figure 4's own `FIG_WIDTH_IN`
variable says 6.5, but its exported PNG is 8.00 x 7.46 in -- the 8.0 is what
matters and what `FIG_SCALE` encodes.

| | Figure 2 | Figure 4 | S8 |
|---|---|---|---|
| built at | 7.1 in | 8.0 in | 8.0 in |
| printed at | 7.1 in | 7.1 in | 7.1 in |
| dpi | 400 | 400 | 600 |
| base font | 7 pt | 9.01 -> **8.0 printed** | 9.01 -> **8.0 printed** |
| axis titles | 6.5 | 8.45 -> 7.5 | 8.45 -> 7.5 |
| tick labels | 6.0 | 7.89 -> 7.0 | 7.89 -> 7.0 |
| panel letters | 8 | 10.14 -> 9.0 | 10.14 -> 9.0 |
| smallest text | — | 7.18 -> 6.4 | 7.18 -> 6.4 |

Figure 4 and S8 therefore run about 1 pt larger in print than
Figure 2, which is a deliberate choice recorded in `fig4_style.R`. Nature
guidance is 5-7 pt for final text; tick labels land exactly at 7.0.

In-panel text uses fig4_style.R's ratios (`LBL_SMALL` 0.28, `LBL_GENE` 0.32,
`LBL_CALL` 0.34, `LBL_STAR` 0.43 of the base). ggplot's `size` is millimetres, so
pt = size x 2.845. `LW` and `GEO` scale linewidths and point sizes with the base.

Canvas is **8.0 x 9.8 in**, printing as 7.1 x 8.7 — a full-page figure, driven by
panel e's 31 gene rows at 8.2 pt.

**If you change the canvas or the base font, re-check four things** that are
hand-tuned to the current geometry: the three stratum labels in panel a (literal
data coordinates; the full "Depleted (5)" forms collide at this width, which is
why the counts were dropped), the ggrepel force and padding for the 16 supertype
labels in a, the manual gene-label coordinates in panel c, and the left gutter
expansion in panel e (`add = c(1.75, 0.55)`, sized for the longest gene name).

## Conventions (all in `_common.R`)## Conventions (all in `_common.R`)

- **Strata palette** purple → lilac → green. Deliberately outside the
  blue/orange direction palette and the magenta/teal module palette, so panel a
  cannot be misread as encoding direction.
- **Direction palette** Figure 2's blue (down) / orange (up), with a light tier
  for the FDR 0.05–0.10 band.
- **Module palette** synaptic (steel), translation (magenta), OxPhos (orange),
  deubiquitination (teal), other (grey).
- **Module definitions** derive from the single `BLOCKS` table. Add a gene set to
  a block there and every panel, module score and supplement follows.
  `module_of()` resolves genes in several modules by `MODULE_PRIORITY`.
- **GSEA** GO BP/CC/MF + Reactome (Hallmark deliberately excluded), size 10–500, `nPermSimple = 10000`,
  ranked on z = estimate/SE. Collections are cached to `.cache/msigdb_sets.rds`.
- **Meta-analysis** `metafor` REML with a DL fallback, k ≥ 5 cohorts, BH-FDR
  within stratum (or within interaction coefficient).

## Notes and gotchas

- **`_common.R` exists because these definitions used to be duplicated.** The
  module gene lists were re-derived in nine scripts and the strata membership
  hard-coded in three. That is how panel c once ended up mislabeled. Import; do
  not re-list.
- **Headline burden numbers are the FDR < 0.10 tier** (130 / 64 / 0 / 55 gene
  sets), matching the labels drawn in panel b. The FDR < 0.05 tier is the dark
  bar underneath.
- **Do not describe our own GSEA results as an integrated stress response.**
  Several top interaction sets have stress-flavoured Reactome names (GCN2 /
  amino-acid deficiency, cellular response to starvation, SRP, NMD) but are
  55–83 % ribosomal-protein genes — the ribosome signal under other labels. EIF2
  belongs in the Discussion as cited literature only.
- **`07_xenium_stratum.R` reconstructs a file that previously had no producer.**
  The recipe (unweighted mean logFC across a stratum's supertypes; per-supertype
  z combined by Stouffer) was recovered by testing candidates against the
  orphaned CSV and reproduces it to 0.0000 on both columns.
- **`archive_ivw/`** holds the superseded inverse-variance-weighted shortcut used
  while waiting for the real pseudobulks. Kept only for the framework-robustness
  supplement, which correlates it against the pseudobulk result. Nothing in the
  pipeline reads it. (An earlier version of that comparison read the top-level
  copies, which by then held pseudobulk output — it was comparing the pseudobulk
  result with itself. Point any such comparison at `archive_ivw/`.)
- `../../../reserve/sst_strata_supp/` holds the supplement and sensitivity scripts, none of them in the paper; see
  `reserve/sst_strata_supp/README.md` for what each one showed and why it was held back.
