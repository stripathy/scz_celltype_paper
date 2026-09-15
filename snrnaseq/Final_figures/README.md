# Final_figures/ — Fig. 1a, Fig. 3, and Supplementary S1, S4, S5, S7

Stage 4. The renderers for every figure this component owns.

Unlike the renderers elsewhere in this repo, these **do not write into
`manuscript/figures/`** — they write to a `Figures/` directory on the cluster
under their own names. So none of the figures below has a rendered output in
this repo, and the output filenames do not carry their figure numbers. Issue 12
in [`../../KNOWN_ISSUES.md`](../../KNOWN_ISSUES.md).

## What builds what

| Figure | Script | Writes |
|---|---|---|
| **Fig. 1a** (UMAPs) | `Code/Figure_1.py` — **this is Python**, see below | `combined_umap_spatial.{png,svg}`, `combined_sst_umap_predicted_id.png` |
| **Fig. 1a** (counts) | `Figure_1a_Barchart.r` | `1a.{png,svg}`, `cohort_donor_counts.png`, `snrnaseq_xenium_cell_counts.png` |
| **Fig. 3** | `Figure_3.r` | `Figure3_composite.{png,svg}` plus each panel standalone |
| **S1** | `Supplemental/Markerplot.py` | `four_compiled_cohort_markers.{png,…}` |
| **S4** | `Supplemental/Supertype_DE.r` | `DE_abundance_and_SST_forest.{png,svg}`, `supertype_DE_vs_abundance.png` |
| **S5** | `Supplemental/NonNeuron_supplement.r` | `Supplement_nonneurons.png` |
| **S7** | `Supplemental/Sensitivity_analysis_Barchart.r` | `Barchart_SCZ_meta_NoSSTDE.png` |

`1_Gene_symbols.r` is a prerequisite, not a figure: it maps Ensembl IDs to gene
symbols in the PsychAD objects (`org.Hs.eg.db`, unmapped IDs kept as-is, made
unique) and writes `*_symbols.rds`. The UMAP and marker-plot scripts read those.

## Figure 3 panel map

The panel variables are **not** named after their final panel letters. The
assembly at the foot of `Figure_3.r` is what decides:

| Panel | Variable | Shows | Main input |
|---|---|---|---|
| a | `p3a` | Abundance change per neuronal supertype, meta-analysed | `crumblr_results_final_meta.csv` |
| b | `p3c` | Sst_25 forest across the seven datasets | `final_results_crumblr_7_cohorts_all.csv` |
| c | `p3b` | Sst_25 per-donor proportion, Control vs SCZ, with Xenium | `neuron_props_7_cohorts` + Xenium proportions |
| d | `p3d` | Spatial rendering of Sst_25 | `Xenium_SCZ_R.rds` |
| e | `p3h` | Abundance change vs cortical depth | Xenium crumblr + depth |
| f | `p3f` | snRNA-seq vs Xenium concordance (Spearman) | Xenium crumblr + 7-dataset meta |

## Cross-platform inputs

Panels c–f read the Xenium results directly from
`/scratch/nendresz/Xenium/…` rather than through this repo's
`spatial/output/crumblr/`. They are meant to be the same results, but the wiring
is a copy on the cluster, not the seam.

Panels d, e and h read Xenium cell metadata through `load_xenium_meta()`. They
use only metadata and the x/y centroids — no expression matrix — so since
2026-09-15 the function prefers **`Xenium_SCZ_R.rds`** when it is present (so
cluster runs are unchanged) and otherwise falls back to the committed
**`Data/xenium_metadata.csv`**, which makes these panels render from a clone.

`Xenium_SCZ_R.rds` itself is a Seurat conversion that no script in this repo
creates and that is not committed. Checked on 2026-09-15: it does derive from the
canonical object, **not** from the superseded one with MGE and CGE transposed. The metadata export beside it (`Data/xenium_metadata.csv`) holds
1,338,922 cells against the canonical 1,339,151 — a gap of exactly the 229
duplicate cell names Seurat de-duplicates on conversion — with all 24 samples,
Br2039 included, and subclass counts agreeing to within 0.02%.
`Data/xenium_crumblr_results_supertype_neuronal.csv` is byte-identical to
`spatial/output/crumblr/crumblr_results_supertype_neuronal.csv`
(md5 `97c546a39016521db5ca58e92afd84ab`). Closed as issue 3 in
[`../../KNOWN_ISSUES.md`](../../KNOWN_ISSUES.md).

## `Code/Figure_1.py` is Python

It opens `from brisc import SingleCell` and runs under `mamba activate brisc`.
It parses as valid Python and fails to parse as R; only the extension is wrong.
Issue 8 in [`../../KNOWN_ISSUES.md`](../../KNOWN_ISSUES.md).

It loads all seven datasets through `brisc`, applies the same donor filters as
the composition analysis plus a **≥ 500 cells per donor** cut and `qc()`, labels
each dataset, and draws the combined UMAP, an Sst-only UMAP coloured by
`predicted.id`, and the spatial panel, into one composite.

`Supplemental/Markerplot.py` is likewise `brisc`-based, and draws a 2×2: marker
expression across the seven datasets and in the SEA-AD reference, at subclass
level (top) and Sst supertype level (bottom), with `Sst_25` highlighted. It also
writes the underlying dot-plot tables as CSVs beside the figure.

## Shared assets

Every R script here reads `cluster_order_and_colors.csv` for cell-type ordering
and the palette — not committed, and load-bearing beyond colour (issue 4).
