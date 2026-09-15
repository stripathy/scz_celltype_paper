# Patch-seq raw data for Figure 4 panels e and f (morphology, traces)

The minimum needed to regenerate those two panels, and nothing else: the
**five exemplar cells**, each with one morphological reconstruction and one
intracellular recording. Held here so panels e and f reproduce without a
checkout of `human_int_patch_seq`, which is otherwise a hard dependency of the
export step.

| cell | supertype | layer | soma depth | files |
|---|---|---|---|---|
| 1079568285 | Sst_20 | L2 | 304 µm | `swc/1079568285_upright.swc` (0.8 MB), `nwb/1079568285.nwb` (31.3 MB) |
| 819770858 | Sst_25 | L2 | 406 µm | `swc/819770858_upright.swc` (1.0 MB), `nwb/819770858.nwb` (19.6 MB) |
| 1037461069 | Sst_3 | L3 | 1109 µm | `swc/1037461069_upright.swc` (1.3 MB), `nwb/1037461069.nwb` (31.0 MB) |
| 758996755 | Sst_5 | L4 | 1460 µm | `swc/758996755_upright.swc` (0.4 MB), `nwb/758996755.nwb` (35.6 MB) |
| 797048104 | Sst_1 | L5 | 1734 µm | `swc/797048104_upright.swc` (1.1 MB), `nwb/797048104.nwb` (30.4 MB) |

SWCs copied byte-identical from `~/Github/human_int_patch_seq/data/morphology/swc/`;
NWBs from `~/Github/human_int_patch_seq/data/nwb_cache/` where cached, otherwise
from DANDI dandiset 000636 (1037461069 = `sub-1036289507_ses-1037460966_icephys.nwb`;
1079568285 = `sub-1078812471_ses-1079568103_icephys.nwb`).

Note 1079568285 is the **only** reconstructed Sst_20 cell in the dataset, which
is why Sst_20 rather than Sst_22 is the upper-layer depleted exemplar in panel
e. The published legend names Sst_20, Sst_25 and Sst_3.

## What each is for

- **SWC** — the reconstruction nodes drawn in panel e. `export_panels_ef.py` parses
  them, flips the y-axis for specimens in `INVERTED_SPECIMEN_IDS`, and writes
  `panel_F_morphology.csv` + `panel_F_meta.csv`.
- **NWB** — panel f picks the single hyperpolarising sweep nearest each cell's
  target current (−100 pA for Sst_25, −30 pA for Sst_5) from the
  `X1PS_SubThresh` / `X3LP_Rheo` / `X4PS_SupraThresh` stimulus sets, and writes
  `panel_G_traces.csv` + `panel_G_meta.csv`. Only that one sweep per cell is
  used; the rest of each NWB is unused but kept so the sweep-selection step
  itself stays reproducible rather than being frozen as a hand-picked extract.

## Storage

The NWBs are git-ignored (binary, 148 MB); the five SWCs listed above are plain
text and are committed, so panel e re-exports from a clone. Either way the *rendered* panels never read this directory —
they read the committed `results/figures/r_panels/panel_[FG]_*.csv`, which is
what makes Figure 4 reproducible from a clone. This directory is what makes the
extraction step reproducible too, and it is what
`r_panels/MANIFEST.tsv` checksums so those panels can no longer go stale
unnoticed.

`export_for_R.py` prefers this directory and falls back to the upstream repo
when it is absent, so nothing breaks on a machine that has one but not the other.

## Still external

The *code* for the extraction — `patchseq_builder.morphology.download.parse_swc`
and `.orientation.{INVERTED_SPECIMEN_IDS, flip_swc_y}` — still comes from
`human_int_patch_seq`. Only the data is vendored here. Re-running the export
therefore still needs that package importable; rendering the figure does not.
