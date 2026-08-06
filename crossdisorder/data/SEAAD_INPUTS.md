# SEA-AD DLPFC (A9) — input data for Figure 4j

The AD axis of Figure 4j is the compositional slope of each supertype along the
SEA-AD **continuous pseudo-progression score (CPS)** in dorsolateral prefrontal
cortex. Two inputs are needed; only the small one is committed.

## 1. Per-nucleus metadata (not committed, ~1.3 GB)

| | |
|---|---|
| **File** | `SEAAD_A9_RNAseq_final-nuclei_metadata.2024-02-13.csv` |
| **Size** | 1.3 GB — 1,395,601 nuclei × 133 columns |
| **Source** | <https://sea-ad-single-cell-profiling.s3.amazonaws.com/index.html#PFC/RNAseq/previous_objects/> |
| **Downloaded** | 2026-08-04 |
| **Citation** | Gabitto et al. 2024, *Nature Neuroscience* ([doi:10.1038/s41593-024-01774-5](https://doi.org/10.1038/s41593-024-01774-5)) |
| **Taxonomy** | `CCN20230505` (class / subclass / supertype) |

Only the metadata CSV is needed — the compositional analysis never touches the
36 GB expression matrix. Place it anywhere and pass `--metadata PATH`, or leave
it in `~/Downloads/` (the script's default).

### Why the 2024-02-13 release and not the 2026 one

The S3 bucket also serves a newer `SEAAD_DFC_RNAseq_final-nuclei.2026-06-22`
release under `PFC/RNAseq/`. **Do not use it for this figure.** It re-annotates
against an expanded taxonomy — 152 supertypes, 15 of them novel `-SEAAD` disease
states — which redistributes cells relative to the 137-supertype taxonomy used
throughout this paper and attenuates the supertypes those cells were drawn from.
The `previous_objects/` 2024-02-13 release shares its taxonomy vintage with the
MTG data used elsewhere here: its `Supertype` column contains 131 labels, all of
which exist in the MTG set (zero novel), including all 16 Sst supertypes.

Measured effect of the choice, on the panel-4j correlation across the 16 Sst
supertypes: **2024-02-13 ρ = 0.87 (94% sign-concordant)** vs 2026-06-22 ρ = 0.79
(81%).

Note the release also ships `Supertype (non-expanded)`. That column is *not* the
right one — it collapses the expansion states and drops several types, leaving
123 labels that match the MTG set less well than `Supertype` does.

## 2. Donor CPS lookup (committed, 2 KB)

`seaad_donor_cps.csv` — `donor, CPS` for the 84 SEA-AD donors that have a score.

CPS is **donor-constant** and derived from quantitative neuropathology measured
in MTG; it is not carried in the A9 per-nucleus metadata, so it is joined by
donor. Extracted from `SEAAD_MTG_RNAseq_obs_lean.csv`
(column `Continuous Pseudo-progression Score`), itself derived from
`SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad` — see `../README.md` §1.
Committed here so the pipeline does not depend on an external checkout.

Range 0.150–0.929. Five MTG donors and three A9 donors have no CPS: these are
the Allen/BICCN neurotypical reference brains (`H18.30.002`, `H19.30.001`,
`H19.30.002` in A9), which are not on the AD trajectory and are dropped.

## Running it

```bash
python3 genetics/scripts/14_seaad_dlpfc_crumblr_input.py
Rscript genetics/scripts/14_seaad_dlpfc_crumblr.R
```

Writes to `genetics/results/intermediates/seaad_dlpfc/`. The second step needs
R with `crumblr`, `dreamlet`, `variancePartition` and `limma`.

## Cohort after filtering

80 donors (83 minus the three reference brains), 1,003,866 neuronal nuclei,
107 of 108 neuronal supertypes passing the ≥50%-donor-presence filter. The nine
severely affected donors are **retained**, following Gabitto et al., who exclude
them only from gene-expression tests and report consistent compositional results
with and without them. Excluding them here gives ρ = 0.77 rather than 0.87.
