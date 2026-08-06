# crossdisorder/ — is the SCZ-depleted Sst population also depleted in Alzheimer's?

Self-contained component, in the same sense as `histology/`: it owns one external
dataset and the analysis on it, and exports a single result that other components
consume.

**The question.** Figure 4j asks whether the upper-layer Sst supertypes we find
compositionally depleted in schizophrenia are the same ones that decline as
Alzheimer's disease advances. If so, they are plausibly an intrinsically
vulnerable population rather than a SCZ-specific casualty.

**The approach.** SEA-AD assigns each of 80 aged DLPFC donors a continuous
pseudo-progression score (CPS) indexing AD neuropathological severity, and uses
the same supertype taxonomy as the rest of this paper. So each supertype's
abundance can be regressed on CPS with the *same* crumblr model used for the SCZ
case/control comparison, and the two slopes compared directly.

## Result

Across the 16 Sst supertypes, SCZ and AD depletion are strongly concordant:

| | |
|---|---|
| Spearman ρ | **0.868** (P = 1.33 × 10⁻⁵) |
| Pearson r | 0.811 (P = 1.39 × 10⁻⁴) |
| Direction agreement | **93.75%** (15/16) |

Reported in the paper as "ρ = 0.87, n = 16" and "94% of supertypes".

## Layout

```
code/
  build_seaad_cps_input.py   per-nucleus metadata + donor CPS -> per-donor supertype counts
  run_seaad_cps_crumblr.R    crumblr ~ CPS (+ covariates) -> per-supertype slope, SE, p, FDR
data/
  seaad_donor_cps.csv        80 donors x CPS. Committed (force-added past the data/ ignore) —
                             small, and not re-derivable from anything else in this repo.
  SEAAD_INPUTS.md            where the large per-nucleus metadata comes from
results/
  crumblr_input_supertype_neurons.csv     per-donor supertype counts
  crumblr_results_supertype_neurons.csv   the AD axis of Figure 4j
```

## Running it

```bash
cd crossdisorder
python3 code/build_seaad_cps_input.py     # needs the metadata below
Rscript  code/run_seaad_cps_crumblr.R
```

The per-nucleus metadata (~1.3 GB) is **not** committed. It is expected at
`shared/seaad/SEAAD_A9_RNAseq_final-nuclei_metadata.2024-02-13.csv`, overridable
with `--metadata` or the `SEAAD_METADATA` environment variable. See
`data/SEAAD_INPUTS.md` for the download URL. Everything else runs from committed
files.

## Where the result goes

`results/crumblr_results_supertype_neurons.csv` is read by
`genetics/scripts/figures/export_fig4_new_panels.py`, which joins it against the
7-cohort SCZ crumblr betas to produce
`genetics/results/figures/r_panels/panel_ad_concordance_{sst,stats}.csv` — the
Figure 4j panel data. The panel CSVs stay under `genetics/` because the Figure 4
renderer lives there; only the analysis lives here.

That join is recorded in `genetics/results/figures/r_panels/MANIFEST.tsv`, so
`shared/verify_provenance.py` will flag it if this component's output changes
without the panel being refreshed.

## History

This analysis previously lived in `genetics/` as `scripts/14_seaad_dlpfc_crumblr*`
with its data under `genetics/data/seaad_dlpfc/`. It was moved out on 2026-08-06:
it is a compositional analysis on snRNAseq data, with no genetics content — it
consumes no GWAS, MAGMA, or fine-mapping input — and it owns its own external
dataset, which is what makes it a component rather than a script in someone
else's. The move preserved the numbers exactly (re-running from the new location
reproduces the committed panel with max |Δ slope| = 0 across all 16 supertypes).
