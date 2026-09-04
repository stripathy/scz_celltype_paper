# Sst interneuron cell-level export

Generated **2026-09-01 14:17:19** on host `tri0354` (SLURM job 2237143); manifest written on `tri-login04`.

Cell-level raw counts for Sst interneuron nuclei, one h5ad per cohort, for the
subsampling sensitivity analysis. 7 cohorts, **101,566 nuclei**, **1,124,476,634 counts**, 0.86 GB total.

## READ FIRST — one expectation did not match, for a definitional reason

Every count-level check reconciles **exactly**. But the donor counts for HBCC
and MSSM2 are lower than the spec's table:

| cohort | spec donors | exported donors | difference |
|---|---|---|---|
| HBCC | 130 | **125** | -5 |
| MSSM2 | 185 | **181** | -4 |

**Cause: those 9 donors have zero Sst nuclei.** They appear in the pseudobulk
grid only because that grid is complete by design — every donor x all 132
supertypes, including all-zero rows (the 'no filtering' instruction for the
pseudobulk stage). A cell-level file cannot contain a donor with no cells.

Verified: the Sst `n_cells` sum for each of these donors is exactly 0 in
`<cohort>_groups.csv`. So the spec's 130 / 185 count donors present in the
*grid*; 125 / 181 is the count of donors with at least one Sst nucleus.
Nuclei and count totals are unaffected and match to the unit.

Donors with zero Sst nuclei:

- **HBCC** (5): `AMPAD_HBCC_0000000004`, `AMPAD_HBCC_0000000056`, `AMPAD_HBCC_0000000084`, `AMPAD_HBCC_0000000152`, `AMPAD_HBCC_0000000398`
- **MSSM2** (4): `AMPAD_MSSM_0000035723`, `AMPAD_MSSM_0000049865`, `AMPAD_MSSM_0000081548`, `AMPAD_MSSM_0000083316`

Nothing has been transferred. Confirm this interpretation before using the export.

## Cell selection

The pseudobulks in `shared/snrnaseq_de/stratum_pseudobulks/` (MANIFEST
2026-08-31) were built by `build_pseudobulk.py`, which sums
`clean/<cohort>.h5ad` over donor x supertype. This export therefore reuses that
exact code path: **the only operation is a row subset of the same
`clean/<cohort>.h5ad` objects** on `obs['supertype']` — the identical column
the pseudobulk export grouped on.

- **No cell QC was applied here.** None. Any QC is inherited from the upstream
  source objects, unchanged since the pseudobulk build.
- **No label transfer was re-run** and no annotation was re-derived.
- **No gene filter** of any kind, including all-zero removal — the full native
  gene space is retained, as the downstream enrichment needs the full ranking.
- Included: the 16 SEA-AD Sst supertypes. `Sst Chodl_1` / `Sst Chodl_2` are
  **excluded** (verified absent from every output).

## Counts slot

For all 7 cohorts the raw integer counts were taken from **`.X`** of the clean
object. This was verified per cohort at export time, not assumed — the script
aborts if `(X[:200].data % 1).sum() != 0`:

| cohort | slot | `X[:200].max()` | `(X[:200].data % 1).sum()` |
|---|---|---|---|
| Batiuk | `.X` | 4,722 | 0.0 |
| HBCC | `.X` | 3,205 | 0.0 |
| McLean | `.X` | 2,405 | 0.0 |
| MSSM1 | `.X` | 1,268 | 0.0 |
| MSSM2 | `.X` | 1,926 | 0.0 |
| Multiome | `.X` | 7,371 | 0.0 |
| OFC | `.X` | 482 | 0.0 |

Note the upstream normalisation this dataset family is prone to *was* present
and was resolved earlier: PsychAD's `.X` is log-normalised, so HBCC and MSSM2
took their counts from its `layers['counts']` when the clean objects were
built. By this stage `.X` is raw everywhere.

## Source objects

| cohort | direct source (subset from) | upstream original | counts slot upstream |
|---|---|---|---|
| Batiuk | `/home/shreejoy/scratch/dataset_compliation/clean/Batiuk.h5ad` | `/scratch/medney/Datasets/Supertypes/Bat_updated.h5ad` | `X` |
| HBCC | `/home/shreejoy/scratch/dataset_compliation/clean/HBCC.h5ad` | `/scratch/medney/Datasets/PsychAD_Data/merged_final_clean.h5ad` | `layers/counts` |
| McLean | `/home/shreejoy/scratch/dataset_compliation/clean/McLean.h5ad` | `/scratch/medney/Datasets/Supertypes/Ruz_updated.h5ad` | `X` |
| MSSM1 | `/home/shreejoy/scratch/dataset_compliation/clean/MSSM1.h5ad` | `/scratch/medney/Datasets/Supertypes/Ruz_updated.h5ad` | `X` |
| MSSM2 | `/home/shreejoy/scratch/dataset_compliation/clean/MSSM2.h5ad` | `/scratch/medney/Datasets/PsychAD_Data/merged_final_clean.h5ad` | `layers/counts` |
| Multiome | `/home/shreejoy/scratch/dataset_compliation/clean/Multiome.h5ad` | `/scratch/medney/Datasets/Brain_scope/Multiome.h5ad` | `X` |
| OFC | `/home/shreejoy/scratch/dataset_compliation/clean/OFC.h5ad` | `/scratch/medney/Datasets/Supertypes/OFC_updated.h5ad` | `X` |

## Schema

`.X` — raw integer counts, cells x genes, CSR, gzip-compressed. Not normalised,
not log-transformed, not scaled.

`.obs` — exactly five columns: `cell_id` (also the index), `donor` (Python
`str` in all cohorts, including numeric-looking Multiome IDs), `supertype`,
`diagnosis` (`Control` / `SCZ` only), `cohort`.

`.var` — the cohort's native gene identifier as index, plus `gene_symbol`
carried through from the existing GENCODE v44 mapping. Nothing was remapped.

| cohort | native ID space | genes | `gene_symbol` non-null |
|---|---|---|---|
| Batiuk | symbol | 60,617 | 43,192 |
| HBCC | Ensembl | 34,890 | 34,863 |
| McLean | symbol | 17,658 | 17,649 |
| MSSM1 | symbol | 17,658 | 17,649 |
| MSSM2 | Ensembl | 34,890 | 34,863 |
| Multiome | symbol | 33,822 | 27,028 |
| OFC | mixed (6,060 Ensembl / 20,131 symbol) | 26,191 | 25,981 |

## Reconciliation against the 2026-08-31 pseudobulk export

For each cohort the exported cells were grouped by donor x supertype, summed
per gene, and compared **gene-by-gene** against the corresponding columns of
`<cohort>_pseudobulk_counts.parquet`, plus per-group `n_cells` and `n_counts`
against `<cohort>_groups.csv`.

| cohort | nuclei | target | counts | target | gene-matrix mismatches | per-group n_cells | per-group n_counts | donors |
|---|---|---|---|---|---|---|---|---|
| Batiuk | 11,576 | 11,576 ✓ | 98,377,040 | 98,377,040 ✓ | **0** | ✓ | ✓ | 15 / 15 |
| HBCC | 26,194 | 26,194 ✓ | 405,055,712 | 405,055,712 ✓ | **0** | ✓ | ✓ | 125 / 130 ⚠ |
| McLean | 11,199 | 11,199 ✓ | 100,614,726 | 100,614,726 ✓ | **0** | ✓ | ✓ | 32 / 32 |
| MSSM1 | 1,573 | 1,573 ✓ | 12,835,163 | 12,835,163 ✓ | **0** | ✓ | ✓ | 35 / 35 |
| MSSM2 | 27,346 | 27,346 ✓ | 396,029,499 | 396,029,499 ✓ | **0** | ✓ | ✓ | 181 / 185 ⚠ |
| Multiome | 2,251 | 2,251 ✓ | 17,519,983 | 17,519,983 ✓ | **0** | ✓ | ✓ | 11 / 11 |
| OFC | 21,427 | 21,427 ✓ | 94,044,511 | 94,044,511 ✓ | **0** | ✓ | ✓ | 61 / 61 |
| **total** | **101,566** | **101,566 ✓** | **1,124,476,634** | **1,124,476,634 ✓** | **0** | | | |

- Gene-by-gene matrix agreement is **exact** in all 7 cohorts: 0 differing
  entries, max absolute difference 0.0. Gene order also matches the parquet.
- Nuclei and total counts match the spec's targets to the unit in all 7.
- Per-group `n_cells` and `n_counts` match `<cohort>_groups.csv` exactly in all 7.
- The only deviation is the HBCC / MSSM2 donor count, explained above.

## Files

| file | nuclei | genes | size | md5 |
|---|---|---|---|---|
| `Batiuk_sst_cells.h5ad` | 11,576 | 60,617 | 89.9 MB | `330e2f69c0c04068a50a7a67047959d4` |
| `HBCC_sst_cells.h5ad` | 26,194 | 34,890 | 275.2 MB | `7fa626e79e37f0a1bb1738f658fec4e3` |
| `McLean_sst_cells.h5ad` | 11,199 | 17,658 | 78.0 MB | `a577c9d3d2adf0221ce1befcdde0cb3e` |
| `MSSM1_sst_cells.h5ad` | 1,573 | 17,658 | 11.6 MB | `d114e7fd5a50a84a20d69af7655c2ff6` |
| `MSSM2_sst_cells.h5ad` | 27,346 | 34,890 | 279.3 MB | `d0ac09dca37b89813b318cc8c4750f43` |
| `Multiome_sst_cells.h5ad` | 2,251 | 33,822 | 14.7 MB | `cc967caebb77f16bb118c9c81d915ce7` |
| `OFC_sst_cells.h5ad` | 21,427 | 26,191 | 113.9 MB | `30f7e5378d236ea056a303047833f192` |
| **total** | **101,566** | | **0.86 GB** | |

Total is 0.86 GB — well under the 5 GB threshold, so no all-zero gene removal
was applied and no expression-based filter was used.

