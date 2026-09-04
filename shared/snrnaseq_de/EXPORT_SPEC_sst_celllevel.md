# Export spec — Sst cell-level counts (for the cell-count confound test)

Paste everything below the line into a Claude session on the cluster that has the
full per-dataset h5ad objects.

---

## Task

Export **cell-level raw counts for Sst interneuron nuclei only**, one file per
cohort, for the seven snRNA-seq datasets in our SCZ meta-analysis. This feeds a
sensitivity analysis that requires subsampling individual nuclei, which cannot be
done from the pseudobulks we already exported.

**Cohorts (7):** Batiuk, HBCC, McLean, MSSM1, MSSM2, Multiome, OFC

**Why.** We stratified 16 Sst supertypes into depleted / intermediate /
non-depleted groups by their abundance change in SCZ, and found translation and
oxidative-phosphorylation programs suppressed in a graded way. But because the
groups are *defined* by cases having fewer cells, SCZ donors contribute ~30%
fewer nuclei and ~36% fewer counts to the depleted-group pseudobulk than
controls, and we have shown that nuclei count predicts these module scores even
among control donors alone. Downsampling reads is not enough — we need to
subsample *cells* from control donors to match the case distribution, rebuild
pseudobulks, and rerun. Hence cell-level data.

## What to produce

One file per cohort: `<cohort>_sst_cells.h5ad`, containing

- **`.X` = raw integer counts**, cells × genes. Not normalised, not log
  transformed, not scaled. If the object's `.X` has been normalised, take the
  counts from `.raw.X` or `.layers['counts']` — whichever holds true integers.
  Verify with `adata.X[:200].max()` (should be a large integer, not ~10) and
  `float(adata.X[:200].data % 1 .sum()) == 0`.
- **`.obs`** with exactly these columns, named exactly this way:
  `cell_id`, `donor`, `supertype`, `diagnosis`, `cohort`
  - `donor` as **string** (Multiome donor IDs look numeric — cast to str or they
    break every join downstream)
  - `diagnosis` with values exactly `Control` / `SCZ`
  - `supertype` = the SEA-AD supertype label, e.g. `Sst_25`
- **`.var`** with the cohort's **native gene identifier as the index**, plus the
  `gene_symbol` column if the object already carries the GENCODE v44 mapping we
  added earlier (HBCC and MSSM2 are Ensembl-native and need it). Do not remap
  anything new; just carry through what is there.

## Which cells

Include a nucleus only if its supertype is one of these **16** (note: **Sst
Chodl is excluded**, it is a separate subclass):

```
Sst_1  Sst_2  Sst_3  Sst_4  Sst_5  Sst_7  Sst_9  Sst_10
Sst_11 Sst_12 Sst_13 Sst_19 Sst_20 Sst_22 Sst_23 Sst_25
```

Apply **exactly the same cell QC and the same supertype labels** that were used
to build `<cohort>_pseudobulk_counts.parquet` in the earlier export
(`shared/snrnaseq_de/stratum_pseudobulks/`, MANIFEST dated 2026-08-31). Do not
re-run label transfer, do not re-filter. If that export was produced by a script
still on the cluster, reuse its cell-selection code path.

## Mandatory reconciliation check — do this before transferring

The extraction is correct only if summing the exported cells reproduces the
pseudobulks we already have. For each cohort:

1. Group the exported cells by `donor` × `supertype` and sum counts per gene.
2. Compare against the corresponding column of `<cohort>_pseudobulk_counts.parquet`
   (column key is `"<donor>|<supertype>"`).
3. Also compare per-group cell counts and total counts against
   `<cohort>_groups.csv` (`n_cells`, `n_counts`).

**These per-cohort totals must match exactly:**

| cohort | donors | Sst nuclei | total counts |
|---|---|---|---|
| Batiuk | 15 | 11,576 | 98,377,040 |
| HBCC | 130 | 26,194 | 405,055,712 |
| McLean | 32 | 11,199 | 100,614,726 |
| MSSM1 | 35 | 1,573 | 12,835,163 |
| MSSM2 | 185 | 27,346 | 396,029,499 |
| Multiome | 11 | 2,251 | 17,519,983 |
| OFC | 61 | 21,427 | 94,044,511 |
| **total** | | **101,566** | **1,124,476,634** |

Report the reconciliation result for every cohort. **If any cohort does not
match exactly, stop and report the discrepancy rather than transferring** — a
silent mismatch means the cell selection differs from what produced our current
results, and every downstream comparison would be invalid.

## Output

Write to a single directory alongside the existing export, with:

- `<cohort>_sst_cells.h5ad` × 7
- `MANIFEST_sst_cells.md` recording: generation date, host, the source object
  path for each cohort, which slot the raw counts came from
  (`.X` / `.raw.X` / `.layers['counts']`), the exact cell-QC applied, the
  reconciliation result per cohort, and per-cohort n_cells / n_genes / file size.

**Size.** Roughly 101,566 nuclei totalling ~1.12 billion counts. Expect a few
hundred MB to a couple of GB compressed. Use gzip compression when writing
(`adata.write_h5ad(path, compression="gzip")`). If the total lands above ~5 GB,
say so before transferring and we will drop genes with zero counts across all
Sst nuclei in that cohort — but **do not** apply any expression-based gene filter
beyond all-zero removal without telling us, because the downstream analysis
re-derives its own filter and needs the full ranking for gene-set enrichment.

## If moving the data is impractical

Tell us the total size first. There is a fallback where the subsampling is done
on the cluster and only the resulting pseudobulks are transferred (roughly 1/1000
the size), but it requires implementing the matching logic there, so we would
send a second, more detailed spec rather than have you improvise it.

## Gotchas seen in this dataset family

- `.X` is normalised in several of these objects; the integer counts live
  elsewhere. Check, do not assume.
- Multiome has no PMI and numeric-looking donor IDs.
- HBCC and MSSM2 are Ensembl-native; the other five are symbol-native.
- Supertype column names differ between objects — confirm you are reading the
  same column the pseudobulk export used, not a different annotation round.
- Some objects contain both `Sst` and `Sst Chodl`; only the 16 listed above go in.
