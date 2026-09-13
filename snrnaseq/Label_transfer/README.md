# Label_transfer/ — SEA-AD supertype labels on all seven datasets

Stage 1 of the pipeline. Each dataset is mapped onto the **SEA-AD neurotypical
supertype taxonomy** by Seurat anchor transfer, and written back out as a Seurat
object carrying a `predicted.id` column plus standardised donor metadata. Every
later stage groups cells by `predicted.id`.

One directory per source dataset. They are independent — run in any order.

## The reference, and the method

All five chains load the same two files:

```r
Counts_ref <- readRDS(".../raw_counts_ref.rds")            # cells x genes
meta_hodge <- readRDS(".../Neurotypical_ref_metadata.rds") # carries $Supertype
```

and then do the same thing:

1. keep genes common to reference and query,
2. `NormalizeData(scale.factor = 1e6)` + `FindVariableFeatures(nfeatures = 3000)` on both,
3. PCA on the reference,
4. `FindTransferAnchors(reference.reduction = "pca", dims = 1:30)`,
5. `TransferData(refdata = ref$Supertype)`,
6. re-attach the transferred metadata to the **unfiltered** counts and save.

Step 6 matters: several chains subset or split the query for the transfer itself
and then rebuild the object on the original count matrix, so no cells are lost
to the labelling step.

Despite the `counts_hodge` / `meta_hodge` variable names this reference is
**SEA-AD** (Gabitto et al. 2024), not Hodge et al. 2019 — see issue 10 in
[`../../KNOWN_ISSUES.md`](../../KNOWN_ISSUES.md).

## The chains

| Directory | Datasets produced | Scripts | Output object(s) |
|---|---|---|---|
| `Batiuk/` | Batiuk | `Label_transfer.r` | `Bat_updated.rds` |
| `Frohlich/` | Fröhlich (`OFC`) | `1_Label_transfer.r`, `2_Compile.r` | `OFC_updated.rds` |
| `Ruzicka/` | McLean + MSSM 1 (`MtSinai`) | `Label_transfer.r` | `Ruz_updated.rds` |
| `Multiome/` | Multiome | `Label_transfer.r` | `Multiome.rds` |
| `PsychAD/` | HBCC + MSSM 2 | `1_py_to_R.py` → `2_Extract_matrices.r` → `3_Label_transfer.r` → `4_compile_objects.r` | `HBCC_{1,2}.rds`, `MSSM_{1..4}.rds` |

`Ruz_updated.rds` holds both Ruzicka datasets; downstream code splits it on
`Cohort == "McLean"` / `"MtSinai"`.

### Why two of them are multi-step

**Fröhlich** and **PsychAD** are too large to anchor in one pass, so both split
the query and transfer in chunks:

- `Frohlich/1_Label_transfer.r` splits into 5 parts, transfers each, saves
  `query_subset_annotated_{1..5}.rds`; `2_Compile.r` rebinds the five metadata
  tables onto the original object and keeps `Classification ∈ {Schizophrenia,
  Control}`.
- `PsychAD/` starts from one large `.h5ad` (`merged_final_clean.h5ad`).
  `1_py_to_R.py` reads it backed and writes 20 Matrix Market chunks — the script
  notes ~316,000 cells per chunk — and `2_Extract_matrices.r` turns those
  into `dgCMatrix` RDS; `3_Label_transfer.r` maps Ensembl IDs to symbols
  (`org.Hs.eg.db`) and transfers per chunk; `4_compile_objects.r` reassembles.

### `PsychAD/4_compile_objects.r` produces more than this paper uses

It splits the 20 chunks by `Source` and `Diagnosis` and writes several
collections. **Only the SCZ objects feed this paper:**

| Output | Used here? |
|---|---|
| `MSSM_{1..4}.rds`, `HBCC_{1,2}.rds` — `Source ∈ {MSSM, HBCC}`, `Diagnosis ∈ {Schizophrenia, Control}` | **yes** |
| `RADC.rds`, `MSSM_AD_part*.rds` — `Diagnosis ∈ {AD, Control}` | no |
| `HBCC_BPD_{1,2,3}.rds` — `Diagnosis ∈ {Bipolar, Control}` | no |

It is a long file (1,475 lines) because the per-chunk metadata joins are written
out one chunk at a time rather than looped.

## Metadata standardisation

Each chain normalises its source's columns to the names every later stage
assumes — `Donor`, `Diagnosis` (`"Schizophrenia"` / `"Control"`), `Sex`, `Age`,
`PMI`, `Cohort`. The mappings differ by source, e.g. Batiuk's `Scz`/`Ctr` and
`Gender` `F`/`M`, Ruzicka's `Phenotype` `SZ`/`CON`, Multiome's `Disorder` and
`Age_death`, PsychAD's PMI in minutes (divided by 60).

Batiuk additionally drops donor `MB8`.

Age and diagnosis **filtering** is not done here — it is applied at load time by
each consuming stage.
