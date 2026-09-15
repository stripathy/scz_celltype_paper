# Stratum pseudobulk export — 7-cohort SCZ meta-analysis

Generated **2026-08-31 13:45:04** on `tri-login02`.
Destination: `scz_celltype_paper/shared/snrnaseq_de/stratum_pseudobulks/`
Purpose: stratified DE/GSEA of Sst interneuron supertypes (Figure S8).

Per cohort: `<cohort>_pseudobulk_counts.parquet` (genes x donor-supertype groups)
and `<cohort>_groups.csv` (one row per group column).

## Read this first

**1. `gene_symbol` provenance — read before using it.** The symbol column was
NOT present in the original upstream cohort objects. It was added earlier in
this pipeline by mapping each cohort's native gene IDs against **GENCODE v44**
(`gencode_v44_gene_map.csv`), and is stored in `var` of the pseudobulk h5ads
this export reads. It is therefore 'already inside the source data' only in
that narrow sense. The `gene` column is always the cohort's untouched native
identifier — use that if you want to avoid the added mapping entirely.

**2. Complete unfiltered grid.** Every donor has a column for all 132
supertypes, including combinations with zero nuclei, which appear as all-zero
columns. This matches the explicit 'no filtering' instruction for the
pseudobulk stage. **Filter on `n_cells` in the groups CSV before modelling** —
17,437 of 61,908 columns have `n_cells == 0`. Including them costs ~2% of
file size because zstd collapses them.

**3. Gene spaces differ between cohorts and are not intersected.** See the
per-cohort table. Cohorts are not directly concatenatable without first
mapping to a common ID space.

## Provenance chain

Each cohort was derived: original object -> donor-subset clean h5ad ->
donor x supertype pseudobulk -> this export. Only the pseudobulk h5ad is
checksummed here; the upstream originals are 3.6-106 GB and live in another
user's scratch (`/scratch/medney/Datasets`), read-only.

- Donor list (defines the analysis set): `/home/shreejoy/scratch/dataset_compliation/7_cohorts_metadata_names.csv`
- Donor metadata (age/sex/diagnosis/PMI) comes from that CSV, not from the
  per-cohort objects, which disagree on vocabulary.

### Taxonomy / label version

SEA-AD cortical taxonomy: **132 supertypes -> 24 subclasses -> 3 classes**.
Supertype calls per cohort:

- **Batiuk** — obs['predicted.id'] in Bat_updated.h5ad (Seurat TransferData)
- **HBCC** — argmax over prediction.score.* in cell_type_bias/data/HBCC_{1,2}.rds, joined on barcodekey
- **MSSM1** — obs['predicted.id'] in Ruz_updated.h5ad (Seurat TransferData)
- **MSSM2** — argmax over prediction.score.* in cell_type_bias/data/MSSM_{1..4}.rds, joined on barcodekey
- **McLean** — obs['predicted.id'] in Ruz_updated.h5ad (Seurat TransferData)
- **Multiome** — obs['predicted.id'] in Multiome.h5ad (Seurat TransferData)
- **OFC** — obs['predicted.id'] in OFC_updated.h5ad (Seurat TransferData)

The SEA-AD *reference build* used for the label transfer is **not recorded**
in any of the source objects; it could not be determined and is not asserted here.

### QC gates applied to nuclei

This pipeline applied **no nucleus-level QC**: cells were subset by donor only,
then summed. Any QC is inherited from the upstream objects.

Those objects carry `uns['QCed'] = False` and `uns['normalized'] = False`, but
the observed distributions contradict the QCed flag — it appears to describe the
RDS->h5ad converter, not the data. Observed evidence of upstream gating:

| source | evidence |
|---|---|
| Ruz (MSSM1, McLean) | `mito.perc` max 10.45%, `nFeature_RNA` min 910 |
| OFC | `pct_counts_mito` max 14.99%, `nFeature_RNA` min 388, `dd_doublet` present |
| Batiuk | `nFeature_RNA` min 364; donor-level `Remaining_nuclei_after_QC` present |
| Multiome | `nFeature_RNA` min 180 |
| PsychAD (HBCC, MSSM2) | no QC columns beyond `n_counts` / `n_genes` |

**Exact thresholds are not recorded anywhere in the objects** — the values above
are observed maxima/minima, not declared parameters.

`uns['normalized'] = False` is consistent with the verified finding that all
values are raw integer counts.

## Per-cohort detail

### Batiuk

- CSV cohort label: `Batiuk`
- Source (checksummed): `/home/shreejoy/scratch/dataset_compliation/pseudobulk/Batiuk_pseudobulk.h5ad`
  - md5 `4220d7ac03fe98dd856dc565dc9c5472`
- Upstream chain: `/scratch/medney/Datasets/Supertypes/Bat_updated.h5ad` (counts from `X`, donor column `Donor`) -> `/home/shreejoy/scratch/dataset_compliation/clean/Batiuk.h5ad`
- Donors: **15** (10 Control, 5 SCZ)
- Supertypes: **132**  |  Genes: **60,617** (symbol)
- Group columns: **1,980** (340 with `n_cells == 0`)
- Raw integer counts: **yes** (non-integer fraction 0.0, negatives: False)
- Supertype normalization: **none required** — all labels already match the
  SEA-AD `<name>_<N>[-SEAAD]` convention (identity mapping)
- Covariates present: age, sex, pmi
- Extra columns carried: `cohort`, `subclass`, `class`, `n_counts`, `n_genes_detected`

### HBCC

- CSV cohort label: `HBCC`
- Source (checksummed): `/home/shreejoy/scratch/dataset_compliation/pseudobulk/HBCC_pseudobulk.h5ad`
  - md5 `4eda7be002e5d09741d83d7dfe4f9714`
- Upstream chain: `/scratch/medney/Datasets/PsychAD_Data/merged_final_clean.h5ad` (counts from `layers/counts`, donor column `SubID_export_synapse`) -> `/home/shreejoy/scratch/dataset_compliation/clean/HBCC.h5ad`
- Donors: **130** (79 Control, 51 SCZ)
- Supertypes: **132**  |  Genes: **34,890** (Ensembl)
- Group columns: **17,160** (4,929 with `n_cells == 0`)
- Raw integer counts: **yes** (non-integer fraction 0.0, negatives: False)
- Supertype normalization: **none required** — all labels already match the
  SEA-AD `<name>_<N>[-SEAAD]` convention (identity mapping)
- Covariates present: age, sex, pmi
- Extra columns carried: `cohort`, `subclass`, `class`, `n_counts`, `n_genes_detected`

### MSSM1

- CSV cohort label: `MtSinai`
- Source (checksummed): `/home/shreejoy/scratch/dataset_compliation/pseudobulk/MSSM1_pseudobulk.h5ad`
  - md5 `71833ece7d4ebdc780ccf90926b2e6e5`
- Upstream chain: `/scratch/medney/Datasets/Supertypes/Ruz_updated.h5ad` (counts from `X`, donor column `Donor`) -> `/home/shreejoy/scratch/dataset_compliation/clean/MSSM1.h5ad`
- Donors: **35** (17 Control, 18 SCZ)
- Supertypes: **132**  |  Genes: **17,658** (symbol)
- Group columns: **4,620** (1,908 with `n_cells == 0`)
- Raw integer counts: **yes** (non-integer fraction 0.0, negatives: False)
- Supertype normalization: **none required** — all labels already match the
  SEA-AD `<name>_<N>[-SEAAD]` convention (identity mapping)
- Covariates present: age, sex, pmi
- Extra columns carried: `cohort`, `subclass`, `class`, `n_counts`, `n_genes_detected`

### MSSM2

- CSV cohort label: `MSSM`
- Source (checksummed): `/home/shreejoy/scratch/dataset_compliation/pseudobulk/MSSM2_pseudobulk.h5ad`
  - md5 `e99b1abe3fa4581a4bf15b2f2dfcde47`
- Upstream chain: `/scratch/medney/Datasets/PsychAD_Data/merged_final_clean.h5ad` (counts from `layers/counts`, donor column `SubID_export_synapse`) -> `/home/shreejoy/scratch/dataset_compliation/clean/MSSM2.h5ad`
- Donors: **185** (140 Control, 45 SCZ)
- Supertypes: **132**  |  Genes: **34,890** (Ensembl)
- Group columns: **24,420** (7,839 with `n_cells == 0`)
- Raw integer counts: **yes** (non-integer fraction 0.0, negatives: False)
- Supertype normalization: **none required** — all labels already match the
  SEA-AD `<name>_<N>[-SEAAD]` convention (identity mapping)
- Covariates present: age, sex, pmi
- Extra columns carried: `cohort`, `subclass`, `class`, `n_counts`, `n_genes_detected`

### McLean

- CSV cohort label: `McLean`
- Source (checksummed): `/home/shreejoy/scratch/dataset_compliation/pseudobulk/McLean_pseudobulk.h5ad`
  - md5 `a35715dfaa0f7a5a6b42d635891d2f33`
- Upstream chain: `/scratch/medney/Datasets/Supertypes/Ruz_updated.h5ad` (counts from `X`, donor column `Donor`) -> `/home/shreejoy/scratch/dataset_compliation/clean/McLean.h5ad`
- Donors: **32** (18 Control, 14 SCZ)
- Supertypes: **132**  |  Genes: **17,658** (symbol)
- Group columns: **4,224** (858 with `n_cells == 0`)
- Raw integer counts: **yes** (non-integer fraction 0.0, negatives: False)
- Supertype normalization: **none required** — all labels already match the
  SEA-AD `<name>_<N>[-SEAAD]` convention (identity mapping)
- Covariates present: age, sex, pmi
- Extra columns carried: `cohort`, `subclass`, `class`, `n_counts`, `n_genes_detected`

### Multiome

- CSV cohort label: `Multiome`
- Source (checksummed): `/home/shreejoy/scratch/dataset_compliation/pseudobulk/Multiome_pseudobulk.h5ad`
  - md5 `36b1d1039862fbe92874f8900c126fbf`
- Upstream chain: `/scratch/medney/Datasets/Brain_scope/Multiome.h5ad` (counts from `X`, donor column `Donor`) -> `/home/shreejoy/scratch/dataset_compliation/clean/Multiome.h5ad`
- Donors: **11** (5 Control, 6 SCZ)
- Supertypes: **132**  |  Genes: **33,822** (symbol)
- Group columns: **1,452** (300 with `n_cells == 0`)
- Raw integer counts: **yes** (non-integer fraction 0.0, negatives: False)
- Supertype normalization: **none required** — all labels already match the
  SEA-AD `<name>_<N>[-SEAAD]` convention (identity mapping)
- Covariates present: age, sex  |  **absent: pmi**
- Extra columns carried: `cohort`, `subclass`, `class`, `n_counts`, `n_genes_detected`

### OFC

- CSV cohort label: `OFC`
- Source (checksummed): `/home/shreejoy/scratch/dataset_compliation/pseudobulk/OFC_pseudobulk.h5ad`
  - md5 `a930687fde86f5183af3403d14c77701`
- Upstream chain: `/scratch/medney/Datasets/Supertypes/OFC_updated.h5ad` (counts from `X`, donor column `Donor`) -> `/home/shreejoy/scratch/dataset_compliation/clean/OFC.h5ad`
- Donors: **61** (29 Control, 32 SCZ)
- Supertypes: **132**  |  Genes: **26,191** (mixed (6,060 Ensembl / 20,131 symbol))
- Group columns: **8,052** (1,263 with `n_cells == 0`)
- Raw integer counts: **yes** (non-integer fraction 0.0, negatives: False)
- Supertype normalization: **none required** — all labels already match the
  SEA-AD `<name>_<N>[-SEAAD]` convention (identity mapping)
- Covariates present: age, sex, pmi
- Extra columns carried: `cohort`, `subclass`, `class`, `n_counts`, `n_genes_detected`

### Covariates NOT available in any cohort

The pseudobulk pipeline carries only donor age / sex / PMI. **No batch, brain
bank, ancestry PCs, PRS or RIN are present in these exports.** Some do exist
in the upstream cell-level objects (e.g. `Ruz_updated.h5ad` has `Batch`,
`EUR/AFR/EAS/AMR/SAS_Ancestry`, `PRS`; `Bat_updated.h5ad` has extensive
clinical and library-prep fields). They were not propagated. Say the word and
they can be added.

PMI is entirely absent for **Multiome** (all 11 donors NA).

## Acceptance checks

All run programmatically at export time; full output in the build log.

### Sst supertype coverage

All 16 expected Sst supertypes (`Sst_1,2,3,4,5,7,9,10,11,12,13,19,20,22,23,25`)
are present **with non-zero cells in all 7 cohorts** — none missing anywhere.
`Sst Chodl_1` and `Sst Chodl_2` are present under their own names in all 7
cohorts and are **not** merged into `Sst_*`.

### Marker sanity (mean counts-per-10k)

| cohort | SST in Sst_* | SST in L2/3 IT_* | ratio | GFAP in Astro_* | AQP4 in Astro_* | GFAP in L2/3 IT_* |
|---|---|---|---|---|---|---|
| Batiuk | 4.34 | 0.14 | **30.9x** | 6.46 | 7.03 | 0.0 |
| HBCC | 2.44 | 0.02 | **134.2x** | 8.59 | 5.68 | 0.05 |
| MSSM1 | 4.58 | 0.01 | **323.6x** | 10.96 | 6.18 | 0.02 |
| MSSM2 | 1.82 | 0.02 | **98.8x** | 9.15 | 6.19 | 0.06 |
| McLean | 2.55 | 0.03 | **74.3x** | 9.03 | 5.69 | 0.02 |
| Multiome | 1.9 | 0.0 | **497.3x** | 5.24 | 7.2 | 0.0 |
| OFC | 6.86 | 0.03 | **271.5x** | 4.99 | 7.11 | 0.02 |

SST is 31-497x enriched in `Sst_*` over excitatory groups in every cohort, and
GFAP/AQP4 are strongly astrocyte-restricted. No column/label misalignment.

### Other checks

- Parquet column names all parse as `<donor>|<supertype>` and the column set
  equals the `key` set in the groups CSV, for all 7 cohorts.
- All values non-negative; non-integer fraction **0.0** in all 7 cohorts.
- Diagnosis: source coding is `Control` / `Schizophrenia` only; mapped to
  `Control` / `SCZ`. No unmappable values in any cohort.

### Donors by diagnosis

| cohort | Control | SCZ | total |
|---|---|---|---|
| Batiuk | 10 | 5 | 15 |
| HBCC | 79 | 51 | 130 |
| MSSM1 | 17 | 18 | 35 |
| MSSM2 | 140 | 45 | 185 |
| McLean | 18 | 14 | 32 |
| Multiome | 5 | 6 | 11 |
| OFC | 29 | 32 | 61 |
| **total** | **298** | **171** | **469** |

MSSM2 is notably case-imbalanced (140 Control / 45 SCZ) because PsychAD is a
general brain bank rather than a matched case-control collection. MSSM1,
McLean, OFC and Multiome are close to balanced; Batiuk is 2:1 control-heavy.

## Flags

1. `gene_symbol` is a GENCODE v44 mapping added by this pipeline, not upstream data (see top).
2. Zero-cell group columns are included by design — filter on `n_cells`.
3. Multiome has no PMI for any donor.
4. No batch / ancestry / RIN covariates in any cohort.
5. SEA-AD reference build for the label transfer is not recorded upstream.
6. Upstream QC thresholds are not recorded; only observed ranges are reported.
7. **Multiome donor IDs are numeric strings** (`1041`, `1320`, ...). `pd.read_csv` infers them as int64, which then fails to string-match the parquet column names. Read with `dtype={'donor': str}`, or just join on the `key` column, which is always a string. Affects Multiome only.

## Output checksums

MANIFEST.md itself is not listed (it cannot contain its own md5); its
checksum is reported with the tarball.

| file | size | md5 |
|---|---|---|
| `Batiuk_groups.csv` | 0.2 MB | `ed2ffdde7a219fce9471763c59630e49` |
| `Batiuk_pseudobulk_counts.parquet` | 45.6 MB | `ab8f7afd2e17c0c0acabcfcfc29b5ba9` |
| `HBCC_groups.csv` | 2.2 MB | `35429ee22f4dc45fff61c5bb93d63ee4` |
| `HBCC_pseudobulk_counts.parquet` | 248.6 MB | `3a35da29b0d40f6df2778ee44fbea3ab` |
| `MSSM1_groups.csv` | 0.4 MB | `e7479de1c281f5db4a65718d4aaebc3a` |
| `MSSM1_pseudobulk_counts.parquet` | 25.5 MB | `de1c2965729cd71815e04493b328ff3b` |
| `MSSM2_groups.csv` | 3.3 MB | `04873213c227675187d8ab1b1ba8fe1e` |
| `MSSM2_pseudobulk_counts.parquet` | 306.3 MB | `e3b147975c853c9756353aa43dd290c9` |
| `McLean_groups.csv` | 0.4 MB | `c90131c5bb771c315427bf30872cff51` |
| `McLean_pseudobulk_counts.parquet` | 47.3 MB | `e887b2ef983fbab6c5f9045469a8f52f` |
| `Multiome_groups.csv` | 0.1 MB | `a4414aa7d9cc4b8e8baa73b66fa64423` |
| `Multiome_pseudobulk_counts.parquet` | 19.0 MB | `c6de6bf8b3c1c744224224a44c3c6bb2` |
| `OFC_groups.csv` | 0.7 MB | `eb0281126f28e29050cc238b64aff4c3` |
| `OFC_pseudobulk_counts.parquet` | 103.3 MB | `591697a75afa4195bb36617060f69d48` |

