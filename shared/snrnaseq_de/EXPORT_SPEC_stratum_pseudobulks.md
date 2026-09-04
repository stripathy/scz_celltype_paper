# Cluster export spec: supertype pseudobulks for the Sst-strata analysis

Paste everything below this line into the Claude session on the cluster.

---

I need you to export per-cohort pseudobulk count files into a standardized format
for transfer to another machine, where they feed a stratified DE/GSEA analysis of
Sst interneuron supertypes (schizophrenia cell-types paper, Figure 5). Work
iteratively: inspect the inputs first, adapt to their actual format, and verify
every claim programmatically — never assume. Do not modify any source file.

## Input situation

Somewhere on this cluster there are pseudobulk files for seven snRNA-seq cohorts:
one file (or set of files) per cohort, containing per-donor × per-supertype
pseudobulk expression for ALL supertypes mixed together in each file. I will tell
you the paths, or you may need to locate them — ask me if the location is unclear.
Before doing anything else, open one file and report its structure: format, axes,
dimensions, how donor and supertype are encoded (column names? separate metadata?),
gene identifier type, and whether values are raw integer counts or normalized.

## Required outputs

Create a new directory `stratum_pseudobulks_export/` (do not write into the source
directories). For EACH cohort produce exactly two files:

### 1. `<cohort>_pseudobulk_counts.parquet`

- Rows = genes. First column named `gene` holding the gene identifier. If the
  source uses gene symbols, use them. If it uses Ensembl IDs, keep them in `gene`
  and ADD a `gene_symbol` column **only if** a symbol mapping already exists inside
  the source data — never invent or download a mapping; if no mapping exists, say
  so and export Ensembl only.
- Remaining columns = one per donor-supertype group, named exactly
  `<donor>|<supertype>` (pipe separator).
- Values = RAW UMI count sums (integers). If the source stores normalized values
  (detect: non-integers), do not silently export them — report it, look for a raw
  layer/version, and only if no raw version exists export the normalized values
  and flag this prominently in the manifest (the downstream model is limma-voom on
  counts, so this matters).
- Include ALL supertypes present, not just Sst (the files are small; downstream
  filters what it needs).
- Supertype naming: normalize to the SEA-AD convention with an underscore before
  the number, e.g. `Sst_2`, `L2/3 IT_5`, `Sst Chodl_1`. If the source uses variants
  (`SST_2`, `Sst 2`, etc.), normalize and record the mapping in the manifest.

### 2. `<cohort>_groups.csv`

One row per column of the parquet (excluding `gene`), with columns:
- `key` — matching the parquet column name `<donor>|<supertype>`
- `donor`, `supertype`
- `n_cells` — number of nuclei aggregated into that group. If the pseudobulk files
  do not store this, recover it from wherever the pipeline recorded it; if it truly
  cannot be recovered, fill NA and flag in the manifest.
- `diagnosis` — standardized to exactly `Control` or `SCZ`. Also keep
  `diagnosis_original` with the source coding. If you find values that do not map
  cleanly to those two (anything other than obvious control/schizophrenia
  labels), STOP and ask me rather than guessing.
- `sex`, `age`, `pmi` — from the donor metadata. Keep source coding; NA if absent
  (flag which cohorts lack which covariates).
- Any additional covariates that the cohort's differential-expression models used
  (batch, brain bank, ancestry PCs, etc.), if they are stored with the pseudobulk
  pipeline — include them as extra columns and list them in the manifest.

### 3. `MANIFEST.md` (one file for the whole export)

For each cohort: the absolute source path(s) and their md5 checksums; the number of
donors by diagnosis; the number of supertypes and genes; whether counts were raw
integers; the supertype-name normalization mapping if any; which covariates are
present; anything flagged. Also record, if discoverable from the pipeline: which
taxonomy/label version produced the supertype calls, and the QC gates applied to
nuclei before pseudobulking. Finish with md5 checksums of every output file.

## Acceptance checks (run and print for every cohort — the export is not done
## until all of these pass or are explicitly flagged)

1. Parquet column names all parse as `<donor>|<supertype>` and match the groups
   CSV `key` column exactly, same order not required but same set required.
2. These 16 Sst supertypes are expected in every cohort:
   Sst_1, Sst_2, Sst_3, Sst_4, Sst_5, Sst_7, Sst_9, Sst_10, Sst_11, Sst_12,
   Sst_13, Sst_19, Sst_20, Sst_22, Sst_23, Sst_25.
   Report which are present/missing per cohort (missing ones are possible in small
   cohorts — flag, don't fail). `Sst Chodl_*` types should be present under their
   own names, not merged into these.
3. Marker sanity check (this catches column/label misalignment): compute mean
   counts-per-10k of the SST gene in Sst_* groups vs in clearly non-Sst groups
   (e.g. Astro_*, L2/3 IT_*). SST must be dramatically higher in the Sst groups.
   Do the same with GFAP or AQP4 for Astro groups if available. Print the numbers.
4. Values are non-negative; report the fraction of non-integer values (should be 0
   for raw counts).
5. Donor counts by diagnosis look plausible for seven SCZ case-control cohorts
   (tens of donors each, roughly balanced); print a per-cohort table.
6. Total export size (should be well under a few GB).

## Packaging

When all cohorts pass, create `stratum_pseudobulks_export_<YYYY-MM-DD>.tar.gz`
containing the parquets, CSVs and MANIFEST.md, print its md5 and size, and tell me
the path. I will transfer it myself; its destination on the receiving machine is
`scz_celltype_paper/shared/snrnaseq_de/stratum_pseudobulks/`.

## Ground rules

- Inspect before acting; print structures and head() of everything you rely on.
- Every number you state must come from code you just ran, not memory.
- If cohorts differ in format, handle each explicitly and note it.
- Do not modify, move, or re-derive the source pseudobulks; read-only.
- If anything is ambiguous (which files are canonical, duplicate versions,
  unclear diagnosis coding), stop and ask rather than choosing silently.
