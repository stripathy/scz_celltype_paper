# Guide to `all_samples_annotated.h5ad`

Structure, columns, and recommended usage of the merged Xenium spatial
transcriptomics object for the Kwon et al. SCZ cohort — the single object every
spatial number in this paper derives from.

For pipeline methodology see [`methods_writeup.md`](methods_writeup.md); for
which object state is canonical and why, see
[`../DATA_FLOW.md`](../DATA_FLOW.md).

---

## File overview

| Property | Value |
|---|---|
| **Path** | `spatial/output/all_samples_annotated.h5ad` — an absolute **symlink** to `~/Github/SCZ_Xenium/output/all_samples_annotated.h5ad` (repo convention: this monorepo tracks code, data is wired in) |
| **Size** | 1.39 GB |
| **State** | the **2026-04-01** pipeline state, reinstated 2026-07-31; md5 `763e6655fc55839f177d57aa98dc5453` |
| **Cells** | 1,339,151 (all 24 samples, no cells removed) |
| **Genes** | 300 (Kwon et al. Xenium panel) |
| **Samples** | 24 (12 Control, 12 SCZ) |
| **`X`** | raw UMI counts, CSR sparse `int32`, 134,482,432 nonzeros |
| **`.obsm['spatial']`** | (1,339,151, 2) — (x, y) in microns |
| **`.layers` / `.uns` / `.varm` / `.obsp`** | empty; `.var` carries gene symbols only |

The 24 per-sample objects (`Br*_annotated.h5ad`) are *not* mirrored into this
repo — they stay in `~/Github/SCZ_Xenium/output/h5ad/`, which is where
`spatial/code/analysis/config.py` (`BASE_DIR`, overridable via `XENIUM_BASE`)
reads them from.

To freeze a physical copy instead of following the upstream file:

```bash
rm spatial/output/all_samples_annotated.h5ad && cp ~/Github/SCZ_Xenium/output/all_samples_annotated.h5ad spatial/output/
```

`spatial/output/` is gitignored, so neither the symlink nor a copy can be
committed by accident.

---

## ⚠️ If you were using `all_samples_annotated_clean.h5ad`

The older collaborator object (`SCZ_Xenium/output/all_samples_annotated_clean.h5ad`,
0.32 GB, built 2026-03-31) is the same cells with columns stripped and renamed.
This object keeps the pipeline's own names, so code written against the clean
object needs this mapping:

| clean object | this object | note |
|---|---|---|
| `qc_pass` | **`corr_qc_pass`** | ⚠️ **`qc_pass` exists in both files and means different things.** See below. |
| `class` | `corr_class` | |
| `subclass` | `corr_subclass` | |
| `supertype` | `corr_supertype` | |
| `subclass_corr` | `corr_subclass_corr` | |
| `subclass_margin` | `corr_subclass_margin` | |
| `supertype_corr` | `corr_supertype_corr` | |
| `diagnosis` = `Control` / `Schizophrenia` | `diagnosis` = `Control` / **`SCZ`** | |
| `rin`, `race` | *absent* | this object's donor metadata is `diagnosis`, `sex`, `age`, `pmi` only, taken from the committed crumblr input (the values the analyses actually used) |
| — | `lieber_*`, `harmony_*`, `*_label`, `nuclear_*`, `banksy_*`, `fail_*`, `hybrid_qc_pass`, `spatial_domain` | present only here |

The clean object also predates the `lieber_*` columns, which is the main reason
to prefer this one.

---

## Loading

```python
import anndata as ad

# backed (lazy) — recommended, the file is 1.39 GB
adata = ad.read_h5ad("spatial/output/all_samples_annotated.h5ad", backed="r")

# obs_names are NOT unique across samples (Xenium cell IDs repeat per sample).
# Use sample_id + obs_name as a compound key, or positional indexing.
adata.obs_names_make_unique()
```

Per-sample analyses should go through the repo loader instead, which applies the
QC gate and renames labels for you:

```python
import sys; sys.path.insert(0, "spatial/code/analysis")
import config
obs = config.load_cells("Br6389", cortical_only=True,
                        extra_obs_columns=["predicted_norm_depth"])
```

---

## QC gates

Three boolean gates are stored; **they are not interchangeable.**

| Column | Cells | % | What it is |
|---|--:|--:|---|
| `qc_pass` | 1,298,687 | 97.0% | step-01 spatial QC only — `~any(fail_*)`: negative-probe / negative-codeword rate, extreme UMI counts, low gene count, unassigned-transcript rate |
| **`corr_qc_pass`** | **1,221,519** | **91.2%** | **the recommended gate.** `qc_pass` ∧ subclass-correlation-margin filter (bottom 5% per sample) ∧ ~`doublet_suspect`. Strictly nested inside `qc_pass`. This is what the clean object calls `qc_pass`. |
| `hybrid_qc_pass` | 1,302,631 | 97.3% | optional variant that re-adjudicates doublets at nuclear resolution and rescues high-UMI cells. Not used for the paper's headline numbers; `config.load_cells(qc_mode="hybrid")` selects it. |

```python
adata_qc = adata[adata.obs["corr_qc_pass"]]      # 1,221,519 cells
```

Three columns are only meaningful on QC-pass cells:

- `corr_class` / `corr_subclass` / `corr_supertype` are `Unassigned` for exactly
  the 40,464 cells that fail `qc_pass`.
- `predicted_norm_depth` is `NaN`, and `banksy_domain` is `''`, for exactly the
  117,632 cells that fail `corr_qc_pass`.
- `layer` and `spatial_domain` are **KNN-imputed** for QC-fail cells (so the
  viewer can render them). Don't treat them as measurements outside
  `corr_qc_pass`.

### The paper's analysis set

```python
obs = adata.obs
gate = obs["qc_pass"] & obs["corr_qc_pass"] & (obs["spatial_domain"] == "Cortical")
neuronal = gate & obs["corr_class"].isin(["GABAergic", "Glutamatergic"])
# 742,103 cortical QC-pass cells; 356,313 of them neuronal  <- the analysis set
```

356,313 neuronal cortical cells is the number the crumblr composition results
were built on. Reproducing it is the provenance check in
[`../DATA_FLOW.md`](../DATA_FLOW.md); `shared/verify_provenance.py` runs it.

---

## `.obs` column reference

### Donor metadata

Sample-level covariates broadcast to every cell, so no separate metadata file is
needed.

| Column | Type | Description |
|---|---|---|
| `sample_id` | category (24) | Brain identifier, e.g. `Br1113` |
| `diagnosis` | category | `Control` or `SCZ` |
| `sex` | category | `F` or `M` |
| `age` | float64 | Age at death (years) |
| `pmi` | float64 | Post-mortem interval (hours) |

### Cell QC

| Column | Type | Description |
|---|---|---|
| `qc_pass`, `corr_qc_pass`, `hybrid_qc_pass` | bool | see [QC gates](#qc-gates) |
| `total_counts`, `n_genes` | int64 | UMIs and genes detected (median 654, IQR 371–1,059; 100 genes — on `corr_qc_pass` cells) |
| `neg_probe_sum`, `neg_codeword_sum`, `unassigned_sum` | int64 | negative-control and unassigned-transcript totals feeding the QC flags |
| `fail_neg_probe`, `fail_neg_codeword`, `fail_total_counts_high`, `fail_total_counts_low`, `fail_n_genes_low`, `fail_unassigned` | bool | the six individual QC failures; `qc_pass` = none of them |
| `corr_subclass_margin` | float32 | gap between the top-2 subclass correlations; the per-sample bottom 5% is what `corr_qc_pass` removes |
| `doublet_suspect` | bool | spatial doublet by marker co-expression — 9,721 cells |
| `doublet_type` | category | `''` / `Glut+GABA` (7,963) / `GABA+GABA` (1,758) |

### Cell-type labels — three parallel sets

**Use the `corr_*` columns.** They are the paper's labels: a two-stage Pearson
correlation classifier (24 subclass centroids, then supertypes within the
assigned subclass) against centroids built from high-confidence HANN exemplars.
The other two sets are kept for comparison only.

| Column | Type | Description |
|---|---|---|
| `corr_class` | category (4) | `Glutamatergic` / `GABAergic` / `Non-neuronal` / `Unassigned` |
| `corr_subclass` | category (25) | 24 SEA-AD subclasses + `Unassigned` |
| `corr_supertype` | category (139) | SEA-AD supertypes, assigned within subclass |
| `corr_subclass_corr`, `corr_supertype_corr` | float32 | best Pearson correlation at each level |
| `class_label`, `subclass_label`, `supertype_label` (+ `*_confidence`) | category / float32 | **MapMyCells HANN** labels from pipeline step 02 (100 bootstrap iterations against the SEA-AD MTG taxonomy). 82.2% identical to `corr_subclass` on QC-pass cells. |
| `harmony_class`, `harmony_subclass`, `harmony_supertype` (+ `*_confidence`) | category / float32 | **legacy** Harmony-integration label transfer (`spatial/code/analysis/validation/02c_run_harmony_transfer.py`). 59.5% identical to `corr_subclass`; it over-assigns Sst and mislabels vascular/glial types against the SEA-AD MERFISH benchmark. Do not use. |

### Nuclear re-segmentation

Populated for all 24 samples from the optional nuclear-resolution step
(removed 2026-09-04; recoverable from tag `pre-prune-2026-09-04`), which re-counts transcripts inside the
nuclear boundary to re-adjudicate whole-cell doublet calls.

| Column | Type | Description |
|---|---|---|
| `nuclear_total_counts`, `nuclear_n_genes` | int64 | counts restricted to the nuclear boundary |
| `nuclear_fraction` | float32 | nuclear / whole-cell counts |
| `nuclear_doublet_status` | category | `clean` 1,326,172 / `resolved` 8,166 / `persistent` 2,490 / `nuclear_only` 2,298 / `insufficient` 25 |
| `nuclear_doublet_suspect`, `nuclear_doublet_type` | bool / category | doublet call recomputed on nuclear counts |

### Spatial position

| Column | Type | Description |
|---|---|---|
| `predicted_norm_depth` | float64 | normalized cortical depth, 0 = pia, 1 = white-matter border (observed range −0.104 to 1.024; `NaN` off `corr_qc_pass`) |
| `layer` | category (7) | final spatially-smoothed layer: `L1` / `L2/3` / `L4` / `L5` / `L6` / `WM` / `Vascular` |
| `spatial_domain` | category (3) | canonical domain: `Cortical` / `WM` / `Vascular` — BANKSY domain, with any cell whose smoothed layer is `WM` forced to `WM`; KNN-imputed for QC-fail cells |
| `banksy_domain` | category | the raw BANKSY call before that reconciliation (`''` off `corr_qc_pass`) |
| `banksy_cluster` | int64 | BANKSY cluster (0–7; −1 = not clustered) |
| `banksy_is_l1` | bool | member of the L1 border cluster — 78,645 cells |

**Which depth model this is:** the 2026-04-01 model (K = 50 spatial neighbours,
trained on all 24 SEA-AD MERFISH donors, unsmoothed). The later K = 100 /
low-CPS / spatially-smoothed model was deployed on 2026-06-11 and then rolled
back with the rest of the June state; it lives in
`SCZ_Xenium/archive/2026-06-11/`. This object therefore has one depth column,
not the `_prev` / `_raw` / smoothed trio the June objects carry.

### Kwon et al. author annotations (`lieber_*`)

The Xenium authors' own per-cell labels, joined onto our cells for validation.
**They are not in the authors' data release** — `LieberInstitute/spatialDLPFC_SCZ_XENIUM`
ships cluster numbers only, its labelled `.RDS` is gitignored, and GEO GSE307404
is raw vendor output. These four columns are the join of their two released
per-cell tables onto our cells **by within-sample raw cell position**:
1,262,965 of 1,339,151 cells matched (**94.31%**); the rest are `unmatched` in
every `lieber_*` column. Built by
[`code/analysis/integrate_lieber_annotations.py`](code/analysis/integrate_lieber_annotations.py);
they drive Supplementary Figure S2a.

| Column | Type | Description |
|---|---|---|
| `lieber_cluster` | category (19) | their BANKSY cluster, `1`–`18`, from `06_cell_type_clustering/banksy_clustering_lambda0.1_res0.7.csv` (λ = 0.1, res = 0.7), or `unmatched` |
| `lieber_celltype` | category (13) | that cluster's coarse cell type (table below), or `unmatched` |
| `lieber_spd` | category (8) | their spatial-domain label `spd01`–`spd07`, from `label_transfer_N24_k50_smoothed_labels.csv`, or `unmatched` |
| `lieber_layer` | category (8) | `lieber_spd` mapped to a layer name, or `unmatched` |

`lieber_cluster` → `lieber_celltype`, with matched cell counts:

| cluster | type | cells | | cluster | type | cells |
|--:|---|--:|---|--:|---|--:|
| 2 | `L2/3 Ex` | 133,500 | | 4 | `Ast` | 106,566 |
| 11 | `L4/5 Ex` | 55,669 | | 16 | `Ast` | 35,898 |
| 18 | `L5 Ex` | 4,825 | | 17 | `Ast` | 22,027 |
| 14 | `L6 Ex` | 45,636 | | 7 | `Mic` | 81,016 |
| **9** | **`MGE`** | **66,513** | | 3 | `Endo` | 109,720 |
| **12** | **`CGE`** | **51,567** | | 5 | `Endo` | 95,784 |
| 1 | `Oligo` | 136,581 | | 13 | `Endo` | 48,430 |
| 6 | `Oligo` | 83,520 | | 10 | `Ambig/Oligo` | 63,996 |
| 8 | `Oligo` | 76,116 | | 15 | `Ambig/In/Endo` | 45,601 |

`lieber_spd` → `lieber_layer`: `spd07` → `L1/M` (162,849), `spd06` → `L2/3`
(294,672), `spd02` → `L3/4` (50,415), `spd05` → `L5` (384,494), `spd03` → `L6`
(117,414), `spd04` → `WM` (222,736), `spd01` → `WMtz` (30,385).

#### ⚠️ The MGE/CGE trap

The authors' cluster→type `case_when` block is copy-pasted into four of their
scripts and **one copy is transposed**:
`06_cell_type_clustering/05_make_updated_cell_type_markers_heatmap_with_annotations.R`
(L47-48) has `9 → CGE`, `12 → MGE`. The final object builder
`07_cell_type_de/04_create_SPE_with_cell_type_info.R` — the one that writes
their labelled object and feeds their published DE — has `9 → MGE`,
`12 → CGE`, as do their other two copies. **Cluster 9 = MGE, cluster 12 = CGE**
is correct; confirmed by S.H. Kwon (email, 2026-07-27) and by markers in the
joined cells:

| our `corr_subclass` | → `MGE` | → `CGE` |
|---|--:|--:|
| Sst | 28,836 | 129 |
| Sst Chodl | 797 | 2 |
| Pvalb | 21,075 | 55 |
| Chandelier | 9,835 | 291 |
| Vip | 150 | 24,504 |
| Lamp5 | 924 | 14,685 |
| Sncg | 202 | 5,017 |
| Pax6 | 182 | 5,635 |

We used their script 05 until 2026-07-27, so **anything derived from the June
objects or from `lieber_celltype` before 2026-07-31 has MGE and CGE swapped** —
which silently inverts the S2a validation claim. The columns in *this* object
are the corrected ones (verified above). The authors' own comment that clusters
10 and 15 are uncertain is justified: by our subclasses cluster 15 is 76% OPC
and cluster 10 is mixed.

---

## Spatial coordinates

```python
xy = adata.obsm["spatial"]          # (n_cells, 2), microns; col 0 = x, col 1 = y
df = adata.obs.copy()
df["x"], df["y"] = xy[:, 0], xy[:, 1]
```

## Layer values

Counts on the `corr_qc_pass` set (these are the numbers in `methods_writeup.md`):

| Layer | Cells | Description |
|---|--:|---|
| `L1` | 72,078 | Lamina 1 (sparse, non-neuronal enriched) |
| `L2/3` | 375,019 | Upper cortical, CUX2+ excitatory |
| `L4` | 77,410 | Granular layer |
| `L5` | 210,910 | Large pyramidal neurons |
| `L6` | 205,809 | Deep excitatory (CT and IT) |
| `WM` | 199,769 | White matter (oligodendrocytes, VLMCs) |
| `Vascular` | 80,524 | Vascular compartment |

Use `predicted_norm_depth` when you need finer resolution than discrete layers.

---

## Cohort

Computed from this object's donor columns:

| | Control | SCZ |
|---|---|---|
| **n samples** | 12 | 12 |
| **Age (mean ± SD)** | 46.6 ± 9.7 | 49.2 ± 8.1 |
| **Sex (F/M)** | 7F / 5M | 6F / 6M |
| **PMI (mean ± SD)** | 27.8 ± 11.2 h | 26.7 ± 9.2 h |

### Br2039 — WM-heavy section

Br2039 (SCZ) has the highest white-matter fraction in the cohort. By the column
the analysis gate actually uses it is **70.1% `spatial_domain == "WM"` against a
cohort median of 19.7%**, the next-highest section being 38.8%; by `layer` the
same section is 41.8% against a median of 14.1%, and by `banksy_domain` 62.1%
against 17.7%. The three disagree because they segment on different signals, so
always say which column a WM fraction came from.

**Br2039 is retained in all analyses**, including the cortical compositional
ones. Restricting to `spatial_domain == "Cortical"` already removes its white
matter cell by cell, so excluding the whole section discards usable cortex; the
decision is recorded in `code/modules/constants.py` (`EXCLUDE_SAMPLES = set()`,
unchanged since the monorepo was scaffolded) and improves concordance with the
snRNA-seq proportions. The paper's Methods state the cohort as all 24 sections.

```python
obs = adata.obs
cortical = obs["corr_qc_pass"] & (obs["spatial_domain"] == "Cortical")
```

> An earlier revision of this guide told readers to drop Br2039 here, which
> contradicted the code and the Methods. Corrected 2026-09-15; see issue 18 in
> [`../KNOWN_ISSUES.md`](../KNOWN_ISSUES.md).

---

## Common analysis recipes

### Subclass proportions (neuronal, within-class)

Neuronal and non-neuronal types are analyzed separately — proportions are
within-class ("Sst as a fraction of all neurons", not of all cells).

```python
obs = adata.obs[adata.obs["corr_qc_pass"]
                & (adata.obs["sample_id"] != "Br2039")
                & (adata.obs["spatial_domain"] == "Cortical")]

for cls in ["Neuronal", "Non-neuronal"]:
    sel = (obs["corr_class"].isin(["GABAergic", "Glutamatergic"])
           if cls == "Neuronal" else obs["corr_class"] == "Non-neuronal")
    n = (obs[sel].groupby(["sample_id", "corr_subclass"], observed=True)
         .size().unstack(fill_value=0))
    props = n.div(n.sum(axis=1), axis=0)
```

### Spatial plot for one sample

```python
import matplotlib.pyplot as plt

mask = (adata.obs["sample_id"] == "Br8667").values
xy = adata.obsm["spatial"][mask]
labels = adata.obs.loc[mask, "corr_subclass"]

fig, ax = plt.subplots(figsize=(10, 10), facecolor="black")
for sc in labels.cat.categories:
    idx = (labels == sc).values
    ax.scatter(xy[idx, 0], xy[idx, 1], s=0.5, label=sc, rasterized=True)
ax.set_aspect("equal"); ax.axis("off")
```

### Depth profile for a cell type

```python
obs = adata.obs[adata.obs["corr_qc_pass"]]
fig, ax = plt.subplots(figsize=(6, 4))
for subclass, color in [("Sst", "steelblue"), ("L2/3 IT", "coral")]:
    d = obs.loc[obs["corr_subclass"] == subclass, "predicted_norm_depth"]
    ax.hist(d, bins=50, alpha=0.6, density=True, label=subclass, color=color)
ax.set_xlabel("Normalized depth (0 = pia, 1 = WM)", fontsize=14)
ax.legend(fontsize=12)
```

### Agreement with the authors' annotations (Supp. Fig. S2a)

```python
obs = adata.obs[adata.obs["lieber_celltype"] != "unmatched"]
ct = pd.crosstab(obs["corr_subclass"], obs["lieber_celltype"])
pct = ct.div(ct.sum(1), axis=0) * 100      # % of each of our subclasses
```

### Expression analysis

`X` is raw UMI counts — normalize before any expression comparison.

```python
import scanpy as sc
a = adata[adata.obs["corr_qc_pass"]].to_memory()
sc.pp.normalize_total(a, target_sum=1e4)
sc.pp.log1p(a)
```

The 300-gene panel has limited within-subclass marker coverage for some
inhibitory types (particularly Sst supertypes). Treat supertype-level
comparisons cautiously — per-supertype classification F1 on this panel, versus
the whole-transcriptome ceiling, is Supplementary Figure S2c/e
(`output/celltyping_supplement/data/resolvability_*.csv`).
