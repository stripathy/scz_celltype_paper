# Data flow & dependency graph

How the analysis components depend on each other. **Every cross-component edge
is realized as a git-ignored symlink** (so the wiring is executable, not just
prose) and listed in the Seams table. The upstream root — the snRNA-seq
pipeline — is not in the repo yet; its outputs are symlinked into
`shared/snrnaseq_de/` from their current scattered locations.

For the figure-to-script map, start from the root [`README.md`](README.md).

## Dependency graph (data flows downward)

```
              snrnaseq/   <- upstream pipeline (Endresz et al.); NOT in repo yet
            (upstream root)
                 |  exports -> shared/snrnaseq_de/
                 |    DE_genes_all_cells_scz.csv         (meta DE)
                 |    meta_results_cohorts_subclass.csv  (per-cohort DE)
                 |    nicole_scz_snrnaseq_betas/         (composition, crumblr)
   +-------------+------------------+-----------------------+
   | DE          | DE + composition | composition           | composition
   v             v                  v                       v
transcriptomic/  spatial/         genetics/            crossdisorder/
   ^   ^         |  Xenium DE      |  MAGMA enrichment       |  AD CPS slopes
   |   +---------+  + crumblr      |  fine-mapping           |
   |                 |             |  patch-seq              |
   |                 +-------------+-------------------------+
   |                               |
   |                               v
   +----------------------> manuscript/figures/
                              main/          Fig 2 (transcriptomic), Fig 4 (genetics)
                              supplementary/ S2, S3 (spatial), S6 (snrnaseq),
                                             S8 (transcriptomic), S9, S10 (genetics)

Figures 1 and 3 are assembled outside this repo; spatial/ supplies the annotated
Xenium objects and crumblr tables they draw on.

reserve/ consumes genetics/ and transcriptomic/ outputs but feeds nothing.

External hub (~/Github/shared_data/, per ~/Github/DATA_LAYOUT.md):
   SEA-AD / MERFISH references ----> genetics/, spatial/  (symlinked, not absorbed)
```

## Per-component inputs → outputs

| component | key inputs (<- from) | produces -> (consumed by) |
|---|---|---|
| `snrnaseq/` *(TBD)* | raw 7-cohort snRNA-seq | DE betas, composition betas -> everyone downstream |
| `genetics/` | composition betas (<- snrnaseq), AD slopes (<- crossdisorder), GWAS set + patch-seq (owns), SEA-AD DLPFC ref (<- shared_data) | Fig 4, S9, S10, T6; the vulnerable-vs-not marker table -> reserve/ |
| `spatial/` | Xenium raw, SEA-AD/MERFISH ref (<- shared_data) | Xenium DE + crumblr composition -> transcriptomic, Figs 1-3; S2, S3; SM1 |
| `transcriptomic/` | DE betas (<- snrnaseq), Xenium DE + depth (<- spatial) | Fig 2, S8 |
| `crossdisorder/` | SEA-AD DLPFC/MTG metadata + CPS (self-contained) | AD crumblr slopes -> genetics (Fig 4i, S10b) |
| `snrnaseq/composition_sensitivity/` | per-donor counts (<- snrnaseq), Xenium crumblr (<- spatial) | S6 |
| `reserve/` | genetics + transcriptomic outputs | nothing in the paper |

## Internal seams (realized symlinks = the edges)

| edge | path in repo | → resolves to (current) |
|---|---|---|
| snrnaseq → transcriptomic | `transcriptomic/data/DE_genes_all_cells_scz.csv` → `shared/snrnaseq_de/` | real file held in `shared/snrnaseq_de/` (was `~/Downloads/`; moved in 2026-07-31) |
| snrnaseq → transcriptomic | `transcriptomic/data/meta_results_cohorts_subclass.csv` → `shared/snrnaseq_de/` | `scz_pathway_enrichment/data/…` |
| snrnaseq → genetics + spatial | `shared/snrnaseq_de/nicole_scz_snrnaseq_betas/` | `SCZ_Xenium/data/nicole_scz_snrnaseq_betas/` |
| ↳ spatial consumes | `spatial/data/nicole_scz_snrnaseq_betas` → `shared/snrnaseq_de/` | (as above) |
| ↳ genetics consumes | `genetics/scripts/figures/build_composition_table.py` reads `shared/snrnaseq_de/…` | (as above) |
| spatial → transcriptomic | `transcriptomic/scripts/09` reads `../spatial/output/de/de_results_subclass.csv` | `SCZ_Xenium/output/de/…` (symlinked) |
| RNAscope data (reserve) | `reserve/histology/coordinates` | `sgACC_cell_depth_analysis/coordinates` |
| Xenium object → spatial | `spatial/output/all_samples_annotated.h5ad` | `SCZ_Xenium/output/all_samples_annotated.h5ad` (schema: [`spatial/all_samples_annotated_guide.md`](spatial/all_samples_annotated_guide.md)) |

## Xenium dataset state (which object the numbers come from)

Everything spatial in this paper derives from **one** annotated object, and there
is more than one on disk. Pin it before trusting any Xenium number.

| | |
|---|---|
| **Canonical** | `SCZ_Xenium/output/all_samples_annotated.h5ad` — the **2026-04-01** state, reinstated 2026-07-31 |
| Identity | md5 `763e6655fc55839f177d57aa98dc5453`, byte-identical to `SCZ_Xenium/archive/2026-04-01/all_samples_annotated_april_augmented.h5ad` |
| In this repo | `spatial/output/all_samples_annotated.h5ad` — symlink to the canonical file; obs schema, QC gates and `lieber_*` semantics in [`spatial/all_samples_annotated_guide.md`](spatial/all_samples_annotated_guide.md) |
| Analysis set | **356,313** neuronal cortical cells (`qc_pass & corr_qc_pass & spatial_domain == 'Cortical'`, neuronal by subclass prefix) |
| Lieber columns | joined by within-sample position, 94.31% matched, using the **corrected** map **cluster 9 = MGE, cluster 12 = CGE** (S.H. Kwon, email 2026-07-27) |
| Override | set `XENIUM_BASE` if the processing repo is not at `~/Github/SCZ_Xenium` |

⚠️ **The archived June object carries the transposed Lieber map** (cluster 9 → CGE,
cluster 12 → MGE) and no donor metadata. It gives 373,514 neuronal cortical cells
instead of 356,313. Anything derived from it has MGE and CGE swapped — which
silently inverts the Supplementary Figure S2a validation claim. Do not reinstall it
without re-deriving the `lieber_*` columns.

The guard against drifting off this state is `shared/verify_provenance.py`: each
component records a `MANIFEST.tsv` naming the upstream files its outputs were built
from, and the script re-checks all of them. Run it before trusting a figure or
committing a result:

```bash
python3 shared/verify_provenance.py
```

## External inputs (documented, not re-wired here)

Per `~/Github/DATA_LAYOUT.md`, `genetics/` and `spatial/` symlink large
SEA-AD / MERFISH references from `~/Github/shared_data/` and
project-canonical locations. These stay external (the hub serves non-SCZ
projects too) and are re-created when each pipeline is set up to run in place.

## Build order (topological)

1. `snrnaseq/` → exports DE + composition betas into `shared/snrnaseq_de/`
2. `spatial/` (Xenium pipeline → DE, crumblr, depth) and `crossdisorder/`
   (AD slopes); both need only their own raw data plus the references
3. `genetics/` (needs composition betas + AD slopes + GWAS + the SEA-AD DLPFC reference)
4. `transcriptomic/` (needs DE betas + Xenium DE and depth)
5. `snrnaseq/composition_sensitivity/` (needs per-donor counts + Xenium crumblr)

Figures are written in place by the renderers at each step, into
`manuscript/figures/`; there is no separate assembly step.

## Current state

Code-complete. **Data is symlinked to the original `~/Github/` repos** (and
`~/Downloads/`) — the accepted setup: the monorepo holds the code, the source
repos hold the data, exactly as data was always external + git-ignored. (Note:
the symlinks resolve on this machine; sharing the repo elsewhere would need the
data transferred separately, same as the source repos always did.) The **one
seam to repoint when `snrnaseq/` is added** is `shared/snrnaseq_de/` → her
pipeline's outputs.
