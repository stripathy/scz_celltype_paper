# Data flow & dependency graph

How the four analysis components depend on each other. **Every cross-component
edge is realized as a git-ignored symlink** (so the wiring is executable, not
just prose) and listed in the Seams table. The upstream root — the student's
snRNA-seq pipeline — is not in the repo yet; its outputs are symlinked into
`shared/snrnaseq_de/` from their current scattered locations.

## Dependency graph (data flows downward)

```
              snrnaseq/   ← student's pipeline (Endresz et al.); NOT in repo yet
            (upstream root)
                 │  exports → shared/snrnaseq_de/
                 │    DE_genes_all_cells_scz.csv         (meta DE)
                 │    meta_results_cohorts_subclass.csv  (per-cohort DE)
                 │    nicole_scz_snrnaseq_betas/         (composition, crumblr)
   ┌─────────────┼──────────────────┬───────────────────────┐
   │ DE          │ DE + composition │ composition           │
   ▼             ▼                  ▼                        │
transcriptomic/  spatial/         genetics/                 │
   ▲   ▲         │  Xenium DE ──►  │  GWAS set ──► transcriptomic/ (script 04)
   │   └─────────┘  (script 09)    │  franken_taxonomy/ (canonical cell-type axis)
   │                               │  enrichment
   └───────────────┬──────────────┘
                   ▼
                figures/  = INDEX of the paper's final figures (each component
                           produces its own, e.g. the composite in
                           transcriptomic/); no separate convergence figure.
                   ▲
histology/ (standalone RNAscope FISH SST density) ──────────┘

External hub (~/Github/shared_data/, per ~/Github/DATA_LAYOUT.md):
   SEA-AD / Siletti / MERFISH references ──► genetics/, spatial/  (symlinked, not absorbed)
```

## Per-component inputs → outputs

| component | key inputs (← from) | produces → (consumed by) |
|---|---|---|
| `snrnaseq/` *(TBD)* | raw 7-cohort snRNA-seq | DE betas, composition betas → everyone downstream |
| `genetics/` | composition betas (← snrnaseq), GWAS set (owns), taxonomy (owns), SEA-AD ref (← shared_data) | cell-type enrichment, `franken_taxonomy/`, GWAS set → transcriptomic, figures |
| `spatial/` | Xenium raw, DE + composition betas (← snrnaseq), SEA-AD/MERFISH ref (← shared_data) | Xenium DE, crumblr composition → transcriptomic, figures |
| `transcriptomic/` | DE betas (← snrnaseq), GWAS set (← genetics), Xenium DE (← spatial) | composite figure, GSEA/pathway results → figures |
| `histology/` | RNAscope FISH counts (self-contained) | SST density results → figures |
| `figures/` | — (index only) | pointers to each component's final paper figures |

## Internal seams (realized symlinks = the edges)

| edge | path in repo | → resolves to (current) |
|---|---|---|
| snrnaseq → transcriptomic | `transcriptomic/data/DE_genes_all_cells_scz.csv` → `shared/snrnaseq_de/` | `~/Downloads/DE_genes_all_cells_scz.csv` |
| snrnaseq → transcriptomic | `transcriptomic/data/meta_results_cohorts_subclass.csv` → `shared/snrnaseq_de/` | `scz_pathway_enrichment/data/…` |
| snrnaseq → genetics + spatial | `shared/snrnaseq_de/nicole_scz_snrnaseq_betas/` | `SCZ_Xenium/data/nicole_scz_snrnaseq_betas/` |
| ↳ spatial consumes | `spatial/data/nicole_scz_snrnaseq_betas` → `shared/snrnaseq_de/` | (as above) |
| ↳ genetics consumes | `genetics/scripts/13` code → `shared/snrnaseq_de/…` (repointed) | (as above) |
| spatial → transcriptomic | `transcriptomic/scripts/09` reads `../spatial/output/de/de_results_subclass.csv` | `SCZ_Xenium/output/de/…` (symlinked) |
| genetics → transcriptomic | `transcriptomic/scripts/04` reads `../genetics/data/gwas/…`; `genetics/data/gwas/*` symlinked | `scz_cell_type_enrichment/data/gwas/…` |
| histology data | `histology/coordinates` | `sgACC_cell_depth_analysis/coordinates` |

## External inputs (documented, not re-wired here)

Per `~/Github/DATA_LAYOUT.md`, `genetics/` and `spatial/` symlink large
SEA-AD / Siletti / MERFISH references from `~/Github/shared_data/` and
project-canonical locations. These stay external (the hub serves non-SCZ
projects too) and are re-created when each pipeline is set up to run in place.

## Build order (topological)

1. `snrnaseq/` → exports DE + composition betas into `shared/snrnaseq_de/`
2. `genetics/`, `spatial/` (need composition + refs; genetics also GWAS + taxonomy)
3. `transcriptomic/` (needs DE + GWAS + Xenium DE)
4. `histology/` (independent — any time)
5. `figures/` — index of the final figures each component already produces

## Current state

Code-complete. **Data is symlinked to the original `~/Github/` repos** (and
`~/Downloads/`) — the accepted setup: the monorepo holds the code, the source
repos hold the data, exactly as data was always external + git-ignored. (Note:
the symlinks resolve on this machine; sharing the repo elsewhere would need the
data transferred separately, same as the source repos always did.) The **one
seam to repoint when `snrnaseq/` is added** is `shared/snrnaseq_de/` → her
pipeline's outputs.
