# Franken taxonomy — canonical artifacts

This directory is the **canonical home** of the combined SEA-AD ⊕ Siletti
("Franken") cell-type taxonomy. If you want to use the Franken taxonomy
in another project, **read this directory** — these two CSVs are the
contract.

For the full methodology (MetaNeighbor reciprocal best hits, AUROC
thresholds, validation against the older 20%-cortical taxonomy), see
[`docs/combined_taxonomy_approach.md`](../docs/combined_taxonomy_approach.md).

## Files

| File | Rows | Description |
|---|---|---|
| `franken_types.csv` | 598 | Provenance table: every candidate cell type (137 SEA-AD + 461 Siletti) with its source, RBH partner if any, and a flag for whether it survived dedup into the combined Franken taxonomy. |
| `franken_rbh.csv` | 95 | The reciprocal-best-hit pairs that drove the dedup. Mirror of `results/tables/metaneighbor_full461_reciprocal_best_hits.csv`. |
| `build_franken_types.py` | — | Builder script that regenerates `franken_types.csv` from the upstream MetaNeighbor outputs. |

## Schemas

### `franken_types.csv`

| Column | Type | Description |
|---|---|---|
| `type` | str | Cell-type name (e.g. `Pvalb_3`, `MGE_259`, `Misc_132`). |
| `source` | str | `seaad` or `siletti`. |
| `rbh_partner` | str | Partner type name in the other dataset, if this row is part of a reciprocal best hit. Empty otherwise. |
| `rbh_mean_auroc` | float | Mean AUROC of the RBH pair (≥0.9 by construction; 1.0 is common). NaN if not RBH. |
| `best_match` | str | (SEA-AD rows only) The closest Siletti cluster by AUROC, even when not reciprocal. |
| `best_match_auroc` | float | (SEA-AD rows only) AUROC of `best_match`. |
| `is_reciprocal` | bool | True iff this type is part of an RBH pair. |
| `in_combined_franken` | bool | True iff this type is one of the 503 types kept in the combined taxonomy. |

The dedup rule is simple:

```
in_combined_franken = (source == "seaad")  OR  (source == "siletti" AND NOT is_reciprocal)
```

Counts:
- 137 SEA-AD types — all kept (`in_combined_franken = True` for all).
- 461 Siletti clusters — 95 dropped (RBH-claimed), 366 kept as novel.
- Combined Franken: 137 + 366 = **503 types brain-wide**.

Downstream consumers may further restrict the 503 (e.g. RSC Xenium
applies an RSC region-of-interest whitelist on the Siletti side, ending
up with 176 types). That's a *consumer*-level decision; this file
defines the brain-wide source of truth.

### `franken_rbh.csv`

| Column | Type | Description |
|---|---|---|
| `seaad_type` | str | The SEA-AD supertype name. |
| `siletti_cluster` | str | The Siletti cluster that is its reciprocal best hit. |
| `auroc_seaad_to_siletti` | float | AUROC for SEA-AD → Siletti direction. |
| `auroc_siletti_to_seaad` | float | AUROC for Siletti → SEA-AD direction. |
| `mean_auroc` | float | Average of the two directions. |

This is a clean copy of `results/tables/metaneighbor_full461_reciprocal_best_hits.csv`,
intended as the stable consumer-facing path.

## How to use these in a downstream project

```python
import pandas as pd

types = pd.read_csv("franken_taxonomy/franken_types.csv")

# All 503 types in the combined Franken taxonomy:
franken = types[types["in_combined_franken"]]

# Just the SEA-AD half (137):
sea = types[(types["source"] == "seaad") & types["in_combined_franken"]]

# Just the novel Siletti clusters (366):
sil_novel = types[(types["source"] == "siletti") & types["in_combined_franken"]]

# RBH pairs (95):
rbh = pd.read_csv("franken_taxonomy/franken_rbh.csv")
```

## How the artifacts are derived

`franken_types.csv` is rebuilt by `build_franken_types.py` from three
upstream MetaNeighbor outputs in `results/tables/`:

1. `metaneighbor_full461_reciprocal_best_hits.csv` — the 95 RBH pairs.
2. `metaneighbor_full461_all_matches.csv` — every SEA-AD type's best
   Siletti match (whether reciprocal or not).
3. `rbh_combined_enrichment.csv` — the 503-type combined enrichment
   table; used as the authoritative list of what's in the combined
   Franken.

Those upstream tables are produced by `scripts/09_metaneighbor_integration.py`
and `scripts/10_combined_taxonomy.py`. Run those first if you ever need
to regenerate from raw inputs (precomputed_stats.h5 files); otherwise
just rerun:

```bash
python3 franken_taxonomy/build_franken_types.py
```

## Consumers

These artifacts are vendored downstream by:

- [`RSC_Xenium`](https://github.com/stripathy/RSC_Xenium) — uses
  `franken_rbh.csv` for cell-typing centroid construction. See
  [`RSC_Xenium/output/reference/franken/`](https://github.com/stripathy/RSC_Xenium/tree/main/output/reference/franken).

If you're vendoring these into a new project, copy the CSVs into your
own repo so the downstream build doesn't depend on cloning this one.
Note the upstream commit you pulled them from (we recommend a `SOURCE.md`
file alongside the vendored copies).

## Citation

If you use the combined Franken taxonomy, please cite both source
datasets:

- **SEA-AD MTG**: Gabitto et al. 2024, *Nature Neuroscience*
  ([doi:10.1038/s41593-024-01774-5](https://doi.org/10.1038/s41593-024-01774-5)).
- **Siletti Whole Human Brain**: Siletti et al. 2023, *Science*
  ([doi:10.1126/science.add7046](https://doi.org/10.1126/science.add7046)).
