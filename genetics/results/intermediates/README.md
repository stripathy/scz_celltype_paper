# Figure 4 normative-expression provenance

**Every normative (non-disease) expression value in Figure 4 comes from a single
dataset: the SEA-AD neurotypical MTG snRNA-seq reference** — 5 neurotypical
donors, 137,303 nuclei × 36,601 genes, the same reference that defines the 137
supertype taxonomy used throughout the paper (Gabitto et al. 2024).

Canonical file (too large to track; kept outside the repo):

    /Users/shreejoy/Github/shared_data/nicole_sea_ad_snrnaseq_reference.h5ad   # raw counts

Normalisation, applied identically everywhere: **counts-per-10,000 + log1p**
(equivalent to Seurat `NormalizeData`).

## Which panel uses what

| Panel | Quantity | Path from the canonical reference |
|---|---|---|
| a, b | MAGMA gene×cell-type specificity | `scripts/01_compute_specificity.py` |
| d, e | per-supertype mean *HCN1* expression | `seaad_sst_supertype_mean_expression.csv` (this directory) |
| h, i | depleted vs not-depleted marker DE, *CALB1* | `scripts/figures/export_fig4_new_panels.py` (loads the h5ad directly; cell-level data is needed, which a per-supertype mean cannot provide) |
| c | GWAS only | — no expression |
| e (sag), f, g | patch-seq physiology / morphology | separate modality, joined by supertype label |
| j | crumblr composition | — no expression |

## Consistency check (re-run if the reference ever changes)

The three routes above are independent code paths, so they are cross-checked on
one value — mean *HCN1* in Sst_25, which must equal **2.4404** in all of:

- `seaad_sst_supertype_mean_expression.csv` (this directory)
- `../figures/r_panels/panel_E_hcn1_vs_sag.csv` → `HCN1_expr`
- a fresh CP10K+log1p normalisation of the canonical h5ad

Verified 2026-07-26: all three agree to 0 decimal places of difference
(max abs diff 0.00e+00). *CALB1* in Sst_25 = 1.1622 by the same routes.

## Files here

- `seaad_sst_supertype_mean_expression.csv` — 36,601 genes × 16 Sst supertypes,
  mean log1p(CP10K) per supertype. A slice of the full 137-supertype matrix
  produced by `scripts/01_compute_specificity.py`; only the Sst columns are
  committed (5.8 MB vs 43 MB) because Figure 4 uses only Sst supertypes.
  The full matrix and the 46 MB specificity matrix remain in
  `~/Github/scz_cell_type_enrichment/results/intermediates/` and are
  regenerable from the canonical h5ad.

## Known remaining inconsistency

Panel **a** reports enrichment computed over the combined **Franken-503**
taxonomy (137 SEA-AD supertypes + 366 non-redundant Siletti clusters; see
Supplementary Methods SM2), so its specificity values are normalised against a
universe that includes subcortical types. Panels b/d/e/h/i are SEA-AD-137 only.
This is deliberate — Supplementary Table T4b holds the cortical-only companion
analysis — but it is the one place where Figure 4 is not on a single reference
universe.
