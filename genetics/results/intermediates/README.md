# Figure 4 normative-expression provenance

**Every normative (non-disease) expression value in Figure 4 comes from the
SEA-AD neurotypical DLPFC (A9) snRNA-seq reference** -- 3 neurotypical donors
(H18.30.002, H19.30.001, H19.30.002; `*_SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad`),
90,579 nuclei, 125 supertypes (Gabitto et al. 2024). This is region-matched to
the frontal-cortex SCZ cohorts. It replaced the SEA-AD MTG reference on
2026-08-12 (panels d/e/h/i) and, for the enrichment panels, on 2026-08-31 when
the combined SEA-AD + Siletti taxonomy was retired on L. Duncan's advice. The
MTG reference (5 donors, 137 supertypes) now appears only as the
reference-region robustness facet of Supplementary Fig. S10.

Canonical files (too large to track; kept outside the repo):

    ~/Downloads/H18.30.002_SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad   # raw counts, one per donor
    ~/Downloads/H19.30.001_SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad
    ~/Downloads/H19.30.002_SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad
    ~/Github/shared_data/nicole_sea_ad_snrnaseq_reference.h5ad             # MTG, S10 robustness only

Normalisation, applied identically everywhere: **counts-per-10,000 + log1p**
(equivalent to Seurat `NormalizeData`).

## Which panel uses what

| Panel | Quantity | Path from the canonical reference |
|---|---|---|
| a, b (+ S9, S10) | MAGMA gene x cell-type specificity across the 125 DLPFC supertypes | `a9_supertype_log1p_mean.csv` (from `scripts/figures/build_dlpfc_specificity.py`) -> `scripts/figures/build_spec_seaad_only.py` -> `spec_dlpfc_a9only.txt` -> `T_a9only_{bigdeli,pgc3}.gsa.out` |
| d | per-supertype mean *HCN1* expression | `a9_sst_supertype_mean_expression_cp10k.csv` (this directory), written by `scripts/figures/export_panels_dehi.py` |
| g, h | depleted vs not-depleted marker DE, *CALB1* | `scripts/figures/export_panels_dehi.py` (loads the three h5ads directly; cell-level data is needed) |
| c | GWAS only | -- no expression |
| d (sag), e, f | patch-seq physiology / morphology | separate modality, joined by supertype label |
| i | crumblr composition | -- no expression |

The S10 MTG facet uses `seaad_supertype_log1p_mean.csv` (MTG, 137 supertypes)
through the same `build_spec_seaad_only.py` recipe -> `spec_mtg_only.txt` ->
`T_mtgonly_{bigdeli,pgc3}.gsa.out`.

## Consistency check (re-run if the reference ever changes)

Panels d/e and h/i are independent code paths over the same h5ads, so they are
cross-checked on one value -- mean *HCN1* in Sst_25, which must equal
**2.3874** (CP10K + log1p, DLPFC) in both:

- `a9_sst_supertype_mean_expression_cp10k.csv` (this directory)
- `../figures/r_panels/panel_E_hcn1_vs_sag.csv` -> `HCN1_expr`

The MTG value this replaced was 2.4404 (`seaad_sst_supertype_mean_expression.csv`,
kept for the record). `export_panels_dehi.py` prints both when it runs.

## Files here

Tracked (whitelisted in `genetics/.gitignore`):

- `T_a9only_bigdeli.gsa.out`, `T_a9only_pgc3.gsa.out`, `namemap_a9only.csv` --
  the MAGMA gene-property runs behind Fig. 4a/b, S9 and the DLPFC facets of S10,
  and the safe-name -> supertype map they need.
- `T_mtgonly_bigdeli.gsa.out`, `T_mtgonly_pgc3.gsa.out`, `namemap_mtgonly.csv` --
  the MTG robustness facets of S10.
- `a9_sst_supertype_mean_expression_cp10k.csv` -- 36,601 genes x 16 Sst
  supertypes, mean log1p(CP10K) per supertype, DLPFC. Feeds panels d/e.
- `seaad_sst_supertype_mean_expression.csv` -- the MTG predecessor of the above
  (5 neurotypical donors). No longer read by Figure 4, but still the input to
  the reserve MTG-vs-DLPFC baseline check in `reserve/sst_strata_supp/`.
The retired combined SEA-AD + Siletti (501-type) runs -- `T_a9rbh_bigdeli.gsa.out`,
`namemap_a9rbh.csv`, `T_a9_pgc3.gsa.out`, `T_mtg_bigdeli.gsa.out`,
`namemap_mtg_sametax.csv`, `franken_rbh_A9.csv` -- were untracked on 2026-09-04
along with their build scripts. They are recoverable from the git tag
`pre-prune-2026-09-04`.

Not tracked (regenerable with `build_dlpfc_specificity.py`, `seaad_supertype_log1p.py`
and `build_spec_seaad_only.py`): the per-supertype mean matrices
(`a9_supertype_log1p_mean.csv` 68 MB, `seaad_supertype_log1p_mean.csv` 74 MB) and
the specificity matrices (`spec_dlpfc_a9only.txt` 38 MB, `spec_mtg_only.txt`, and the
retired `spec_dlpfc_a9rbh.txt` / `spec_mtg_sametax.txt`).

## Reference universes

Since 2026-08-31 every expression-derived panel of Figure 4 is on the same
reference: the 125-supertype SEA-AD DLPFC taxonomy. Panel a's specificity is
normalised across those 125 types, and S9 shows all 125, so the figure and the
multiple-testing universe coincide. The earlier note that panel a alone sat on a
501-type universe including subcortical types no longer applies.
