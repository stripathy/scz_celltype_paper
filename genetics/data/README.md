# Inputs for Figure 4, S9, S10 and Table T6

Most of this directory is git-ignored: GWAS summary statistics, the MAGMA
binary and reference panel, and the patch-seq recordings are far too large to
track. What **is** committed is the small set of files the exporters cannot
regenerate on their own, so a clone can rebuild every panel CSV once the large
inputs below are in place.

## Committed here

| Path | What it is |
|---|---|
| `seaad_supertype_colors.json` | SEA-AD supertype -> hex colour, 139 entries. The palette for every Figure 4 scatter and for S9. |
| `patchseq/swc/*.swc` (5 files, 4.3 MB) | Morphological reconstructions of the five exemplar Sst cells drawn in Figure 4e. See `patchseq/README.md`. |

## Not committed — fetch or point at these

| Path | Size | Where it comes from |
|---|---|---|
| `gwas/bigdeli_eur_scz_sum_stats.gz` | 329 MB | Bigdeli et al. 2026, European-ancestry SCZ GWAS. **The paper's primary GWAS.** |
| `gwas/magma_bigdeli/bigdeli.step2.genes.{raw,out}` | 1.5 MB | MAGMA gene analysis of the above (step 2). |
| `gwas/magma_pgc3/pgc3.step2.genes.raw` | | Same for PGC3 (Trubetskoy et al. 2022); robustness only, S10. |
| `gwas/ncbiRefSeq_hg38.txt.gz` | 7.0 MB | UCSC RefSeq track, hg38. Gene track under the HCN1 locus zoom (Fig. 4c). `wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz` |
| `magma/magma_mac/magma` | | MAGMA v1.10 binary, <https://cncr.nl/research/magma/> |
| `patchseq/nwb/*.nwb` (5 files) | 148 MB | Intracellular recordings for Figure 4f. DANDI dandiset 000636, or the `human_int_patch_seq` cache. See `patchseq/README.md`. |

The gene-location file with the extended MHC excluded
(`NCBI37.3.gene.loc.extendedMHCexcluded`) is read from a sibling checkout of
`scz_cell_type_enrichment`; it ships with the MAGMA reference bundle.

## The expression reference

Specificity and the marker panels come from the **SEA-AD neurotypical DLPFC
(A9)** snRNA-seq reference: three donors, 90,579 nuclei, 125 supertypes
(Gabitto et al. 2024). Those three h5ads are expected at
`~/Downloads/*_SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad` and can be
downloaded from
<https://sea-ad-single-cell-profiling.s3.amazonaws.com/index.html#DLPFC/RNAseq/>.

The MTG reference (5 donors, 137 supertypes,
`~/Github/shared_data/nicole_sea_ad_snrnaseq_reference.h5ad`) is used **only**
for the reference-region robustness facets of Supplementary Fig. S10.

See `../results/intermediates/README.md` for which panel draws on which, and
`../README.md` for the full build order.
