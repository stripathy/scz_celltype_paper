# Figure 5 + S10 — Methods text

Two pieces, per the agreed split. Every number verified against the pipeline
outputs on 2026-09-01 (`transcriptomic/scripts/fig5/09_verify.R` and the
supplement scripts). Per your decisions: no gene-permutation caveat, no
linear-trend model, no statement about when the grouping was fixed, and no
per-section code pointer (single code-availability statement elsewhere).

---

## (A) Insert into `Cell type-specific differential gene expression`

Place after the existing sentence ending "…only genes with ≥ 1 count in at least
80% of retained donors were tested." — i.e. as the third aggregation level,
before the TMM/voom sentence, or as its own short paragraph at the end of the
**snRNA-seq** block.

> **Depletion-group aggregation.** The same pipeline was applied at a third
> aggregation level, in which nuclei were pooled across all supertypes of a given
> Sst depletion group (defined in *Gene-set enrichment and depletion-group
> interaction analyses*) to give one pseudobulk profile per donor per dataset.
> Donor and gene retention, the model (1), the meta-analysis and the BH-FDR
> correction were as above, with the correction applied within each group. This
> retained 395, 411 and 342 donors and 44,744, 42,814 and 13,157 nuclei for the
> depleted, intermediate and non-depleted groups respectively; the pooled "All
> Sst" comparison applied the same procedure across all 16 Sst supertypes. Gene
> identifiers were harmonised to symbols before pooling, with the two
> Ensembl-native datasets (HBCC, MSSM2) mapped using GENCODE v44 and counts for
> duplicate symbols summed.

*(~120 words)*

---

## (B) New section, placed after `Patch-seq electrophysiology and morphology`

> ## Gene-set enrichment and depletion-group interaction analyses
>
> **Depletion groups.** The depleted group is the five-supertype set defined
> above (Molecular characterisation of the depleted Sst supertypes): Sst_2,
> Sst_3, Sst_20, Sst_22 and Sst_25, reduced in SCZ at FDR < 0.20. The eleven
> not-depleted supertypes were further split by the sign of their compositional
> estimate into an intermediate group, the six with a negative but
> non-significant estimate, and a non-depleted group, the five with an estimate
> ≥ 0. The FDR < 0.20 threshold was used rather than 0.10 so as to include
> Sst_20 (β = −0.17, FDR = 0.19), which shares the laminar position and
> genetic-risk profile of the other four. The three groups order from upper to
> deep cortical layers (Fig. 3f).
>
> **Gene-set enrichment.** Gene-set enrichment analysis was performed with fgsea
> [(Korotkevich et al. 2021)] on genes ranked by the meta-analytic z statistic
> (estimate/standard error), separately for each depletion group and for the
> pooled Sst subclass. Gene sets were taken from the Molecular Signatures
> Database (MSigDB) release 2025.1.Hs [(Liberzon et al. 2015)], accessed through
> the msigdbr R package (v25.1.1), and comprised the Gene Ontology biological
> process, cellular component and molecular function collections as distributed
> in that release, together with Reactome; the Hallmark collection was excluded
> because its sets are coarse signatures assembled across pathways and are
> difficult to interpret alongside named pathways. Sets were
> restricted to 10–500 genes, giving 6,024–6,839 sets tested per group, and
> significance was assessed at BH-FDR < 0.10 computed within each group. The gene
> sets shown in Fig. S8d and S8e are representative members of four blocks; the
> module assignments used to colour Fig. S8c are the union of the leading-edge
> genes of each block in the depleted group.
>
> **Software.** Analyses used R 4.5.1 with limma 3.64.3, edgeR 4.6.3, metafor
> 4.8.0, fgsea 1.34.2 and msigdbr 25.1.1.

*(~560 words)*

---

## References to add

- Korotkevich G, Sukhov V, Budin N, et al. Fast gene set enrichment analysis.
  bioRxiv 2021. doi:10.1101/060012  *(fgsea)*
- Liberzon A, Birger C, Thorvaldsdóttir H, et al. The Molecular Signatures
  Database hallmark gene set collection. Cell Syst 2015;1:417–425.  *(MSigDB /
  msigdbr — check whether the paper already cites MSigDB elsewhere)*
- limma, edgeR and metafor are already cited in the DE section.

## Notes

  numbered equation intervenes. Check against the Doc before pasting.
- The within-donor interaction model and the nuclei-count control were removed from the
  Methods on 2026-09-01 together with the controls supplement (former S10). The analyses
  stand in the repo (`05_interaction.R`, `05b_interaction_gsea.R`, `supp/subsample_cells.py`,
  `supp/sensitivity_cellsubsample.R`, `supp/figS_strata_controls.R`) and the text can be
  restored from git history if a reviewer asks.
- The GSEA paragraph is the paper's only description of gene-set enrichment,
  since no other analysis in the manuscript uses it.
- The depth correlation previously quoted here (ρ = 0.74, P = 0.0016, median Xenium depth) was removed on 2026-09-01: Fig. 3f already reports the same relationship with mean depth (ρ = 0.76), and two ρs for one relationship invites confusion. Decide whether Fig. S8a should adopt Fig. 3f's statistic.
