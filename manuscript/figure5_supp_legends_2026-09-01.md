# Figure 5 supplementary figure legend

One combined supplement (was two). Numbered S10 per the current supplement
ordering.

---

**Fig. S10 | Cell-count control and diagnosis × depletion-group interaction for Sst interneurons.**
(**a**) Case–control difference in nuclei per donor within each depletion group,
before and after cell matching (mean ± 95% CI, from a within-cohort linear model
on log nuclei count). Because the groups are defined by reduced abundance in SCZ,
cases contribute progressively fewer nuclei across the three groups as observed
(−29.5%, −17.9% and −4.0%); subsampling control nuclei to the case distribution
removes this (−0.4%, +1.8% and +4.9%). (**b**) Gene sets significant at
FDR < 0.10 in each group. Grey bars, the analysis as reported in Fig. 5b; red
points, five independent cell-matched draws. Matching is performed within cohort
and group by assigning each control donor a target nuclei count drawn from the
case distribution at that donor's percentile, subsampling their nuclei without
replacement, and rebuilding the pseudobulks from cell-level counts before
repeating the differential expression and gene-set enrichment of Fig. 5 unchanged
(Methods). (**c**) Interaction normalized enrichment scores for the gene sets
shown in Fig. 5d, together with the two sets named in the text (diamonds).
Negative values indicate that the effect of diagnosis is more negative in the
depleted group than in the non-depleted group. Each donor contributes one
pseudobulk per depletion group, so the contrast is taken within donor and
donor-level factors such as medication exposure, age and postmortem interval
cannot generate it (Methods). Significance throughout: \*\*\*FDR < 0.01,
\*\*FDR < 0.05, \*FDR < 0.10, +FDR 0.10–0.20.

---

## Text changes needed

1. The downsampling sentence already cites `(Fig. SXXX; Methods)` — point it at
   this figure.
2. The interaction sentence cites no figure. Add one:

   > …and oxidative phosphorylation (electron transport chain, FDR = 0.053)
   > **(Fig. S10c)**.

3. **The SRP FDR in the text is stale.** It reads 4.4 × 10⁻¹⁰; the current build
   gives **7.9 × 10⁻¹⁰**. The value changed when Hallmark was dropped from the
   gene-set universe, which altered the multiple-testing correction. The
   electron-transport-chain value (0.053) is unchanged.

   Consider citing **ribosomal subunit (NES = −2.14, FDR = 3.4 × 10⁻⁶)** instead
   of, or alongside, SRP: it is plotted in Fig. 5d, so a reader can find it,
   whereas SRP appears only in the supplement.
