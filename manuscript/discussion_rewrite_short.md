# Discussion — half-length version (2026-08-10)

**Companion to `discussion_rewrite_draft.md` (the full version).** Same claims and locked
framing, ~40% shorter. Hedges live only in Limitations; no references to our own figures;
statistics cut to the load-bearing minimum (all still verbatim from the Doc Results).
**2026-08-11:** ¶1 replaced with Shreejoy's own rewrite (one grammar fix: ", a view" inserted
so "based largely on in situ evidence" attaches to the view, not the transcripts); new ¶2
("how the data were combined") added per his plan — cite-free, carries the ~26%-invisible-at-
subclass payoff and the pooling/replication logic; L6b paragraph now opens with the finding,
since ¶1 no longer introduces L6b.

---

Whether cortical inhibitory interneurons in schizophrenia (SCZ) are lost or persist in an altered molecular state has been contested for two decades. The dominant consensus view is that interneurons are altered in their transcriptional state without detectable cell loss, including lowered SST and GABA-synthetic transcripts, a view based largely on in situ evidence (Volk et al., 2000; Hashimoto et al., 2003; Dienel et al., 2022, 2023). In contrast, a competing view is that such cells are indeed lost in SCZ, supported by histological counting and large-scale transcriptomic studies (Benes et al., 1991; Beasley et al., 2002; Toker et al., 2018; Batiuk et al., 2022; Kiss et al., 2026). Our data, comprising a seven-dataset snRNAseq meta-analysis and spatial transcriptomics-based replication, indicate that both accounts are correct, in partly different cells. Namely, we find a subclass-wide reduction of SST mRNA within Sst neurons that co-occurs with the depletion (i.e., reduced abundance) of specific upper-layer Sst supertypes.

These findings were resolvable because of how the data were combined. Mapping every dataset onto a common fine-grained taxonomy separated closely related Sst subtypes rather than averaging them: a ~26% reduction (95% CI 14–36%) of a single supertype (Sst_25) is nearly invisible at subclass resolution, explaining how coarser studies could disagree. Meta-analytic pooling supplied power that no individual dataset had — the Sst_25 depletion reached nominal significance in just one of the seven datasets on its own — and spatial replication guarded the result against platform- and dissociation-specific artifacts that single-modality studies cannot exclude.

The depletion is tiered and laminar: Sst_2, Sst_3 and Sst_25 were FDR-significant in the meta-analysis, replicated in Xenium, and survived re-annotation excluding SCZ-DE genes, whereas Sst_22 and Sst_20 (the latter nominal-only) did not replicate spatially and we treat them as provisional. More superficial Sst supertypes were more depleted, and an orthogonal RNAscope re-analysis of subgenual anterior cingulate (Arbabi et al., 2025) found a trend-level SST density reduction, likewise larger in upper layers.

The SST reduction held across both platforms; PVALB was unchanged in snRNAseq but nominally reduced in Xenium, consistent with a partly cytoplasmic transcript that whole-cell Xenium captures and nucleus-restricted snRNAseq under-samples [TODO: cite PVALB mRNA localization]. Because SST transcription is activity- and BDNF-dependent [TODO: cite Tripp et al.; Guilloux et al.], the pan-Sst reduction, though measured within Sst cells, may index reduced excitatory drive onto them rather than a cell-autonomous defect.

Reduced density of calbindin-immunoreactive interneurons was reported in SCZ prefrontal cortex more than twenty years ago (Beasley et al., 2002); our data localize that observation to named transcriptomic supertypes. The depleted supertypes share a molecular signature that includes elevated CALB1, and in primates CALB1-expressing interneurons include the double-bouquet cell, a dendrite-targeting morphology with no clear rodent counterpart to which the depleted supertypes partly correspond (Ballesteros-Yáñez et al., 2005; Raghanti et al., 2010). Because rodents lack this specialization, mouse models are blind to this axis of pathology; mechanistic follow-up will require human or primate systems.

Common-variant genetics places the depleted population closer to cause than consequence: germline risk is fixed before disease onset, and across the 16 Sst supertypes genetic enrichment tracked depletion (see Limitations for GWAS-version sensitivity). Decomposing that signal nominated HCN1 — fine-mapped as the likely causal gene at its locus, specifically expressed in Sst_25, with expression tracking both depletion and HCN-dependent voltage sag across supertypes. HCN1 loss-of-function causes an epileptic channelopathy [TODO: cite Nava et al., 2014] and I_h shapes dendritic integration [TODO: cite Magee], converging with double-bouquet inhibition on the dendritic compartments of upper-layer pyramidal cells.

SCZ cases also showed increased abundance of deep-layer L6b excitatory supertypes, replicated spatially. The increase parallels the long-standing excess of interstitial white-matter neurons in SCZ (Akbarian et al., 1993, 1996), but lies within gray-matter L6b — our spatial analysis excluded white matter — consistent with a subplate-lineage abnormality spanning the L6b/white-matter boundary [TODO: cite subplate development, e.g., Kanold & Luhmann, 2010] and, like the germline anchoring of the Sst depletion, with early-established rather than degenerative processes.

The depleted Sst supertypes are not vulnerable only in schizophrenia. Sst supertypes are among the earliest populations lost as Alzheimer's disease (AD) advances (Gabitto et al., 2024); whether the cells affected in SCZ are the same ones had not been tested. Because SEA-AD supplies the taxonomy used throughout our study, we compared exactly, supertype to supertype: SCZ depletion and decline along the AD pseudo-progression trajectory agreed in direction for 94% of the 16 Sst supertypes (Spearman ρ = 0.87). Cortical somatostatin deficits were among the first neurochemical findings in AD (Davies et al., 1980); the present correspondence sharpens that link to named supertypes. Depletion of one population by a germline-anchored, presumably early process in SCZ and by a late, pathology-graded process in AD is what a vulnerability intrinsic to the cells, rather than to either disease, predicts (Saxena and Caroni, 2011). Whether the same population also declines in normal aging is untested here.

Limitations. Because SCZ alters the transcripts that define cell identity, apparent depletion could reflect either neuronal loss or erosion of transcriptomic identity below confident classification; cross-platform replication and the DE-gene-excluded re-annotation narrow but do not eliminate that ambiguity; the abundance results therefore identify which cells are most affected, not whether they die or dedifferentiate. Cross-sectional sampling cannot order molecular and compositional change, nor separate disease from exposures collinear with diagnosis, notably antipsychotic medication — although germline risk fixed before treatment argues against a purely iatrogenic account, as would the AD-trajectory decline if SEA-AD donors were antipsychotic-naive [TODO: verify SEA-AD donor medication status]. The genetics–depletion convergence attenuates under the larger Bigdeli et al. (2026) GWAS. The SCZ–AD comparison crosses designs and demographics — a case-control contrast versus a pathology-graded progression score in aged donors — and identifies shared cellular targets, not a shared mechanism or time-course. The 300-gene Xenium panel limits supertype resolution, so subclass-level spatial conclusions are firmer than individual supertype allocations. Our data are predominantly DLPFC; generalization requires additional cortical regions.

Pairing fine-grained taxonomies with cross-dataset meta-analysis is a generalizable route to cell-type-resolved neuropathology where individual post-mortem datasets are small and heterogeneous. Next steps follow directly: expression-independent labels (protein, morphology, or lineage) to settle loss versus erosion; aging and longitudinal designs to test how far the vulnerability generalizes; and pharmacology aimed at the circuit these cells define — α5-subunit-containing GABA-A receptors, which mediate dendrite-targeting Sst inhibition [TODO: cite α5-PAM work, e.g., Prevot et al.], and HCN1/I_h itself — toward alternative drug targets for treatment-resistant SCZ. Upper-layer, CALB1-expressing, HCN1-high Sst interneurons, depleted in schizophrenia and declining along the AD trajectory, constitute an intrinsically vulnerable population — and now a named, testable one.

---

## What was cut relative to the full version

- **Moved into Limitations** (hedges allowed there per your rule): the standalone
  loss-vs-erosion paragraph; the SCZ–AD design/demographics caveats; the GWAS-version caveat
  (genetics paragraph now just points there).
- **Deleted hedges:** shared-taxonomy-cuts-both-ways; sub-threshold-power and
  conservative-bias circularity arguments; mixed calbindin histology (Daviss & Lewis);
  L6b label-drift defusal; "remains to be tested" for HCN1 cause-vs-marker; "We read SST
  as..." verdict sentence; the Hashimoto PV-mRNA reconciliation clause.
- **Deleted figure/supplement pointers:** all of them (Figs. 2j/3f/4a/c/f/h/i/j, S2/S13/S14,
  SM1/SM3).
- **Deleted statistics:** SST/PVALB Xenium P values, DE concordance r = 0.73/76%, depth
  ρ = 0.76, genetics ρ = 0.56 (also removes the PGC3-sync problem from this section entirely),
  HCN1 ρ = 0.79/0.65, 234-gene signature, L6b re-annotation βs. Retained: ~26% (CI 14–36%),
  94% / ρ = 0.87 / 16 supertypes, 300-gene panel.
- **Deleted entirely:** Fröhlich-OFC aside; meta-analysis-privileges-reproducible-effects
  limitation; Erraji-Benchekroun aging citation (sentence kept, citation dropped).

**If you need it shorter still**, in cut order: the Davies 1980 sentence; the Saxena &
Caroni sentence (keep the claim, lose the anchor); the sgACC sentence; the BDNF/activity
interpretation sentence.
