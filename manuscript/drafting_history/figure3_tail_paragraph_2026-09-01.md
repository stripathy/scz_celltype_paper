# Closing paragraph for the Figure 3 section — Sst depletion groups (2026-09-01)

Placement: end of "Upper-layer Sst supertypes are less abundant and L6b excitatory
supertypes more abundant in schizophrenia", after the paragraph ending "...per-supertype
abundance changes were positively, if modestly, correlated (Fig. 3e; Spearman ρ = 0.40,
P < 0.001)" and before the Fig. 3 image. Every number verified against the pipeline
2026-09-01. S8 is the former Figure 5, unchanged; the controls supplement (former S10)
is not in this version. If Fig. 3f is recolored by group, cite "Fig. 3f" in place of
"Fig. S8a".

---

**Working text (Shreejoy's drop-in of 2026-09-01, second pass; rodent-literature sentence
placed last, per the ordering decision, so the interpretive sentence follows the data it
interprets and the paragraph hands off to the Discussion).**

> We then asked whether the transcriptional state of Sst interneurons differs between
> the supertypes that are depleted and those that persist. We grouped the 16 Sst
> supertypes by their depletion status in Fig. 3a into three broad groups: depleted (FDR
> < 0.20; n = 5), intermediate (negative but non-significant; n = 6) and non-depleted (n
> = 5). We repeated the differential expression analysis of Fig. 2 followed by gene set
> enrichment analysis within each group (Fig. S8, Methods). We found that *SST* and
> *VGF* mRNA and synaptic gene sets were reduced in all three groups. In contrast, gene
> sets for cytosolic translation were downregulated only in the depleted and
> intermediate groups, and gene sets for oxidative phosphorylation were downregulated in
> the depleted group alone. Although these changes at the level of individual genes were
> modest (median log₂FC −0.11 for the leading-edge genes of both programs), these
> results suggest that while many aspects of transcriptional dysregulation are broadly
> shared across Sst interneurons, the suppression of protein synthesis and oxidative
> phosphorylation is largely confined to the supertypes undergoing depletion. This
> pattern is consistent with mechanistic work in rodents in which chronic stress
> suppresses protein translation and induces endoplasmic reticulum stress selectively in
> Sst interneurons (Lin and Sibille 2015; Tomoda et al. 2022).

<details><summary>Previous draft (superseded)</summary>

> We then asked whether the transcriptional disease state of Sst interneurons differs
> between the supertypes that are depleted and those that persist. We grouped the 16 Sst
> supertypes by their depletion in Fig. 3a into depleted (FDR < 0.20; n = 5),
> intermediate (negative but non-significant; n = 6) and non-depleted (n = 5) groups,
> which broadly track cortical depth (Fig. S8a), pooled each donor's nuclei within each
> group, and repeated the differential expression analysis of Fig. 2 followed by gene
> set enrichment analysis (Methods). *SST* and *VGF* mRNA and synaptic gene sets were
> reduced in all three groups. In contrast, gene sets for cytosolic translation were
> suppressed in the depleted and intermediate groups, and gene sets for oxidative
> phosphorylation in the depleted group alone, whereas no gene set reached significance
> in the non-depleted group (Fig. S8b–e). The underlying changes in individual genes
> were small (log₂FC −0.16 to −0.23 for the exemplar genes in Fig. S8e) and became
> detectable only as coordinated programs. Together, these results suggest that while
> transcriptional dysregulation is broadly shared across Sst interneurons, the
> suppression of protein synthesis and oxidative phosphorylation is largely confined to
> the supertypes undergoing depletion.

</details>

---

## What was deliberately left out, and why

- **The pooling-dilution sentence** (130 sets in the depleted group versus 55 pooled,
  "despite greater power"). True at the gene-set level, false at the gene level
  (19 DE genes in the depleted group versus 342 pooled). It is the most attackable
  claim in the current section and is not needed for the paragraph's point.
- **Gene-set counts** (130 / 64 / 0). The 130 sets collapse to roughly three
  programs; counts invite the redundancy objection. The counts remain as bar labels
  in the S8 legend, where they are descriptive.
- **Citations.** Two now sit in the paragraph at Shreejoy's request, both to add to Zotero:
  Lin LC, Sibille E. Somatostatin, neuronal vulnerability and behavioral emotionality.
  Mol Psychiatry 2015;20:377–387 (PMID 25600109); Tomoda T, Sumitomo A, Newton D, Sibille E.
  Molecular origin of somatostatin-positive neuron vulnerability. Mol Psychiatry
  2022;27:2304–2314 (doi 10.1038/s41380-022-01463-4). The Discussion sentences below
  can then drop their Lin and Sibille citation or keep it; either is fine.
- **NES ranges and exemplar gene names.** Moved to the S8 legend; the paragraph names
  programs, not sets.
- **The within-donor interaction test and the nuclei-count control** (former S10).
  Dropped with the controls supplement on 2026-09-01 to keep the section short; the
  analyses stay in the repo and can be revived if a reviewer raises the confound.

## Numbers behind each clause

| Clause | Value |
|---|---|
| Group sizes | 5 / 6 / 5 |
| Translation sets, depleted | NES −2.28 to −2.72, FDR ≤ 6.8 × 10⁻⁸ |
| Translation sets, intermediate | NES −1.92 to −2.72, FDR ≤ 5 × 10⁻⁴ |
| Translation sets, non-depleted | NES −0.58 to −0.72, FDR = 1.0 |
| OxPhos sets, depleted | NES −1.93 / −2.02, FDR 0.0032 / 0.00048 |
| OxPhos sets, intermediate / non-depleted | FDR 0.64–0.87 |
| Trans-synaptic signaling | NES −1.52 / −1.48 / −1.38 (FDR 0.023 / 0.078 / 0.46) |
| *SST* log₂FC by group | −0.43 / −0.51 / −0.33 |
| *VGF* log₂FC by group | −0.60 / −0.77 / −0.65 |
| Sets at FDR < 0.10, non-depleted | 0 |
| Exemplars RPL36, EIF3G, NDUFS8, COX5B (depleted) | log₂FC −0.16 to −0.23; FDR 0.27–0.58 |
| Leading-edge genes, translation block (depleted; n = 130) | log₂FC −0.04 to −0.31, median −0.11 |
| Leading-edge genes, oxphos block (depleted; n = 66) | log₂FC −0.03 to −0.25, median −0.11 |
| All 16 exemplars drawn in S8e (depleted) | log₂FC −0.05 to −0.27, median −0.15 |

## Companion edits this paragraph implies

1. **Discussion, intrinsic-vulnerability paragraph** (after "...raising the question of
   whether that specialization contributes to their vulnerability."), two sentences:

   > The depleted supertypes are also the Sst neurons whose protein-synthesis and
   > oxidative-phosphorylation programs are most suppressed in SCZ, consistent with
   > earlier evidence that cortical Sst neurons reduce translation under chronic
   > stress [(Lin and Sibille 2015)]. Whether this state precedes their loss, follows
   > it, or reflects which cells survive cannot be resolved in post-mortem tissue, but
   > the coordinated withdrawal of a neuron's costliest activities in the very cells
   > declining in number is the signature expected of a population under sustained
   > strain.

2. **Supplementary figures.** The current Fig. 5 image and legend move to the
   supplement as **Figure S8**, legend text unchanged apart from the title prefix
   (below). The Fig. S10 controls legend and image are removed. The genetics
   supplements shift by one: legend headers "Figure S8" → "Figure S9" (enrichment
   landscape) and "Figure S9" → "Figure S10" (robustness). The three in-text
   citations in the genetics section already read S9 and S10 and become correct
   without edits.
3. **Results section 5** (the Fig. 5 text, legend and image) is removed.
4. **Methods.** Insert (A) and the trimmed section (B) in
   `manuscript/figure5_methods_2026-09-01.md`; the interaction and nuclei-count
   subsections are gone, so no new equation is needed.
5. **Abstract and Introduction** need no change in this configuration.

### Figure S8 legend (the Fig. 5 legend as it stands in the Doc, retitled)

> **Figure S8 | Shared and depletion-graded transcriptional dysregulation of Sst
> interneurons in schizophrenia.** (**a**) Definition of the three Sst supertype
> groups based on depletion status. Per-supertype abundance change in SCZ (crumblr
> estimate ± SE, from Fig. 3a) against median cortical depth measured in Xenium
> (0 = pia), for the 16 Sst supertypes. Depleted, FDR < 0.20 (n = 5); intermediate,
> negative but non-significant (n = 6); non-depleted, estimate ≥ 0 (n = 5). Cortical
> layer boundaries based on Xenium sections. (**b**) Gene sets significantly enriched
> in each depletion group; bar labels denote significant gene set counts at
> FDR < 0.10. "All Sst" denotes the pooled Sst subclass-level analysis of Fig. 2.
> (**c**) Per-gene meta-analytic z (SCZ versus control) in the depleted versus the
> non-depleted group, for the shared genes tested in both. Points are coloured by
> membership of the leading-edge union of the gene sets in d; dashed line, identity;
> ρ, Spearman. Labelled genes are shown in bold in e. (**d, e**) Normalized
> enrichment scores for representative gene sets (**d**) and meta-analytic z for
> exemplar leading-edge genes (**e**), grouped into blocks shared across the subclass
> (synaptic, ubiquitin) and blocks graded by depletion (cytosolic translation,
> oxidative phosphorylation). Significance throughout: \*\*\*FDR < 0.01,
> \*\*FDR < 0.05, \*FDR < 0.10.

## Wording note on the closing sentence

"Largely confined to the supertypes undergoing depletion" is chosen over "specific to the
depleted supertypes" because translation suppression is as strong in the intermediate group
(negative, non-significant depletion) as in the depleted group; only oxidative phosphorylation
is depleted-only. "Undergoing depletion" covers both groups with negative estimates and
excludes the non-depleted group (estimate ≥ 0), which shows neither program. "Largely"
acknowledges that module scores show a weaker version of both programs in the non-depleted
group, even though no gene set reaches significance there.
