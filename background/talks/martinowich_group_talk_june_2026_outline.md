# Group talk — "Resolving cellular pathology in schizophrenia through single-cell meta-analysis"

Lightweight, greppable text outline of `martinowich_group_talk_june_2026.pptx`
(Shreejoy Tripathy, CAN meeting, May 19 2026; 60 slides). The full 67 MB deck sits
alongside this file; this outline is the durable, committable reference (the deck
is the fullest narrative of the paper's argument and figure emphasis).

## Main talk (slides 1–34)
1. Title — Resolving cellular pathology in schizophrenia through single-cell meta-analysis
3. SCZ affects 0.5–1% adults; lifelong treatment; ~15 yr reduced life expectancy; not gross neurodegeneration
4. Key questions: which brain cells most changed? how — dying/less abundant vs changing molecular state? does pathology track symptoms?
5. **The debate:** GABAergic interneuron deficits are a hallmark; some studies find reduced *numbers* (esp. SST, PVALB), others find altered *gene expression* without cell loss
6. SEA-AD taxonomy (Gabitto 2024): 3 classes / 24 subclasses / 127 supertypes
7–8. Patch-seq integrates atlases with multimodal features (Lee & Dalley 2023; 16 Sst supertypes)
9. Study goals: can fine-grained taxonomies reconcile competing hypotheses? shared vulnerability across disorders? mitigate heterogeneity via meta-analysis?
10. Meta-analysis design across post-mortem cohorts + modalities (Nicole Endresz)
11–15. Xenium cell typing (two-pass correlation label transfer; depth labels; panel resolves subclasses but mixed/limited at supertype)
17. **snRNAseq meta-analysis → cell-type-wide SCZ DE**
18. **★ FIG 2 headline: "Somatostatin cells show cell type-INTRINSIC mRNA dysregulation in snRNAseq and Xenium"** — replicates Dienel et al. 2023 & Lewis-lab RNAscope work
19. **★ FIG 2 headline: "Xenium shows EXCELLENT replication of DE effects from snRNAseq meta-analysis"** — genes DE in meta (FDR<0.10) show concordant direction in Xenium
21. AD comparison (Gabitto 2024): upper-layer (L2/3) cells most impacted in AD; SST subtypes lost earliest
22. **FIG 3:** meta-analysis of 7 datasets, crumblr (Hoffman & Roussos 2025) → vulnerable neocortical subtypes
23. **FIG 3 headline:** SST subtypes LESS abundant, L6b subtypes MORE abundant; vulnerable Sst = Sst_25, Sst_2, Sst_22, Sst_3, Sst_20 (all also lost in AD); L6b_1/L6b_4 increased (NOT shared with AD)
24. **Sst_25: ~26% reduction in SCZ (95% CI 14–36%)** across 7 cohorts; abundance effects underpowered per-dataset → illustrates meta-analysis power
25. Spatial transcriptomics replicates decreased vulnerable-SST abundance in superficial layers
26. Xenium corroborates altered abundances
27. **Synthesis:** global downregulation of Sst mRNA in Sst cells + Sst-subtype-specific reduced abundance = bona fide loss OR drastic erosion of cell-type identity
28–31. Vulnerable SST express CALB1 → calbindin+ double-bouquet cells (Beasley 2002 calbindin reduction; primate-specialized, Yanez 2005; vertical column inhibition, Raghanti 2010)
32. Summary: fine-grained taxonomies = meaningful divisions; SCZ alters subtype numbers
33. **Limitation:** cell-intrinsic mRNA dysregulation may confound cell typing at fine subtype level — can't fully exclude that cells are "there" but transcriptionally unrecognizable
34. Acknowledgements (Nicole Endresz, Etienne Sibille, TripLab/KCNI)

## Bonus slides (35–60)
37. RNAscope corroborates reduced SST density in sgACC (Arbabi & Newton 2025; Sibille lab, Pitt Brain Bank)
38/50. Increased deep-layer/white-matter neuron density in SCZ reported prior (Kubo 2020)
39/53. Dataset citations (Batiuk, Lee/PsychAD, Ruzicka, Emani, Fröhlich, Kwon, Arbabi)
40–41,46. **FIG 4:** vulnerable Sst enriched for SCZ GWAS risk (Duncan et al. 2025); GWAS implicates genes specifically expressed in upper-layer SST
45. Patch-seq (Cadwell 2016; Lipovsek 2021)
48–49,51. L6b: proportions, human morphologies (Dalley 2026), cortico-thalamic role (Zolnik 2026)
52. Inclusion criteria: 10x v3/v3.1; age 19–70; frontal/PFC; ~2.3M cells
54. **Shared SCZ↔early-AD vulnerability:** 5 of 8 SST supertypes vulnerable in early AD also affected in SCZ
58. Marker genes of Sst supertypes

## Key emphasis takeaways for the manuscript
- **Fig 2 (DE):** headline is *cell-type-INTRINSIC* SST mRNA dysregulation (slide 18; replicates Dienel/Lewis), **co-headlined with** Xenium cross-platform replication of DE (slide 19). Talk is **SST-primary** — PVALB is not co-headlined in the talk. DE-burden-tracks-abundance is not a talk headline.
- **Fig 3 (composition):** SST subtypes less abundant + L6b more abundant; Sst_25 ~26% reduction (95% CI 14–36%) — fills the manuscript's `[TODO: Sst_25 effect size]`.
