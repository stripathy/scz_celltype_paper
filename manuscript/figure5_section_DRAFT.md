# Figure 5 results section — DRAFT v2 (2026-09-01)

**Status.** Full rewrite around the 7-cohort **pseudobulk** meta-DE (Nicole's
framework replicated; the earlier IVW draft is in git history). **Every number
below is asserted programmatically** by `transcriptomic/scripts/fig5/09_verify.R`
against the pipeline outputs — run it before any renumber and before pasting into
the Doc; it exits non-zero on any mismatch.

Headline burden numbers are the **FDR < 0.10** tier, matching the labels drawn in
panel b (changed from FDR < 0.05 on 2026-09-01).

Figure 5 is built by `transcriptomic/scripts/fig5/08_figure5.R` →
`figure5_sst_strata.(png|pdf)`: **a** strata definition (depletion × depth, layer
gutter), **b** gene-set burden (two FDR tiers), **c** gene-level z, depleted vs
non-depleted, colored by module, **d** NES heatmap (4 groups × 8 sets),
**e** exemplar-gene z heatmap. See `transcriptomic/scripts/fig5/README.md`.

Remaining gates before Doc paste: co-author sign-off on FDR < 0.20 strata
(Sst_20 in), Fig. 5 as last main figure, supplementary figure numbers (S# below),
and the [VERIFY] flags at the bottom.

**VERSION 0 (added 2026-09-01) is the working base**: combined from the three
versions below, citation-free per Shreejoy's call (literature context lives in
the Discussion paragraph). Versions 1–3 kept for salvage while editing.

---

## Section header (declarative, house style)

Option 1 (recommended): **Suppression of translation and oxidative-phosphorylation
programs in Sst interneurons is graded by depletion status**

Option 2: **Depleted Sst supertypes carry a graded suppression of protein-synthesis
and energy-metabolism programs superimposed on a shared synaptic deficit**

---

## CURRENT WORKING TEXT (2026-09-01, post-Hallmark-removal) — Shreejoy's edit + verified numbers

Gene sets restricted to GO (BP/CC/MF) + Reactome; Hallmark excluded, KEGG not
added. Every number below passes `fig5/09_verify.R`.

Finally, we asked whether the transcriptional disease state of Sst interneurons
differs between the supertypes being depleted and those that persist. While
prior cell-resolved studies of Sst neurons have reported suppressed protein
translation under chronic stress (Lin and Sibille 2015), to our knowledge, no
prior study has asked how such dysregulation might differ between Sst subtypes.
To answer this, we grouped the 16 Sst supertypes into 3 broad groups based on
their depletion status in Fig. 3a: depleted (FDR < 0.20; n = 5 supertypes),
intermediate (negative, non-significant; n = 6) and non-depleted (n = 5),
broadly tracking upper-, middle-, and deep-layers in cortex. We pooled nuclei
annotated to each depletion group per donor and repeated the
differential-expression analysis of Fig. 2 followed by gene set enrichment
analysis (GSEA) (Fig. 5a; Methods).

Sst interneuron gene-set dysregulation was strongly graded by depletion status
(Fig. 5b). We identified 130 gene sets altered in the depleted group at
FDR < 0.10 (88% downregulated), 64 in the intermediate group, and none in the
non-depleted group. Notably, pooling all Sst cells together yielded only 55
significant gene sets, fewer than the depleted group alone, despite the pooled
analysis drawing on more than twice as many nuclei per donor and therefore
having greater power to detect dysregulation. Pooling across supertypes
therefore dilutes signal concentrated in the depleted supertypes.

Two broad transcriptional patterns organize the dysregulation landscape of Sst
interneurons (Fig. 5c–e). The first pattern was shared across Sst interneurons
regardless of depletion status. Synaptic and neuropeptide programs were
downregulated to a similar degree in all three groups (regulation of
trans-synaptic signaling, normalized enrichment score, NES −1.4 to −1.5 in every
group, exemplar genes NMU, APBA2 and PICK1), as were the individual genes SST
and VGF. Similarly, protein-deubiquitination machinery was upregulated in all
three groups, reaching significance in the depleted group (K63-linked
deubiquitination, NES = +2.15, FDR = 0.007, exemplar genes USP8, STAMBPL1 and
VCP).

The second transcriptional component was largely restricted to the depleted
group and, to a lesser extent, the intermediate group. Gene sets for cytosolic
translation and ribosome components were suppressed in the depleted and
intermediate groups (NES −1.9 to −2.6, all FDR < 10⁻³, exemplar genes RPL36,
RPL10 and EIF3G) but showed no change in the non-depleted group (NES ≥ −0.7,
FDR = 1.0). Gene sets for oxidative phosphorylation, the mitochondrial pathway
that supplies neuronal energy, were suppressed only in the depleted group
(NES −1.9 to −2.0, both FDR ≤ 0.0033, exemplar genes NDUFS8, COX5B and CYCS).
Formal models testing the interaction between depletion group status and SCZ
diagnosis supported these findings for protein synthesis (SRP-dependent
cotranslational protein targeting, FDR = 7.9 × 10⁻¹⁰). Oxidative phosphorylation
showed the same direction at trend level in this gene-set analysis (electron
transport chain, FDR = 0.053) and reached significance when the module was
tested as a single measurement per donor (P = 0.032). Our analysis therefore separates two components of the Sst
disease state, one shared across the subclass and including the reduction of SST
itself, which corresponds to what earlier work has described, and a second
confined to the depleted supertypes that becomes visible only when Sst neurons
are resolved by their layer-graded depletion.

---

## VERSION 0 — COMBINED (~640 words; citation-free; THE BASE FOR MANUAL EDITING)

Rationale framed from within the paper (Figs 2–4), no literature citations —
the Discussion paragraph carries the Sibille-corpus context. All numbers
identical to the verified set below.

**Cell-count confound: tested and excluded (2026-09-01).** The groups are defined
by cases having fewer cells, so SCZ donors contribute ~30% fewer nuclei to the
depleted-group pseudobulk. Subsampling control donors' nuclei down to the case
distribution and rebuilding the pseudobulks from cell-level data removes that
imbalance entirely (-30% becomes -0.4%, P = 0.96) and leaves the result intact
across 5 draws: translation significant in 20/20 tests in the depleted group and
20/20 in intermediate, 0/20 in non-depleted; oxidative phosphorylation 10/10 in
depleted, 0/10 elsewhere. A sentence to this effect belongs in the results or
methods. Detail in `FIGURE5_FINALIZATION_PLAN.md`; do NOT cite the
log(nuclei)-covariate version, which over-adjusts for a mediator.

**Opening sentence — three hypothesis framings (2026-09-01).** Shreejoy's edit
adds an explicit hypothesis; the choice of *which* hypothesis changes what the
results have to deliver.

- **(a) Quantity** — "…hypothesizing that depleted Sst cells may have a greater
  degree of coordinated transcriptional dysregulation." Simplest, and directly
  motivates panel b. Cost: "greater degree" invites the power objection in the
  very next breath (the depleted stratum could simply be better powered), so the
  nuclei-count sentence has to do defensive work immediately.
- **(b) Kind, ADOPTED IN THE TEXT ABOVE** — "…hypothesizing that the depleted
  supertypes carry dysregulation beyond the changes shared across the subclass."
  Still directional and still predicts panel b, but "beyond … shared" pre-loads
  the two-component decomposition that the results actually deliver, so the
  payoff lands as confirmation rather than as a surprise. Also less exposed to
  the power objection, because it is a claim about an additional component
  rather than about a bigger count.
- **(c) Vulnerability-linked** — "…reasoning that if these supertypes are
  intrinsically vulnerable (Fig. 4), their disease state should differ from that
  of the Sst neurons that persist." Ties Fig. 5 back to the genetics figure most
  explicitly. Use if the section should read as continuous with Fig. 4 rather
  than as a new question.

Small fixes folded in: "between … compared to" → "between … and" (mixed
construction); doubled "their their"; "followed up by" → "followed by"; "per
strata" → "per stratum"; serial "and" before postmortem interval. The Spearman
ρ was kept as a parenthetical, since without it "corresponding approximately to
upper-, middle- and deep-layer" is an assertion the reader cannot check. The
per-donor pseudobulk and nuclei counts moved out to Methods.

Finally, we asked whether the transcriptional disease state of Sst interneurons
differs between the supertypes that are depleted and those that persist,
hypothesizing that the depleted supertypes carry dysregulation beyond the
changes shared across the subclass. To answer this, we first grouped the 16 Sst
supertypes into three broad strata based on their abundance change in SCZ
(Fig. 3a): depleted (FDR < 0.20: Sst_2, Sst_3, Sst_20, Sst_22, Sst_25),
intermediate (reduced but non-significant compositional change; n = 6
supertypes) and non-depleted (estimate ≥ 0; n = 5 supertypes) (Fig. 5a),
corresponding approximately to upper-, middle- and deep-layer cortical Sst
populations (Spearman ρ = 0.74 between abundance change and median cortical
depth, *P* = 0.0016). Within each stratum we performed differential expression
as in Fig. 2, followed by gene-set enrichment (GSEA) to identify dysregulated
pathways per stratum (Methods). Because the three strata are measured within the
same donors, donor-level factors, including medication exposure, age, sex and
postmortem interval, are shared across strata by construction.

Gene-set dysregulation was strongly graded by depletion status (Fig. 5b). At
FDR < 0.10, 145 gene sets were altered in the depleted stratum (88%
downregulated), 61 in the intermediate stratum, and none in the non-depleted
stratum; the pooled Sst subclass analysis recovered 56, indicating that pooling
across subtypes dilutes signal concentrated in the affected types. Individually
significant genes were few and did not track this gradient (19, 99 and 35 genes
at FDR < 0.10, versus 342 for the pooled subclass), indicating that the stratum
differences arise from coordinated shifts across many genes rather than from
large single-gene effects.
The non-depleted stratum contains fewer nuclei (13,157, versus 44,744 and
42,814), reducing its power; power cannot, however, explain the gradient between
the two well-sampled strata, where the intermediate stratum has slightly more
donors and comparable nuclei counts yet carries roughly a third of the depleted
stratum's burden.

The burden resolved into distinct components (Fig. 5c–e). A shared synaptic and
neuropeptide program was downregulated with comparable effect size in all three
strata: regulation of trans-synaptic signaling (NES −1.4 to −1.5 in every
stratum), *SST* itself (meta-analytic z = −3.4, −3.6 and −2.1 in the depleted,
intermediate and non-depleted strata), *VGF* (z = −3.1, −2.9, −2.9) and *NMU*
(z = −2.8, −2.7, −3.1; FDR = 0.018 in the pooled subclass); treating each
donor's mean expression of the synaptic module as a single measure confirmed
downregulation in every stratum, including the non-depleted one (−0.23, −0.18
and −0.10 SD; *P* = 7.0 × 10⁻⁸, 1.6 × 10⁻⁶ and 3.9 × 10⁻³). The few upregulated
gene sets (18 of 145 in the depleted stratum) prominently included
protein-deubiquitination machinery (protein K63-linked deubiquitination,
NES = +2.15, FDR = 0.006; NES positive in all strata; leading edge *USP8*,
*STAMBPL1*, *ATXN3*, *VCP*). By contrast, cytosolic translation and ribosomal
gene sets were suppressed in the depleted (NES = −2.6 and −2.3,
FDR = 2.3 × 10⁻¹² and 7.3 × 10⁻⁸) and intermediate strata (NES = −2.2 and −1.9)
but showed no signal in the non-depleted stratum (NES ≥ −0.7, FDR = 1.0), and
oxidative-phosphorylation gene sets were suppressed only in the depleted stratum
(NES = −1.9 to −2.5, all FDR ≤ 0.0023; intermediate smallest FDR = 0.20).
Exemplar leading-edge genes show the same graded pattern (Fig. 5e): *RPL36*
(z = −2.5, −1.7, +0.5 across strata), *EIF3G* (−1.8, −0.4, +1.1), *NDUFS8*
(−2.9, −0.1, +0.5) and *COX5B* (−2.4, −0.7, +0.9).

A formal diagnosis × stratum interaction analysis — each donor contributing one
pseudobulk per stratum, modeled within donor (Methods) — confirmed the
dissociation. Interaction effects were concentrated in translation
(SRP-dependent cotranslational protein targeting, NES = −2.7,
FDR = 4.4 × 10⁻¹⁰; eukaryotic translation elongation, FDR = 2.3 × 10⁻⁹) and
oxidative phosphorylation (NES = −2.0, FDR = 1.6 × 10⁻⁴, with a monotonic trend
across strata, FDR = 7.5 × 10⁻⁴), whereas synaptic gene sets (FDR ≥ 0.21), *SST*
and *VGF* (interaction z = −0.8 and −0.6) showed none: their reduction is
stratum-independent. Few individual genes carried significant interactions (7 at
FDR < 0.05, e.g., *CELF2*, *GRIK4*, *QKI* — crossover patterns, not members of
the graded modules), indicating stratum dependence carried by coordinated small
shifts across many genes. Where the 300-gene Xenium panel permitted a
cross-platform check it agreed — *SST* (Xenium stratum-level z = −3.2, −2.2,
−2.5) and *VGF* (−4.6, −3.3, −2.3) were reduced in all three strata on both
platforms (Fig. S#) — but the translation and oxidative-phosphorylation modules
are essentially absent from the panel (2 of 130 and 2 of 114 leading-edge genes)
and remain snRNA-seq observations. Unlike the germline genetic convergence of
Fig. 4, these are cross-sectional state differences and carry no causal
direction; the shared component is by construction immune to compositional
artifacts, whereas depletion-graded signal could partly reflect altered subtype
composition among the surviving cells of the depleted stratum (Discussion).

---

## VERSION 1 — LONG (~1,000 words, full references; the everything-on-the-table draft)

Having identified which Sst supertypes are depleted (Fig. 3) and shown that
genetic risk converges on them (Fig. 4), we finally asked how the transcriptional
disease state differs between the Sst subtypes being depleted and those that
persist. Prior cell-resolved work supplies specific expectations. In
laser-captured cortical SST neurons from chronically stressed mice, protein
translation through EIF2 signaling is the most suppressed pathway (Lin and
Sibille 2015), a suppression attributed to endoplasmic-reticulum stress arising
from the processing load of the SST precursor peptide itself (Tomoda et al.
2022); ribosomal-protein genes are likewise downregulated in bulk prefrontal
cortex in human depression and mouse chronic stress (Zhang et al. 2023). A
mitochondrial and oxidative liability of SST neurons was proposed over a decade
ago (Lin and Sibille 2013) but has not been observed at cell-type resolution: in
stressed mouse cortex the oxidative-phosphorylation decrement localizes to
pyramidal neurons (Newton et al. 2022), and in human subgenual cingulate cortex
across psychiatric disorders, SST neurons in SCZ show reduced ATPase-activity
gene expression alongside elevated ER-stress genes (Arbabi et al. 2025).
Critically, none of these studies resolved SST-neuron subtypes, and laminar in
situ work concluded that the SST deficit reflects a general vulnerability of SST
neurons independent of cell subtype (Seney et al. 2015) — a conclusion our
compositional data now allow us to test directly.

We grouped the 16 Sst supertypes into three strata by their compositional change
in SCZ (Fig. 3a): depleted (crumblr FDR < 0.20: Sst_2, Sst_3, Sst_20, Sst_22,
Sst_25), intermediate (negative but non-significant; n = 6), and non-depleted
(estimate ≥ 0; n = 5) (Fig. 5a). Because depletion tracks laminar position
(Spearman ρ = 0.74 between crumblr estimate and median cortical depth,
*P* = 0.0016; stratum mean depths 0.31, 0.49 and 0.72), the strata correspond
approximately to upper-layer, mid-depth and deep Sst populations. For each
stratum we pooled counts across its supertypes into one pseudobulk per donor per
dataset (donors with ≥ 10 nuclei of the stratum: n = 395, 411 and 342 across the
seven datasets; median 85, 90 and 31 nuclei per donor), repeated the
differential-expression meta-analysis using the identical framework as Fig. 2,
and performed preranked gene-set enrichment (GSEA) on the meta-analytic
statistics (Methods). Because the three strata are measured within the same
donors, donor-level factors — medication exposure, age, sex, postmortem
interval — are shared across strata by construction.

Gene-set dysregulation was strongly graded by depletion status (Fig. 5b). At
FDR < 0.10, 145 gene sets were altered in the depleted stratum (88%
downregulated), 61 in the intermediate stratum, and none in the non-depleted
stratum; the pooled Sst subclass analysis recovered 56, indicating that pooling
across subtypes dilutes signal concentrated in the affected types. Individually
significant genes were few and did not track this gradient (19, 99 and 35 genes
at FDR < 0.10 in the depleted, intermediate and non-depleted strata, versus 342
for the pooled subclass), indicating that the stratum differences arise from
coordinated shifts across many genes rather than from large single-gene effects. The non-depleted stratum contains
fewer nuclei (13,157, versus 44,744 and 42,814), which reduces its power; power
differences cannot, however, explain the gradient between the two well-sampled
strata, where the intermediate stratum has slightly more donors and comparable
nuclei counts yet carries roughly a third of the depleted stratum's gene-set
burden.

The burden resolved into distinct components (Fig. 5c–e). First, a shared
synaptic and neuropeptide program was downregulated with comparable effect size
in all three strata: regulation of trans-synaptic signaling (NES −1.4 to −1.5 in
every stratum), *SST* itself (meta-analytic z = −3.4, −3.6 and −2.1 in the
depleted, intermediate and non-depleted strata), *VGF* (z = −3.1, −2.9, −2.9)
and *NMU* (z = −2.8, −2.7, −3.1; FDR = 0.018 in the pooled subclass). Treating
each donor's mean expression of the 186-gene synaptic module as a single measure
confirmed downregulation in every stratum, including the non-depleted one
(−0.23, −0.18 and −0.10 SD; *P* = 7.0 × 10⁻⁸, 1.6 × 10⁻⁶ and 3.9 × 10⁻³;
Fig. S#). Second, the few upregulated gene sets (18 of the 145 in the depleted
stratum) prominently included protein-deubiquitination machinery: protein
K63-linked deubiquitination (NES = +2.15, FDR = 0.006 in the depleted stratum;
NES positive in all strata), with leading-edge deubiquitinases *USP8*,
*STAMBPL1* and *ATXN3* and the proteostasis ATPase *VCP* (Fig. 5d,e).

Third — and in sharpest contrast to the shared program — cytosolic translation
and oxidative-phosphorylation programs were suppressed in a depletion-graded
manner (Fig. 5c–e). Ribosomal-subunit and cytoplasmic-translation gene sets were
strongly downregulated in the depleted (NES = −2.6 and −2.3, FDR = 2.3 × 10⁻¹²
and 7.3 × 10⁻⁸) and intermediate strata (NES = −2.2 and −1.9, FDR = 6.3 × 10⁻⁷
and 5.1 × 10⁻⁴) but showed no signal in the non-depleted stratum (NES = −0.7 and
−0.6, FDR = 1.0). Oxidative-phosphorylation and respiratory-electron-transport
sets were suppressed only in the depleted stratum (NES = −1.9 to −2.5, all
FDR ≤ 0.0023; not significant in the intermediate stratum, smallest FDR = 0.20).
Exemplar leading-edge genes show the same pattern (Fig. 5e): *RPL36* (z = −2.5,
−1.7, +0.5 across the three strata), *EIF3G* (−1.8, −0.4, +1.1), *NDUFS8*
(−2.9, −0.1, +0.5) and *COX5B* (−2.4, −0.7, +0.9).

A formal diagnosis × stratum interaction analysis — each donor contributing one
pseudobulk per stratum, modeled within donor (Methods) — confirmed the
dissociation. Gene sets enriched among depleted-versus-non-depleted interaction
effects were overwhelmingly translational (SRP-dependent cotranslational protein
targeting, NES = −2.7, FDR = 4.4 × 10⁻¹⁰; eukaryotic translation elongation,
NES = −2.7, FDR = 2.3 × 10⁻⁹; cytoplasmic translation, FDR = 1.6 × 10⁻⁵), with
oxidative phosphorylation close behind (Hallmark set, NES = −2.0,
FDR = 1.6 × 10⁻⁴) and additionally showing a monotonic trend across the three
strata (linear-trend GSEA FDR = 7.5 × 10⁻⁴). Synaptic gene sets showed no
interaction (FDR ≥ 0.21), and neither did *SST* or *VGF* individually
(interaction z = −0.8 and −0.6): their reduction is stratum-independent. A
complementary module-score test detected a modest within-donor
depleted-versus-non-depleted difference for all three modules (translation
*P* = 0.028, oxidative phosphorylation *P* = 0.023, synaptic *P* = 1.3 × 10⁻⁴),
consistent with a graded global component to the disease state; only the
translation and oxidative-phosphorylation modules, however, exceeded that
transcriptome-wide background in the interaction GSEA. Few individual genes
carried formally significant interactions (7 at FDR < 0.05, e.g., *CELF2*,
*GRIK4*, *QKI* — crossover patterns rather than members of the graded modules;
the nearest module gene was *COX5B*, interaction z = −3.8, FDR = 0.13),
indicating that the stratum dependence is carried by coordinated small shifts
across many genes.

Where the 300-gene Xenium panel permitted a cross-platform check, it agreed:
*SST* (Xenium stratum-level z = −3.2, −2.2 and −2.5) and *VGF* (−4.6, −3.3,
−2.3) were reduced in all three strata on both platforms (Fig. S#). The
translation and oxidative-phosphorylation modules, however, are essentially
absent from the panel (2 of 130 and 2 of 114 leading-edge genes) and remain
snRNA-seq observations. We emphasize that, unlike the germline genetic
convergence of Fig. 4, these are cross-sectional state differences and carry no
causal direction; and whereas the shared component is by construction immune to
compositional artifacts, depletion-graded signal could partly reflect altered
subtype composition among the surviving cells of the depleted stratum
(Discussion).

---

## VERSION 2 — PARED (~430 words; rationale-led, light references; paper-adjacent)

Finally, we asked whether the transcriptional disease state of Sst interneurons
differs between the supertypes being depleted and those that persist.
Cell-resolved studies of SST neurons have reported suppressed EIF2-mediated
protein translation under chronic stress (Lin and Sibille 2015) and reduced
ATPase-activity gene expression in SCZ (Arbabi et al. 2025), but none resolved
SST subtypes, and laminar in situ work concluded that the SST deficit is
subtype-independent (Seney et al. 2015). We therefore grouped the 16 Sst
supertypes into depleted (FDR < 0.20 in Fig. 3a; n = 5), intermediate (negative,
non-significant; n = 6) and non-depleted (n = 5) strata — an axis that also
tracks laminar depth (Spearman ρ = 0.74, *P* = 0.0016) — pooled each stratum
into one pseudobulk per donor per dataset, and repeated the
differential-expression meta-analysis and GSEA of Fig. 2 (Fig. 5a; Methods).
Because the strata are measured within the same donors, donor-level factors such
as medication exposure cannot generate stratum differences.

Gene-set dysregulation was strongly graded by depletion status (Fig. 5b): 145
gene sets were altered in the depleted stratum at FDR < 0.10 (88%
downregulated), 61 in the intermediate stratum, and none in the non-depleted
stratum, whereas the pooled Sst subclass analysis recovered 56 — pooling across
subtypes dilutes signal concentrated in the affected types. This burden resolved
into two components (Fig. 5c–e). A shared synaptic and neuropeptide program —
including *SST* itself (meta-analytic z = −3.4, −3.6 and −2.1 across the
depleted, intermediate and non-depleted strata) and *VGF* (−3.1, −2.9, −2.9) —
was downregulated with comparable effect size in all three strata, accompanied
by upregulated protein-deubiquitination machinery (K63-linked deubiquitination:
NES = +2.15, FDR = 0.006; leading edge *USP8*, *STAMBPL1*, *VCP*). By contrast,
cytosolic translation and ribosomal gene sets were suppressed in the depleted
and intermediate strata (NES −1.9 to −2.6, all FDR < 10⁻³) but not in the
non-depleted stratum (NES ≥ −0.7, FDR = 1.0), and oxidative-phosphorylation gene
sets were suppressed only in the depleted stratum (NES −1.9 to −2.5,
FDR ≤ 0.0023) — a graded pattern echoed by exemplar genes such as *RPL36*
(z = −2.5, −1.7, +0.5) and *NDUFS8* (−2.9, −0.1, +0.5). A formal
diagnosis × stratum interaction analysis, modeled within donors (Methods),
confirmed this dissociation: interaction effects concentrated in translation
(SRP-dependent cotranslational targeting, FDR = 4.4 × 10⁻¹⁰) and oxidative
phosphorylation (FDR = 1.6 × 10⁻⁴), whereas synaptic gene sets, *SST* and *VGF*
showed none. Where the Xenium panel allowed a cross-platform check it agreed —
*SST* and *VGF* were reduced in all three strata on both platforms (Fig. S#) —
but the translation and oxidative-phosphorylation modules are essentially absent
from the 300-gene panel and remain snRNA-seq observations. Unlike the germline
genetic convergence above, these are cross-sectional state differences and carry
no causal direction (Discussion).

---

## VERSION 3 — SHORT (~290 words; tightest, matches Doc section density; paste candidate)

Finally, we asked whether the transcriptional disease state of Sst interneurons
differs between the supertypes being depleted and those that persist — prior
studies of SST-neuron pathology could not resolve subtypes and had concluded the
deficit was subtype-independent (Seney et al. 2015; Lin and Sibille 2015). We
grouped the 16 Sst supertypes into depleted (FDR < 0.20 in Fig. 3a; n = 5),
intermediate (n = 6) and non-depleted (n = 5) strata, an axis that also tracks
laminar depth (Spearman ρ = 0.74; Fig. 5a), pooled each stratum into one
pseudobulk per donor per dataset, and repeated the differential-expression
meta-analysis and GSEA of Fig. 2 (Methods). Because strata are compared within
the same donors, donor-level factors such as medication cannot generate stratum
differences.

Gene-set dysregulation was strongly graded (Fig. 5b): 145 gene sets at
FDR < 0.10 in the depleted stratum (88% downregulated), 61 in the intermediate,
none in the non-depleted, and 56 in the pooled subclass analysis. The burden
resolved into two components (Fig. 5c–e). A shared synaptic and neuropeptide
program — including *SST* (z = −3.4, −3.6, −2.1 across strata) and *VGF* (−3.1,
−2.9, −2.9) — was downregulated comparably everywhere. By contrast, cytosolic
translation and ribosomal gene sets were suppressed in the depleted and
intermediate strata only (NES −1.9 to −2.6, FDR < 10⁻³; non-depleted FDR = 1.0),
and oxidative-phosphorylation sets only in the depleted stratum (NES −1.9 to
−2.5, FDR ≤ 0.0023). A within-donor diagnosis × stratum interaction analysis
confirmed the dissociation: interaction effects concentrated in translation
(FDR = 4.4 × 10⁻¹⁰) and oxidative phosphorylation (FDR = 1.6 × 10⁻⁴), whereas
synaptic gene sets, *SST* and *VGF* showed none. *SST* and *VGF* were reduced in
all three strata on the Xenium platform as well (Fig. S#); the translation and
oxidative-phosphorylation modules are essentially absent from the 300-gene panel
and remain snRNA-seq observations. Unlike the genetic convergence above, these
are cross-sectional state differences and carry no causal direction
(Discussion).

---

## DISCUSSION PARAGRAPH (~250 words; no figure refs, per Doc convention)

**Placement: between the AD/intrinsic-vulnerability paragraph ("Intriguingly,
the Sst supertypes we find depleted in SCZ are also…") and the L6b paragraph
("Beyond our focus on GABAergic interneurons…").** It extends the vulnerability
question (what does vulnerability look like transcriptionally?) and feeds the
two-part therapeutic close, which already reasons from it ("restoring their
inhibitory output will not be enough"). **Required companion edit: the Fig. 4
results section currently opens "Finally, we asked…" — change to "Next, we
asked…" once Fig. 5 becomes the last figure.** Optional light touch, not
required: ¶1 of the Discussion could gain a clause noting the state itself
decomposes into shared and depletion-graded parts.

> Stratifying the transcriptional disease state by depletion revealed two
> superimposed programs, a distinction that speaks to both literatures
> reconciled above. The component the in situ era measured — reduced *SST*,
> *VGF* and synaptic gene expression — is present in depleted and spared
> supertypes alike, consistent with earlier laminar analyses that concluded the
> SST deficit is subtype-independent (Seney et al. 2015), and is most
> parsimoniously read as a circuit-level state imposed on all Sst neurons, for
> example through reduced trophic support: SST-neuron gene programs depend on
> dendritically targeted BDNF (Tripp et al. 2012; Oh et al. 2019), and BDNF is
> itself reduced in pyramidal populations in our meta-analysis [VERIFY breadth
> against Table T2; the L2/3 IT reduction is already stated in the Fig. 2
> section]. The depletion-graded component — coordinated suppression of
> cytosolic translation and oxidative phosphorylation — instead marks the
> vulnerable subtypes specifically, and connects to a mechanistic literature on
> SST-neuron fragility: suppressed EIF2-mediated translation in stressed SST
> neurons (Lin and Sibille 2015), attributed to endoplasmic-reticulum stress
> from preproSST processing load (Tomoda et al. 2022), and ribosomal-protein
> downregulation in bulk cortex in depression (Zhang et al. 2023). Our findings
> extend this translational suppression to schizophrenia, localize it to the
> depleted subtypes, and provide the first subtype-resolved support for the
> mitochondrial liability of SST neurons proposed over a decade ago (Lin and
> Sibille 2013) — a deficit cell-resolved studies had instead located in
> pyramidal neurons (Newton et al. 2022), with reduced ATPase-activity genes in
> SST cells the nearest antecedent (Arbabi et al. 2025). These are state
> differences, not causes: suppressed protein synthesis and energy metabolism in
> the depleted types could precede and promote their depletion, could reflect
> the response of cells to an ongoing disease process, or could partly reflect
> the composition of surviving cells. Whichever holds, protein synthesis and
> oxidative ATP production are among the costliest activities of a neuron, and
> their coordinated withdrawal in the very cells declining in number is the
> signature expected of a population under sustained strain.

---

## Figure 5 legend — DRAFT (updated to v15 panels; descriptive noun-phrase titles)

**Fig. 5 | Transcriptional disease state of Sst neurons across depletion
strata.** (**a**) Stratum definitions: per-supertype SCZ abundance change
(snRNA-seq crumblr estimate ± SE, Fig. 3a) versus median cortical depth from
Xenium (0 = pia). Depleted, crumblr FDR < 0.20 (n = 5 supertypes); intermediate,
negative but non-significant (n = 6); non-depleted, estimate ≥ 0 (n = 5). Right
gutter, manually drawn layer boundaries from the Xenium sections. (**b**) Gene
sets significant per group in preranked GSEA on the stratum-level meta-analytic
differential expression, by direction and FDR tier; "All Sst" is the pooled
subclass analysis (Fig. 2). (**c**) Per-gene meta-analytic z (SCZ vs control) in
the depleted versus non-depleted stratum; genes colored by module membership
(leading-edge union of the Fig. 5d gene sets); labeled genes are highlighted in
e; dashed line, identity; r, Pearson. (**d**) Normalized enrichment scores for
representative gene sets, grouped into shared (synaptic; ubiquitin) and
depletion-graded (cytosolic translation; oxidative phosphorylation) blocks;
\*\* FDR < 0.01, \* FDR < 0.05, + FDR < 0.10. (**e**) Stratum-level
meta-analytic z for exemplar leading-edge genes of each block (bold, labeled in
c); *SST* shown for reference. Stratum-level differential expression: one
pseudobulk per donor per stratum per dataset (donors with ≥ 10 nuclei),
meta-analyzed across the seven datasets exactly as in Fig. 2 (Methods).

---

## Methods stub (to expand, ~150 of the eventual ~400 words)

- Strata: crumblr 7-dataset meta (Fig. 3a); FDR < 0.20 & estimate < 0 =
  depleted; estimate < 0 otherwise = intermediate; estimate ≥ 0 = non-depleted.
- Stratum DE: per donor per stratum pseudobulk (sum of counts across the
  stratum's supertypes); donors ≥ 10 nuclei; genes with ≥ 1 count in ≥ 80% of
  retained donors; TMM → voom → limma (~ diagnosis + age + sex + PMI; PMI
  omitted for Multiome); metafor REML random-effects across datasets, k ≥ 5;
  BH-FDR within stratum. Identical framework to the subclass DE (Fig. 2).
- GSEA: fgsea preranked on z = estimate/SE; MSigDB Hallmark + GO BP/CC/MF +
  Reactome; minSize 10, maxSize 500; BH within stratum.
- Interaction: donor-blocked (duplicateCorrelation) limma,
  ~ diagnosis × stratum + covariates, reference = non-depleted; also a linear
  stratum-score term; GSEA on interaction z.
- Module scores: per donor-stratum mean of within-dataset z-scored log-CPM over
  each module's leading-edge union (translation 130 genes, OxPhos 114, synaptic
  186); score ~ diagnosis + covariates per dataset → REML meta; paired test on
  within-donor (depleted − non-depleted) scores.
- Xenium check: stratum-level pseudobulk DE on the 24 sections, panel genes.

## Citations (PMIDs verified in background/sibille_sst_literature.md)

- Lin & Sibille 2013, Front Pharmacol — PMID 24058344 (mitochondrial/oxidative liability proposal)
- Lin & Sibille 2015, Mol Psychiatry — PMID 25600109 (EIF2 translation suppression, stressed SST neurons)
- Seney et al. 2015, Neurobiol Dis — PMID 25315685 ("general vulnerability … independent of specific cell type")
- Tomoda et al. 2022, Mol Psychiatry — PMID 35145229 (PERK/eIF2α, preproSST processing load)
- Zhang et al. 2023, PNAS Nexus — PMID 37822767 (ribosomal-protein genes down, MDD + chronic stress)
- Newton et al. 2022, Biol Psychiatry — PMID 34861977 (OxPhos decrement in pyramidal, not SST)
- Arbabi et al. 2025, Mol Psychiatry — PMID 39237723 (sgACC LCM-seq; ATPase down in SST in SCZ)
- Tripp et al. 2012, Am J Psychiatry — PMID 23128924 (NTRK2/GAD2/SST down; BDNF dependency)
- Oh et al. 2019, Biol Psychiatry — PMID 30449530 (dendritic long-3'UTR BDNF → Sst)

## Flags / open items

1. **[VERIFY] BDNF breadth**: "reduced in pyramidal populations" is verified only
   for L2/3 IT (stated in the Doc's Fig. 2 section). Check Table T2 for other
   excitatory subclasses before the Discussion sentence stands as written.
2. **Do NOT use integrated-stress-response / eIF2α language for our own GSEA
   results**: the stress-named Reactome interaction sets (SRP, GCN2/amino-acid
   deficiency, starvation, NMD) are 55–83% ribosomal-protein genes — they are the
   ribosome signal under other names. EIF2 appears only as *cited literature*.
3. Supp figure numbers: Doc currently ends at S10; the Fig. 5 companions
   (IVW-vs-pseudobulk concordance S8-internal, interaction/volcanoes S11-internal,
   module scores S13-internal, Xenium stratum check) need final S numbers.
4. "Finally → Next" edit in the Fig. 4 section opener when Fig. 5 lands last.
5. Held in reserve (repo-only, not drafted in): Jens-cohort replication;
   SEA-AD A9 CPS contrast (translation declines in *spared* strata along CPS —
   argues the graded modules are not a generic dying-neuron signature).
6. Panel-d strip wording: "Shared: ubiquitin" block is FDR-significant only in
   the depleted stratum (NES positive everywhere); "shared" there means
   direction, not significance — legend wording covers it, keep an eye on it.
7. Sensitivities not yet run: S1 (FDR < 0.05 strata), S2 (downsample),
   S3 (leave-one-cohort-out), S4 (Kiss age), S5 (leave-Sst_22-out),
   S9 (depth dichotomy), S10 (depleted-vs-rest).
