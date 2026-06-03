# Literature context — what's known vs novel

A targeted search of the SCZ literature for the OxPhos / mitochondrial findings in inhibitory neurons (especially Sst), and the cholesterol / SREBF2 findings in glia and excitatory neurons. Goal: position our findings against the field and identify what's genuinely novel.

Search performed 2026-05-27. Caveats at the bottom.

---

## 1. The mitochondrial-interneuron story — convergent with Ni et al. 2019, extending it in three ways

### The most directly relevant prior work: Ni et al. 2019

**[Ni, Wang, Brennand et al. 2019, *Mol Psychiatry*](https://pmc.ncbi.nlm.nih.gov/articles/PMC6813882/)** is the foundational prior result that our findings should be positioned against. They used iPSC-derived developing cortical interneurons (CINs, MGE-lineage, mixed PV+SST) from 15 HC vs. 15 SCZ lines and showed:

- **Transcriptional**: dysregulated OxPhos genes, particularly Complex I subunits ND2 and ND4L (both *mtDNA-encoded*); also NSUN3 (mito tRNA methyltransferase), SLC25A15 (mito ornithine transporter), AGO3, FOSL2 — all DOWN in SCZ
- **Functional**: reduced maximal respiration, reserve capacity, NAD levels in SCZ CINs; elevated oxidative stress (DCF assay)
- **Cell-type-specific**: deficit NOT observed in iPSC-derived glutamatergic neurons → interneuron-specific
- **Therapeutically reversible**: ALA/ALC supplementation rescued OxPhos function and arborization defects

What Ni et al. did *not* do: separate SST from PV; investigate protein import / ETC assembly / mitochondrial quality control specifically; discuss PGC-1α / NRF1 / TFAM transcriptional axis.

### Our work is convergent and extends Ni et al. in three ways

1. **Postmortem adult brain (vs. developing iPSC-derived cells)** — confirms the OxPhos deficit survives differentiation, maturation, and ongoing disease processes. Removes the developmental-only / iPSC-artifact alternative interpretations.
2. **Cell-type-resolved within interneurons**: SST > Vip > Lamp5 > Pvalb. Ni et al. couldn't address this with their mixed CIN population.
3. **Mechanism refinement**: our data point specifically to **mitochondrial protein import + ETC complex assembly + quality control + mitoribosome failure**, not just "ETC dysfunction broadly." Key gene-level candidates: AFG3L2, TOMM40, UQCC2, PINK1 (none discussed in Ni et al.).

### A striking mechanistic consistency between the two studies

**Ni et al. 2019**: ND2 and ND4L (mtDNA-encoded Complex I subunits) are DOWN in SCZ CINs.
**Our work**: Mitoribosome (`GOMF_STRUCTURAL_CONSTITUENT_OF_RIBOSOME`) is DOWN in Sst (NES = -2.00, padj = 4e-4).

These connect mechanistically: mtDNA-encoded ETC subunits are translated by mitoribosomes. If mitoribosomes are downregulated (as we see at the transcript level), mtDNA-encoded subunits can't be made — which Ni et al. observed directly at the mtDNA transcript level. 10x snRNA-seq typically misses mtDNA reads themselves, but the *upstream machinery* is captured in our data and is consistent with their direct measurement.

This is the kind of cross-platform convergence that significantly strengthens both studies. Worth foregrounding if writing up.

### What the field's other relevant prior work shows

- **[Chung, Volk, Arion, Lewis et al. 2018, *Mol Psychiatry*](https://www.nature.com/articles/mp2017216)** — laser-capture microarray on DLPFC layer-3 PV neurons (36 SCZ-control pairs); >800 DEGs, mitochondrial function the top enriched pathway. **The** canonical "PV neurons have a mitochondrial deficit" reference. Our work adds that SST is *more* affected than PV at the single-cell-resolved level.
- **[Steullet, Cabungcal, Do et al. 2017, *Mol Psychiatry*](https://www.nature.com/articles/mp201747)** and **[Cabungcal et al. 2013, *PNAS*](https://www.pnas.org/doi/10.1073/pnas.1300454110)** (Do lab) — established the PV oxidative-stress / glutathione / perineuronal-net hypothesis. The standard "PV is the mitochondrially-vulnerable cell" framing.
- **[Glausier, Enwright & Lewis 2020, *Am J Psychiatry*](https://pmc.ncbi.nlm.nih.gov/articles/PMC8195258/)** — bulk DLPFC: 41% of mito genes DE, 83% down. Centred on layer-3/5 pyramidal neurons and PV. **SST not flagged.**

### The SST loss literature treats SST mRNA reduction and mitochondrial dysfunction as separate

- **[Dienel, Dowling, Barile et al. 2023, *JAMA Psychiatry*](https://pubmed.ncbi.nlm.nih.gov/37647039/)** — established SST mRNA loss is *more severe* than PV in DLPFC and SCZ-specific. **Did not invoke mitochondria.**
- **[Morris, Hashimoto & Lewis 2008, *Cereb Cortex*](https://pmc.ncbi.nlm.nih.gov/articles/PMC2888087/)** — foundational SST mRNA reduction paper.

### Recent snRNA-seq studies that should have caught the SST OxPhos pattern but didn't

- **[Ruzicka et al. 2024, *Science*](https://www.science.org/doi/10.1126/science.adg5136)** — largest SCZ snRNA-seq atlas; headlined excitatory upper-layer neurons / neurodevelopment / synapse, **not interneuron OxPhos**.
- **[Bast, Hjerling-Leffler, Sullivan et al. 2025, medRxiv](https://www.medrxiv.org/content/10.1101/2025.03.14.25323827v1.full)** — concluded "mitochondrial dysregulation" is central, but assigned it primarily to **layer 2-3 excitatory neurons**, not SST.

### Revised verdict on novelty

The bare statement "interneurons have OxPhos deficits in SCZ" is **not** novel — Ni et al. 2019 established this with both transcriptional and functional evidence. What our work adds is:

| Contribution | Status |
|---|---|
| OxPhos deficit in cortical interneurons in SCZ | **Established** (Ni et al. 2019) |
| OxPhos deficit confirmed in postmortem adult human brain | **Convergent novel evidence** |
| Cell-type resolution within interneurons (SST > Vip > Lamp5 > Pvalb) | **Novel** |
| Mechanistic refinement to import/assembly/QC (vs broad OxPhos) | **Novel framing** |
| Mitoribosome down — mechanistic bridge to Ni's mtDNA-encoded subunit finding | **Novel mechanistic link** |
| Specific gene candidates: AFG3L2, TOMM40, UQCC2, PINK1, GOT2 | **Novel candidates** |
| SST being more affected than PV directly challenges PV-vulnerability paradigm | **Novel cell-type claim** |

A translational thread worth pulling: **Ni et al. showed ALA/ALC supplementation rescues the iPSC OxPhos phenotype.** If the postmortem signal is real and originates from the same mechanism, this could motivate biomarker/treatment-response follow-up.

---

## 2. Specific candidate genes — appear to be novel SCZ leads

| Gene | Established disease context | Prior SCZ literature link? |
|---|---|---|
| **AFG3L2** (mito AAA protease) | SCA28 ([Almajan et al. 2012, *JCI*](https://www.jci.org/articles/view/64604); [Patron et al. 2019, *EMBO J*](https://pmc.ncbi.nlm.nih.gov/articles/PMC6618114/)) | None found |
| **TOMM40** (outer-membrane translocase) | AD / APOE locus ([Zeitlow et al. 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC8226536/)) | None found in SCZ snRNA-seq |
| **PINK1** (mitophagy initiator) | Early-onset Parkinson's | None retrievable in SCZ context |
| **UQCC2** (Complex III assembly factor) | — | None found |

The **collective framing** — mitochondrial protein import + ETC complex assembly + quality control failure, rather than PGC-1α-axis biogenesis suppression — does not appear to be articulated this way in the published SCZ literature. This is a distinct mechanistic hypothesis from the standard "decreased OxPhos transcription" story and predicts different therapeutic targets (chaperone/import-machinery modulators rather than biogenesis activators).

---

## 3. Cholesterol / SREBF2 story — partly known, partly novel

### Already-established components
- **SREBF2 / cholesterol biosynthesis association with SCZ**: **[Le Hellard et al. 2008, *Mol Psychiatry*](https://pubmed.ncbi.nlm.nih.gov/18936756/)** (Scandinavian/German association); **[Steen et al. 2017, *Eur Neuropsychopharmacol*](https://www.sciencedirect.com/science/article/abs/pii/S0924977X16301213)** (SREBP-system genetic evidence).
- **Lipid biosynthesis replication**: **[Hubler et al. 2018, *Sci Rep*](https://www.nature.com/articles/s41598-018-25280-4)**.
- **Oligodendrocyte mechanism**: **[Zhou et al. 2021, *eLife*](https://elifesciences.org/articles/60467)** — Qki/SREBP2-dependent cholesterol biosynthesis controls myelination. Directly relevant to our Oligo cholesterol-down finding.

### Novel contribution
- The *cell-type-resolved, coordinated suppression of the entire mevalonate pathway via SREBF2 in excitatory neurons + Astro + Oligo, with OPCs going the opposite direction*, has not been reported with this resolution. Provides an explicit transcriptional readout consistent with the long-standing genetic association.

### To verify before manuscript claim
- Whether **SREBF2** appears in PGC3 SCZ prioritized credible-set gene list (Trubetskoy 2022 Supplementary Table 12). The search did not retrieve this directly — we should check our copy of the supplement.

---

## 4. Field-level reframing context (2023–2026)

The field is actively shifting from a pure synaptic frame to an integrated mito ↔ synaptic frame:

- **[Henkel et al. 2024, *Psychiatry Res*](https://www.sciencedirect.com/science/article/pii/S0165178124005055)** — "SCZ as impaired dynamic metabolic flexibility" — reframes the disease bioenergetically.
- **[Frydecka et al. 2025, *Int J Mol Sci*](https://www.mdpi.com/1422-0067/26/9/4415)** — review of OxPhos dysfunction in SCZ.
- **[Howes & Onwordi 2023, *Mol Psychiatry*](https://pmc.ncbi.nlm.nih.gov/articles/PMC10575788/)** — "synaptic hypothesis of SCZ v3" still privileges synaptic mechanisms (PV-pyramidal microcircuits) as the master mechanism.
- **[Bitanihirwe et al. 2025, *Schizophrenia*](https://www.nature.com/articles/s41537-025-00638-6)** — recent cell-type GABAergic gene-expression reductions in cingulate.

Our SST-specific OxPhos finding lands directly in this debate as evidence for the metabolic side — with the additional twist that it's the *less metabolically demanding* interneuron, which directly challenges the standard PV-vulnerability argument.

---

## 5. Summary verdict per finding

| Finding | Status | Prior reference(s) |
|---|---|---|
| Pan-cell-type synaptic suppression in SCZ | Known | multiple |
| Mito dysfunction in SCZ at pathway level | Known | Glausier 2020, Bast 2025 |
| Mito dysfunction in PV interneurons | Known | Chung 2018 |
| **OxPhos deficit in cortical interneurons broadly** | **Established (iPSC)** | Ni et al. 2019 |
| **Postmortem-adult confirmation of interneuron OxPhos deficit** | **Convergent novel evidence** | extends Ni 2019 |
| **Cell-type resolution within interneurons (SST > PV)** | **Novel** | — (contradicts PV-vulnerability model) |
| SST mRNA loss in SCZ | Known | Morris 2008, Dienel 2023 |
| **Linking SST mRNA loss mechanistically to OxPhos** | **Novel bridge** | connects Dienel 2023 & Ni 2019 |
| **Mitoribosome down — links to Ni's mtDNA-encoded subunit finding** | **Novel mechanistic link** | bridges to Ni 2019 |
| **AFG3L2 / TOMM40 / PINK1 / UQCC2 in SCZ** | **Novel candidate signal** | — |
| **Mito import + assembly + QC framing (vs PGC-1α biogenesis)** | **Novel framing** | — |
| ALA/ALC rescue of OxPhos deficits as treatment lead | Known (iPSC) | Ni et al. 2019 — translational thread to revisit |
| SREBF2 / cholesterol genetic association | Known | Le Hellard 2008, Steen 2017 |
| **Coordinated SREBF2-axis suppression at cell-type resolution** | **Novel resolution** | — |
| Cholesterol-myelin axis in oligos | Known | Zhou 2021, Hubler 2018 |

## Search caveats — honesty notes

1. **Paywall blocks**: full text of some Nature / AJP papers couldn't be retrieved. Claims about those rest on abstracts + PMC summaries — should be verified against full text before citing in a manuscript.
2. **Trubetskoy 2022 SREBF2 prioritization** — claim that SREBF2 is in PGC3 prioritized credible-set list should be verified against Supplementary Table 12 before being asserted.
3. **"No prior SCZ link" claims** for AFG3L2 / TOMM40 / PINK1 / UQCC2 are based on keyword searches and the absence of hits. A literature gap of this scale should be double-checked with a more targeted PubMed search (controlled vocabulary, MeSH terms) before being asserted as definitive.

## Implications for write-up framing

If we were to write this up as a manuscript or preprint, **Ni et al. 2019 should be cited as the most directly relevant prior work** and our findings positioned as "convergent and extending" rather than "novel." The angles with the most genuine novelty:

1. **Cross-platform convergence with Ni 2019** as the foundation: iPSC (Ni) and postmortem (us) both show interneuron-specific OxPhos deficits. The mitoribosome-down (us) ↔ mtDNA-encoded-Complex-I-subunits-down (Ni) connection is the kind of mechanistic bridge that strengthens both papers.
2. **Cell-type resolution within interneurons** (SST > Vip > Lamp5 > Pvalb) — directly challenges the PV-vulnerability paradigm. Ni et al. couldn't address this; the SST literature (Dienel et al.) didn't connect it to mitochondria.
3. **Mechanistic refinement** to import + assembly + QC machinery (AFG3L2, TOMM40, UQCC2, PINK1) — distinct from the broad-OxPhos framing and predicts different therapeutic targets.
4. **The SREBF2-axis collapse** as a secondary story — known SCZ genetic association, novel cell-type-resolved transcriptional readout.

The cholesterol story is biologically rich but less novel as the headline; better positioned as a "second hit" alongside the OxPhos story to argue for cell-class-selective metabolic dysfunction as a unifying theme.

**Translational connection to follow up**: Ni et al.'s ALA/ALC rescue result in iPSC interneurons suggests the SST OxPhos phenotype, if it shares the same mechanism, may be biomarker-relevant for treatment response. Worth a directed literature pull on ALA/ALC trials in SCZ and on antioxidant trials more broadly.
