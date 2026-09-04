# Sibille-lab SST interneuron literature: molecular / pathway claims

Focused review compiled 2026-08-29 to contextualize the 7-cohort SCZ snRNA-seq
Sst-supertype strata analysis. Scope: ~2005–2026, Etienne Sibille senior or first
author, plus close trainees (Lin, Seney, Guilloux, Douillard-Guilloux, Tripp,
Oh, Newton, Shukla, Fee, Prevot, Tomoda, Arbabi, Kiss).

Every citation, p-value, effect size and quoted phrase below was pulled
programmatically from PubMed E-utilities, Europe PMC full-text XML, or the
bioRxiv/medRxiv full text on 2026-08-29 — none written from memory.
Where a claim could not be verified against a retrieved source it is flagged.

**Method caveat for the "absence" claims.** Gene-absence statements below
(VGF, RASGRF2, TOMM40) come from PubMed `Sibille E[au] AND <gene>` searches,
which index **title/abstract/MeSH only**, not full text. Read them as "never
promoted to a headline finding", not "never mentioned anywhere". Counts as run:
VGF 0, RASGRF2 0, TOMM40 0, calbindin 0, PINK1 1 (Glorioso 2011), CALB1 1
(Pabba 2017), ribosomal 2 (Zhang 2023 + one unrelated 2003 record). Full-text
keyword checks were run separately and are noted inline where they apply.

**Legend for module relevance** — the three modules from our analysis:
- **[TRANS]** cytosolic translation / ribosome
- **[OXPHOS]** oxidative phosphorylation / mitochondria
- **[SYN-NT]** synaptic + neurotrophic program (VGF, NTRK2/BDNF, RASGRF2, GAD2)
- **[SST]** SST mRNA itself
- **[CALB1]** upper-layer / CALB1+ SST subtype identity

---

## 1. The human postmortem SST-reduction series (depression, aging, corticolimbic)

### 1.1 Erraji-Benchekroun L, Underwood MD, Arango V, Galfalvy H, Pavlidis P, Smyrniotopoulos P, Mann JJ, Sibille E (2005). *Molecular aging in human prefrontal cortex is selective and continuous throughout adult life.* Biol Psychiatry 57(5):549-58. PMID 15737671.

- **Tissue**: two postmortem human PFC areas, 39 subjects, ages 13–79.
- **Claims**: ≥540 genes change robustly with age; changes progressive across
  adult life. Age-**up** transcripts "mostly of glial origin and related to
  inflammation and cellular defenses"; age-**down** transcripts "mostly
  neuron-enriched … relating to cellular communication and signaling."
  Concludes the restricted scope "suggests cellular populations or functions
  that are selectively vulnerable during aging."
- **Relevance**: the founding "selective cellular vulnerability" frame. No
  ribosome or OxPhos claim. Sets up the aging axis that later becomes the
  SST-vulnerability line. **[SYN-NT]** loosely (signaling/communication down).

### 1.2 Sibille E, Wang Y, Joeyen-Waldorf J, Gaiteri C, Surget A, Oh S, Belzung C, Tseng GC, Lewis DA (2009). *A molecular signature of depression in the amygdala.* Am J Psychiatry 166(9):1011-24. PMID 19605536.

- **Tissue/model**: human amygdala + ACC (N=14–16 matched male MDD/control
  pairs), cross-referenced against mouse UCMS arrays.
- **Claims**: a phylogenetically conserved MDD "molecular signature" in
  amygdala, reversed by antidepressants in mice, resolving into "two distinct
  oligodendrocyte and neuronal phenotypes" embedded in cohesive coexpression
  networks.
- **Relevance**: establishes the cross-species design the lab reuses. Not
  SST-specific and not module-specific; useful only as a lineage citation.

### 1.3 Sibille E, Morris HM, Kota RS, Lewis DA (2011). *GABA-related transcripts in the dorsolateral prefrontal cortex in mood disorders.* Int J Neuropsychopharmacol 14(6):721-34. PMID 21226980.

- **Tissue**: human DLPFC (BA9), triads of control/BPD/MDD, n=19 each.
- **Claims**: MDD showed **reduced SST mRNA** with *no* change in GAD67, GAD65
  or CR; BPD showed **reduced PV** with only trend-level SST reduction.
  Confirmed at the SST protein-precursor level. Notably: "The characteristic
  age-related decline in SST expression was not observed in MDD, as low
  expression was detected across age."
- **Relevance**: **[SST]** direct. Important negative for **[SYN-NT]** — in
  *this* region/diagnosis GAD1/GAD2 were **not** co-reduced with SST, unlike
  the sgACC result below (§1.6). The GAD2 co-reduction we see is therefore
  region- and diagnosis-contingent in their data.

### 1.4 Tripp A, Kota RS, Lewis DA, Sibille E (2011). *Reduced somatostatin in subgenual anterior cingulate cortex in major depression.* Neurobiol Dis 42(1):116-24. PMID 21232602.

- **Tissue**: human sgACC, 26 male + 25 female MDD vs 51 matched controls;
  Western blot subset n=42 pairs.
- **Claims**: SST mRNA reduced **38% in females, 27% in males**; prepro-SST
  protein reduced 19% in both. Age-related SST decline present in controls
  (Pearson R = −0.357, p = 0.005) but **absent in MDD** (R = −0.104, p = 0.234).
- **Relevance**: **[SST]** direct, and the "low SST across ages" pattern is the
  lab's core observation. Their conclusion — SST loss reflects "SST signaling
  and/or SST-bearing GABA neurons" — deliberately leaves cell-number vs
  per-cell ambiguous, which §1.6/§1.7 later resolve.

### 1.5 Guilloux JP, Douillard-Guilloux G, Kota R, Wang X, Gardier AM, Martinowich K, Tseng GC, Lewis DA, Sibille E (2012). *Molecular evidence for BDNF- and GABA-related dysfunctions in the amygdala of female subjects with major depression.* Mol Psychiatry 17(11):1130-42. PMID 21912391.

- **Tissue/model**: human lateral/basolateral/basomedian amygdala, n=21 MDD/
  control pairs; plus BDNF-heterozygous and BDNF exon-IV-KO mice.
- **Claims**: "Among the most robust findings were downregulated transcripts
  for genes coding for γ-aminobutyric acid (GABA) interneuron-related peptides,
  including somatostatin (SST), tachykinin, neuropeptide Y (NPY) and
  cortistatin." BDNF itself down at RNA **and** protein. The core MDD gene
  profile is named as **SST, NPY, TAC1, RGS4, CORT**, and is recapitulated by
  constitutive or activity-dependent BDNF reduction in mice, "with a common
  effect on SST and NPY". Frames the result as "linking the neurotrophic and
  GABA hypotheses of depression."
- **Relevance**: **[SYN-NT]** strong precedent — a *co-regulated neuropeptide/
  neurotrophic module* falling together with SST, driven by BDNF. This is the
  closest published analogue of our shared synaptic/neurotrophic program, minus
  VGF and RASGRF2 (neither is named here). **[SST]** direct.

### 1.6 Tripp A, Oh H, Guilloux JP, Martinowich K, Lewis DA, Sibille E (2012). *Brain-derived neurotrophic factor signaling and subgenual anterior cingulate cortex dysfunction in major depressive disorder.* Am J Psychiatry 169(11):1194-1202. PMID 23128924.

- **Tissue/model**: human sgACC, N=102 (51 MDD + 51 controls, 49% women); mouse
  Bdnf(+/−) and Bdnf(KIV) for the "BDNF-dependency" classification of 15 genes.
- **Claims**: In human MDD sgACC, **TRKB (NTRK2) expression was reduced but BDNF
  itself was not**. Depressed subjects showed downregulation of high- and
  intermediate-BDNF-dependency genes: "markers of dendritic targeting
  interneurons (SST, NPY, and CORT) and a GABA synthesizing enzyme (GAD2)."
  Changes extended to **BDNF-independent** genes PVALB and GAD1. Effects larger
  in men here (opposite to amygdala).
- **Relevance**: **[SYN-NT]** the single most directly citable Sibille hook for
  our shared module. They name **NTRK2 down** and **GAD2 down** alongside SST,
  in cingulate cortex, and explicitly organize the genes by BDNF/TrkB
  dependency. Our shared module (VGF, NTRK2/BDNF signaling, RASGRF2, GAD2, SST)
  is essentially this gene set observed at single-cell resolution.

### 1.7 Seney ML, Tripp A, McCune S, Lewis DA, Sibille E (2015). *Laminar and cellular analyses of reduced somatostatin gene expression in the subgenual anterior cingulate cortex in major depression.* Neurobiol Dis 73:213-9. PMID 25315685.

- **Tissue**: human sgACC, in situ hybridization in cohorts with known
  tissue-level SST reduction.
- **Claims**: "SST mRNA levels were lower **across all cortical layers**…
  Expression levels **per cell** were also lower, but the **density of labeled
  neurons did not differ**." More robust in females. Conclusion: the pattern
  suggests "a **general vulnerability of SST neurons independent of specific
  cell type**."
- **Relevance**: **CRITICAL POINT OF CONTRAST.** This is the lab's explicit
  claim that the SST deficit is *laminarly uniform* and *not subtype-specific*.
  Our stratum-specific, upper-layer/CALB1+-graded result directly challenges it
  — at higher resolution, in a different disorder, and with abundance rather
  than per-cell mRNA as the depletion axis. Cite this when framing the strata
  analysis as an advance. **[CALB1]**, **[SST]**.

### 1.8 Douillard-Guilloux G, Lewis D, Seney ML, Sibille E (2017). *Decrease in somatostatin-positive cell density in the amygdala of females with major depression.* Depress Anxiety 34(1):68-78. PMID 27557481.

- **Tissue**: human amygdala (lateral, basolateral, basomedial nuclei), ISH,
  N=10/group female.
- **Claims**: significant **reduction in density of SST-labeled neurons** across
  nuclei; SST mRNA per neuron unchanged in lateral/basolateral but lower in
  basomedial; **total cell density unchanged**. Conclusion: "these results
  suggest the possibility of a change in SST cell **phenotype rather than cell
  death**."
- **Relevance**: **[SST]** + directly relevant to how we interpret Sst supertype
  *depletion*. Their preferred reading of apparent SST-cell loss is
  dedifferentiation / phenotype loss, not death — the same interpretive fork we
  face for CALB1+ upper-layer Sst depletion in snRNA-seq.

### 1.9 McKinney BC, Lin CW, Oh H, Tseng GC, Lewis DA, Sibille E (2015). *Hypermethylation of BDNF and SST genes in the orbital frontal cortex of older individuals: a putative mechanism for declining gene expression with age.* Neuropsychopharmacology 40(11):2604-13. PMID 25881116.

- **Tissue**: human OFC, 22 younger (<42y) vs 22 older (>60y).
- **Claims**: 10/26 BDNF CpGs and **8/9 SST CpGs hypermethylated** in older
  individuals. DNAm in SST 5'UTR and first exon/intron negatively correlated
  with SST expression (r = −0.48, p<0.01; r = −0.63, p<0.001). "An expanded set
  of BDNF- and GABA-related genes exhibited similar age-related changes."
- **Relevance**: **[SST]**, **[SYN-NT]**. Offers an epigenetic mechanism for the
  *shared* (all-subtype) SST/BDNF-module downregulation, distinct from whatever
  drives the graded stratum-specific effects.

---

## 2. The SST-neuron vulnerability hypothesis (mechanism papers)

### 2.1 Lin LC, Sibille E (2013). *Reduced brain somatostatin in mood disorders: a common pathophysiological substrate and drug target?* Front Pharmacol 4:110. PMID 24058344. **[open access]**

- **Type**: review; the lab's explicit statement of *why* SST neurons should be
  intrinsically vulnerable.
- **Named candidate mechanisms** (verbatim from the abstract): "nitric oxide
  induced **oxidative stress**, **mitochondrial dysfunction**, high
  inflammatory response, **high demand for neurotrophic environment**, and
  overall aging processes."
- **Body text, oxidative/mitochondrial section** (retrieved from PMC3766825
  full-text XML): "Depressed states in mood disorders are associated with
  **decreased brain energy generation**"; "high baseline oxidative stress could
  be an intrinsic characteristic of vulnerable neuronal populations"; and the
  mechanistic anchor — neuronal nitric oxide synthase (nNOS) and NADPH
  diaphorase "are extensively and almost exclusively co-localized with
  somatostatin and neuropeptide Y … hence providing a neurochemical basis for
  high susceptibility of somatostatin-expressing neurons to generate oxidative
  stress."
- **Neurotrophic section**: BDNF–TrkB "is one of the key mediators for
  maintaining normal somatostatin gene expression"; TrkB signaling is itself
  vulnerable to inflammation and glucocorticoids, and "mild oxidative stress
  inhibits tyrosine phosphatases activity … potentially leading to impaired
  TrkB downstream signaling." CORT and NPY reductions are "partly downstream"
  of BDNF.
- **Relevance**: **[OXPHOS] — the key hypothesis-level precedent.** The lab
  predicted, in 2013, that SST neurons carry an intrinsic oxidative/
  mitochondrial liability. They never measured it in SST cells. Our
  stratum-graded OxPhos/mitochondrial downregulation (PINK1, TOMM40) is the
  first cell-resolved test of that prediction. Also **[SYN-NT]**.
  Note: **no translation/ribosome mechanism is proposed in this review** —
  that appears only in 2015 (§2.2).

### 2.2 Lin LC, Sibille E (2015). *Somatostatin, neuronal vulnerability and behavioral emotionality.* Mol Psychiatry 20(3):377-87. PMID 25600109.

- **Tissue/model**: Sst-KO and Sst-HZ mice; laser capture microdissection of
  cortical SST-positive interneurons vs pyramidal neurons after chronic stress.
- **Claims**:
  - Sst-KO mice show elevated behavioral emotionality, high basal plasma
    corticosterone, and "reduced gene expression of **Bdnf, Cortistatin and
    Gad67**."
  - "cortical SST-positive interneurons display significantly **greater
    transcriptome deregulations after chronic stress compared with pyramidal
    neurons**."
  - **The central mechanistic sentence**: "**Protein translation through
    eukaryotic initiation factor 2 (EIF2) signaling**, a pathway previously
    implicated in neurodegenerative diseases, **was most affected and
    suppressed in stress-exposed SST neurons**."
  - Pharmacologically activating EIF2 signaling (EIF2 kinase inhibition)
    mitigated stress-induced emotionality.
  - Explicit conclusion: "deregulated EIF2-mediated protein translation may
    represent a mechanism for vulnerability of SST neurons."
- **Relevance**: **[TRANS] — THE precedent.** This is the answer to the
  question of whether the translation finding has *any* prior basis in the
  Sibille corpus: yes, and it is central, cell-type-specific, and
  directionally concordant (suppressed translation in SST neurons). Caveats
  for our framing: it is (i) mouse, (ii) chronic stress not SCZ, (iii)
  *initiation-factor* signaling rather than ribosomal-protein/cytosolic-
  translation gene sets, and (iv) pan-SST, with no subtype stratification.
  Also **[SYN-NT]** via Bdnf/Gad67/CORT co-reduction in Sst-KO.

### 2.3 Lin LC, Sibille E (2015). *Transcriptome changes induced by chronic psychosocial/environmental or neuroendocrine stressors reveal a selective cellular vulnerability of cortical somatostatin (SST) neurons, compared with pyramidal (PYR) neurons.* Mol Psychiatry 20(3):285. PMID 25754192.

- Companion/summary item to §2.2 (title-only record in PubMed). Useful as a
  compact citation for the "SST > PYR vulnerability" claim.

### 2.4 Tomoda T, Sumitomo A, Newton D, Sibille E (2022). *Molecular origin of somatostatin-positive neuron vulnerability.* Mol Psychiatry 27(4):2304-2314. PMID 35145229. (Preprint: bioRxiv 2021.02.16.431515.)

- **Tissue/model**: mouse prefrontal cortex, chronic psychosocial stress;
  cell-type-restricted genetic manipulation in SST+ vs pyramidal neurons.
- **Claims**:
  - "a selective vulnerability of SST+ neurons through **exacerbated unfolded
    protein response (UPR) of the endoplasmic reticulum (ER), or ER stress**,
    in the prefrontal cortex."
  - Of the three UPR arms, **SST+ neurons are enriched in the PERK (EIF2AK3)
    pathway** over IRE1 and ATF6 (per the preprint's GSEA). PERK phosphorylates
    eIF2α and suppresses translation — mechanistically continuous with §2.2.
  - Genetically suppressing ER stress **in SST+ neurons but not pyramidal
    neurons** normalized behavioral emotionality.
  - **Cause identified as the SST peptide itself**: forced expression of the
    precursor **preproSST** in SST+ neurons induces ER stress, whereas mature
    SST or processing-incompetent preproSST does not; "psychosocial stress
    induces **SST protein aggregation** under elevated ER stress conditions."
  - Analogy drawn to preproinsulin/β-cell dysfunction in diabetes — "a
    universal mechanism for proteinopathy … induced by excess processing of
    native endogenous proteins."
  - Human corroboration cited (not generated): increased ER chaperones
    GRP78/HSPA5/BiP, GRP94 and calreticulin in postmortem temporal cortex of
    MDD suicides.
- **Explicitly absent** (checked the preprint full text): no ribosome /
  ribosomal-protein claim, no mTOR, **no mitochondrial, OxPhos, electron
  transport, PINK1 or TOMM claim**.
- **Relevance**: **[TRANS]** — the mechanistic bridge. If our translation module
  is down in the depleted upper-layer strata, the Sibille prediction is that
  it reflects PERK/eIF2α-mediated translational attenuation driven by
  proteostatic overload from SST peptide processing. Note the tension: their
  model says ER stress is **up** and translation down as a *consequence*.
  **[SST]** — and it inverts the usual causal arrow: SST peptide load is the
  *cause* of vulnerability, not merely a marker of it.

---

## 3. Cell-type-resolved profiling (LCM-seq) — the direct methodological ancestors

### 3.1 Shukla R, Prevot TD, French L, Isserlin R, Rocco BR, Banasr M, Bader GD, Sibille E (2019). *The relative contributions of cell-dependent cortical microcircuit aging to cognition and anxiety.* Biol Psychiatry 85(3):257-267. PMID 30446205.

- **Tissue/model**: mouse frontal cortex, young (2mo, n=9) vs old (22mo, n=12);
  LCM + RNA-seq of PYR, Vip, Sst, Pvalb cells; FISH validation.
- **Claims**: the four cell types show distinct age-related transcriptomes
  "affecting **metabolic** and cell signaling pathways, and selective markers of
  neuronal vulnerability (**Ryr3**), resilience (**Oxr1**), and **mitochondrial
  dynamics (Opa1)**", implying high age-related vulnerability of PYC and
  "variable degree of adaptation" in GABAergic neurons. Behavioral correlations
  implicate "age-independent decreases in **synaptic and signaling pathways**,
  notably in PYCs and **somatostatin neurons**."
- **Relevance**: **[OXPHOS]** partial — mitochondrial *dynamics* and oxidative
  *resilience* genes are named, but the direction in Sst cells is described as
  adaptive rather than depleted (per the secondary summary, Opa1 was
  downregulated in pyramidal cells and **up**regulated in interneurons; this
  specific direction was retrieved from a secondary source, not the primary
  text — **flag before citing the direction**). **[SYN-NT]** for the synaptic
  decrement in Sst cells.

### 3.2 Newton DF, Oh H, Shukla R, Misquitta K, Fee C, Banasr M, Sibille E (2022). *Chronic stress induces coordinated cortical microcircuit cell-type transcriptomic changes consistent with altered information processing.* Biol Psychiatry 91(9):798-809. PMID 34861977. (Preprint: bioRxiv 2020.08.18.249995.)

- **Tissue/model**: mouse mPFC (cingulate 24a/24b/32), 5 weeks UCMS; LCM of 130
  cells per type per mouse for PYR, SST, PV, VIP; RNA-seq + GSEA + WGCNA.
- **DE counts** (preprint): 217 DE genes in PYR, **371 in SST**, 406 in PV, 545
  in VIP. Enriched pathways: 400 PYR, **411 SST (256 up / 155 down)**, 389 PV,
  345 VIP.
- **SST-cell result** (verbatim): "SST-cells were characterized by **reduced
  growth factor signalling** and up-regulation of pre-synaptic CAMs and
  **response to ER-stress** after UCMS." Down: insulin, EGF, FGF, PDGF
  signaling (4 clusters, 41 pathways). Up: proteostasis/ER-stress (4 clusters,
  31 pathways), driven by chaperones **Hspa5, Hspa14, Tor1a** and mRNA
  degradation genes **Eri1, Pde5a, Cnot7**. Up: post-synaptic cadherin
  adhesion (Smad4, Celsr1, Clstn2, Ctnnd1).
- **Their own translation gloss on SST cells**: "**ER-stress induces adaptive
  responses of reduced translation** and increased chaperone, ROS generation,
  protein quality control, and mRNA and protein degradation. SST-cells showed
  evidence of all such changes."
- **Bioenergetics, by cell type** — this is the key contrast:
  - **PYR**: down-regulation of oxidative phosphorylation, "primarily
    cytochrome C subunits (**Cox5a, Cox7b, mt-Co3**)"; also down translational
    machinery and proteasome.
  - **PV**: **up**-regulated bioenergetics/biosynthesis (45 pathways) driven by
    glycolysis (Gpd1), oxidative phosphorylation (**mt-Cytb, mt-Atp6**),
    nucleotide biosynthesis (Prps1); plus "increased mitochondrial translation
    and oxidative phosphorylation."
  - **VIP**: increased oxidative stress response (Sod2), increased ER-stress
    (Bcap31, Hspa9), increased Map3k5/ASK1 apoptotic signaling.
  - **SST**: OxPhos is *not* among the highlighted SST clusters.
- **Marker note**: "**Calb1** (calbindin) in SST and PV-cells, and Calb2
  (calretinin) in SST and VIP-cells" — CALB1 was used only as a validation
  marker; SST cells were **not** stratified by it.
- **Microcircuit-level**: WGCNA identified a network "enriched in **synaptic,
  bioenergetic, and oxidative stress response** genes that correlated with
  UCMS-induced behaviors"; UCMS increased PYR–SST (p=0.029) and PYR–PV
  (p=1.0×10⁻⁴) module co-expression.
- **Relevance**: **[OXPHOS]** — the most informative comparison in the whole
  corpus, and a partial *contradiction*: in their mouse chronic-stress LCM data
  the OxPhos decrement sits in **pyramidal** cells, while **PV** cells go up and
  SST cells are characterized by ER stress + growth-factor loss rather than
  bioenergetic collapse. **[TRANS]** — translational-machinery downregulation is
  reported for PYR explicitly, and inferred for SST as an ER-stress consequence.
  **[SYN-NT]** — reduced neurotrophic signaling in SST cells, though via
  EGF/FGF/PDGF/insulin rather than BDNF/NTRK2 (they note "these other growth
  factor pathways may represent additional neurotrophic deficits").
  **[CALB1]** — the negative: no subtype stratification attempted.

### 3.3 Arbabi K\*, Newton DF\*, Oh H, Davie MC, Lewis DA, Wainberg M, Tripathy SJ, Sibille E (2025). *Transcriptomic pathology of neocortical microcircuit cell types across psychiatric disorders.* Mol Psychiatry 30(3):1057-1068. PMID 39237723. (Preprint: medRxiv 2023.10.26.23297640.) *\*equal contribution.*

- **Tissue**: human **subgenual ACC**, Pittsburgh Brain Tissue Donation Program,
  76 subjects evenly split MDD / BD / SCZ / control; LCM-seq of 130 pooled cells
  per neuronal subtype (VIP, SST, PVALB, superficial and deep PYR) → 380 bulk
  transcriptomes from ~50,000 neurons.
- **DE burden**: 87.3% of DE genes in interneurons. PVALB largest (239 total;
  **SCZ: 91 = 25 up / 66 down**). SST: **SCZ 38 DE genes (31 up / 7 down)**, BD
  40, MDD 12. VIP: SCZ 24.
- **SST-cell pathway results** (Figure 2B, verbatim): "For SST cells, this
  included the **downregulation of genes involved in ATPase activity in SCZ
  (p = 5.5 × 10⁻³)**, glial-cell derived neurotrophic factor receptor signaling
  in SCZ and BD (p = 6.0 × 10⁻³, p = 2.3 × 10⁻³), lysosomal protein catabolic
  process (p = 2.3 × 10⁻³), and retinoic acid receptor signaling in SCZ and MDD
  (p = 0.025, p = 1.8 × 10⁻³). Notably, there was an **upregulation of genes
  involved in response to endoplasmic reticulum (ER) stress in SST cells for all
  three disorders (SCZ: p = 0.034, BD: p = 0.02, MDD: p = 0.023)**,
  corresponding to cell type-specific findings in a chronic psychosocial stress
  model of depression in rodents."
- **Contrast cells**: VIP showed "increased negative regulation of genes
  involved in … **translation in SCZ and BD** (p = 0.011, p = 5.0 × 10⁻⁴)";
  PVALB showed "downregulation of genes involved in **oxidoreductase activity**
  in BD and MDD (p = 4.8 × 10⁻⁴, p = 5.6 × 10⁻³)."
- **Transdiagnostic shared genes** include **RBFOX1** (listed among neurogenesis
  genes: AUST2, CNTNAP2, ERBB4, FGF14, LSAMP, NRG3, NPAS3, RBFOX1, ZNF536).
  Transdiagnostic shared processes: "glial neurotrophic factor signaling,
  response to endoplasmic reticulum stress, dendrite development, and
  **regulation of translation**."
- **Genetics**: "Genetic risk prominently manifested in **PVALB cells in SCZ and
  MDD, less so in SST and VIP cells in SCZ**."
- **Relevance**: this is the closest existing human comparison to our analysis
  and the most citable. **[OXPHOS]** — *partial precedent*: **ATPase activity
  genes down in SST cells specifically in SCZ**, which is bioenergetics-adjacent
  and directionally concordant with our OxPhos finding, though not framed as
  mitochondrial/ETC and not tied to PINK1/TOMM40. **[TRANS]** — translation
  down-regulation is reported in **VIP**, not SST, in SCZ; "regulation of
  translation" is named as transdiagnostically disrupted. **[SYN-NT]** —
  GDNF-receptor signaling down in SST cells in SCZ/BD; RBFOX1 shared. **[SST]**
  — SST cells in SCZ were mostly **up** (31/38), a point worth reconciling with
  our per-cell SST downregulation (theirs is pooled-cell LCM with the cells
  selected *by* SST immunoreactivity, which biases against detecting per-cell
  SST loss). **[CALB1]** — no subtype resolution.

---

## 4. Bioenergetic / oxidative / ribosomal claims from the lab

### 4.1 Glorioso C, Oh S, Douillard GG, Sibille E (2011). *Brain molecular aging, promotion of neurological disease and modulation by sirtuin 5 longevity gene polymorphism.* Neurobiol Dis 41(2):279-90. PMID 20887790.

- **Tissue**: cross-cohort microarray, four human brain areas.
- **Claims**: neurological disease pathways largely overlap molecular aging;
  a low-expressing SIRT5 promoter polymorphism (SIRT5prom2) associates with
  accelerated molecular aging in cingulate (+9 years, p=0.004), driven by a core
  transcript set (+24 years, p=0.0004) "many of which were **mitochondrial,
  including Parkinson's disease genes, PINK-1 and DJ-1/PARK7**."
- **Relevance**: **[OXPHOS]** — the only place in the Sibille corpus where
  **PINK1** is named. Not SST-specific (bulk tissue), but it is the lab's own
  precedent for a mitochondrial-aging axis in cingulate cortex, and it ties
  PINK1 to accelerated molecular aging. Worth citing alongside our PINK1/TOMM40
  result as a same-region, same-gene antecedent.

### 4.2 Lin LC, Lewis DA, Sibille E (2011). *A human-mouse conserved sex bias in amygdala gene expression related to circadian clock and energy metabolism.* Mol Brain 4:18. PMID 21542937.

- **Claims**: "**mitochondrial-related gene groups** were identified as the top
  biological pathways associated with sexual dimorphism in both species";
  proposes "baseline differences in amygdalar circadian regulation of cellular
  metabolism" as a substrate for female MDD vulnerability.
- **Relevance**: **[OXPHOS]** background only — establishes the lab's interest
  in mitochondrial gene sets, and flags **sex** as a confound for any
  mitochondrial contrast.

### 4.3 Shukla R, Newton DF, Sumitomo A, Zare H, McCullumsmith R, Lewis DA, Tomoda T, Sibille E (2022). *Molecular characterization of depression trait and state.* Mol Psychiatry 27(2):1083-1094. PMID 34686766.

- **Tissue**: human sgACC across five clinical states (first episode n=20,
  remission after single episode n=15, recurrent episode n=20, remission after
  recurrence n=15, control n=20).
- **Claims**: "MDD-**trait** was associated with genes involved in inflammation,
  immune activation, and **reduced bioenergetics** (q<0.05) whereas MDD-**states**
  were associated with altered neuronal structure and reduced neurotransmission
  (q<0.05)." Deconvolution showed significant change in density of CRH-, SST-,
  or VIP-positive GABA interneurons (p<3×10⁻³). Bayesian network analysis put
  **oxidative stress** (q<2.05×10⁻³) among the putative causal pathways across
  phases.
- **Relevance**: **[OXPHOS]** — "reduced bioenergetics" is assigned to the
  *trait* (stable, illness-liability) axis, and oxidative stress is nominated as
  causal. Useful for arguing that a bioenergetic decrement is a persistent
  cortical feature rather than an episode artifact. **[SST]** via deconvolved
  SST-cell density change.

### 4.4 Zhang X, Eladawi MA, Ryan WG, Fan X, Prevoznik S, Devale T, Ramnani B, Malathi K, **Sibille E**, McCullumsmith R, **Tomoda T**, Shukla R (2023). *Ribosomal dysregulation: a conserved pathophysiological mechanism in human depression and mouse chronic stress.* PNAS Nexus 2(10):pgad299. PMID 37822767. **[open access]**

- **Tissue/model**: human postmortem DLPFC, OFC and ACC (MDD vs control) and
  mouse PFC after 4 weeks chronic variable stress (CVS); NAc as negative
  control; DEX-treated primary PFC neurons in vitro.
- **Claims**:
  - "**Ribosomal protein genes (RPGs) were down-regulated, and associated
    ribosomal protein (RP) pseudogenes were up-regulated in both conditions**"
    — 15 RPGs shared between human DLPFC MDD and mouse CVS PFC. Confirmed in
    independent DLPFC and ACC datasets. **Not** enriched in NAc in either
    species (regional specificity).
  - Seeded coexpression: "down-regulated RPGs **homeostatically regulated the
    synaptic changes** in both groups through an RP-pseudogene-driven
    mechanism"; RPGs were **negatively correlated** with synaptic-infrastructure
    themes (vesicles, axon, dendrite, presynapse, postsynapse).
  - "the RPG dysregulation was a **glucocorticoid-driven endocrine response to
    stress**" (reproduced by dexamethasone in primary neurons, blocked by the
    GR antagonist RU-486).
  - Reversed in remission; responds selectively to ketamine but not imipramine.
- **Limitation for our purposes**: bulk tissue, **not cell-type resolved**; SST
  neurons are not examined (verified by keyword search of the full text —
  "somatostatin" does not appear in the results).
- **Relevance**: **[TRANS] — the second and much more literal precedent.** This
  is the lab-affiliated paper that reports *ribosomal protein gene
  downregulation* per se, in human cingulate and prefrontal cortex, conserved
  across species, glucocorticoid-driven, and mechanistically coupled to synaptic
  gene changes. Our finding that the translation/ribosome module and the
  synaptic/neurotrophic program move together in Sst cells has a direct
  conceptual antecedent here — with the important difference that they read the
  coupling as *homeostatic* (RPGs down to compensate synaptic change), which is
  a specific alternative hypothesis we could test or acknowledge.

---

## 5. The BDNF/TrkB-dependent SST reduction line

### 5.1 Oh H, Piantadosi SC, Rocco BR, Lewis DA, Watkins SC, Sibille E (2019). *The role of dendritic brain-derived neurotrophic factor transcripts on altered inhibitory circuitry in depression.* Biol Psychiatry 85(6):517-526. PMID 30449530.

- **Tissue/model**: human PFC (n=19/group MDD vs control) + C57BL/6J chronic
  stress (n=12/group) + shRNA knockdown of long-3'UTR Bdnf transcripts.
- **Claims**: BDNF mRNAs bearing long 3'UTRs — the pool that traffics to distal
  pyramidal dendrites — are **selectively reduced** in MDD PFC, and their
  expression is **highly correlated with SST expression**. Same downregulation
  after chronic stress in mice. Bdnf L-3'UTR knockdown alone is sufficient to
  produce dendritic shrinkage, "cell-specific MDD-like gene changes (**including
  Sst downregulation**)", and depressive/anxiety-like behavior.
- **Relevance**: **[SYN-NT]** — provides the *spatial* mechanism for the shared
  module: SST-cell health depends on locally translated dendritic BDNF at the
  very compartment SST cells innervate. Directly supports interpreting a shared
  NTRK2/BDNF-signaling decrement across all Sst subtypes as an
  environment-driven (rather than subtype-intrinsic) effect.

### 5.2 Adjacent, non-Sibille but foundational: Glorioso C, Sabatini M, Unger T, Hashimoto T, Monteggia LM, Lewis DA, **Mirnics K** (2006). *Specificity and timing of neocortical transcriptome changes in response to BDNF gene ablation during embryogenesis or adulthood.* Mol Psychiatry 11(7):633-48. PMID 16702976.

- **Claim**: "BDNF appeared to be required to maintain gene expression in the
  **SST-NPY-TAC1 subclass of GABA neurons**, although the absence of BDNF did
  not alter their general phenotype as inhibitory neurons"; also altered IEGs
  (ARC, EGR1, EGR2, FOS, DUSP1, DUSP6) and RGS4.
- **Note**: Sibille is **not** an author; this is Mirnics/Lewis. It is the paper
  the Sibille reviews cite for "BDNF maintains SST expression", so cite it
  directly rather than attributing it to the Sibille lab.

---

## 6. The therapeutic arc: α5-GABA-A positive allosteric modulators

The consistent design logic: SST cells inhibit pyramidal **dendrites** via
**α5-subunit-containing GABA-A receptors**. If SST cells are hypofunctional, you
cannot easily restore the cell — but you can potentiate its postsynaptic
receptor and bypass the deficit.

### 6.1 Prevot TD, Li G, Vidojevic A, Misquitta KA, Fee C, Santrac A, Knutson DE, Stephen MR, Kodali R, Zahn NM, Arnold LA, Scholze P, Fisher JL, Marković BD, Banasr M, Cook JM, Savic M, Sibille E (2019). *Novel benzodiazepine-like ligands with various anxiolytic, antidepressant, or pro-cognitive profiles.* Mol Neuropsychiatry 5(2):84-97. PMID 31192221.
The introduction of **GL-II-73**; reversed stress-induced and age-related
working memory deficits where diazepam did not.

### 6.2 Prevot TD, Sumitomo A, Tomoda T, Knutson DE, Li G, Mondal P, Banasr M, Cook JM, Sibille E (2021). *Reversal of age-related neuronal atrophy by α5-GABAA receptor positive allosteric modulation.* Cereb Cortex 31(2):1395-1408. PMID 33068001.
GL-II-73 increases dendritic branching and spine number in vitro; 3-month
chronic treatment in old mice "significantly reverses age-related dendritic
shrinkage and spine loss in frontal cortex and hippocampus." Morphological
benefit outlasted the behavioral benefit after a 1-week washout.

### 6.3 Fee C, Prevot TD, Misquitta K, Knutson DE, Li G, Mondal P, Cook JM, Banasr M, Sibille E (2021). *Behavioral deficits induced by somatostatin-positive GABA neuron silencing are rescued by alpha 5 GABA-A receptor potentiation.* Int J Neuropsychopharmacol 24(6):505-518. PMID 33438026.
Chemogenetic **brain-wide SST+ cell silencing** produced elevated neuronal
activity, raised corticosterone, anxiety- and anhedonia-like behavior, and
short-term memory impairment; GL-II-73 rescued these. The cleanest causal
statement that the α5-PAM "bypasses low SST+ cell function."

### 6.4 Prevot T, Sibille E (2021). *Altered GABA-mediated information processing and cognitive dysfunctions in depression and other brain disorders.* Mol Psychiatry 26(1):151-167. PMID 32346158.
The framing review. "reduced signaling of the **SST+ neuron/α5-GABA-A receptor
pathway** contributes to cognitive dysfunctions, and … represents a novel
therapeutic target." Also Fee C, Banasr M, Sibille E (2017), *Somatostatin-
positive GABA interneuron deficits in depression: cortical microcircuit and
therapeutic perspectives*, Biol Psychiatry 82(8):549-559, PMID 28697889 — the
standard review citation for "selective vulnerability of GABAergic interneurons
that coexpress the neuropeptide somatostatin."

- **Relevance to our modules**: the deficit these compounds compensate is
  **functional output** (dendritic inhibitory tone), which is the downstream
  consequence of exactly the shared SST/GAD2/synaptic module we see down across
  all Sst subtypes. If our stratum-specific translation/OxPhos decrements
  identify *which* Sst cells are failing hardest, that is a plausible
  patient-stratification angle for this therapeutic program (and Sibille is
  co-founder/CSO of Damona Pharmaceutical, per the Arbabi 2025 COI statement).

---

## 7. Recent, directly adjacent (SCZ / aging / genetics)

### 7.1 Kiss D, Zhou X, Endresz N, Arbabi K, Segura AG, Felsky D, Diaconescu AO, Sibille E, Tripathy SJ (2026). *Cortical GABAergic neuron dysregulation in schizophrenia is age dependent.* Biol Psychiatry Glob Open Sci 6(1):100606. PMID 41159096.
14 bulk and cell-type-specific RNA-seq datasets, 1408 individuals (672 SCZ, 736
control), 3 neocortical regions. "Younger SCZ cases (age < 70 years) showed
**reduced PVALB and SST cell proportions**, while older cases showed unchanged
or increased proportions." Earlier onset → greater reductions. "robust evidence
for **reduced per-cell SST and vasoactive intestinal peptide mRNA** among
younger cases."
**Relevance**: **[SST]** — direct SCZ precedent for *both* Sst proportion
depletion *and* per-cell SST mRNA reduction, and establishes **age at death /
age of onset as a first-order moderator** we should check as a confound in any
strata analysis. (PDF already in `background/relevant papers/`.)

### 7.2 Chen Y, Hunter E, Arbabi K, Guet-McCreight A, Consens M, Felsky D, Sibille E, Tripathy SJ (2023). *Robust differences in cortical cell type proportions across healthy human aging inferred through cross-dataset transcriptome analyses.* Neurobiol Aging 125:49-61. PMID 36841202.
1142 subjects / 1429 samples, ages 15–90: "the largest changes reflecting
**fewer somatostatin- and vasoactive intestinal peptide-expressing
interneurons**, more astrocytes."
**Relevance**: the normative aging baseline against which SCZ Sst depletion must
be read.

### 7.3 Pabba M, Scifo E, Kapadia F, Nikolova YS, Ma T, Mechawar N, Tseng GC, Sibille E (2017). *Resilient protein co-expression network in male orbitofrontal cortex layer 2/3 during human aging.* Neurobiol Aging 58:180-190. PMID 28750307.
Mass-spec proteomics of human OFC **layer 2/3**, 15 young (15–43y) vs 18 old
(62–88y) males; 4193 proteins, 127 DE (65 up / 62 down), "e.g., GFAP, **CALB1**."
Hallmark categorization gave altered cell-cell communication (54%), deregulated
nutrient sensing (39%), and **loss of proteostasis (35%)**. Co-expression modules
themselves were conserved across age.
**Relevance**: **[CALB1]** — this is the **only** appearance of CALB1 in the
Sibille corpus, and notably it is in *upper-layer* cortex, as a
differentially-expressed aging protein. Thin, but it is the one place the lab
touches the upper-layer/CALB1 axis our strata analysis is built on. The
proteostasis theme also recurs.

### 7.4 Dos Santos FC, Zhou X, Clifford KP, Segura AG, Syeda AS, Felsky D, Lenze EJ, Mulsant BH, Tripathy S, Sibille E, Nikolova YS (2026). *A novel polygenic risk score indexing somatostatin-expressing inhibitory neurons predicts somatostatin-expressing cell proportions and severity of symptoms in late-life depression.* Biol Psychiatry Glob Open Sci 6(3):100685. PMID 41783159.
SST-PRS built from cis-eQTLs of SST-coexpressed genes; associated with reduced
dlPFC SST+ neuron proportion (R = 0.5, p = 1.783 × 10⁻⁴) and with higher MADRS
in LLD (t₃₆₉ = 2.267, p = .024).
**Relevance**: the lab's move toward genetic indexing of SST-cell burden; a
natural comparator for any genetics arm of our paper.

---

## Synthesis

**Translation / ribosome — strong precedent, and it is the lab's flagship
mechanism.** This is the module where the Sibille corpus is *least* silent.
Lin & Sibille 2015 (Mol Psychiatry 20:377) reports that in laser-captured
cortical SST neurons after chronic stress, "protein translation through
eukaryotic initiation factor 2 (EIF2) signaling … was most affected and
suppressed" — cell-type-specific, directionally concordant with our finding, and
causally probed (EIF2 kinase inhibition rescued behavior). Tomoda 2022 supplies
the upstream driver: SST+ neurons are selectively enriched for the **PERK/eIF2α**
arm of the UPR, and the trigger is proteostatic overload from preproSST
processing, with SST peptide aggregation. Zhang 2023 (PNAS Nexus, with Sibille
and Tomoda) is the most literal antecedent: **ribosomal protein genes are
downregulated in human MDD DLPFC/ACC and mouse CVS PFC**, glucocorticoid-driven,
and inversely coupled to synaptic gene sets. So our translation/ribosome result
is *not* new as a phenomenon. What is new is (a) that it is **graded across Sst
strata** rather than uniform, (b) that it is in **schizophrenia**, and (c) that
it appears as **cytosolic ribosome/translation gene sets in the depleted
upper-layer CALB1+ types specifically** — none of the above work resolved SST
subtypes at all.

**OxPhos / mitochondria — hypothesized by the lab in 2013, never demonstrated in
SST cells.** Lin & Sibille 2013 explicitly nominated "nitric oxide induced
oxidative stress, mitochondrial dysfunction" as the intrinsic vulnerability of
SST neurons, grounding it in nNOS/NADPH-diaphorase co-localization with SST.
But every cell-resolved test since has landed elsewhere: Newton 2022 found the
OxPhos decrement in **pyramidal** cells (Cox5a, Cox7b, mt-Co3) and an OxPhos
**increase** in PV cells, with SST cells characterized by ER stress and lost
growth-factor signaling. Arbabi 2025 is the nearest hit — **ATPase-activity genes
down in SST cells specifically in SCZ (p = 5.5 × 10⁻³)** — bioenergetics-adjacent
and same direction, but not framed as mitochondrial and not gene-level. PINK1
appears in the corpus only in Glorioso 2011 (SIRT5, bulk cingulate molecular
aging). **Our OxPhos/mitochondrial finding is therefore the first cell-resolved
confirmation of a 13-year-old Sibille prediction**, and that is the strongest
framing available.

**Synaptic / neurotrophic — extensive precedent, near-identical gene set.**
Tripp 2012 (AJP) reports **NTRK2 down, GAD2 down, SST/NPY/CORT down** in human
sgACC MDD, organized explicitly by BDNF dependency; Guilloux 2012 gives the
amygdala version (SST, NPY, TAC1, RGS4, CORT); Oh 2019 shows dendritic
long-3'UTR BDNF loss is sufficient to drive Sst downregulation. Our *shared*
module is essentially this, at single-cell resolution. **VGF and RASGRF2 never
appear as headline findings anywhere in the Sibille corpus** (zero title/abstract
hits each) — those two are ours.

**Directly citable hooks.** For translation: Lin & Sibille 2015 (PMID 25600109),
Tomoda 2022 (35145229), Zhang 2023 (37822767). For the OxPhos prediction:
Lin & Sibille 2013 (24058344) and, for the contrast, Newton 2022 (34861977).
For the human sgACC cell-resolved comparison in SCZ: Arbabi 2025 (39237723) —
same cortical region, same disorder, SST-cell ER stress up and ATPase down.
For the shared synaptic/neurotrophic module: Tripp 2012 (23128924) and Oh 2019
(30449530). For SST-cell depletion vs phenotype loss: Douillard-Guilloux 2017
(27557481). And the sharpest contrast to cite when introducing the strata:
**Seney 2015 (25315685)**, which concluded from laminar ISH that SST reduction
implies "a general vulnerability of SST neurons independent of specific cell
type" — precisely the claim our graded, CALB1-linked result revises.

**One caveat to carry.** Kiss 2026 (41159096) shows SCZ-associated SST
proportion loss is strongly **age-dependent** (present below ~70, absent or
reversed above). Any stratum-level depletion or module claim should be checked
against age at death before it is presented as a stable disease feature.
