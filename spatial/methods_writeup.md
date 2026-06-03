# Methods: Cell Typing, Depth Inference, and Platform Validation

## 1. Overview

We mapped cell types in 24 human DLPFC Xenium sections (12 SCZ, 12 control; 1.3M cells; 300-gene panel) to the SEA-AD MTG taxonomy (24 subclasses, 137 supertypes). The core challenge is cross-platform label transfer: Xenium's 300 genes capture only ~1% of the transcriptome measured by the snRNA-seq reference (137,303 cells × 36,601 genes), and platform-specific artifacts (optical detection limits, probe specificity) make standard integration approaches unreliable.

Our solution is a self-referencing correlation classifier that operates entirely within the Xenium feature space, combined with a neighborhood composition-based depth model trained on the SEA-AD MERFISH reference. Outputs are validated against the independent MERFISH spatial atlas (341,595 cortical cells, 27 neurotypical donors, 180-gene panel, manual layer annotations). The pipeline produces 1,225,037 QC-pass cells with cell type assignments, continuous cortical depth, and spatially-smoothed layer labels.

---

## 2. Cell Type Classification

### 2.1 Pipeline summary

```
Step 00: Create h5ad from raw Xenium data
Step 01: Spatial QC (negative probes, gene counts, total counts)
    → 1,298,687 cells pass basic QC (24 samples)

Step 02: MapMyCells HANN label transfer
    → class/subclass/supertype labels + confidence scores
    → Used as input for exemplar selection, not as final labels

Step 02b: Two-stage correlation classifier
    → Build centroids from top-100 HANN exemplars per type (from Xenium data)
    → Stage 1: Subclass assignment via Pearson correlation (24 types)
    → Stage 2: Supertype assignment within subclass
    → QC: Flag bottom 5% margin per sample + spatial doublets
    → Save centroids to disk for optional nuclear resolution (see code/nuclear_resolution/)

Step 03: Export transcript coordinates (for viewer)

Step 04: Cortical depth model (trained on MERFISH)
Step 05: Spatial domain annotation + layer assignment

Step 06: Viewer export (per-sample JSON + standalone HTML)
Step 07: Cell + nucleus boundary polygon export

Optional: Nuclear doublet resolution → code/nuclear_resolution/
    → See nuclear_resolution/README.md for details (not part of numbered sequence)

Analysis: Cortical cells → crumblr compositional regression (SCZ vs Control)
```

| QC Step | Cells | Lost | % Lost |
|---------|-------|------|--------|
| Raw (all cells) | 1,339,151 | — | — |
| Step 01: spatial QC (`qc_pass`) | 1,298,687 | 40,464 | 3.0% |
| **Step 02b: `corr_qc_pass` (default gate)** | **1,225,037** | **114,114** | **8.5%** |

*Cell counts reflect all 24 samples. `corr_qc_pass` combines spatial QC (step 01), 5th-percentile margin filter, and doublet suspect exclusion. An optional nuclear doublet resolution step (in `code/nuclear_resolution/`, not part of the numbered pipeline) can produce `hybrid_qc_pass` but was found to have negligible impact on downstream biology. Br2039 is excluded from downstream disease comparisons due to high white matter content (65%).*

### 2.2 Centroid construction from HANN exemplars

Allen Institute's MapMyCells HANN classifier provides initial cell type labels with bootstrap-based confidence scores (100 iterations, 0.5 subsampling). Rather than using these labels directly — HANN confidence thresholds aggressive enough to be useful removed ~25% of cells, disproportionately deep excitatory neurons and rare interneurons — we use the top 100 highest-confidence HANN exemplars per type (pooled across all Xenium samples) to build per-type centroid expression vectors. Expression is normalized (counts per 10k, log1p) and averaged. Because centroids are built from Xenium data itself, they inherently capture platform-specific expression characteristics, eliminating cross-platform normalization issues.

### 2.3 Two-stage hierarchical Pearson correlation

**Stage 1 (Subclass):** Each cell's normalized expression is correlated against all 24 subclass centroids. The winning subclass is assigned along with the correlation value and *margin* (best minus second-best correlation). Both cell and centroid vectors are z-scored before correlation.

**Stage 2 (Supertype):** Within the assigned subclass, each cell is correlated against only the supertype centroids belonging to that subclass. This prevents biologically impossible assignments (e.g., an L2/3 IT cell being called an Sst supertype).

The 300-gene panel provides sufficient marker resolution to validate assignments at the single-molecule level:

![Exemplar cell transcript-level classification](output/presentation/exemplar_transcript_classification.png)
*Figure 1. Transcript-level cell type identity for six exemplar subclasses. Each panel shows a single cell with individual marker transcript molecules (e.g., CUX2 for L2/3 IT, PVALB+GAD1+GAD2 for Pvalb). Scale bar = 10 um.*

### 2.4 QC: margin filtering and spatial doublet detection

**Margin filtering:** The bottom 5th percentile of subclass correlation margins per sample are flagged as low-confidence. Per-sample thresholds account for variation in data quality across sections. This threshold was calibrated against SEA-AD MERFISH ground-truth labels, where cells below the 5th percentile have ~80% subclass accuracy — see `docs/pipeline_qc_audit.md` for the full calibration analysis.

**Spatial doublet detection** identifies cells with biologically implausible marker co-expression:
- *Glut+GABA doublets:* Cells expressing 4+ of 7 GABAergic markers (GAD1, GAD2, SLC32A1, SST, PVALB, VIP, LAMP5) while also expressing glutamatergic markers. False-positive rate validated at 0.098% in snRNA-seq.
- *GABA+GABA doublets:* GABAergic cells co-expressing SST + PVALB + LAMP5 simultaneously (<0.01% in snRNA-seq).

Flagged cells show the expected doublet signatures: mixed marker expression from both neuronal classes and ~1.7-1.9x elevated UMI counts.

![Marker detection rates in doublets vs normal cells](output/presentation/doublet_marker_barplot.png)
*Figure 2. Marker detection rates comparing normal cells and doublets. Glut+GABA doublets express markers from both neuronal classes — a pattern expected from two co-captured cells.*

Transcript-level visualization directly confirms that doublet cells contain interleaved GABAergic and glutamatergic marker molecules within a single segmented cell body:

![Transcript molecules in doublet vs normal cells](output/presentation/doublet_transcript_examples.png)
*Figure 3. Individual transcript molecules within cell boundaries. Red = GABAergic markers; blue = glutamatergic markers. Normal cells show transcripts from one class only; doublets contain both.*

### 2.5 Why not Harmony?

We benchmarked against Harmony + kNN label transfer (the standard cross-dataset integration approach). Harmony showed systematic failure modes with 300-gene Xenium data: Sst was inflated to 12.1% of cells (vs 2.5% expected), driven by non-neuronal cells (VLMC, Astrocyte) being misassigned into GABAergic space. Overall subclass agreement with both the correlation classifier and HANN was only 69%. The fundamental issue is that Harmony corrects batch effects between same-technology datasets, but cross-modality integration (snRNA-seq vs Xenium) introduces systematic differences — detection budgets, probe specificity, and PCA on 300 curated marker genes — that go beyond batch effects. Self-referencing centroids built from Xenium data avoid these problems entirely.

![Subclass proportions by method](output/presentation/subclass_proportions_by_method.png)
*Figure 4. Subclass proportions across classification methods. Note Sst inflation and non-neuronal distortion in Harmony results.*

### 2.6 Nuclear doublet resolution (optional)

> **Note:** This step was empirically shown to have negligible impact on downstream compositional analysis (see `docs/pipeline_qc_audit.md`). The simplified QC pipeline (spatial QC + 5th-percentile margin filter + doublet exclusion) is now the default. This section is retained for scientific interest.

**Motivation:** Inspection of doublet cells revealed that mixed-type marker transcripts concentrate in the cytoplasmic compartment, suggesting most doublets arise from mRNA spillover between neighboring cells during segmentation rather than true co-expression. Since nuclei are spatially more isolated, nuclear-restricted transcripts should be less affected by spillover.

**Method:** For each sample, the pipeline builds a nuclear-only count matrix by intersecting per-molecule transcript coordinates with nucleus boundary polygons (Shapely STRtree spatial indexing). Each doublet-flagged cell is reclassified using nuclear-only counts against the pre-built correlation centroids:
- **Resolved:** Nuclear classification agrees with the non-doublet class identity (cytoplasmic spillover confirmed)
- **Persistent:** Nuclear classification still shows mixed-class identity (true doublet)
- **Nuclear-only:** Not flagged by whole-cell, but nuclear classification reveals a different class
- **Insufficient:** <50 nuclear UMIs

**Results across 24 samples:**

| Category | Count | % of Flagged |
|----------|-------|-------------|
| Whole-cell doublets (QC-pass) | 10,505 | 100% |
| Resolved (cytoplasmic spillover) | 7,941 | 75.6% |
| Persistent (true doublets) | 2,539 | 24.2% |
| Insufficient evidence | 25 | 0.2% |
| Nuclear-only (new detections) | 2,128 | — |

The overall resolution rate is **75.8%** (range: 68-83% across samples). Resolved doublets show high marker scores in the whole-cell compartment but low scores in the nuclear compartment, directly confirming cytoplasmic spillover as the mechanism.

![Nuclear doublet resolution summary](output/presentation/nuclear_doublet_resolution_summary.png)
*Figure 5. Nuclear doublet resolution outcomes across all 24 samples. The majority of whole-cell doublets (75.8%) are resolved by nuclear-only classification.*

![Nuclear doublet marker evidence](output/presentation/nuclear_doublet_marker_evidence.png)
*Figure 6. Marker expression evidence for doublet resolution. Resolved doublets show high whole-cell but low nuclear marker scores. Persistent doublets maintain high scores in both compartments.*

**Impact on downstream biology:** The nuclear resolution machinery does not meaningfully change compositional SCZ vs Control results. Switching from `hybrid_qc_pass` to `corr_qc_pass` changes the snRNAseq meta-analysis correlation from r=0.405 to r=0.397, with identical top FDR-significant cell types. See `docs/pipeline_qc_audit.md` for the full comparison.

---

## 3. Cortical Depth Prediction

### 3.1 Feature construction

Rather than predicting depth from a cell's own gene expression (which would be sensitive to individual cell misclassification), the model uses **local neighborhood composition** as features. For each cell, the K=50 nearest spatial neighbors are identified (ball tree algorithm), and the fraction of each subclass among those neighbors is computed. This produces a feature vector of length 2 × n_subclass: the first half encodes neighbor composition fractions, the second half is a one-hot encoding of the cell's own subclass. The key insight is that neighborhood composition captures local tissue context — a region dominated by L2/3 IT neurons is likely superficial cortex regardless of any individual cell's label — making the model robust to sporadic cell misclassification.

### 3.2 MERFISH reference training

The model is trained on SEA-AD MERFISH data (cells with manual "Normalized depth from pia" annotations). Training uses a donor-level split: the 3 donors with the fewest depth-annotated cells are held out for testing; the remaining donors form the training set. This ensures the model generalizes across individuals rather than memorizing section-specific patterns.

### 3.3 Model and performance

The depth model is a `GradientBoostingRegressor` (n_estimators=300, max_depth=5, learning_rate=0.1, subsample=0.8, min_samples_leaf=20). Performance on the held-out MERFISH test donors:

| Metric | Train | Test |
|--------|-------|------|
| R² | 0.93 | 0.89 |
| MAE | 0.050 | 0.069 |
| Pearson r | 0.96 | 0.95 |

Predictions are deliberately **not clamped** to [0, 1]. Cells in white matter receive depth > 1 and cells above the pia receive depth < 0, providing natural tissue boundary detection without requiring hard cutoffs.

![Depth model diagnostics](output/depth_model_normalized_diagnostics.png)
*Figure 7. Depth model training and validation. Predicted vs actual normalized depth for train and test (held-out donor) sets, with R² and MAE metrics.*

### 3.4 Depth coordinate system

| Normalized depth | Cortical region |
|-----------------|----------------|
| < 0.00 | Above pia (meninges) |
| 0.00 – 0.12 | Layer 1 |
| 0.12 – 0.47 | Layer 2/3 |
| 0.47 – 0.54 | Layer 4 |
| 0.54 – 0.71 | Layer 5 |
| 0.71 – 0.93 | Layer 6 |
| > 0.93 | White matter |

*Boundaries derived from pairwise excitatory neuron marker crossovers in SEA-AD MERFISH, validated against Xenium Control samples. See `code/analysis/derive_layer_boundaries.py` for derivation details.*

![Depth comparison: MERFISH vs Xenium](output/presentation/slide_depth_comparison.png)
*Figure 8. Spatial depth maps comparing MERFISH manual annotations and Xenium model predictions. Viridis colormap encodes normalized depth (dark = superficial, bright = deep).*

---

## 4. Spatial Domain Classification

### 4.1 Motivation

Not all tissue on a Xenium section is cortex. Sections may include pia/meningeal tissue (cell-sparse, dominated by astrocytes and microglia), vascular clusters (concentrated Endothelial and VLMC cells), and white matter. The depth model is trained on cortical tissue from MERFISH, so its predictions for non-cortical regions are unreliable. Spatial domain classification identifies these regions so they can be handled appropriately.

### 4.2 BANKSY spatial clustering

We use BANKSY (Singhal et al., Nature Genetics 2024) for spatial domain classification. BANKSY augments each cell's gene expression with spatial neighbor expression and expression gradients, then performs dimensionality reduction and Leiden clustering in this augmented feature space. This produces spatially coherent clusters by construction — neighboring cells with similar expression are grouped together, naturally recovering tissue domains without requiring separate neighborhood feature engineering.

**BANKSY parameters:** λ=0.8 (spatial weighting), Leiden resolution=0.3, k_geom=15 spatial neighbors, 20 PCA dimensions. Preprocessing: library-size normalization (target 10,000), log1p, z-scoring.

Each resulting cluster is classified by its cell type composition and mean predicted depth (from step 04):

| Domain | Rule |
|--------|------|
| **Vascular** | >50% Endothelial + VLMC |
| **White Matter** | >40% Oligodendrocyte AND mean depth > 0.80 |
| **L1 Border** | >50% non-neuronal AND mean depth < 0.20 → classified as **Cortical** with `banksy_is_l1=True` flag |
| **Neuronal Cortex** | Neuronal fraction > 20% AND 0 ≤ mean depth ≤ 0.90 |
| **Deep WM** | Mean depth > 0.80 (fallback) |
| **Cortical** | Default |

Key advantages:
1. **L1 border detection**: Shallow non-neuronal-dominated BANKSY clusters are correctly identified as L1 cortex rather than pia/meninges. Validated by MERFISH: L1 has ~81% non-neuronal composition (astrocytes, microglia, endothelial), yet is cortical tissue.
2. **White matter detection**: BANKSY clusters with high oligodendrocyte fraction (>40%) and deep mean depth (>0.80) are classified as WM, providing explicit white matter identification.
3. **Lower vascular threshold**: BANKSY clusters are spatially coherent by construction, so a 50% Endo+VLMC threshold reliably identifies vascular regions without false positives from scattered cortical vascular cells.

### 4.3 Aggregate domain breakdown

Across all 24 samples (1,225,037 QC-pass cells):

| Domain | Cells | % |
|--------|-------|---|
| Cortical (including L1 border) | 744,878 | 60.8% |
| Vascular | 219,030 | 17.9% |
| White Matter | 261,129 | 21.3% |

Note: The Vascular domain fraction (17.9%) is larger than the final Vascular layer fraction (6.8%) because BANKSY captures spatially coherent vascular-associated tissue including border cells. Spatial smoothing (Section 5) reassigns most border Vascular cells to cortical layers.

![Layer composition by sample](output/presentation/slide_layer_stacked_bar.png)
*Figure 9. Per-sample layer composition showing the proportion of cells in each cortical layer, white matter, and vascular domains.*

---

## 5. Layer Assignment

### 5.1 Depth binning

Discrete layer labels are assigned by binning the continuous depth predictions:

| Layer | Depth range |
|-------|------------|
| L1 | < 0.12 |
| L2/3 | 0.12 – 0.47 |
| L4 | 0.47 – 0.54 |
| L5 | 0.54 – 0.71 |
| L6 | 0.71 – 0.93 |
| WM | > 0.93 |

Vascular-domain cells are overridden to "Vascular" regardless of predicted depth, since depth predictions are unreliable for spatially isolated vascular clusters. L1 border cells (identified by BANKSY with `banksy_is_l1=True`) retain their depth-bin layer, which is typically L1 given their shallow position.

### 5.2 Spatial layer smoothing

Raw depth-bin layers produce noisy boundaries: individual cells may receive incorrect layer assignments due to local depth prediction noise, and border Vascular cells may be classified as Vascular despite being surrounded by cortical tissue. A 3-step spatial smoothing pipeline (`smooth_layers_spatial()` in `depth_model.py`) addresses these issues:

**Step 1: Within-domain majority vote (k=30, 2 rounds).** For each cell, the layer labels of its k=30 spatial nearest neighbors *within the same BANKSY domain* are tallied, and the cell is reassigned to the majority layer. This smooths noisy cortical layer boundaries without allowing reassignments to cross BANKSY domain borders (e.g., a cortical cell cannot be voted into Vascular). Two rounds of voting are applied for convergence.

**Step 2: Vascular border trim.** Border Vascular cells — those with >33% of spatial neighbors in cortical layers (L2/3, L4, L5, L6) — are reassigned to the most common non-Vascular layer among their neighbors. A secondary rule reassigns Vascular cells with >66% of neighbors in any non-Vascular layer (including WM and L1). This trims the Vascular domain from 17.9% (BANKSY domain) to 6.8% (smoothed layer).

**Step 3: BANKSY-anchored L1 contiguity.** Two sub-steps refine L1 assignment: (a) *Promotion*: cells flagged as `banksy_is_l1` with predicted depth < 0.20 and at least 5% L1 neighbors are promoted to L1, recovering cells that depth binning alone missed. (b) *Removal*: isolated L1 cells (non-BANKSY L1 cells with <20% L1 neighbors, or BANKSY L1 cells with <5% L1 neighbors) are reassigned to their neighbors' majority layer.

### 5.3 Aggregate layer distribution (spatially smoothed)

| Layer | Cells | % |
|-------|-------|---|
| L1 | 72,078 | 5.9% |
| L2/3 | 375,019 | 30.7% |
| L4 | 77,410 | 6.3% |
| L5 | 210,910 | 17.3% |
| L6 | 205,809 | 16.8% |
| WM | 199,769 | 16.4% |
| Vascular | 80,524 | 6.6% |

---

## 6. Cross-Platform Validation

All cell type assignments and depth predictions are validated against the SEA-AD MERFISH spatial atlas — an independent dataset from a different spatial transcriptomics platform (MERFISH, 180 genes) measuring the same tissue type (human temporal cortex).

### 6.1 Subclass-level proportion concordance

Median subclass proportions correlate at **r = 0.84** between Xenium Control samples (n=12) and MERFISH (n=27). Most subclasses fall near the identity line, with systematic deviations limited to a few types:

- **Oligodendrocyte**: slightly underrepresented in Xenium relative to MERFISH (~15% vs ~25%), likely reflecting the DLPFC vs MTG tissue difference and possibly lower detection efficiency for oligodendrocyte markers in the 300-gene panel
- **Astrocyte**: slightly overrepresented in Xenium (~19% vs ~13%)
- **Microglia-PVM**: overrepresented in Xenium (~8% vs ~4%), consistent with known Xenium detection biases for immune markers

![Subclass proportions: MERFISH vs Xenium](output/presentation/slide_proportion_scatter.png)
*Figure 10: Subclass proportions — Xenium Control median vs SEA-AD MERFISH median. r = 0.84, n = 24 subclasses. Dashed line = perfect agreement.*

The per-layer concordance is also strong. Within each cortical layer (L2/3, L4, L5, L6), Xenium and MERFISH subclass proportions correlate at r = 0.70–0.92, confirming that the depth model and layer assignments produce biologically accurate laminar distributions.

The correlation classifier tracks the MERFISH reference more closely than Harmony at both resolution levels:

| Level | Correlation Classifier r | Harmony r | n types |
|-------|-------------------------|-----------|---------|
| Subclass (Pearson, log-scale) | **0.80** | 0.73 | 24 |
| Supertype (Pearson, log-scale) | **0.73** | 0.60 | 110 |

![Subclass proportion benchmark](output/presentation/cell_typing_benchmark_subclass_proportions.png)
*Figure 11. Subclass proportions: Xenium vs MERFISH reference. Left: Correlation Classifier (r = 0.80). Right: Harmony (r = 0.73). Dashed line = perfect agreement.*

### 6.2 Supertype-level concordance

At the supertype level (131 types), the concordance is necessarily weaker — the 300-gene Xenium panel lacks discriminating markers for many within-subclass distinctions. Overall accuracy against MERFISH ground truth is **84.9%** at the subclass level when the correlation classifier is applied to MERFISH data directly. Most subclasses achieve >80% accuracy (Astrocyte 95%, Oligodendrocyte 96%, L2/3 IT 90%, Pvalb 90%). Notable exceptions include **L6 IT Car3** (69%), **Lamp5** (55%), **Lamp5 Lhx6** (55%), and **Sst** (69%) — types where the 300-gene panel lacks sufficient within-subclass markers to resolve supertypes cleanly.

### 6.3 Depth profile concordance

Median cortical depth per subclass correlates at **r = 0.96** (subclass) and **r = 0.95** (supertype) between Xenium predicted depth and MERFISH manual annotations, confirming the depth model is well-calibrated:

![Median depth comparison](output/presentation/slide_median_depth_by_celltype.png)
*Figure 12: Median depth from pia — Xenium vs MERFISH. Left: subclass level (r = 0.96, n = 23). Right: supertype level (r = 0.95, n = 131). Near-perfect agreement confirms spatially coherent laminar assignments.*

At finer resolution, supertype depth distributions within each subclass show close agreement between MERFISH manual depth annotations and Xenium predicted depth across the full laminar gradient:

![Supertype depth distributions — Glutamatergic](output/presentation/supertype_depth_violins_glutamatergic.png)
*Figure 13. Supertype depth distributions for glutamatergic subclasses. Green = MERFISH manual depth, orange = Xenium predicted depth. Supertypes ordered by median MERFISH depth.*

![Supertype depth distributions — GABAergic](output/presentation/supertype_depth_violins_gabaergic.png)
*Figure 14. Supertype depth distributions for GABAergic subclasses. Interneuron supertypes show broader depth distributions than excitatory types, consistent with their wider laminar spread.*

![Supertype depth distributions — Non-neuronal](output/presentation/supertype_depth_violins_nonneuronal.png)
*Figure 15. Supertype depth distributions for non-neuronal subclasses. Non-neuronal types are distributed across all cortical depths, with astrocytes and microglia spanning the full column and oligodendrocytes concentrated in deep cortex/white matter.*

Full depth profile comparisons (proportion vs depth curves for all subclasses) show strong agreement between Xenium Control samples and MERFISH reference, with minor divergences only at the WM boundary. See Figures 9–10 in the [Depth-Stratified Analysis Report](output/depth_proportions/DEPTH_STRATIFIED_ANALYSIS_REPORT.md).

### 6.4 Classifier validation summary

The correlation classifier (step 02b) achieves r = 0.81 against MERFISH proportions (Controls only), compared to r = 0.73 for Harmony-based integration. Critically, Harmony misclassified non-neuronal types into GABAergic categories (e.g., VLMC classified as OPC 82% of the time) and inflated Sst proportions to 12.1% vs the expected 2.5%.

---

## 7. Robustness

The top disease signals are robust across QC configurations. At the supertype level with the default pipeline (corr QC, 5th-percentile margin), 5 supertypes reach FDR < 0.05 (L5/6 NP_3, L2/3 IT_8, L6b_5, Pax6_3, Pvalb_13) and 8 reach FDR < 0.10. Effect sizes correlate with Nicole's snRNAseq meta-analysis at r = 0.432 (all supertypes) and r = 0.466 (neuronal only).

The same top hits persist regardless of:
- **QC gate**: Switching between `corr_qc_pass` and the legacy `hybrid_qc_pass` (which includes nuclear doublet resolution) changes median |delta logFC| by only 0.0065 across all 106 shared supertypes, with no cell type changing direction of effect (see `docs/pipeline_qc_audit.md`).
- **Margin threshold**: Varying the margin filter from 1st to 10th percentile sharpens the signal but does not change the top hits.
- **Classifier hierarchy**: Flat vs hierarchical correlation classifier yields the same disease signals.
- **Doublet handling**: Including or excluding resolved doublets has negligible impact.

---

## 8. Hierarchy of Evidence

| Level | What | Concordance | Confidence |
|-------|------|-------------|------------|
| **Subclass proportions** | Xenium vs MERFISH | r = 0.84 | High — both are spatial platforms measuring same tissue types |
| **Subclass depth** | Xenium vs MERFISH | r = 0.96 | High — near-perfect laminar agreement |
| **Per-layer proportions** | Xenium vs MERFISH | r = 0.70–0.92 | High within layers, lower for sparse types |
| **Supertype proportions** | Xenium vs MERFISH | r = 0.35 (per-donor) | Low — insufficient markers for many supertype distinctions |

---

## 9. Related Documents

| Document | Focus |
|----------|-------|
| [SCZ Compositional Findings](scz_compositional_findings.md) | Disease effects, snRNAseq concordance, confidence tiers |
| [Depth-Stratified Analysis Report](output/depth_proportions/DEPTH_STRATIFIED_ANALYSIS_REPORT.md) | Per-layer and CLR depth × diagnosis results with figures |
| [Supertype Classification Confidence](output/marker_analysis/SUPERTYPE_CLASSIFICATION_CONFIDENCE_REPORT.md) | Per-supertype confidence ratings and Sst fragility analysis |
| [Panel Design & Supertype Classification](output/marker_analysis/XENIUM_PANEL_DESIGN_AND_SUPERTYPE_CLASSIFICATION.md) | Cross-platform marker adequacy, add-on gene recommendations |
