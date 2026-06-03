# Cross-Platform Validation

## Purpose

This document assesses how well the Xenium spatial transcriptomics measurements agree with independent reference datasets — the SEA-AD MERFISH spatial atlas and a snRNA-seq schizophrenia meta-analysis. It establishes the measurement credibility that underpins all downstream disease analyses. For SCZ-specific findings and their interpretation, see [SCZ Compositional Findings](scz_compositional_findings.md).

---

## 1. Xenium vs SEA-AD MERFISH: Cell Type Proportions

The SEA-AD MERFISH dataset (341,595 cortical cells, 27 neurotypical donors, 180-gene panel, manual layer annotations) provides an independent spatial reference for validating Xenium cell type assignments. Because both platforms measure spatially resolved single cells in human temporal/prefrontal cortex, proportion agreement directly tests whether the correlation classifier produces biologically accurate cell type labels.

### 1.1 Subclass-level concordance

Median subclass proportions correlate at **r = 0.84** between Xenium Control samples (n=12) and MERFISH (n=27). Most subclasses fall near the identity line, with systematic deviations limited to a few types:

- **Oligodendrocyte**: slightly underrepresented in Xenium relative to MERFISH (~15% vs ~25%), likely reflecting the DLPFC vs MTG tissue difference and possibly lower detection efficiency for oligodendrocyte markers in the 300-gene panel
- **Astrocyte**: slightly overrepresented in Xenium (~19% vs ~13%)
- **Microglia-PVM**: overrepresented in Xenium (~8% vs ~4%), consistent with known Xenium detection biases for immune markers

![Subclass proportions: MERFISH vs Xenium](output/presentation/slide_proportion_scatter.png)
*Figure 1: Subclass proportions — Xenium Control median vs SEA-AD MERFISH median. r = 0.84, n = 24 subclasses. Dashed line = perfect agreement.*

The per-layer concordance is also strong. Within each cortical layer (L2/3, L4, L5, L6), Xenium and MERFISH subclass proportions correlate at r = 0.70–0.92, confirming that the depth model and layer assignments produce biologically accurate laminar distributions.

### 1.2 Supertype-level concordance

At the supertype level (131 types), the concordance is necessarily weaker — the 300-gene Xenium panel lacks discriminating markers for many within-subclass distinctions. Overall accuracy against MERFISH ground truth is **84.9%** at the subclass level when the correlation classifier is applied to MERFISH data directly. Most subclasses achieve >80% accuracy (Astrocyte 95%, Oligodendrocyte 96%, L2/3 IT 90%, Pvalb 90%). Notable exceptions include **L6 IT Car3** (69%), **Lamp5** (55%), **Lamp5 Lhx6** (55%), and **Sst** (69%) — types where the 300-gene panel lacks sufficient within-subclass markers to resolve supertypes cleanly.

### 1.3 Depth profile concordance

Median cortical depth per subclass correlates at **r = 0.96** (subclass) and **r = 0.95** (supertype) between Xenium predicted depth and MERFISH manual annotations, confirming the depth model is well-calibrated:

![Median depth comparison](output/presentation/slide_median_depth_by_celltype.png)
*Figure 2: Median depth from pia — Xenium vs MERFISH. Left: subclass level (r = 0.96, n = 23). Right: supertype level (r = 0.95, n = 131). Near-perfect agreement confirms spatially coherent laminar assignments.*

Full depth profile comparisons (proportion vs depth curves for all subclasses) show strong agreement between Xenium Control samples and MERFISH reference, with minor divergences only at the WM boundary. See Figures 9–10 in the [Depth-Stratified Analysis Report](output/depth_proportions/DEPTH_STRATIFIED_ANALYSIS_REPORT.md).

### 1.4 Classifier validation summary

The correlation classifier (step 02b) achieves r = 0.81 against MERFISH proportions (Controls only), compared to r = 0.73 for Harmony-based integration. Critically, Harmony misclassified non-neuronal types into GABAergic categories (e.g., VLMC classified as OPC 82% of the time) and inflated Sst proportions to 12.1% vs the expected 2.5%. See [Methods Writeup](methods_writeup.md) for the full comparison.

---

## 2. Hierarchy of Evidence

| Level | What | Concordance | Confidence |
|-------|------|-------------|------------|
| **Subclass proportions** | Xenium vs MERFISH | r = 0.84 | High — both are spatial platforms measuring same tissue types |
| **Subclass depth** | Xenium vs MERFISH | r = 0.96 | High — near-perfect laminar agreement |
| **Per-layer proportions** | Xenium vs MERFISH | r = 0.70–0.92 | High within layers, lower for sparse types |
| **Supertype proportions** | Xenium vs MERFISH | r = 0.35 (per-donor) | Low — insufficient markers for many supertype distinctions |

---

## 3. Related Documents

| Document | Focus |
|----------|-------|
| [SCZ Compositional Findings](scz_compositional_findings.md) | Disease effects, snRNAseq concordance, confidence tiers |
| [Methods: Cell Typing, Depth Inference, and Validation](methods_writeup.md) | Pipeline methods, classifier validation, depth model, QC calibration |
| [Depth-Stratified Analysis Report](output/depth_proportions/DEPTH_STRATIFIED_ANALYSIS_REPORT.md) | Per-layer and CLR depth × diagnosis results with figures |
| [Supertype Classification Confidence](output/marker_analysis/SUPERTYPE_CLASSIFICATION_CONFIDENCE_REPORT.md) | Per-supertype confidence ratings and Sst fragility analysis |
| [Panel Design & Supertype Classification](output/marker_analysis/XENIUM_PANEL_DESIGN_AND_SUPERTYPE_CLASSIFICATION.md) | Cross-platform marker adequacy, add-on gene recommendations |
