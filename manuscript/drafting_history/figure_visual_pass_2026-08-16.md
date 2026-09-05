# Figures 1–4: legend-blind visual pass (2026-08-16)

Process: each figure was inspected **without** reading its legend, the Methods, or any code, and described as a naive reader would see it. The legend was then rewritten to match that reading (minimal jargon), and the Results text was checked sentence-by-sentence against what the figure actually shows.

Sources inspected: the four high-resolution PNGs supplied 2026-08-16, plus the versions embedded in the Google Doc (`1cO5Z…`). **Figs 2 and 4 in the Doc are identical to the high-res files. Figs 1 and 3 in the Doc are older versions** (Fig 1 lacks the subclass key in panel e; Fig 3 has different Xenium numbers, see §3.4).

Conventions applied to the proposed text: descriptive legend titles, US spelling, `snRNA-seq`, no em-dashes/colons inside sentences, lowercase panel letters, ~250–300 words per legend, statistics stated once, color encodings stated once. Anything I could not determine from the figure alone is in `[square brackets]` for you to fill or confirm.

---

## 0. Cross-figure discrepancies to resolve first

These are places where the figure, the legend, and the Results text disagree with each other. I have not tried to decide which is right; they need to be re-derived from the current source outputs and then made identical in all three places.

| # | Where | Figure shows | Legend says | Results text says |
|---|---|---|---|---|
| 1 | Fig 1a flowchart | **398 donors** | — | 469 donors (298 + 171; the donor bars in 1a also sum to 469) |
| 2 | Fig 2b | Multiome estimate **positive**, MSSM 2 ≈ 0 | — | "directionally consistent across all seven datasets" |
| 3 | Fig 2i inset | ρ = **0.84** | — | Spearman ρ = 0.86, P = 1.6e-7, n = 23 |
| 4 | Fig 2i inset | ~23 points | "22 subclasses with ≥ 1 DE gene" | n = 23; "23 of 24 tested subclasses" |
| 5 | Fig 2j | **ρ = 0.75, 73 % concordant** | Pearson r = 0.73; 76 % | Pearson r = 0.73; 76 % (126/166); sign-test P = 7.2e-12 |
| 6 | Fig 3b (new) | Xenium **\*** | (Doc figure had \*\*) | — |
| 7 | Fig 3c (new) | Xenium Sst_25 **p = 0.0083** | — | Sst_25 Xenium P = 0.0024 (Doc figure also had 0.0024) |
| 8 | Fig 3e (new) | ρ = 0.37; L6b_1 ≈ +0.65 | (Doc figure had ρ = 0.4, L6b_1 ≈ +0.8) | ρ = 0.37; L6b_1 = +0.80, P = 0.0068 |
| 9 | Fig 3a (new) | trend marker is "**·**" | "+ FDR 0.10–0.20" (Doc figure used "+") | — |
| 10 | Fig 3b legend | — | "+ FDR 0.10–0.20, ** unadjusted p < 0.05" (missing glyph, should be "•") | — |
| 11 | Results ¶ on Xenium composition | — | — | cites "Fig. 3d" for the Xenium crumblr and "Fig 3f" for the cross-platform ρ = 0.37; both should be **Fig. 3e** |
| 12 | Results Fig 3 ¶ | — | — | "~[1.0]%" and "~[0.75]%" placeholders still present |
| 13 | Results Fig 4 ¶ | — | — | "Supp. Fig. XX", "(CITE …)" placeholders still present |

Items 6–8 mean the *new* Fig 3 was rebuilt from a Xenium composition run that differs from the one the Results paragraph was written against. Every Xenium composition number in that paragraph (Sst_25, Sst_3, Sst_2, L6b_1, L6b_4, and the ρ) needs to be re-extracted from whatever produced the new figure.

Also worth a sentence somewhere: Fig 1a says **137 supertypes**, but Fig 3a tests **109 neuronal** and Supp Fig Sx tests **19 non-neuronal** (= 128). A reader will notice the missing nine.

---

## 1. Figure 1

### 1.1 What the figure shows (legend-blind)

**a.** A three-tier flowchart on the left. Two inputs at the top, "snRNA-seq meta-analysis: seven datasets, 398 donors, 2.4 million cells" and "Spatial transcriptomics: Xenium, 24 donors, 1.2 million cells". Both feed one box, "Harmonized cell types across datasets and technologies: 24 subclasses → 137 supertypes", which feeds two outputs, "Cell-type-specific differential expression" and "Cell-type composition changes". Beside it, a bar chart of number of cells: the snRNA-seq bar (~2.4 M) is stacked in seven colors, one per dataset, with two datasets making up most of it; the Xenium bar (~1.2 M) is one color. To the right, horizontal bars of donors per dataset, largest to smallest (MSSM 2 140/45, HBCC 79/51, Fröhlich 29/32, MSSM 1 17/18, McLean 18/14, Batiuk 10/5, Multiome 5/6), each split into a solid part (control) and a hatched part (schizophrenia), with a separate "Spatial" row for Xenium (12/12). At the far right, a brain icon, a pie of dataset colors labeled DLPFC, and a single dot in the Fröhlich color labeled OFC.
*Take-away:* seven snRNA-seq datasets of very unequal size plus one Xenium dataset, all frontal cortex, all mapped to one cell-type taxonomy.

**b.** A UMAP of all nuclei in which each island is one labeled subclass (glia to the left, excitatory types top-center, interneurons lower-center). The Sst island is boxed, with an arrow to a second UMAP of only the Sst nuclei, colored and labeled by 16 supertypes. Some supertypes form their own islands (Sst_9/10/13, Sst_1/12, Sst_4/5/7); the rest merge into one large connected mass (Sst_2, 3, 11, 19, 20, 22, 23, 25).
*Take-away:* subclasses separate cleanly; Sst supertypes separate only partly.

**c, d.** The same UMAP recolored by dataset (seven colors) and by diagnosis (control blue, schizophrenia red). Both look mostly gray, i.e., thoroughly mixed inside every island; a few glial islands in c show single-dataset patches (an orange Fröhlich patch in astrocytes/oligodendrocytes).
*Take-away:* nuclei group by cell type, not by dataset or diagnosis.

**e, f.** Two Xenium tissue sections (Control, Schizophrenia), each cell a colored dot. In e cells are colored by subclass (24-entry key at right): a fine multicolor pepper with a yellow-green L2/3 IT band near the surface and a gray-green oligodendrocyte mass at the bottom. In f the same cells are colored by layer (L1 red, L2/3 light blue, L4 green, L5 yellow, L6 orange, WM pink, vascular gray) and form concentric ordered bands from pia to white matter.
*Take-away:* every Xenium cell has a subclass and a layer, and the layers look anatomically right.

### 1.2 Things a figure-only reader cannot decode / figure fixes

- "398 donors" contradicts the donor bars in the same panel (sum = 469). Fix in figure.
- Panel c key is titled "Cohort"; the paper's unit is "dataset". Retitle "Dataset".
- The DLPFC pie has equal wedges, so it reads as a proportion chart but only means "these datasets sampled DLPFC". Seven colored dots would say the same thing without implying proportions.
- "1.2 million cells" for Xenium; the Results say 1.34 M profiled and 1.22 M after QC. Say which in the legend (I assumed post-QC).
- Figure text uses "snRNA-seq", which is now the convention for prose too (decided 2026-08-16); the Doc still has ~40 unhyphenated "snRNAseq" occurrences to sweep.

### 1.3 Proposed legend

**Fig. 1 | Study design, cross-dataset snRNA-seq harmonization, and Xenium spatial annotation.**
(**a**) Overview of the study. Left, workflow. Seven snRNA-seq datasets (469 donors, ~2.4 million nuclei) and one Xenium spatial transcriptomics dataset (24 donors, ~1.2 million cells after quality control) were annotated to a shared taxonomy of 24 subclasses and 137 supertypes and analyzed for cell-type-specific differential expression and cell-type composition changes. Middle, number of cells per platform, with the snRNA-seq bar divided by dataset (colors as in c). Right, donors per dataset, split into controls (solid) and schizophrenia cases (hatched), with control/schizophrenia counts at the bar ends. Far right, brain region sampled by each dataset (all DLPFC except Fröhlich, OFC). (**b**) UMAP embedding of all snRNA-seq nuclei, colored and labeled by subclass (left). Right, the Sst nuclei alone, colored and labeled by their 16 supertypes. (**c, d**) The embedding in b colored by dataset (c) and by diagnosis (d); nuclei from different datasets and diagnoses largely overlap within every cluster. (**e, f**) One control and one schizophrenia Xenium section, with each cell colored by its assigned subclass (e; key at right) or by its inferred cortical layer (f; L1 to L6, white matter, and vascular). DLPFC, dorsolateral prefrontal cortex; OFC, orbitofrontal cortex; WM, white matter. Source data are provided as a Source Data file.

### 1.4 Results text alignment (Results ¶1–2)

The two paragraphs already match the figure well. Suggested edits only:

- ¶1: "including 469 donors in total (298 neurotypical controls and 171 SCZ cases; approximately 2.4 million nuclei)" ✓ matches the figure once the flowchart is fixed.
- ¶2: "…comprising 1.34 million cells profiled …" then "After quality control, this yielded 1,221,519 cells…". Fine, but the figure only shows the 1.2 M number. Consider "(1.22 million cells after quality control; Fig. 1a)" so the reader can find the number in the figure.
- ¶2 cites "Fig. 1e, f" for the 742,103 cortical cells; the figure shows sections, not the count. Cite Fig. 1e,f for the sections and leave the count uncited (or cite Table 1 / SM1).

---

## 2. Figure 2

### 2.1 What the figure shows (legend-blind)

The top two rows are the same four panels repeated for two genes: SST in Sst cells (a–d) and PVALB in Pvalb cells (e–h).

**a.** Volcano plot for Sst cells: x is schizophrenia log2 fold change, y is −log10 P. Points above a horizontal dashed line are blue (down) or orange (up), darker for a stricter tier; points below are gray. SST (bold) sits on the down side (~−0.45) just above the line. Also labeled: NAT16 (most significant, down), DRD3, AFG3L2 (down); SMAD1 (most significant, up), SLC9A9, STAC, KCTD4 (up). Blue and orange clouds are of similar size.

**b.** Forest plot of the SST effect per dataset (gray circles of varying size, horizontal 95 % CI bars): five datasets left of zero (Fröhlich, HBCC, McLean, MSSM 1, Batiuk), MSSM 2 essentially at zero with a tight interval, Multiome to the right of zero with a very wide interval. Small dots next to Fröhlich and McLean, a star next to HBCC. Below, a large black diamond (meta-analysis, \*\*) at ~−0.45 and, under a dotted separator, a green triangle (Xenium) at ~−0.3 with no marker.

**c.** Per-donor SST expression in Xenium Sst cells, control (blue, 12) vs SCZ (orange, 12), boxes with points; SCZ median lower and much more spread, with several donors far below any control; p = 0.052.

**d.** Two drawn cells (gray outline, dashed nucleus) with red dots for SST transcripts: 35 in the control cell, 28 in the SCZ cell.

**e–h.** Same layout for PVALB. In e, PVALB (bold) sits at the bottom center of the volcano (essentially unchanged); labeled genes ANXA2, NAT16, VGF (down), SMAD1, SCN3A, TCAF2, CIRBP, FGF10 (up). In f, the seven dataset estimates straddle zero, the diamond is "n.s.", and the Xenium triangle sits at ~−0.2 with a small dot. In g, PVALB in Xenium Pvalb cells is slightly lower in SCZ, p = 0.044. In h, 10 vs 9 transcripts; 5 µm scale bar.

**i.** Horizontal bars of DE-gene counts for 23 subclasses, down (blue, left) and up (orange, right), dark = stricter tier, light = looser, counts at the bar ends. Astro leads (619 down / 324 up), then L5 IT, L2/3 IT, Vip, L4 IT, Micro-PVM, Oligo…, down to L5 ET (2/1). Down exceeds up in most rows; notable exceptions where up wins are L6b (176 vs 80), Vip, L4 IT, Lamp5, Sncg. Inset: DE-gene count vs cell proportion on a log axis rises steadily, ρ = 0.84 (L5 IT, Astro, L6b, Vip labeled); points are green/purple/yellow (plus one gray).

**j.** Scatter of Xenium log2FC (y) against snRNA-seq meta-analysis log2FC (x) for gene–cell-type pairs, colored green/purple/yellow, in two point sizes, with a dashed identity line. Most points fall in the lower-left and upper-right quadrants along the diagonal. Ten ringed, labeled examples: SERPING1 (Astro), FKBP5 (OPC), CALB1 (L6 IT), ATP2B4 (Sst), SMAD1 (Pvalb) up in both; FGFR3 (Astro), SST (Sst), CX3CR1 (Micro-PVM), VGF (Chandelier), BDNF (L2/3 IT) down in both. Text: ρ = 0.75, p < 0.001, 73 % concordant.

*Take-away of the figure:* SST mRNA is down in Sst cells in the meta-analysis and (nearly) in Xenium; PVALB is not down in Pvalb cells in snRNA-seq but is slightly down in Xenium; DE-gene counts scale with subclass abundance; DE effects agree in direction across the two platforms.

### 2.2 Things a figure-only reader cannot decode / figure fixes

- The green/purple/yellow (and one gray) class colors in i-inset and j have no key. Either add a small key or spell out the mapping in the legend (I inferred green = excitatory, purple = inhibitory, yellow = non-neuronal from the labeled points; the gray inset point is unexplained).
- Circle sizes in b/f vary (MSSM 2 largest) with no explanation. Say what size encodes.
- Panel titles "SST / Sst" and "PVALB / Pvalb" are cryptic; "SST in Sst cells" reads instantly.
- The red numbers in d/h ("35", "28") need a word ("SST transcripts") or a legend definition.
- Subclass names in i use underscores (L2_3 IT, L5_6 NP, Lamp5_Lhx6) while Figs 1/3 use "L2/3 IT". Harmonize.
- The statistics printed on j (ρ = 0.75, 73 %) must be synced with the legend/text (Discrepancy #5).
- Legend leftovers to delete: "(provenance table)", "data-driven axes", the two donor IDs, "DerSimonian–Laird", "edgeR quasi-likelihood" (Methods material).

### 2.3 Proposed legend

**Fig. 2 | Cell-type-specific differential expression in schizophrenia across snRNA-seq datasets and Xenium.**
(**a–d**) SST expression in Sst interneurons. (**a**) Differential expression in Sst cells in the seven-dataset snRNA-seq meta-analysis. Each point is a gene, plotted by schizophrenia log2 fold change (x) and −log10 P (y). Genes below the dashed line (FDR ≥ 0.10) are gray; genes above it are blue (down in schizophrenia) or orange (up), with the darker shade marking FDR < 0.05. Selected genes are labeled, SST in bold. (**b**) SST log2 fold change in Sst cells in each snRNA-seq dataset (gray circles, sized by [meta-analytic weight]; bars, 95 % CI), in the pooled random-effects meta-analysis (black diamond), and in the independent Xenium dataset (green triangle). (**c**) Per-donor SST expression in Xenium Sst cells (counts per 1,000 transcripts, CP1K), 12 control versus 12 schizophrenia donors; boxes, median and interquartile range; P from the Xenium pseudobulk test. (**d**) One representative Xenium Sst cell per group, showing the cell outline (gray), nucleus (dashed), and SST transcripts (red dots; count in red). Scale bar, 5 µm (shown in h). (**e–h**) PVALB expression in Pvalb interneurons, shown as in a–d. (**i**) Number of differentially expressed genes per subclass in the meta-analysis, down-regulated (blue, left) and up-regulated (orange, right) at FDR < 0.10 (light) and FDR < 0.05 (dark); counts at the bar ends; the 23 subclasses with at least one DE gene are shown. Inset, DE-gene count per subclass against the subclass's mean proportion of cells per donor [in snRNA-seq / in Xenium, confirm] (log axis); Spearman ρ. (**j**) Xenium log2 fold change against snRNA-seq meta-analysis log2 fold change for the 166 gene × cell-type pairs with meta-analysis FDR < 0.10 that are measurable on the Xenium panel; larger points, meta-analysis FDR < 0.05; dashed line, identity; ten pairs labeled. Points in the inset of i and in j are colored by cell class (green, excitatory; purple, inhibitory; yellow, non-neuronal). Significance markers throughout, \*\*\*FDR < 0.01, \*\*FDR < 0.05, \*FDR < 0.10; •, nominal P < 0.05; n.s., not significant. Source data are provided as a Source Data file.

### 2.4 Results text alignment (section "Somatostatin interneurons show cell-intrinsic reduction…")

Sentence-level suggestions (current text → proposed), with the reason:

1. "The reduction was directionally consistent across all seven datasets and reached nominal significance in three (Fig. 2b)."
   → "SST was reduced in five of the seven datasets and reached nominal significance in three (Fröhlich, HBCC, McLean), while the largest dataset (MSSM 2) and the smallest (Multiome) showed no reduction (Fig. 2b)."
   *Reason:* the forest plot shows Multiome on the positive side and MSSM 2 at zero. **This is a substantive change to the claim; please check against the per-dataset estimates.**

2. "(log₂FC = −0.22, P = 0.044; Fig. 2f,g)" → "Fig. 2f–h" so panel h is cited.

3. "DE-gene count scaled with mean per-donor proportion (Spearman ρ = 0.86, P = 1.6 x 10⁻⁷, n = 23; Fig. 2i inset)" → the figure prints ρ = 0.84. Regenerate and make figure = text (Discrepancy #3/#4).

4. "positively correlated (Pearson r = 0.73) and agreed in direction for 76 % of pairs (126/166; sign-test P = 7.2 x 10⁻¹²; Fig. 2j)" → figure prints Spearman ρ = 0.75 and 73 %. Decide whether the figure reports Pearson or Spearman, regenerate, and make figure = legend = text (Discrepancy #5).

5. "CX3CR1 in Micro-PVM cells (… Xenium −0.34 P = 0.50 x 10-3)" → formatting typo; presumably 5.0 × 10⁻⁴.

6. "across 23 of 24 tested subclasses" → name the subclass without DE genes (it is the one missing from the 23 rows of Fig. 2i) so the reader can reconcile the bar chart with the sentence.

7. The concordant examples listed (VGF/Chandelier, BDNF/L2/3 IT, CX3CR1/Micro-PVM, ATP2B4/Sst, SERPING1/Astro) are all ringed and labeled in j ✓. Values quoted match the point positions to the eye ✓.

8. Volcano genes: text lists AFG3L2, NAT16, SLC9A9, STAC, SMAD1 (Sst) and ANXA2, NAT16, VGF, SCN3A, SMAD1, TCAF2 (Pvalb) ✓ all labeled in a/e. DRD3, KCTD4, CIRBP, FGF10 are labeled in the figure but not mentioned; either mention or unlabel.

---

## 3. Figure 3

### 3.1 What the figure shows (legend-blind)

**a.** A wide bar chart across 109 neuronal supertypes, colored by subclass and mostly (not strictly) grouped by subclass along the x-axis (Lamp5 Lhx6, Lamp5, Pax6, Sncg, Vip, Sst Chodl, Sst, Pvalb, Chandelier, L2/3 IT, L6 IT, L4 IT, L5 IT, L6 CT, L5 ET, L6b, L6 IT Car3, L5/6 NP). y is "SCZ abundance change (β ± SE)", annotated "Increased proportion in SCZ" above zero and "Decreased proportion in SCZ" below. Almost every bar is short with an error bar crossing zero. Two things stand out: a cluster of five Sst bars pointing down (Sst_3, Sst_20, Sst_22, Sst_25, Sst_2; red labels, \*\*/\*\*\*/· marks) and two tall L6b bars pointing up (L6b_1 \*\*\*, L6b_4 \*\*). Four more labels are italic red with a small dot (Pvalb_14, L2/3 IT_7, L6 CT_1, L5/6 NP_4), all pointing up.
*Take-away:* composition changes are rare and concentrated in a few Sst supertypes (down) and L6b (up).

**b.** Forest plot for Sst_25: every dataset estimate is at or left of zero (MSSM 2 ~0, Batiuk most negative but widest), a small dot next to Fröhlich, meta-analysis diamond (\*\*) at ~−0.26, and under a dotted separator a green Xenium triangle (\*) at ~−0.66.

**c.** Boxplots with donor points of Sst_25 as % of all neurons, CON (blue) vs SCZ (red), for the seven datasets and Xenium, with a p-value over each pair (0.13, 0.16, 0.019, 0.22, 0.16, 0.73, 0.75, 0.0083). Controls sit around 1 %; SCZ boxes are lower in every dataset except MSSM 2, but the distributions overlap heavily.
*Take-away of b–c:* the Sst_25 drop is small and noisy in any one dataset and only clear when pooled; Xenium shows it more sharply.

**d.** Two vertical strips of Xenium tissue (Control, Schizophrenia) from L1 at the top to WM at the bottom, layer labels at left. All cells gray; Sst cells red ("Depleted") or navy ("Not depleted"). Red cells crowd L1–L2/3, navy cells L4–L6; the smoothed depth curves beside each strip show the same (red peak superficial, navy peak deep). The SCZ strip has visibly fewer red cells up top.
*Take-away:* the depleted Sst supertypes are the superficial ones.

**e.** Scatter of Xenium abundance change (y) vs meta-analysis abundance change (x) for ~100 neuronal supertypes, purple vs green points, with a dashed line through the origin; the significant/trend supertypes from a are enlarged and labeled. Sst_25, Sst_3, Sst_2 sit lower-left; L6b_1, L6b_4 upper-right; Sst_20 and Sst_22 sit near zero on the Xenium axis. ρ = 0.37, p < 0.001.
*Take-away:* the snRNA-seq composition changes broadly reproduce in Xenium, including the headline Sst and L6b hits, but not Sst_20/Sst_22.

**f.** Scatter of the 16 Sst supertypes, mean cortical depth (y; 0.2 at top to 0.8 at bottom, so superficial is up) vs meta-analysis abundance change (x), points in a brown gradient, dashed fitted line, ρ = 0.76, p < 0.001, plus a horizontal dotted line at ~0.53 and a vertical one at 0. The five depleted supertypes are all in the upper-left (superficial and depleted); Sst_23 and Sst_11 are superficial but not depleted; deep supertypes (Sst_7, 1, 4, 5, 10, 12) cluster near zero change.
*Take-away:* the more superficial the Sst supertype, the more depleted.

### 3.2 Things a figure-only reader cannot decode / figure fixes

- Trend-level marker: the new figure uses "·", the legend says "+". Pick one and match (Discrepancy #9).
- a: the x-axis order is mostly by subclass but a few supertypes sit inside another subclass's block (L2/3 IT_3 and L4 IT_1 among the L5 IT bars; Sst Chodl_2 among the Sst bars), presumably reference-taxonomy order. Either group strictly by subclass or say what the order is. (I counted the bars, 109 ✓.)
- b: what does the Xenium star mean? A single dataset has no meta-analytic FDR; the star must be defined as nominal P or as an FDR across the Xenium supertypes.
- b: circle sizes vary (as in Fig 2b) with no definition.
- e: the dashed line is undefined (identity vs fit; it looks like identity through the origin).
- f: the horizontal dotted line at depth ≈ 0.53 is undefined (define or remove); the reversed y-axis should be flagged in the legend ("superficial at top") or with an arrow.
- f: the brown point shading is undefined here and in Fig 4 (see §4.2).
- Sync the Xenium numbers on the figure with the text (Discrepancies #6–8).

### 3.3 Proposed legend

**Fig. 3 | Neuronal cell-type composition in schizophrenia across snRNA-seq datasets and Xenium.**
(**a**) Change in the abundance of each neuronal supertype in schizophrenia in the seven-dataset snRNA-seq meta-analysis (109 supertypes; 298 control and 171 schizophrenia donors; β ± SE, adjusted for age, sex, and post-mortem interval). Bars above zero indicate a higher proportion in schizophrenia and bars below zero a lower proportion; bars are colored by subclass [and ordered following the reference taxonomy]. Supertype labels in bold red reach FDR < 0.10 and in italic red FDR < 0.20. Markers, \*\*\*FDR < 0.01, \*\*FDR < 0.05, \*FDR < 0.10, [•/+] FDR 0.10–0.20. (**b**) The Sst_25 abundance change in each snRNA-seq dataset (gray circles [sized by …]; bars, 95 % CI), in the meta-analysis (black diamond), and in Xenium (green triangle); markers as in a; •, nominal P < 0.05 within one dataset [confirm the Xenium star definition]. (**c**) Sst_25 as a percentage of all neurons per donor in each snRNA-seq dataset and in Xenium, control (CON, blue) versus schizophrenia (SCZ, red); P values are unadjusted. (**d**) A control and a schizophrenia Xenium section shown from layer 1 (top) to white matter (bottom); all cells gray, with Sst cells colored by whether their supertype is depleted in schizophrenia (red; Sst_2, Sst_3, Sst_20, Sst_22, Sst_25) or not (navy). Curves at the right of each section show the depth distribution of the two groups. (**e**) Xenium abundance change (y) against meta-analysis abundance change (x) for the 106 neuronal supertypes present in both (purple, inhibitory; green, excitatory); supertypes reaching FDR < 0.20 in a are enlarged and labeled; dashed line, [identity]; Spearman ρ. (**f**) Mean cortical depth of each Sst supertype in Xenium (0 = pial surface, 1 = white matter; superficial at top) against its meta-analysis abundance change; dashed line, linear fit; Spearman ρ [; horizontal dotted line, …]. Source data are provided as a Source Data file.

### 3.4 Results text alignment (section "Upper-layer Sst supertypes are less abundant…")

1. Sst_2/22/25/3/20 and L6b_1/L6b_4 effect sizes and FDR tiers ✓ consistent with the marks in a (Sst_2 \*\*\*, Sst_22/25/3 \*\*, Sst_20 trend; L6b_1 \*\*\*, L6b_4 \*\*).

2. "Interestingly, we did not observe abundance changes among the PVALB supertypes."
   → "No Pvalb supertype changed at FDR < 0.10 (Pvalb_14 showed a trend-level increase, FDR < 0.20), and three excitatory supertypes showed trend-level increases (L2/3 IT_7, L6 CT_1, L5/6 NP_4; Fig. 3a)."
   *Reason:* the figure flags Pvalb_14 in italic red with a marker; the sentence as written reads as a contradiction.

3. "…going from an average of ~[1.0]% of neurons in controls to ~[0.75]% in SCZ." → fill from data (panel c shows control medians ≈ 1 % and lower SCZ boxes but the exact means are not readable).

4. "…applying the same crumblr analysis to a 24-section Xenium DLPFC dataset (Fig. 3d)." → **Fig. 3e** (3d is the section image).

5. "Three of the five vulnerable Sst supertypes were significantly depleted in Xenium, including Sst_25 (logFC = −0.66, P = 0.0024; Fig. 3e), Sst_3 (−0.56, P = 0.0019) and Sst_2 (−0.31, P = 0.028). We also saw that L6b supertypes were also more abundant, including L6b_1 (+0.80, P = 0.0068), and L6b_4 (+1.05, P = 0.0022)."
   → All of these must be re-extracted from the run behind the new figure (panel c now prints Sst_25 p = 0.0083 and L6b_1 sits near +0.65 in e). Also "were also … also" (drop one "also").

6. "…showed a good degree of cross-platform agreement (Spearman  = 0.37, P < 0.001, Fig 3f)." → "(Spearman ρ = 0.37, P < 0.001; Fig. 3e)". Missing ρ glyph and wrong panel.

7. "…across all 16 Sst supertypes the more superficial a supertype's mean cortical depth, the more it tended to be depleted in SCZ (ρ = 0.76, P < 0.001; Fig. 3f)." ✓ matches f. Consider adding the visible exceptions in one clause: "(Sst_23 and Sst_11 are superficial but not depleted)". Optional.

8. "Two Sst supertypes failed to replicate, Sst_22 and Sst_20…" ✓ visible in e (both near zero on the Xenium axis).

9. "…the depleted Sst cells were concentrated in superficial cortical layers whereas the not depleted Sst supertypes lay deeper (examples in Fig. 3d)" ✓.

---

## 4. Figure 4

### 4.1 What the figure shows (legend-blind)

All Sst-supertype scatter panels (a, d, i) share a look: 16 labeled points in a brown gradient, five with thick black outlines (Sst_2, 3, 20, 22, 25), a dashed fitted line with a gray band, and a Spearman ρ.

**a.** x = SCZ GWAS enrichment (−log10 P), y = cell depletion in SCZ (−β). ρ = 0.66, p = 0.00685. The outlined five sit at the top; three of them (Sst_2, Sst_3, Sst_20) are also the most enriched (x ≈ 5–6). Two exceptions are obvious: Sst_23 is the most enriched of all (x ≈ 7) but barely depleted, and Sst_22/Sst_25 are strongly depleted with modest enrichment (x ≈ 2–3). Sst_1, 4, 5, 7, 10 are low on both.
*Take-away:* genetically loaded Sst supertypes tend to be the depleted ones, with exceptions on both sides.

**b.** A cloud of all genes: x = how specific the gene is to Sst_3 (log scale), y = −log10 SCZ MAGMA gene p. Genes right of a vertical dashed line (~0.003) and above a lower dashed line are red; everything else gray. HCN1 (bold, large brown dot) sits at y ≈ 11, among the strongest associations in the red set; a few red genes (LACC1, EYS, GSX2) cross the purple genome-wide line at ~7.3; RBFA, ZNF664, TMED6, SLC35F6, DLX2, CDH13, OPRD1 are also labeled.
*Take-away:* HCN1 is both Sst_3-specific and strongly SCZ-associated.

**c.** A regional GWAS plot across the HCN1 locus on chr5 (~45.75 to ~45.05 Mb, hg38, coordinates decreasing left to right), gene body shaded, gene model with exons and arrows below. Gray variants form the floor, pink variants sit above the purple genome-wide line across the gene body, and four large markers (a red diamond and three circles colored yellow→red by PIP) sit at the top right, at one end of the gene and just beyond it.
*Take-away:* the locus is genome-wide significant and the fine-mapped credible set falls at HCN1.

**d.** x = HCN1 expression (log2 CP10K+1), y = sag ratio (Patch-seq), 16 supertypes; ρ = 0.65, p = 0.00786. Sst_25 is highest on both; the outlined five are mostly upper-right; Sst_12 is a low-sag outlier at high HCN1.

**e.** Five reconstructed neurons placed by soma depth against pia/L1/L2/L3 lines: Sst_25, Sst_22, Sst_3 ("Depleted in SCZ") in L2–L3, the first two with narrow, vertically bundled dendrites; Sst_5 and Sst_1 ("Not depleted") deeper, with spreading multipolar arbors. A "100 µm" scale bar that renders as a tiny square.

**f.** Voltage traces of the same five cells during a 1-s hyperpolarizing step: Sst_25 has the biggest sag (circle at the early minimum, square at steady state, arrow labeled "Sag") and a rebound after the step; Sst_3/Sst_22 modest sag; Sst_5 and Sst_1 hyperpolarize furthest with almost no sag.
*Take-away of d–f:* depleted supertypes are superficial, HCN1-high, and show sag; non-depleted are deep and flat.

**g.** Volcano of depleted vs not-depleted Sst supertypes: brown points (higher in depleted) to the right of +0.25, light-orange points (higher in not-depleted) to the left of −0.25, gray between; y capped at 300 (triangles); "FDR = 0.05" line near the bottom. CALB1 and HCN1 in bold on the depleted side; SLIT2, TRHDE, ZMAT4, GALNTL6 and NELL2, PCDH9, ADAMTS9-AS2, ZNF536 also labeled.

**h.** Violins of per-cell CALB1: not-depleted cells are almost all at zero; depleted cells spread up to ~2.3 with a box; FDR < 0.001.

**i.** x = cell depletion in SCZ (−β), y = cell depletion in AD (−β per SD CPS); ρ = 0.87, p < 10⁻⁴. The outlined five are all upper-right (depleted in both); Sst_23 and Sst_11 are depleted in AD but only slightly in SCZ; the deep supertypes are lower-left; Sst_13 is the one slightly discordant point (15/16 same-sign).
*Take-away:* the SCZ-depleted supertypes are the same ones lost with AD progression.

### 4.2 Things a figure-only reader cannot decode / figure fixes

- **The brown gradient on points, traces, and reconstructions is undefined** (a, d, e, f, i, and Fig 3f). It is not depletion (Sst_23 is dark but not depleted; Sst_2 is light but depleted) and not depth (Sst_2 is superficial but light). If it is a fixed per-supertype palette, say so once ("colored by supertype, palette shared with Fig. 3f"); if it encodes a variable, add a key; otherwise drop it, since the thick outline already carries the grouping.
- e: the 100 µm scale bar renders as a small square; the lowest dashed layer line is unlabeled (L3/L4?).
- c: the diamond is the only unlabeled key element; a tiny "rs10035564" next to it would let the panel stand alone. "SuSiE-R PIP" in the key is jargon; "PIP (fine-mapping)" plus a legend definition is enough.
- b: the axis title "Gene specificity in Sst_3" needs a one-phrase definition in the legend.
- g: state the comparison universe in the legend (which cells, how many).
- Panel letters are now lowercase ✓ (previously uppercase per memory).

### 4.3 Proposed legend

**Fig. 4 | Genetic risk, HCN1-linked physiology, and cross-disorder vulnerability of the depleted upper-layer Sst supertypes.**
Panels a, d, and i plot the 16 Sst supertypes; thick outlines mark the five depleted in schizophrenia (Sst_2, Sst_3, Sst_20, Sst_22, Sst_25; FDR < 0.20 in Fig. 3a); dashed lines are linear fits with 95 % confidence bands; ρ, Spearman correlation [; points colored by supertype as in Fig. 3f]. (**a**) Enrichment of schizophrenia common-variant risk among the genes specific to each supertype (−log10 P, x) against that supertype's depletion in schizophrenia (−β from Fig. 3a, y). (**b**) All genes plotted by how specific their expression is to Sst_3 (x, log scale, [fraction of total expression contributed by Sst_3]) and by their gene-level schizophrenia GWAS association (y). Genes right of the vertical dashed line (top 10 % most Sst_3-specific) and above the lower dashed line (FDR < 0.05) are red; purple dotted line, genome-wide significance; HCN1 highlighted. (**c**) Schizophrenia GWAS association for variants across the HCN1 locus (chr5, hg38; gene body shaded, gene model below; purple dashed line, genome-wide significance). Variants in the fine-mapped credible set are enlarged and colored by their posterior probability of being causal (PIP); diamond, lead variant rs10035564 (PIP = 0.56; four-variant credible set, cumulative PIP = 0.985). (**d**) Mean HCN1 expression per supertype against mean sag ratio measured by patch-seq in an independent set of human Sst interneurons. (**e**) Reconstructed morphology of one patch-seq cell from each of five supertypes, positioned at its recorded soma depth (black dot); dashed lines, layer boundaries; scale bar, 100 µm. (**f**) Voltage responses of the cells in e to a hyperpolarizing current step (0–1,000 ms); sag is the difference between the early minimum (circle) and the later steady state (square), marked on the Sst_25 trace. (**g**) Genes differentially expressed between depleted and not-depleted Sst supertypes in the SEA-AD reference snRNA-seq data (positive, higher in depleted); brown, higher in depleted; light orange, higher in not-depleted; dashed lines, |log2 fold change| = 0.25 and FDR = 0.05; triangles, P values beyond the axis; HCN1 and CALB1 highlighted. (**h**) CALB1 expression per cell in not-depleted versus depleted Sst supertypes (violins with median and interquartile box); FDR from a cell-level Wilcoxon test. (**i**) Depletion in schizophrenia (x, as in a) against depletion along the Alzheimer's disease continuum in the SEA-AD DLPFC data (y, −β per SD of the continuous pseudo-progression score, CPS). Source data are provided as a Source Data file.

### 4.4 Results text alignment (section "Schizophrenia genetic risk converges…")

1. "(Fig. 4a; Spearman ρ = 0.66, p = 0.0069 …)" ✓ (figure 0.00685). Suggest one honest clause on the visible exceptions after this sentence: "The correspondence was not perfect. Sst_23 carried the strongest enrichment but was not depleted, and Sst_22 and Sst_25 were depleted despite modest enrichment (Fig. 4a)."

2. "HCN1 … is specifically expressed in Sst_3, carries a strong gene-level SCZ association (MAGMA −log₁₀p = 10.9)" ✓ matches b (~10.9). Consider "is among the 10 % of genes most specific to Sst_3" so the reader sees why the vertical line is where it is.

3. "…lies under a fine-mapped SCZ locus that nominates HCN1 as the likely causal gene (Fig. 4c)" ✓.

4. "(Fig. 4d; ρ = 0.65, p = 0.008, example voltage traces in Fig. 4f)" ✓ (figure 0.00786).

5. "…identified 580 genes differentially expressed … (FDR < 0.05, |log₂FC| > 0.25; Fig. 4g)" ✓ thresholds match the dashed lines. "HCN1 … (log₂FC = +0.43)" and "CALB1 (log₂FC = +0.55; FDR < 0.001, Fig. 4h)" ✓ consistent with point positions to the eye.

6. "…the reconstructed exemplars of depleted and not-depleted supertypes differ systematically in both laminar position and voltage sag (Fig. 4e,f)" ✓; consider naming what is visible: "the depleted exemplars (Sst_25, Sst_22, Sst_3) sit in L2–L3 with narrow vertically oriented dendritic bundles and show pronounced sag, whereas the not-depleted exemplars (Sst_5, Sst_1) sit deeper, are multipolar, and show little sag".

7. "…agreeing in direction for 94% of supertypes … (Spearman ρ = 0.87; n = 16 supertypes; Fig. 4i …)" ✓ (15/16 in the figure; Sst_13 is the exception, could be named).

8. Placeholders remain: "(CITE)", "(CITE Bigdeli, Trubetskoy)", "Supp. Fig. XX" (twice).

---

## 5. Suggested order of operations

1. Regenerate and freeze the numbers behind Discrepancies #2–#8 (Fig 2b direction claim; Fig 2i inset ρ; Fig 2j r/ρ and %; all Xenium composition numbers in the Fig 3 paragraph).
2. Apply the figure fixes (§1.2, §2.2, §3.2, §4.2), then re-export to the Doc so the embedded Figs 1 and 3 match the high-res files.
3. Paste the proposed legends, filling the `[…]` items.
4. Apply the sentence-level Results edits.
