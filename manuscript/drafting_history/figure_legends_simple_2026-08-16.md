# Figures 1–4: simplified legends (full rewrite, 2026-08-16)

Written for a reader who has the figure and nothing else. Each panel says what a point, bar, or dot *is* before saying what colors and lines mean. Numbers appear only where the reader needs them to read the panel. Conventions: descriptive titles, US spelling, `snRNA-seq`, no em-dashes/colons inside sentences, lowercase panel letters. `[…]` marks the few things only you can fill (see the list at the end).

Companion file with the panel-by-panel visual descriptions and Results edits: `figure_visual_pass_2026-08-16.md`.

---

**Fig. 1 | Study design, harmonized cell types across seven snRNA-seq datasets, and Xenium spatial annotation.**
(**a**) Overview. Seven published single-nucleus RNA sequencing (snRNA-seq) datasets of frontal cortex (469 donors, about 2.4 million nuclei) and one Xenium spatial transcriptomics dataset (24 donors, about 1.2 million cells) were annotated to a shared taxonomy of 24 subclasses and 137 finer supertypes, and each cell type was then tested for changes in gene expression and in abundance in schizophrenia. The bar chart gives the number of cells per platform, with the snRNA-seq bar split by dataset. The horizontal bars give donors per dataset, controls in solid color and schizophrenia cases hatched, with control/schizophrenia counts beside each bar. Far right, the brain region each dataset sampled (all DLPFC except Fröhlich, OFC). (**b**) A two-dimensional map (UMAP) of all snRNA-seq nuclei, in which nuclei with similar gene expression sit near each other. Left, colored and labeled by subclass. Right, the Sst nuclei alone, colored and labeled by their 16 supertypes. (**c**, **d**) The same map colored by dataset (c) and by diagnosis (d). The colors mix within every cluster, so nuclei group by cell type rather than by dataset or diagnosis. (**e**, **f**) One control and one schizophrenia Xenium tissue section, each cell drawn as a dot. In e, dots are colored by the subclass assigned to the cell (key at right). In f, the same cells are colored by the cortical layer assigned to them, from layer 1 at the surface (L1) through L6 to white matter (WM), with blood vessels in gray. DLPFC, dorsolateral prefrontal cortex; OFC, orbitofrontal cortex. Source data are provided as a Source Data file.

*(~250 words. Note the figure currently prints "398 donors"; the legend uses 469 to match the donor bars and the text.)*

---

**Fig. 2 | Cell-type-specific gene expression changes in schizophrenia in snRNA-seq and Xenium.**
Throughout, *SST* and *PVALB* (italic) are genes and Sst and Pvalb are the interneuron types named after them. (**a**–**d**) *SST* expression in Sst cells. (**a**) Every gene tested in Sst cells in the seven-dataset snRNA-seq meta-analysis, plotted by its change in schizophrenia (log2 fold change; negative means lower in schizophrenia) against the strength of evidence (−log10 P). Genes above the dashed line pass FDR < 0.10 and are blue if lower in schizophrenia and orange if higher, in a darker shade when FDR < 0.05; genes below the line are gray. Selected genes are labeled, *SST* in bold. (**b**) The *SST* change in Sst cells in each dataset on its own (gray circles with 95 % confidence intervals[; circle size, weight in the meta-analysis]), pooled across the seven datasets (black diamond), and in the independent Xenium dataset (green triangle). (**c**) *SST* expression in the Sst cells of each Xenium donor, 12 controls versus 12 schizophrenia cases (*SST* transcripts per 1,000 transcripts counted, CP1K). (**d**) One typical Xenium Sst cell from each group, with each red dot marking one detected *SST* transcript (count in red); gray outline, cell boundary; dashed outline, nucleus; scale bar, 5 µm (drawn in h). (**e**–**h**) The same four panels for *PVALB* expression in Pvalb cells. (**i**) Number of genes changed in schizophrenia in each subclass in the meta-analysis, lower in schizophrenia to the left (blue) and higher to the right (orange); dark shade, FDR < 0.05; light shade, additional genes at FDR < 0.10; totals at the bar ends. Inset, the number of changed genes per subclass against how common that subclass is (mean proportion of cells per donor, log axis). (**j**) For every gene and cell type that changed in the meta-analysis and could be measured on the Xenium gene panel (166 pairs), the change in Xenium (y) against the change in the meta-analysis (x). Points near the dashed diagonal changed similarly on both platforms; ten examples are labeled. In the inset of i and in j, points are colored by cell class (green, excitatory neurons; purple, inhibitory neurons; yellow, non-neuronal cells), and in j larger points reached FDR < 0.05 in the meta-analysis. Marks, \*\*\* FDR < 0.01, \*\* FDR < 0.05, \* FDR < 0.10; •, unadjusted P < 0.05; n.s., not significant. FDR, false discovery rate. Source data are provided as a Source Data file.

*(~370 words for ten panels. Trim candidates if needed: the CP1K parenthetical, "ten examples are labeled".)*

---

**Fig. 3 | Cell-type abundance changes in schizophrenia in snRNA-seq and Xenium.**
(**a**) Change in the abundance of each of 109 neuronal supertypes in schizophrenia, pooled across the seven snRNA-seq datasets (298 control and 171 schizophrenia donors). Bars show the fitted change (β) with its standard error, adjusted for donor age, sex, and post-mortem interval. Bars above zero mean the supertype makes up a larger share of neurons in schizophrenia, and bars below zero a smaller share. Bars are colored by subclass. Supertype names in bold red changed at FDR < 0.10 and in italic red at FDR < 0.20; marks, \*\*\* FDR < 0.01, \*\* FDR < 0.05, \* FDR < 0.10, [•] FDR 0.10–0.20. (**b**) The Sst_25 result in detail. The Sst_25 abundance change in each snRNA-seq dataset on its own (gray circles with 95 % confidence intervals), pooled across datasets (black diamond), and in Xenium (green triangle)[; •, unadjusted P < 0.05 within one dataset; the Xenium star denotes …]. (**c**) The same result as raw proportions. Each dot is one donor, plotted by the share of that donor's neurons that are Sst_25, controls (CON, blue) versus schizophrenia (SCZ, red), for each dataset and for Xenium; unadjusted P values above each pair. (**d**) A control and a schizophrenia Xenium section, oriented with layer 1 at the top and white matter (WM) at the bottom. All cells are gray except Sst cells, which are red if they belong to one of the five supertypes depleted in a (Sst_2, Sst_3, Sst_20, Sst_22, Sst_25) and navy otherwise. The curves beside each section show how the two groups are distributed across cortical depth. (**e**) The abundance change of each neuronal supertype in Xenium (y) against its change in the snRNA-seq meta-analysis (x), for the 106 supertypes found in both datasets (purple, inhibitory; green, excitatory). Supertypes that changed in a are enlarged and labeled; dashed line, [identity]; ρ, Spearman correlation. (**f**) For each of the 16 Sst supertypes, its average cortical depth in Xenium (0 = pial surface, 1 = white matter; surface at the top of the plot) against its abundance change in the meta-analysis. Dashed line, linear fit; ρ, Spearman correlation[; horizontal dotted line, …]. Source data are provided as a Source Data file.

*(~330 words.)*

---

**Fig. 4 | Genetic risk, HCN1 physiology, and Alzheimer's disease overlap of the Sst supertypes depleted in schizophrenia.**
Panels a, d, and i show the 16 Sst supertypes, one point each. Thick outlines mark the five supertypes depleted in schizophrenia in Fig. 3a (Sst_2, Sst_3, Sst_20, Sst_22, Sst_25); dashed lines are linear fits with 95 % confidence bands; ρ is the Spearman correlation[; point shading, …]. (**a**) How strongly schizophrenia genetic risk is concentrated in the genes specific to each supertype (x, −log10 P) against how depleted that supertype is in schizophrenia (y, the negative of β from Fig. 3a; higher means more depleted). (**b**) All genes, plotted by how specific their expression is to Sst_3 (x, log scale[; the share of the gene's expression that comes from Sst_3]) against how strongly the gene is associated with schizophrenia in GWAS (y, gene-level −log10 P). Genes right of the vertical dashed line are among the 10 % most Sst_3-specific, genes above the lower dashed line are associated at FDR < 0.05, and genes meeting both are red. Purple dotted line, genome-wide significance. *HCN1* is highlighted. (**c**) Schizophrenia GWAS association of common variants across the *HCN1* gene on chromosome 5 (gene body shaded, exon structure below; purple dashed line, genome-wide significance). Variants that statistical fine-mapping flags as likely causal are enlarged and colored by that probability (PIP). The diamond is the lead variant rs10035564 (PIP 0.56), one of four variants whose probabilities sum to 0.985. (**d**) Average *HCN1* expression in each supertype against the average sag ratio of its cells in an independent human patch-seq dataset. Sag is a voltage response produced by HCN channels (see f). (**e**) One reconstructed patch-seq neuron from each of five supertypes, drawn at the depth its cell body was recorded (black dot); dashed lines, layer boundaries; scale bar, 100 µm. (**f**) Voltage responses of the same five cells to a one-second hyperpolarizing current step. Sag is the rebound from the early minimum (circle) to the later plateau (square), marked on the Sst_25 trace. (**g**) Genes expressed differently between the depleted and the not-depleted Sst supertypes in reference snRNA-seq data from neurotypical donors (SEA-AD). x, log2 fold change (positive means higher in the depleted supertypes); y, −log10 P (values beyond 300 drawn as triangles). Brown, higher in depleted; light orange, higher in not-depleted; dashed lines, fold change of ±0.25 and FDR = 0.05. *HCN1* and *CALB1* are highlighted. (**h**) *CALB1* expression in single cells from the not-depleted versus the depleted Sst supertypes (violin, distribution across cells; box, median and interquartile range). (**i**) Depletion in schizophrenia (x, as in a) against depletion as Alzheimer's disease pathology advances in the SEA-AD DLPFC cohort (y, decline per standard deviation of the pseudo-progression score; higher means more depleted). Source data are provided as a Source Data file.

*(~400 words for nine panels. Trim candidates: the "(see f)" cross-reference, "one of four variants…" clause.)*

---

## Jargon replaced

| Was | Now |
|---|---|
| crumblr β / compositional depletion (crumblr −β) | change in abundance (β) / how depleted the supertype is |
| MAGMA gene-property −log10 P | how strongly schizophrenia genetic risk is concentrated in the genes specific to each supertype |
| SuSiE-R PIP; credible set; cumulative PIP | probability that the variant is causal, from statistical fine-mapping; four variants whose probabilities sum to 0.985 |
| DerSimonian–Laird random-effects; edgeR quasi-likelihood; NB mixed model | pooled across the seven datasets; (test named in Methods only) |
| kernel-density profiles | curves showing how the two groups are distributed across depth |
| library-normalised expression (CP1K) | transcripts per 1,000 transcripts counted |
| direction × FDR tier | blue if lower, orange if higher, darker if FDR < 0.05 |
| SEA-AD neurotypical DLPFC reference (3 donors, 36,601 genes) | reference snRNA-seq data from neurotypical donors (SEA-AD) |
| continuous pseudoprogression score (CPS), −β per SD | decline per standard deviation of the pseudo-progression score, as Alzheimer's pathology advances |
| voltage-sag ratio, HCN-dependent | sag, a voltage response produced by HCN channels; rebound from early minimum to later plateau |
| provenance table; data-driven axes; donor IDs Br6432/Br5973 | removed |
| SST vs Sst unstated | one sentence up front in Fig. 2 defining gene (italic) vs cell type |

## Left for you to fill (`[…]`)

1. Fig 2b / 3b: what circle size encodes (or make circles equal size in the figure).
2. Fig 3a: final trend-level glyph (new figure "•", old legend "+").
3. Fig 3b: what the Xenium star means for a single dataset (unadjusted P, or FDR across the Xenium supertypes).
4. Fig 3e: dashed line, identity or fit.
5. Fig 3f: the horizontal dotted line at depth ≈ 0.53 (define or remove).
6. Fig 3f / Fig 4: the brown point shading (define once, e.g., "colored by supertype", or drop).
7. Fig 4b: exact definition of "specificity" (share of the gene's total expression from Sst_3?).
8. Fig 2i inset: whether the proportion is from snRNA-seq or Xenium (old legend said Xenium).
