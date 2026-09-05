# Representative Microscopy Images — Lipofuscin Suppression

This directory contains the pipeline and outputs for generating representative
RNAscope microscopy images for the sgACC SST/VIP interneuron analysis.

## The problem

Postmortem human brain tissue contains **lipofuscin**, an autofluorescent
pigment that accumulates in aged neurons. Lipofuscin is a well-known nuisance
for fluorescence microscopy of aged tissue because it is **broadband** — it
fluoresces across multiple excitation/emission wavelengths simultaneously,
contaminating the 488 nm (SST) and 568 nm (VIP) channels in our RNAscope FISH
data. In raw composite images, lipofuscin appears as bright yellow/white
granular clumps that visually dominate the true RNAscope puncta, making it
difficult for a reader to distinguish genuine signal from artifact.

This is a universal issue for postmortem human cortex FISH, and essentially
every major publication in the field addresses it in some form. See
"Literature context" below.

## The strategy

We use **channel-based spectral subtraction**: because lipofuscin emits in
both the 488 nm and 568 nm channels, we estimate the lipofuscin component
of each pixel as the minimum of the red and green channel intensities, and
subtract it from each channel.

```
lipo_estimate = min(red, green)                      # broadband = bright in both
SST_clean     = green - 1.2 × lipo - 0.15 × red     # retain green-specific signal
VIP_clean     = red   - 1.2 × lipo - 0.15 × green   # retain red-specific signal
```

The **LIPO_FACTOR = 1.2** slightly over-subtracts the minimum to account for
cases where lipofuscin's intensity differs slightly between the two channels.
The **CROSS_CH_FACTOR = 0.15** adds a small cross-channel penalty to suppress
residual bleedthrough. Any pixels below **threshold = 0.03** after subtraction
are zeroed out to remove dim noise.

Note that this processing is applied **only for display** (representative
image figures). All quantitative analyses use the original cell counts from
Dwight F. Newton's manual SlideBook segmentation of the uncleaned images.

## Literature context

Computational lipofuscin suppression for representative images is standard
practice in postmortem human brain FISH publications. Common approaches:

1. **Spectral unmixing using a dedicated lipofuscin channel** (the gold
   standard, used by the Lewis lab). A separate far-red or empty channel
   is imaged specifically to capture lipofuscin, and subtracted.
2. **Channel subtraction of broadband autofluorescence** (what we do). Used
   in Sibille lab work and others when a dedicated lipofuscin channel is
   not available.
3. **Chemical quenching** (Sudan Black B, TrueBlack, CuSO₄/NH₄OAc). Applied
   at the bench during tissue prep; not computational.
4. **Display-only manipulations** (brightness/contrast adjustment).
   Acceptable when applied uniformly across all conditions and disclosed.

References (Lewis lab and Sibille lab examples):
- Fish, Rocco, Lewis et al. (2018) *Cerebral Cortex* — SST/PV in schizophrenia DLPFC
- Hoftman et al. (2018) *Biological Psychiatry*
- Joshi et al. (2015) — lipofuscin channel masking for RNAscope quantification
- Arbabi et al. (2025) *Molecular Psychiatry* — original source of this dataset

Our approach is (2). It is acceptable because it is (i) applied identically
to all representative images shown, (ii) disclosed in the methods, and
(iii) does not affect any quantitative analyses.

## Workflow

### 1. Load raw images
Raw RNAscope composite TIFFs (1024×1024 RGB, 8-bit) from
`data/dwight microscopy data representative images/`. These are named
`{subject}-R-{site}.tif` for R-section images (SST + VIP staining).

### 2. Crop to the 800×800 counting frame
Per the stereology protocol (`documentation/Protocol-Cell densitometry*.docx`),
cells were quantified within an 800×800 px inset from the full 1024×1024
image (frame at [112:912, 112:912]). We crop to this same region for the
figure images so that the displayed region matches what was counted.

### 3. Apply lipofuscin suppression (see "strategy" above)

### 4. Manually annotate representative cells in Inkscape
We use a hybrid approach:
- An initial algorithmic detection pass (based on the cleaned channels)
  places stars near DAPI-centered clusters of bright puncta.
- The user then opens the SVG in Inkscape and manually adjusts each star
  to sit on visually convincing cells, adding/removing markers as needed
  based on their own judgment.

The annotation is **illustrative**: stars mark *examples* of labeled cells,
not an exhaustive count. Reported cell counts come from Dwight's original
quantification, not from the annotated stars.

### 5. Export the annotations and build the manuscript panel
`marker_coordinates.csv` holds the coordinates of every marker in the curated
SVG (panel, type, pixel x/y in the 800-px counting frame). The manuscript
figure panel is built outside this directory: `code/extract_micrograph_panels.py`
pulls the four cleaned composites embedded in the SVG into `panels/`, and
`code/plot_figS_rnascope.R` draws them with the markers, labels, scale bar and
an on-image channel key (Fig. S panel b). The stars are **illustrative**:
they mark examples of labeled cells, not an exhaustive count; reported cell
counts come from the original quantification.

## Files

| File | Description |
|------|-------------|
| `lipofuscin_removal.py` | Loads raw composites (needs the raw TIFFs, not in the repo), applies the display-only lipofuscin suppression, and exports the SVG marker coordinates. Its legacy figure functions (stars-on-raw renders) are superseded by the R panel above. |
| `RNAscope_fig_representative_inkscape_SJT.svg` | Curated annotation SVG (edit in Inkscape); embeds the four cleaned composites. |
| `marker_coordinates.csv` | Exported marker coordinates (panel 0-3 = Control L2/3, SCZ L2/3, Control L5/6, SCZ L5/6). |
| `panels/panel{0..3}_*.png` | The four cleaned composites extracted from the SVG, display orientation (input to the R panel). |
| `fig_original_vs_cleaned_L23.png` | Original vs lipofuscin-suppressed comparison for two L2/3 frames (documents the display processing). |

## Methods-section template

Suggested text for the manuscript methods section:

> *"For representative image figures, RGB composite TIFFs (488 nm: SST;
> 568 nm: VIP; 405 nm: DAPI) were cropped to the 800×800 px counting frame
> used by the original stereological quantification (Arbabi et al. 2025),
> normalized per channel by 1st–99.5th percentile stretch (DAPI) or
> fixed-range (SST, VIP), and processed to suppress lipofuscin
> autofluorescence by spectral subtraction. Lipofuscin was estimated per
> pixel as the minimum of the SST (488 nm) and VIP (568 nm) channels and
> subtracted from each signal channel (multiplicative factor 1.2), with
> an additional cross-channel correction of 0.15 and an intensity floor
> of 0.03. This processing was applied uniformly to all displayed images
> and was used for visualization only; all quantitative cell density
> analyses reported in Results used the original SlideBook mask-based
> cell counts from Arbabi et al. (2025). Representative cells indicated
> by markers were identified by visual inspection; marker counts are
> illustrative and may differ from the original quantitative counts."*

## Parameter tuning

The default parameters were selected to balance lipofuscin suppression
against preservation of genuine RNAscope puncta:

- `LIPO_FACTOR = 1.2` — slight over-subtraction handles minor channel
  imbalance in lipofuscin emission.
- `CROSS_CH_FACTOR = 0.15` — catches residual bleedthrough not captured by
  the min() estimate.
- `THRESHOLD = 0.03` — removes dim pixelwise noise in the cleaned channel.

To adjust (e.g., for a different imaging setup), edit the constants near the
top of `lipofuscin_removal.py` or pass custom values to `load_clean()`.
