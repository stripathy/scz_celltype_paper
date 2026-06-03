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

### 5. Extract annotations and overlay on raw images
To produce the final publication figures, annotations from the Inkscape-edited
SVG are detected by rendering the SVG to a PNG and finding the colored markers
visually (bypasses matplotlib SVG coordinate-math gotchas). The markers are
then plotted on the **raw** (non-cleaned) composites, so the reader sees the
original microscopy data with the annotator's cell calls overlaid.

## Files

### Scripts

| File | Description |
|------|-------------|
| `lipofuscin_removal.py` | All processing: loads raw images, applies cleaning, reads SVG annotations, regenerates all figures |

Run as: `python lipofuscin_removal.py` (from the repo root, or from this
directory). Will regenerate all `.png`/`.svg` outputs below.

### User-editable annotation SVG

| File | Description |
|------|-------------|
| `RNAscope_fig_representative_inkscape_SJT.svg` | **User's manual annotations.** Open in Inkscape to add/move/delete star markers. Regenerate the final figures after editing by running `lipofuscin_removal.py`. |

### Data

| File | Description |
|------|-------------|
| `marker_coordinates.csv` | Exported coordinates of every manual annotation marker (panel, type, pixel x/y). Regenerated automatically by the script. |

### Figures

| File | Description |
|------|-------------|
| `fig_original_vs_cleaned_L23.png` | Side-by-side comparison showing the visual effect of lipofuscin suppression on two L2/3 panels |
| `fig_L23_manual_stars_on_raw.png/.svg` | **Main figure:** L2/3 Control vs SCHIZ with manual annotations on raw (non-cleaned) composites |
| `fig_all_manual_stars_on_raw.png/.svg` | Full 4-panel version (includes L5/6 panels) |

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
