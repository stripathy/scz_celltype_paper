# Nuclear Doublet Resolution (Optional — Not Part of the Main Pipeline)

> **This is an optional side investigation**, not part of the numbered pipeline sequence (steps 00-07). The main pipeline runs without this step. Nuclear doublet resolution was found to have negligible impact on downstream biology (see below). It is retained here for scientific interest and reproducibility.

Side investigation into using nuclear-only transcript counts to arbitrate doublet calls in Xenium spatial transcriptomics data.

## Motivation

In Xenium data, cells flagged as doublets based on whole-cell marker co-expression (e.g., Glut+GABA) may actually be single cells with cytoplasmic RNA spillover from adjacent neurons. The nuclear compartment should be cleaner — if a cell's nucleus contains markers from only one type, the whole-cell doublet signal is likely spillover.

This approach was found to be **scientifically interesting but not necessary for downstream analysis**. The simplified pipeline (spatial QC + 5th-percentile margin filter + doublet suspect exclusion) produces equivalent biological results without the rescue machinery. See `docs/pipeline_qc_audit.md` for the empirical comparison.

## How It Works

### Nuclear count building (`nuclear_counts.py`)

1. Load per-cell nucleus boundary polygons from Xenium cell boundaries CSV
2. Build an STRtree spatial index over all nucleus polygons
3. For each transcript, determine which nucleus (if any) it falls within using point-in-polygon queries
4. Aggregate to per-cell nuclear count matrices (nuclear UMI, nuclear gene count, nuclear fraction)

### Doublet resolution (`hybrid_qc.py`, `04_run_nuclear_doublet_resolution.py`)

For each cell flagged as a whole-cell doublet suspect:

1. Check if nucleus has sufficient UMI (>= 50, configurable via `NUCLEAR_MIN_UMI`)
2. Apply the same doublet marker rules to nuclear-only counts
3. Classify the outcome:
   - **resolved** — whole-cell doublet=True, nuclear doublet=False → cytoplasmic spillover, rescue the cell
   - **persistent** — whole-cell doublet=True, nuclear doublet=True → likely real doublet, exclude
   - **insufficient** — nuclear UMI below threshold, can't arbitrate
   - **nuclear_only** — whole-cell clean, nuclear shows doublet markers (excluded in original pipeline)
   - **clean** — not doublet in either assay

### Key findings

- ~76% of whole-cell doublet suspects are resolved as spillover (nuclear compartment is clean)
- The rescue/exclusion has negligible impact on compositional SCZ analysis (crumblr logFC correlation with snRNAseq changes from r=0.405 to r=0.397)
- Top FDR-significant cell types are identical with or without the rescue

## Files

### Scripts
| File | Description |
|------|-------------|
| `04_run_nuclear_doublet_resolution.py` | Main pipeline step: builds nuclear counts, runs resolution, writes `hybrid_qc_pass` |
| `nuclear_counts.py` | Nuclear count matrix builder (point-in-polygon on nucleus boundaries) |
| `hybrid_qc.py` | Doublet resolution logic and `hybrid_qc_pass` computation |
| `plot_nuclear_doublet_validation.py` | Multi-panel validation figures (resolution rates, marker evidence) |
| `plot_doublet_examples.py` | Spatial zoom panels of example doublet/resolved cells |
| `plot_doublet_transcripts.py` | Transcript-level visualization of nuclear vs cytoplasmic markers |
| `archive/build_nuclear_counts.py` | Standalone nuclear count builder (precursor to integrated pipeline step) |
| `archive/compare_nuclear_vs_wholecell.py` | Nuclear vs whole-cell concordance analysis |

### Key figures (in `output/presentation/`)
- `nuclear_doublet_resolution_summary.png` — Resolution category breakdown across all samples
- `nuclear_doublet_marker_evidence.png` — Nuclear marker evidence for resolved vs persistent doublets
- `nuclear_fraction_by_subclass.png` — Nuclear fraction distributions by cell type
- `nuclear_vs_wholecell_concordance.png` — Concordance between nuclear and whole-cell classifications
- `nuclear_hybrid_qc_impact.png` — Impact of hybrid QC on cell retention
- `doublet_spatial_zoom.png` — Spatial examples of doublet suspects and their resolution
- `doublet_transcript_examples.png` — Transcript-level evidence for spillover vs real doublets

### Data outputs (in `output/presentation/`)
- `nuclear_doublet_resolution_all_samples.csv` — Per-cell resolution results
- `nuclear_counting_stats.csv` — Nuclear counting statistics per sample
- `nuclear_vs_wholecell_comparison.csv` — Nuclear vs whole-cell classification comparison

## Running

These scripts expect the main pipeline to have been run through step 03 (transcript export). This optional step is independent of the main pipeline sequence (steps 04-07 do not depend on it). The `hybrid_qc_pass` column remains in the h5ad files for reference but is no longer the default QC gate.

```bash
# Run nuclear doublet resolution on all samples
python3 -u 04_run_nuclear_doublet_resolution.py

# Run on a specific sample
python3 -u 04_run_nuclear_doublet_resolution.py Br8667

# Generate validation figures
python3 plot_nuclear_doublet_validation.py
python3 plot_doublet_examples.py
python3 plot_doublet_transcripts.py
```

## Dependencies

Same as the main pipeline, plus:
- `shapely` (for point-in-polygon nuclear boundary queries; optional, not in `environment.yml` — install with `pip install shapely`)
