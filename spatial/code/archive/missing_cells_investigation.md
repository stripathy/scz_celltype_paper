# Missing / Unsegmented Cells: Investigation Notes

**Status: Parked** — Background research completed; not actively pursued.

---

## 1. Problem Statement

In Xenium brain tissue data, a large fraction of detected transcripts fall outside segmented cell boundaries. 10x Genomics reports that across 22 human and mouse brain sections, only 32–60% of transcripts are assigned to cells (mean ~48%). The unassigned transcripts arise from three sources:

1. **Truly missing cells** — cell bodies where DAPI was out of the focal plane or too faint for nuclear detection
2. **Cytoplasmic/neuritic transcripts** — RNA in cell processes extending beyond the nucleus expansion zone (5 μm in XOA v2)
3. **Diffused transcripts** — molecules that migrated from source cells onto the glass surface (artifact)

For the SCZ study, systematic under-segmentation could bias cell type proportions if certain morphologies (e.g., small glia, cells with faint DAPI) are preferentially missed.

---

## 2. Our Data

Each sample has 100–282 million transcripts across 300 genes. Per-gene transcript coordinates are already exported as quantized JSON files in `output/viewer/transcripts/{sample_id}/` (0.2 μm resolution, uint16 coordinates). Cell and nucleus boundary polygons are in `data/raw/` as gzipped CSVs.

Key sizes for feasibility:
- Smallest sample (Br2039): ~101M transcripts, 1.0 GB zarr
- Largest sample (Br8667): ~282M transcripts, 2.6 GB zarr
- Mitochondrial genes account for ~67% of all transcripts — filtering these dramatically reduces working set size
- Laptop: 24 GB RAM, 524 GB disk available

---

## 3. Published Tools

### Proseg (Nature Methods, 2025)
- Rust-based cellular Potts model simulation. Initializes from DAPI, expands boundaries to fit transcript spatial patterns.
- Uniquely models **transcript diffusion** (repositions leaked molecules).
- ~3 GB peak memory for 5M transcripts; scales to hundreds of thousands of cells. Full 300M-transcript samples may take hours/days.
- CLI: `proseg --xenium transcripts.parquet`
- GitHub: [dcjones/proseg](https://github.com/dcjones/proseg)

### Baysor (Nature Biotechnology, 2021)
- Julia-based Bayesian/MRF probabilistic model. Best accuracy when combined with Cellpose nuclear prior.
- **Major scalability issue**: ~53 GB peak memory for 5M transcripts; 250M transcripts would need ~5–6 TB. Requires spatial tiling (e.g., via Sopa or Nextflow).
- GitHub: [kharchenkolab/Baysor](https://github.com/kharchenkolab/Baysor)

### FastReseg (Scientific Reports, 2025)
- R-based refinement tool (not de novo segmentation). Corrects existing segmentation boundaries.
- Processed 79M transcripts in ~45 minutes with 45 GB memory — more scalable than Baysor/Proseg at full Xenium scale.

### Segger (bioRxiv, 2025)
- GNN-based (PyTorch Geometric). Frames segmentation as transcript-to-cell link prediction.
- Supports transfer learning across sections. Requires GPU.
- GitHub: [gerstung-lab/segger](https://github.com/gerstung-lab/segger)

### Sopa (Nature Communications, 2024)
- Python pipeline wrapper that tiles large datasets, runs Baysor/Proseg/Cellpose per tile, stitches results.
- Best option for production-scale processing of all 24 samples.
- GitHub: [gustaveroussy/sopa](https://github.com/gustaveroussy/sopa)

### Xenium Ranger (10x Genomics)
- `import-segmentation` module can import Proseg/Baysor results into Xenium-compatible format.
- `resegment` module can re-run with larger expansion distance (up to 100 μm, though >10 μm is mostly noise).

---

## 4. Proposed Approach (If Pursued)

### Phase 1 — Diagnostic (laptop-feasible)
1. For one sample (Br2039), compute what fraction of transcripts are unassigned to any cell boundary
2. Run KDE on unassigned transcript coordinates to find spatial hotspots of "missing" transcripts
3. Compare transcript density maps vs cell density maps — find regions with high transcript density but few cells
4. Correlate hotspot locations with cortical depth and spatial domain
5. **Key question**: Are unassigned transcripts randomly scattered (diffusion noise) or spatially concentrated (missed cells)?

### Phase 2 — Exemplar identification (laptop-feasible)
1. DBSCAN on unassigned transcript coordinates (eps ~5–10 μm, min_samples ~20–50) to find candidate missing-cell clusters
2. Cross-reference with nucleus boundary polygons to confirm clusters are outside all segmented cells
3. Characterize candidate clusters: gene composition, spatial context, proximity to nearest segmented cell
4. Visualize exemplars showing transcript clouds with no overlapping boundaries

### Phase 3 — Resegmentation (likely needs server)
1. Run Proseg on one sample; compare cell counts, types, and transcript capture rates vs default segmentation
2. If successful, scale via Sopa pipeline to all 24 samples
3. Validate: do cell type proportions change systematically? Are certain types preferentially recovered?

### Computational notes
- Phase 1–2 are feasible on laptop if mitochondrial transcripts are filtered (reduces 282M → ~90M transcripts)
- Phase 3 (full Proseg/Baysor) likely needs a machine with 64+ GB RAM for the largest samples, or spatial tiling
- The existing per-gene JSON transcript export and `code/nuclear_resolution/nuclear_counts.py` point-in-polygon infrastructure can be reused for Phase 1–2

---

## 5. Key References

- 10x Genomics KB: [Brain fraction of transcripts within cells is low](https://kb.10xgenomics.com/hc/en-us/articles/29913019227405)
- 10x Genomics KB: [Rescue Xenium transcripts outside cells](https://kb.10xgenomics.com/hc/en-us/articles/29912501365901)
- Optimizing Xenium data utility, Nature Methods 2025 — transcripts >10.7 μm from centroid correlate with background, not cell signal
- Proseg (dcjones/proseg), Nature Methods 2025 — cellular Potts model resegmentation
- Baysor (kharchenkolab/Baysor), Nature Biotechnology 2021 — Bayesian transcript segmentation
- Sopa (gustaveroussy/sopa), Nature Communications 2024 — scalable spatial pipeline
