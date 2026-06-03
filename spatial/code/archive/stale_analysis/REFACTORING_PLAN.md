# Refactoring Plan: Probe Design & Cross-Platform Analysis Codebase

## Current State

Seven scripts in `code/analysis/` handle probe design, cross-platform gene validation, and panel comparison. They were written incrementally over ~2 months and have significant code duplication, inconsistent import patterns, and no shared utilities beyond `config.py` (which most of them don't use).

### Scripts in scope (execution order)

| # | Script | Lines | Purpose |
|---|--------|-------|---------|
| 1 | `merscope_panel_assessment.py` | 511 | Load MERSCOPE 4K data, calibrate detection efficiency, compare panels |
| 2 | `compare_snrnaseq_merscope_expression.py` | 595 | Pseudobulk snRNAseq vs MERSCOPE 4K correlations |
| 3 | `characterize_poor_genes.py` | 685 | Gene properties predicting poor cross-platform correlation |
| 4 | `cross_platform_gene_corr.py` | 668 | Extend gene correlations to MERFISH + Xenium platforms |
| 5 | `supertype_markers_panel_overlap.py` | ~300 | Wilcoxon 1-vs-rest at supertype level, panel coverage |
| 6 | `nsforest_supertype_markers.py` | ~250 | NS-Forest combinatorial marker discovery |
| 7 | `hierarchical_probe_selection.py` | 1219 | Final probe panel design (hierarchical markers + dropout scoring) |
| 8 | `visualize_marker_quality.py` | ~370 | Dot plots / matrix plots of selected markers (depth-ordered) |

### Problems

1. **Duplicated constants**: `SUBCLASS_TO_CLASS` defined in 4 separate files (constants.py, compare_snrnaseq_merscope_expression.py, cross_platform_gene_corr.py, characterize_poor_genes.py)
2. **Duplicated functions**: `compute_pseudobulk_mean()` in 2 files; gene biotype classification logic in 3 files; panel gene loading in 4 files
3. **Inconsistent imports**: Scripts 1-6 hardcode all paths; only ~35 visualization scripts use `config.py`
4. **No shared utilities**: Each script re-implements reference loading, normalization, subsampling
5. **Flat directory**: 48 scripts in one folder with no subdirectories or logical grouping
6. **Output scattered**: Results split between `output/marker_analysis/` and `output/presentation/`

---

## Proposed Structure

```
code/
├── modules/
│   ├── constants.py              (existing - unchanged)
│   ├── panel_utils.py            (NEW - panel loading, gene quality, detection efficiency)
│   ├── reference_utils.py        (NEW - snRNAseq reference loading, normalization, subsampling)
│   ├── pseudobulk.py             (NEW - pseudobulk computation, used by 3+ scripts)
│   ├── gene_properties.py        (NEW - biotype classification, gene quality filtering)
│   └── ...                       (existing modules unchanged)
│
├── analysis/
│   ├── config.py                 (existing - add panel paths + new module imports)
│   │
│   ├── # ── Cross-platform validation (run in order) ──
│   ├── compare_snrnaseq_merscope_expression.py   (refactored: use shared modules)
│   ├── merscope_panel_assessment.py               (refactored: use shared modules)
│   ├── characterize_poor_genes.py                 (refactored: use shared modules)
│   ├── cross_platform_gene_corr.py                (refactored: use shared modules)
│   │
│   ├── # ── Probe design ──
│   ├── supertype_markers_panel_overlap.py          (refactored: use shared modules)
│   ├── nsforest_supertype_markers.py               (refactored: use shared modules)
│   ├── hierarchical_probe_selection.py             (refactored: use shared modules)
│   ├── visualize_marker_quality.py                 (NEW: dot plots / matrix plots, uses shared modules)
│   │
│   ├── # ── Visualization scripts (unchanged) ──
│   ├── plot_*.py                  (30+ scripts, no changes needed)
│   └── ...
```

---

## Step-by-Step Refactoring

### Step 1: Create `code/modules/panel_utils.py`

Extract panel-related utilities used across 4+ scripts:

```python
"""Panel gene list loading and detection efficiency lookup."""

def load_xenium_panels():
    """Load Xenium 5K Prime and v1 Brain gene lists.
    Returns: dict of {'xenium_5k': set, 'xenium_v1': set}
    """

def load_detection_efficiency(path):
    """Load snRNAseq vs MERSCOPE 4K detection efficiency table.
    Returns: dict of gene -> efficiency
    """

def load_gene_quality(path):
    """Load gene properties vs correlation table.
    Returns: DataFrame with gene quality metrics
    """

def load_spatial_validation(xenium_corr_path, merfish_corr_path, r_threshold=0.7):
    """Load cross-platform gene correlations for spatial validation.
    Returns: set of validated genes, dict of gene -> best_r
    """
```

### Step 2: Create `code/modules/reference_utils.py`

Extract reference data loading pattern used in 5+ scripts:

```python
"""snRNAseq reference loading, normalization, and subsampling."""

def load_snrnaseq_reference(path, normalize=True, min_cells=10):
    """Load SEA-AD snRNAseq reference with optional normalization.
    Applies: normalize_total(1e4), log1p, filter_genes(min_cells).
    """

def subsample_by_group(adata, groupby, max_cells=500, min_cells=20, seed=42):
    """Subsample to max_cells per group, dropping groups below min_cells."""
```

**Note**: `config.py` already has `load_snrnaseq_reference()` but it does NOT normalize. The new version should have a `normalize` flag to handle both use cases.

### Step 3: Create `code/modules/pseudobulk.py`

Extract pseudobulk computation from `compare_snrnaseq_merscope_expression.py` and `cross_platform_gene_corr.py`:

```python
"""Pseudobulk computation utilities."""

def compute_pseudobulk_mean(adata, groupby_col, gene_subset=None):
    """Compute mean expression per group (pseudobulk).
    Returns: DataFrame (groups x genes)
    """
```

### Step 4: Create `code/modules/gene_properties.py`

Extract gene biotype classification from `characterize_poor_genes.py`:

```python
"""Gene property classification utilities."""

def classify_gene_biotype(gene_name):
    """Classify gene biotype based on naming conventions.
    Returns: str (protein_coding, lncRNA, antisense, etc.)
    """

def filter_eligible_genes(gene_quality_df, min_pearson_r_quantile=0.2,
                          exclude_biotypes=None, min_detection_rate=0.005):
    """Build set of eligible genes for probe selection."""
```

### Step 5: Update `config.py` with panel paths

Add panel paths and new module imports to `config.py`:

```python
# Panel gene lists
PANEL_5K_PATH = os.path.expanduser(
    "~/Downloads/XeniumPrimeHuman5Kpan_tissue_pathways_metadata.csv")
PANEL_V1_PATH = os.path.expanduser(
    "~/Downloads/Xenium_hBrain_v1_metadata.csv")

# Cross-platform outputs
MARKER_ANALYSIS_DIR = os.path.join(BASE_DIR, "output", "marker_analysis")
```

### Step 6: Refactor each analysis script

For each of the 7 scripts:
1. Replace hardcoded paths with imports from `config.py`
2. Replace duplicated `SUBCLASS_TO_CLASS` with import from `constants.py`
3. Replace inline panel loading with `panel_utils.load_xenium_panels()`
4. Replace inline reference loading with `reference_utils.load_snrnaseq_reference()`
5. Replace inline pseudobulk with `pseudobulk.compute_pseudobulk_mean()`
6. Replace inline biotype classification with `gene_properties.classify_gene_biotype()`

**Priority order** (by duplication severity):
1. `compare_snrnaseq_merscope_expression.py` — has SUBCLASS_TO_CLASS + pseudobulk + reference loading
2. `cross_platform_gene_corr.py` — has SUBCLASS_TO_CLASS + pseudobulk + reference loading
3. `characterize_poor_genes.py` — has SUBCLASS_TO_CLASS + biotype classifier
4. `supertype_markers_panel_overlap.py` — has panel loading + reference loading
5. `merscope_panel_assessment.py` — has panel loading + reference loading
6. `nsforest_supertype_markers.py` — has panel loading + reference loading
7. `hierarchical_probe_selection.py` — most complex but already well-structured; just needs shared imports
8. `visualize_marker_quality.py` — NEW script, built with shared modules from the start (reference_utils)

### Step 7: Add execution order documentation

Add a brief comment block at the top of each script indicating its position in the pipeline:

```python
# Pipeline position: 1 of 7 (cross-platform validation)
# Upstream: none
# Downstream: characterize_poor_genes.py (reads snrnaseq_vs_merscope_gene_corr.csv)
```

---

## What NOT to refactor

1. **Don't create a monolithic pipeline script** — each script should remain independently runnable. The PI runs them manually and iterates on results.
2. **Don't move visualization scripts** — the 30+ `plot_*.py` scripts work fine in the flat structure and already use `config.py`.
3. **Don't restructure output directories** — downstream analyses and the PI reference these paths. Just document the data flow.
4. **Don't merge scripts** — `hierarchical_probe_selection.py` (1219 lines) is already at the upper limit. Merging scripts would make them harder to iterate on.

---

## Estimated effort

| Step | New lines | Files touched | Risk |
|------|-----------|---------------|------|
| 1. panel_utils.py | ~80 | 1 new | Low |
| 2. reference_utils.py | ~60 | 1 new | Low |
| 3. pseudobulk.py | ~40 | 1 new | Low |
| 4. gene_properties.py | ~100 | 1 new | Low |
| 5. config.py paths | ~10 | 1 existing | Low |
| 6. Refactor 7 scripts | ~-200 net | 7 existing | Medium (test each after) |
| 7. visualize_marker_quality.py | ~370 | 1 new | Low (uses shared modules) |
| 8. Documentation | ~50 | 8 existing | Low |

**Total**: ~4 new module files + 1 new analysis script, ~710 new lines, ~200 lines removed from duplication, 8 files modified.

**Verification**: After each script refactoring, re-run it and diff the outputs against current versions to ensure identical results.
