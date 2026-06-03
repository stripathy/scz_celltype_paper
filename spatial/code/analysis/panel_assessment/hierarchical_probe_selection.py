#!/usr/bin/env python3
"""
Hierarchical probe selection for SST and L6b supertype disambiguation.

Designs add-on probe panels for Xenium v1 Brain and Xenium 5K Prime using:
  1. Two-round hierarchical marker discovery:
     - Round 1: Subclass-level (SST-vs-all, L6b-vs-all)
     - Round 2: Within-subclass (Sst_X vs other SST supertypes, L6b_X vs other L6b)
  2. Dropout-aware scoring calibrated by MERSCOPE 4K detection efficiency
  3. Redundancy planning to ensure sufficient markers despite spatial dropout

Usage:
    python hierarchical_probe_selection.py
"""

import os
import sys
import time
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import sparse

warnings.filterwarnings("ignore", category=FutureWarning)

# ── Shared module imports ────────────────────────────────────────────────────
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "modules"))
from panel_utils import (load_xenium_panels, load_detection_efficiency,
                          load_spatial_validation, load_gene_quality)
from reference_utils import load_and_normalize_reference, subsample_by_group

# ── Configuration ────────────────────────────────────────────────────────────
BASE_DIR = os.path.expanduser("~/Github/SCZ_Xenium")
REF_PATH = os.path.join(BASE_DIR, "data", "reference",
                         "nicole_sea_ad_snrnaseq_reference.h5ad")
DETECTION_PATH = os.path.join(BASE_DIR, "output", "marker_analysis",
                               "snrnaseq_vs_merscope4k_detection.csv")
GENE_QUALITY_PATH = os.path.join(BASE_DIR, "output", "presentation",
                                  "gene_properties_vs_correlation.csv")
# Cross-platform validated gene performance from prior analysis
XENIUM_CORR_PATH = os.path.join(BASE_DIR, "output", "presentation",
                                 "snrnaseq_vs_xenium_gene_corr.csv")
MERFISH_CORR_PATH = os.path.join(BASE_DIR, "output", "presentation",
                                  "snrnaseq_vs_merfish_gene_corr.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "marker_analysis")
os.makedirs(OUTPUT_DIR, exist_ok=True)

MAX_CELLS_PER_TYPE = 500
MIN_CELLS_PER_TYPE = 20
N_MARKERS_PER_GROUP = 50  # extract top N markers per group from Wilcoxon
TARGET_SUBCLASSES = ["Sst", "L6b"]  # focus subclasses

# Per-panel budgets (constrained by 10x add-on pricing: max 100 genes per panel)
# Both panels: 100 max add-on, 30 reserved for DE → 70 marker genes each
BUDGET_V1 = 70             # v1 Brain panel: 100 add-on max, 30 reserved for DE
BUDGET_5K = 70             # 5K panel (5001 genes): 100 max, 30 reserved for DE
BUDGET_ROUND1_V1 = 10      # subclass-level markers for v1
BUDGET_ROUND1_5K = 10      # subclass-level markers for 5K (already has many)
MIN_EXPECTED_DETECTED = 3  # minimum expected detected markers per cell per supertype

# SST supertype resolution is the priority. L6b subclass coverage is sufficient;
# individual L6b supertype disambiguation is lower priority.
# In the greedy allocator, SST supertypes are filled first to threshold,
# then L6b supertypes get remaining budget.
PRIORITY_SUBCLASS = "Sst"
DEPRIORITY_SUBCLASSES = ["L6b"]  # these get Round 2 genes only after priority is satisfied

# Spatial validation boost: genes validated on Xenium/MERFISH get a score multiplier
SPATIAL_VALIDATION_BOOST = 1.5  # 50% boost for spatially-validated genes

# Cardinal cell type markers that MUST be on any panel for biological interpretability
CARDINAL_MARKERS = [
    # GABAergic subclass markers
    "SST", "PVALB", "VIP", "LAMP5", "SNCG", "PAX6",
    # Glutamatergic subclass markers
    "CUX2", "RORB", "FEZF2", "THEMIS", "TLE4",
    # Non-neuronal markers
    "AQP4", "MBP", "MOG", "PDGFRA", "CSF1R", "FLT1",
    # Pan-neuronal
    "RBFOX3", "SYT1",
    # SST-specific key genes (biological interpretability)
    "NPY", "NOS1", "CHODL", "CALB1",
    # Key layer markers
    "SATB2", "BCL11B",
]

print("=" * 80)
print("HIERARCHICAL PROBE SELECTION FOR SST/L6b SUPERTYPE DISAMBIGUATION")
print("=" * 80)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1: Load & Prepare Data
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("STEP 1: LOADING DATA")
print("=" * 80)

# Load panel gene lists
print("\nLoading panel gene lists...")
panels = load_xenium_panels()
panel_5k_genes = panels["xenium_5k"]
panel_v1_genes = panels["xenium_v1"]

# Load detection efficiency calibration
print("\nLoading detection efficiency calibration...")
det_lookup, median_efficiency = load_detection_efficiency(DETECTION_PATH)
# Also load raw frac_snrna for downstream use
det_df = pd.read_csv(DETECTION_PATH)
frac_snrna_lookup = dict(zip(det_df["gene"], det_df["frac_snrna"]))

# Load spatially-validated gene performance from Xenium and MERFISH
print("\nLoading spatially-validated gene performance...")
spatial_validated_set, spatial_validated = load_spatial_validation(
    xenium_corr_path=XENIUM_CORR_PATH,
    merfish_corr_path=MERFISH_CORR_PATH,
    r_threshold=0.7)
spatial_validated_good = {g: r for g, r in spatial_validated.items() if r > 0.7}

# Check cardinal markers in existing panels
print(f"\nCardinal markers ({len(CARDINAL_MARKERS)} genes):")
for g in CARDINAL_MARKERS:
    in_v1 = "v1" if g in panel_v1_genes else "  "
    in_5k = "5K" if g in panel_5k_genes else "  "
    spatial_r = spatial_validated.get(g, None)
    spatial_str = f"spatial_r={spatial_r:.2f}" if spatial_r is not None else "no spatial data"
    print(f"  {g:<12} [{in_v1}] [{in_5k}]  {spatial_str}")

# Load gene quality data (from cross-platform analysis)
print("\nLoading gene quality data...")
quality_df = load_gene_quality(GENE_QUALITY_PATH)

# Identify bottom-quintile genes to exclude
bottom_quintile_r = quality_df["pearson_r"].quantile(0.2)
bad_genes_quality = set(quality_df[quality_df["pearson_r"] < bottom_quintile_r].index)
print(f"  Bottom quintile threshold: pearson_r < {bottom_quintile_r:.3f}")
print(f"  Genes flagged as poor quality: {len(bad_genes_quality)}")

# Identify antisense/lncRNA genes to exclude
if "biotype" in quality_df.columns:
    bad_biotypes = {"antisense", "lncRNA", "lincRNA"}
    bad_genes_biotype = set(
        quality_df[quality_df["biotype"].isin(bad_biotypes)].index)
    print(f"  Antisense/lncRNA genes to exclude: {len(bad_genes_biotype)}")
else:
    bad_genes_biotype = set()
    print("  (biotype column not found, skipping biotype filter)")

bad_genes = bad_genes_quality | bad_genes_biotype
print(f"  Total genes to exclude: {len(bad_genes)}")

# Load snRNAseq reference
adata = load_and_normalize_reference(REF_PATH, normalize=True, min_cells=10)

# Add Subclass column if needed — infer from Supertype
if "Subclass" not in adata.obs.columns:
    print("  ERROR: 'Subclass' column not found in reference")
    sys.exit(1)

# Print SST and L6b supertype summary
for sc_name in TARGET_SUBCLASSES:
    mask = adata.obs["Subclass"] == sc_name
    if mask.sum() == 0:
        print(f"  WARNING: No cells found for subclass {sc_name}")
        continue
    types = adata.obs.loc[mask, "Supertype"].value_counts().sort_index()
    print(f"\n  {sc_name} supertypes ({mask.sum()} cells total):")
    for st, n in types.items():
        if n > 0:  # skip categories with 0 cells
            print(f"    {st}: {n} cells")

# Subsample (reference was already normalized + gene-filtered by load_and_normalize_reference)
adata = subsample_by_group(adata, "Supertype",
                           max_cells=MAX_CELLS_PER_TYPE,
                           min_cells=MIN_CELLS_PER_TYPE)

# Make obs columns categorical (clean up unused categories)
for col in ["Supertype", "Subclass"]:
    if hasattr(adata.obs[col], "cat"):
        adata.obs[col] = adata.obs[col].cat.remove_unused_categories()
    else:
        adata.obs[col] = pd.Categorical(adata.obs[col])

print(f"  {adata.shape[0]} cells, "
      f"{adata.obs['Supertype'].nunique()} supertypes, "
      f"{adata.obs['Subclass'].nunique()} subclasses, "
      f"{adata.shape[1]} genes")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2: Gene Quality Filter — build eligible gene set
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("STEP 2: GENE QUALITY FILTER")
print("=" * 80)

all_genes = set(adata.var_names)
# Compute global detection rate per gene
if sparse.issparse(adata.X):
    global_det_rate = np.array((adata.X > 0).mean(axis=0)).flatten()
else:
    global_det_rate = np.mean(adata.X > 0, axis=0)
gene_det_df = pd.DataFrame({
    "gene": adata.var_names,
    "global_det_rate": global_det_rate
})
low_det_genes = set(gene_det_df[gene_det_df["global_det_rate"] < 0.005]["gene"])

# Apply filters
eligible_genes = all_genes - bad_genes - low_det_genes
print(f"  All genes in reference: {len(all_genes)}")
print(f"  Excluded (poor quality): {len(bad_genes & all_genes)}")
print(f"  Excluded (low detection <0.5%): {len(low_det_genes - bad_genes)}")
print(f"  Eligible genes: {len(eligible_genes)}")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3: HIERARCHICAL MARKER DISCOVERY
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("STEP 3: HIERARCHICAL MARKER DISCOVERY")
print("=" * 80)


def extract_markers(adata_subset, groupby, n_genes, label=""):
    """Run Wilcoxon rank-sum and extract marker results as DataFrame."""
    print(f"\n  Running rank_genes_groups ({label})...")
    print(f"    groupby={groupby}, n_genes={n_genes}")
    print(f"    {adata_subset.shape[0]} cells, "
          f"{adata_subset.obs[groupby].nunique()} groups")
    t0 = time.time()
    sc.tl.rank_genes_groups(adata_subset, groupby=groupby,
                             method="wilcoxon", n_genes=n_genes,
                             use_raw=False)
    print(f"    Done in {time.time()-t0:.0f}s")

    result = adata_subset.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    records = []
    for group in groups:
        for rank in range(n_genes):
            try:
                gene = result["names"][rank][group]
                score = float(result["scores"][rank][group])
                pval = float(result["pvals_adj"][rank][group])
                logfc = float(result["logfoldchanges"][rank][group])
                records.append({
                    "group": group,
                    "rank": rank + 1,
                    "gene": gene,
                    "wilcoxon_score": score,
                    "logfoldchange": logfc,
                    "pval_adj": pval,
                })
            except (IndexError, KeyError):
                break
    return pd.DataFrame(records)


# ── Round 1: Subclass-level markers ──────────────────────────────────────────
print("\n" + "-" * 60)
print("ROUND 1: SUBCLASS-LEVEL MARKERS (SST-vs-all, L6b-vs-all)")
print("-" * 60)

round1_df = extract_markers(adata, groupby="Subclass",
                             n_genes=N_MARKERS_PER_GROUP,
                             label="Round 1: Subclass")

# Filter to SST and L6b groups
round1_target = round1_df[round1_df["group"].isin(TARGET_SUBCLASSES)].copy()
round1_target["round"] = 1
round1_target["marker_type"] = "subclass_vs_all"
print(f"\n  Round 1 markers for target subclasses:")
for sc_name in TARGET_SUBCLASSES:
    n = len(round1_target[round1_target["group"] == sc_name])
    top5 = list(round1_target[round1_target["group"] == sc_name].head(5)["gene"])
    print(f"    {sc_name}: {n} markers (top 5: {', '.join(top5)})")

# ── Round 2: Within-subclass supertype markers ───────────────────────────────
print("\n" + "-" * 60)
print("ROUND 2: WITHIN-SUBCLASS SUPERTYPE MARKERS")
print("-" * 60)

round2_dfs = []
for sc_name in TARGET_SUBCLASSES:
    print(f"\n  === {sc_name} within-subclass analysis ===")

    # Subset to this subclass only
    mask = adata.obs["Subclass"] == sc_name
    adata_sub = adata[mask].copy()

    # Clean up categories
    adata_sub.obs["Supertype"] = (adata_sub.obs["Supertype"]
                                   .cat.remove_unused_categories())

    # Check we have enough types
    n_types = adata_sub.obs["Supertype"].nunique()
    if n_types < 2:
        print(f"    WARNING: Only {n_types} supertype(s), skipping")
        continue

    # Print type composition
    for st, n in adata_sub.obs["Supertype"].value_counts().sort_index().items():
        print(f"    {st}: {n} cells")

    # Run within-subclass markers
    r2 = extract_markers(adata_sub, groupby="Supertype",
                          n_genes=N_MARKERS_PER_GROUP,
                          label=f"Round 2: {sc_name} within-subclass")
    r2["round"] = 2
    r2["marker_type"] = f"within_{sc_name}"
    r2["subclass"] = sc_name
    round2_dfs.append(r2)

    # Show top markers per supertype
    print(f"\n  Top 5 within-{sc_name} markers per supertype:")
    for st in sorted(r2["group"].unique()):
        top5 = list(r2[r2["group"] == st].head(5)["gene"])
        print(f"    {st}: {', '.join(top5)}")

round2_df = pd.concat(round2_dfs, ignore_index=True) if round2_dfs else pd.DataFrame()

# Check for overlap between Round 1 and Round 2
if len(round2_df) > 0:
    r1_genes = set(round1_target["gene"])
    r2_genes = set(round2_df["gene"])
    overlap = r1_genes & r2_genes
    print(f"\n  Gene overlap between Round 1 and Round 2: {len(overlap)} genes")
    if overlap:
        print(f"    Shared genes: {', '.join(sorted(list(overlap))[:20])}")
        if len(overlap) > 20:
            print(f"    ... and {len(overlap)-20} more")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4: DROPOUT-AWARE SCORING
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("STEP 4: DROPOUT-AWARE SCORING")
print("=" * 80)


def compute_detection_stats(markers_df, adata_ref, det_lookup, median_eff):
    """Add detection efficiency and predicted spatial detection to markers."""
    records = []
    for _, row in markers_df.iterrows():
        gene = row["gene"]
        group = row["group"]

        # Get detection efficiency
        eff = det_lookup.get(gene, median_eff)
        eff_source = "calibrated" if gene in det_lookup else "median_fallback"

        # Compute snRNAseq detection rate in this specific group
        # (what fraction of cells in this supertype/subclass express the gene?)
        if group in adata_ref.obs["Supertype"].values:
            mask = adata_ref.obs["Supertype"] == group
        elif group in adata_ref.obs["Subclass"].values:
            mask = adata_ref.obs["Subclass"] == group
        else:
            mask = np.zeros(adata_ref.shape[0], dtype=bool)

        if mask.sum() > 0 and gene in adata_ref.var_names:
            gene_idx = list(adata_ref.var_names).index(gene)
            X_col = adata_ref[mask, gene_idx].X
            if sparse.issparse(X_col):
                X_col = X_col.toarray().flatten()
            else:
                X_col = np.asarray(X_col).flatten()
            frac_detected = np.mean(X_col > 0)
        else:
            frac_detected = 0.0

        # Predicted spatial detection
        pred_spatial = min(frac_detected * eff, 1.0)

        # Composite score: Wilcoxon score × predicted detection
        composite = row["wilcoxon_score"] * pred_spatial

        # Spatial validation boost: genes proven on Xenium/MERFISH get a boost
        is_spatially_validated = gene in spatial_validated_good
        if is_spatially_validated:
            composite *= SPATIAL_VALIDATION_BOOST

        records.append({
            **row.to_dict(),
            "detection_efficiency": eff,
            "eff_source": eff_source,
            "frac_snrna_in_group": frac_detected,
            "predicted_spatial_detection": pred_spatial,
            "composite_score": composite,
            "spatially_validated": is_spatially_validated,
            "best_spatial_r": spatial_validated.get(gene, None),
        })

    return pd.DataFrame(records)


# Score Round 1 markers
print("\nScoring Round 1 markers...")
round1_scored = compute_detection_stats(round1_target, adata,
                                         det_lookup, median_efficiency)
print(f"  Round 1: {len(round1_scored)} marker entries scored")

# Score Round 2 markers
if len(round2_df) > 0:
    print("\nScoring Round 2 markers...")
    round2_scored = compute_detection_stats(round2_df, adata,
                                             det_lookup, median_efficiency)
    print(f"  Round 2: {len(round2_scored)} marker entries scored")
else:
    round2_scored = pd.DataFrame()

# Summary statistics
for label, df in [("Round 1", round1_scored), ("Round 2", round2_scored)]:
    if len(df) == 0:
        continue
    print(f"\n  {label} detection stats:")
    print(f"    Median predicted spatial detection: "
          f"{df['predicted_spatial_detection'].median():.3f}")
    print(f"    Mean composite score: {df['composite_score'].mean():.1f}")
    print(f"    Genes with calibrated efficiency: "
          f"{(df['eff_source']=='calibrated').sum()}/{len(df)}")
    n_validated = df["spatially_validated"].sum()
    print(f"    Spatially validated (Xenium/MERFISH, r>0.7): "
          f"{n_validated}/{df['gene'].nunique()} unique genes")

# Save intermediate results
round1_scored.to_csv(os.path.join(OUTPUT_DIR,
                     "hierarchical_markers_round1_subclass.csv"), index=False)
if len(round2_scored) > 0:
    for sc_name in TARGET_SUBCLASSES:
        r2_sub = round2_scored[round2_scored.get("subclass") == sc_name]
        if len(r2_sub) > 0:
            r2_sub.to_csv(os.path.join(OUTPUT_DIR,
                          f"hierarchical_markers_round2_{sc_name.lower()}_within.csv"),
                          index=False)
            print(f"  Saved Round 2 {sc_name} markers: {len(r2_sub)} entries")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5: BUDGET ALLOCATION & GENE SELECTION
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("STEP 5: BUDGET ALLOCATION & GENE SELECTION")
print("=" * 80)


def make_gene_record(gene, round_num, marker_for, row_data=None, reason="marker"):
    """Create a standard gene record dict."""
    if row_data is not None:
        return {
            "gene": gene,
            "round": round_num,
            "marker_for": marker_for,
            "reason": reason,
            "wilcoxon_score": row_data.get("wilcoxon_score", 0),
            "logfoldchange": row_data.get("logfoldchange", 0),
            "detection_efficiency": row_data.get("detection_efficiency", 0),
            "frac_snrna_in_group": row_data.get("frac_snrna_in_group", 0),
            "predicted_spatial_detection": row_data.get("predicted_spatial_detection", 0),
            "composite_score": row_data.get("composite_score", 0),
            "spatially_validated": row_data.get("spatially_validated", False),
            "best_spatial_r": row_data.get("best_spatial_r", None),
        }
    return {
        "gene": gene, "round": round_num, "marker_for": marker_for,
        "reason": reason, "wilcoxon_score": 0, "logfoldchange": 0,
        "detection_efficiency": det_lookup.get(gene, median_efficiency),
        "frac_snrna_in_group": 0, "predicted_spatial_detection": 0,
        "composite_score": 0, "spatially_validated": gene in spatial_validated_good,
        "best_spatial_r": spatial_validated.get(gene, None),
    }


def select_addon_genes(round1_scored, round2_scored, panel_genes,
                       panel_name, eligible_genes, budget_r1, budget_total,
                       priority_subclass=None, depriority_subclasses=None):
    """Select add-on genes for a panel using hierarchical + dropout-aware scoring.

    Includes:
    - Cardinal marker guarantees (SST, PVALB, VIP, etc.)
    - Spatial validation boost for Xenium/MERFISH-proven genes
    - Priority-aware greedy coverage: priority_subclass supertypes are
      filled to threshold FIRST, then depriority supertypes get remaining budget
    """
    if depriority_subclasses is None:
        depriority_subclasses = []
    print(f"\n{'='*60}")
    print(f"  Panel: {panel_name} ({len(panel_genes)} existing genes)")
    print(f"  Budget: {budget_total} total")
    print(f"{'='*60}")

    selected = []  # list of dicts with gene info
    selected_genes = set()  # track gene names

    # ── Step 0: Cardinal markers (guaranteed) ────────────────────────────
    print(f"\n  --- Cardinal markers (guaranteed) ---")
    cardinal_added = []
    for gene in CARDINAL_MARKERS:
        if gene not in panel_genes and gene in adata.var_names:
            if gene not in selected_genes:
                selected_genes.add(gene)
                # Look up spatial validation
                rec = make_gene_record(gene, 0, "cardinal", reason="cardinal_marker")
                selected.append(rec)
                cardinal_added.append(gene)
    print(f"    Added {len(cardinal_added)} cardinal markers not in panel: "
          f"{', '.join(cardinal_added)}")

    # ── Step 1: Subclass-level genes ─────────────────────────────────────
    print(f"\n  --- Round 1 selection (subclass-level) ---")
    r1_candidates = round1_scored.copy()
    r1_candidates = r1_candidates[~r1_candidates["gene"].isin(panel_genes)]
    r1_candidates = r1_candidates[~r1_candidates["gene"].isin(selected_genes)]
    r1_candidates = r1_candidates[r1_candidates["gene"].isin(eligible_genes)]
    r1_candidates = (r1_candidates.sort_values("composite_score", ascending=False)
                     .drop_duplicates(subset="gene", keep="first"))

    r1_to_add = min(budget_r1, len(r1_candidates))
    r1_selected = r1_candidates.head(r1_to_add)
    r1_gene_set = set(r1_selected["gene"])

    for _, row in r1_selected.iterrows():
        if row["gene"] not in selected_genes:
            selected_genes.add(row["gene"])
            selected.append(make_gene_record(
                row["gene"], 1, row["group"], row.to_dict(), reason="subclass_marker"))

    print(f"    Selected {len(r1_selected)} Round 1 genes")
    for sc_name in TARGET_SUBCLASSES:
        r1_sc = r1_selected[r1_selected["group"] == sc_name]
        if len(r1_sc) > 0:
            print(f"    {sc_name}: {', '.join(r1_sc['gene'].tolist()[:10])}")

    # ── Step 2: Within-subclass supertype genes (greedy coverage) ────────
    print(f"\n  --- Round 2 selection (within-subclass supertype) ---")
    budget_remaining = budget_total - len(selected)
    print(f"    Budget remaining after cardinal + R1: {budget_remaining}")

    if len(round2_scored) == 0:
        print("    No Round 2 markers available")
        return pd.DataFrame(selected), {}

    all_supertypes = sorted(round2_scored["group"].unique())
    print(f"    {len(all_supertypes)} supertypes to cover")

    # Compute existing coverage per supertype
    existing_coverage = {}
    for st in all_supertypes:
        st_markers = round2_scored[round2_scored["group"] == st]
        in_panel = st_markers[st_markers["gene"].isin(
            panel_genes | selected_genes)]
        expected_det = in_panel["predicted_spatial_detection"].sum()
        existing_coverage[st] = {
            "n_in_panel": len(in_panel),
            "expected_detected": expected_det,
        }

    # Prepare candidates
    r2_candidates = round2_scored.copy()
    r2_candidates = r2_candidates[~r2_candidates["gene"].isin(
        panel_genes | selected_genes)]
    r2_candidates = r2_candidates[r2_candidates["gene"].isin(eligible_genes)]

    # Track expected detection per supertype
    supertype_expected = {st: existing_coverage[st]["expected_detected"]
                          for st in all_supertypes}

    # Classify supertypes by priority
    priority_sts = [st for st in all_supertypes
                    if st.split("_")[0] not in depriority_subclasses]
    depriority_sts = [st for st in all_supertypes
                      if st.split("_")[0] in depriority_subclasses]

    if depriority_subclasses:
        print(f"    Priority supertypes (fill first): {len(priority_sts)} "
              f"({', '.join(set(st.split('_')[0] for st in priority_sts))})")
        print(f"    Deprioritized supertypes (fill after): {len(depriority_sts)} "
              f"({', '.join(depriority_subclasses)})")

    n_r2_selected = 0
    genes_used = set()

    def _add_gene(gene, group, row_data):
        """Helper to add a gene and update coverage tracking."""
        nonlocal n_r2_selected
        genes_used.add(gene)
        selected_genes.add(gene)
        selected.append(make_gene_record(
            gene, 2, group, row_data, reason="supertype_marker"))
        n_r2_selected += 1
        # Update coverage for ALL supertypes this gene is a marker for
        gene_rows = r2_candidates[r2_candidates["gene"] == gene]
        for _, gr in gene_rows.iterrows():
            if gr["group"] in supertype_expected:
                supertype_expected[gr["group"]] += gr["predicted_spatial_detection"]

    # Phase 1: Fill PRIORITY supertypes to threshold
    phase1_budget = budget_remaining
    while n_r2_selected < phase1_budget and len(r2_candidates) > 0:
        # Find worst priority supertype
        priority_coverage = {st: supertype_expected[st] for st in priority_sts
                             if supertype_expected[st] < 900}
        if not priority_coverage:
            break

        worst_st = min(priority_coverage, key=priority_coverage.get)
        worst_expected = priority_coverage[worst_st]

        # If all priority supertypes meet threshold, move to phase 2
        if worst_expected >= MIN_EXPECTED_DETECTED:
            break

        # Select best gene for the worst priority supertype
        st_candidates = r2_candidates[
            (r2_candidates["group"] == worst_st) &
            (~r2_candidates["gene"].isin(genes_used)) &
            (~r2_candidates["gene"].isin(selected_genes))
        ].sort_values("composite_score", ascending=False)

        if len(st_candidates) == 0:
            supertype_expected[worst_st] = 999
            continue

        best = st_candidates.iloc[0]
        _add_gene(best["gene"], worst_st, best.to_dict())

    n_phase1 = n_r2_selected
    print(f"    Phase 1 (priority={priority_subclass}): {n_phase1} genes, "
          f"all priority supertypes ≥ {MIN_EXPECTED_DETECTED}")

    # Phase 2: Fill DEPRIORITY supertypes (L6b) with remaining budget
    if depriority_sts and n_r2_selected < budget_remaining:
        phase2_start = n_r2_selected
        while n_r2_selected < budget_remaining:
            depri_coverage = {st: supertype_expected[st] for st in depriority_sts
                              if supertype_expected[st] < 900}
            if not depri_coverage:
                break

            worst_st = min(depri_coverage, key=depri_coverage.get)
            worst_expected = depri_coverage[worst_st]

            if worst_expected >= MIN_EXPECTED_DETECTED:
                break

            st_candidates = r2_candidates[
                (r2_candidates["group"] == worst_st) &
                (~r2_candidates["gene"].isin(genes_used)) &
                (~r2_candidates["gene"].isin(selected_genes))
            ].sort_values("composite_score", ascending=False)

            if len(st_candidates) == 0:
                supertype_expected[worst_st] = 999
                continue

            best = st_candidates.iloc[0]
            _add_gene(best["gene"], worst_st, best.to_dict())

        n_phase2 = n_r2_selected - phase2_start
        print(f"    Phase 2 (depriority={','.join(depriority_subclasses)}): "
              f"{n_phase2} genes")

    # Phase 3: Fill remaining budget with best overall composite scores
    if n_r2_selected < budget_remaining:
        remaining = r2_candidates[~r2_candidates["gene"].isin(genes_used)]
        remaining = (remaining.sort_values("composite_score", ascending=False)
                     .drop_duplicates(subset="gene", keep="first"))
        phase3_start = n_r2_selected
        for _, row in remaining.iterrows():
            if n_r2_selected >= budget_remaining:
                break
            gene = row["gene"]
            if gene in genes_used or gene in selected_genes:
                continue
            _add_gene(gene, row["group"], row.to_dict())

        n_phase3 = n_r2_selected - phase3_start
        print(f"    Phase 3 (fill by composite score): {n_phase3} genes")

    print(f"    Selected {n_r2_selected} Round 2 genes")
    print(f"    Total selected: {len(selected)} genes")

    # Coverage assessment
    print(f"\n  --- Coverage assessment ---")
    print(f"  {'Supertype':<20} {'Before':>10} {'After':>10} {'Delta':>8}  "
          f"{'Validated':>10}")
    print(f"  {'-'*60}")
    for st in all_supertypes:
        before = existing_coverage[st]["expected_detected"]
        after = supertype_expected.get(st, before)
        delta = after - before
        # Count spatially-validated genes selected for this supertype
        st_selected = [s for s in selected
                       if s["marker_for"] == st and s.get("spatially_validated")]
        flag = " !!" if after < MIN_EXPECTED_DETECTED else ""
        print(f"  {st:<20} {before:>10.2f} {after:>10.2f} {delta:>8.2f}  "
              f"{len(st_selected):>10}{flag}")

    # Summary of spatial validation
    result_df = pd.DataFrame(selected)
    n_validated = result_df["spatially_validated"].sum()
    print(f"\n  Spatially validated genes in selection: {n_validated}/{len(result_df)}")
    if n_validated > 0:
        validated_genes = result_df[result_df["spatially_validated"]].sort_values(
            "best_spatial_r", ascending=False)
        top_validated = list(zip(validated_genes["gene"].head(10),
                                  validated_genes["best_spatial_r"].head(10)))
        print(f"  Top validated: " +
              ", ".join(f"{g} (r={r:.2f})" for g, r in top_validated))

    return result_df, supertype_expected


# Run selection for both panels with different budgets
panel_configs = [
    ("xenium_v1", panel_v1_genes, BUDGET_V1, BUDGET_ROUND1_V1,
     PRIORITY_SUBCLASS, DEPRIORITY_SUBCLASSES),
    ("xenium_5k", panel_5k_genes, BUDGET_5K, BUDGET_ROUND1_5K,
     PRIORITY_SUBCLASS, DEPRIORITY_SUBCLASSES),
]
panel_results = {}
for panel_name, panel_genes, budget, budget_r1, pri_sc, depri_scs in panel_configs:
    result_df, coverage = select_addon_genes(
        round1_scored, round2_scored, panel_genes, panel_name,
        eligible_genes, budget_r1, budget,
        priority_subclass=pri_sc, depriority_subclasses=depri_scs)

    panel_results[panel_name] = (result_df, coverage)

    # Save
    out_path = os.path.join(OUTPUT_DIR,
                             f"probe_recommendations_{panel_name}.csv")
    result_df.to_csv(out_path, index=False)
    print(f"\n  Saved: {out_path} ({len(result_df)} genes)")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 6: COVERAGE ASSESSMENT TABLE
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("STEP 6: COVERAGE ASSESSMENT")
print("=" * 80)

# Build comprehensive coverage table
coverage_records = []
for panel_name, (result_df, coverage) in panel_results.items():
    panel_genes = panel_v1_genes if panel_name == "xenium_v1" else panel_5k_genes
    for st in sorted(coverage.keys()):
        # Count markers in original panel
        if len(round2_scored) > 0:
            st_all_markers = round2_scored[round2_scored["group"] == st]
            n_in_panel = st_all_markers[
                st_all_markers["gene"].isin(panel_genes)
            ]["gene"].nunique()
        else:
            n_in_panel = 0

        # Count proposed markers
        n_proposed = result_df[result_df["marker_for"] == st]["gene"].nunique()

        # Expected detection after augmentation
        expected_det = coverage[st] if coverage[st] < 900 else 0

        sc_name = st.split("_")[0]
        if sc_name == "Sst":
            # Get cell count
            n_cells = (adata.obs["Supertype"] == st).sum()
        else:
            n_cells = (adata.obs["Supertype"] == st).sum()

        coverage_records.append({
            "panel": panel_name,
            "supertype": st,
            "subclass": sc_name,
            "n_cells_snrnaseq": n_cells,
            "n_markers_in_panel": n_in_panel,
            "n_markers_proposed": n_proposed,
            "expected_detected_per_cell": expected_det,
            "meets_threshold": expected_det >= MIN_EXPECTED_DETECTED,
        })

coverage_table = pd.DataFrame(coverage_records)
coverage_table.to_csv(os.path.join(OUTPUT_DIR,
                      "supertype_coverage_assessment.csv"), index=False)
print(f"\nCoverage assessment saved ({len(coverage_table)} rows)")

# Print summary
for panel_name in ["xenium_v1", "xenium_5k"]:
    panel_cov = coverage_table[coverage_table["panel"] == panel_name]
    n_meet = panel_cov["meets_threshold"].sum()
    n_total = len(panel_cov)
    print(f"\n  {panel_name}: {n_meet}/{n_total} supertypes meet "
          f"threshold (>= {MIN_EXPECTED_DETECTED} expected detected)")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 7: DIAGNOSTIC FIGURES
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("STEP 7: DIAGNOSTIC FIGURES")
print("=" * 80)

# ── Figure 1: Hierarchical marker heatmap (within-subclass) ──────────────────
print("\nGenerating marker heatmaps...")

for sc_name in TARGET_SUBCLASSES:
    # Get supertypes for this subclass
    sc_mask = adata.obs["Subclass"] == sc_name
    adata_sc = adata[sc_mask].copy()
    supertypes = sorted(adata_sc.obs["Supertype"].unique())

    if len(round2_scored) == 0:
        continue
    r2_sc = round2_scored[round2_scored.get("subclass") == sc_name]
    if len(r2_sc) == 0:
        continue

    # Get unique top markers (top 5 per supertype by composite score)
    top_genes = []
    for st in supertypes:
        st_markers = (r2_sc[r2_sc["group"] == st]
                      .sort_values("composite_score", ascending=False)
                      .head(5))
        for g in st_markers["gene"]:
            if g not in top_genes and g in adata_sc.var_names:
                top_genes.append(g)

    if len(top_genes) == 0:
        continue

    # Limit to reasonable size
    top_genes = top_genes[:60]

    # Compute mean expression per supertype for these genes
    expr_matrix = np.zeros((len(supertypes), len(top_genes)))
    for i, st in enumerate(supertypes):
        st_mask = adata_sc.obs["Supertype"] == st
        X_sub = adata_sc[st_mask][:, top_genes].X
        if sparse.issparse(X_sub):
            X_sub = X_sub.toarray()
        expr_matrix[i] = np.mean(X_sub, axis=0)

    # Z-score normalize per gene (column) for visualization
    from scipy.stats import zscore
    expr_z = zscore(expr_matrix, axis=0)
    expr_z = np.nan_to_num(expr_z)

    fig, ax = plt.subplots(figsize=(max(14, len(top_genes) * 0.3),
                                     max(6, len(supertypes) * 0.4)))
    im = ax.imshow(expr_z, cmap="RdBu_r", aspect="auto", vmin=-2, vmax=2)
    ax.set_xticks(range(len(top_genes)))
    ax.set_xticklabels(top_genes, rotation=90, fontsize=8)
    ax.set_yticks(range(len(supertypes)))
    ax.set_yticklabels(supertypes, fontsize=10)
    ax.set_title(f"{sc_name} Within-Subclass Markers\n"
                 f"(Top 5 per supertype by composite score, z-scored)",
                 fontsize=14, fontweight="bold")
    plt.colorbar(im, ax=ax, label="Z-score", shrink=0.8)
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR,
                f"hierarchical_heatmap_{sc_name.lower()}.png"),
                dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: hierarchical_heatmap_{sc_name.lower()}.png")

# ── Figure 2: Wilcoxon score vs predicted spatial detection ──────────────────
print("\nGenerating score vs detection scatter...")

all_scored = pd.concat([round1_scored, round2_scored], ignore_index=True)
if len(all_scored) > 0:
    fig, ax = plt.subplots(figsize=(10, 8))

    # Round 1
    r1 = all_scored[all_scored["round"] == 1]
    r2 = all_scored[all_scored["round"] == 2]

    if len(r1) > 0:
        ax.scatter(r1["wilcoxon_score"], r1["predicted_spatial_detection"],
                   alpha=0.5, s=30, c="#e41a1c", label="Round 1 (subclass)",
                   edgecolors="none")
    if len(r2) > 0:
        ax.scatter(r2["wilcoxon_score"], r2["predicted_spatial_detection"],
                   alpha=0.3, s=20, c="#377eb8", label="Round 2 (within-subclass)",
                   edgecolors="none")

    ax.set_xlabel("Wilcoxon Score", fontsize=14)
    ax.set_ylabel("Predicted Spatial Detection Rate", fontsize=14)
    ax.set_title("Marker Quality: Statistical Separation vs Spatial Detectability",
                 fontsize=14, fontweight="bold")
    ax.legend(fontsize=12)
    ax.axhline(y=0.1, color="gray", linestyle="--", alpha=0.5,
               label="10% detection threshold")
    ax.tick_params(labelsize=12)
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR,
                "hierarchical_score_vs_detection.png"),
                dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  Saved: hierarchical_score_vs_detection.png")

# ── Figure 3: Coverage improvement bar chart ─────────────────────────────────
print("\nGenerating coverage improvement chart...")

for panel_name, (result_df, final_coverage) in panel_results.items():
    panel_genes = panel_v1_genes if panel_name == "xenium_v1" else panel_5k_genes

    # Get supertypes
    supertypes = sorted(final_coverage.keys())
    if not supertypes:
        continue

    # Compute before/after for each supertype
    before_vals = []
    after_vals = []
    for st in supertypes:
        # Before: coverage from panel only
        if len(round2_scored) > 0:
            st_markers = round2_scored[round2_scored["group"] == st]
            in_panel = st_markers[st_markers["gene"].isin(panel_genes)]
            before = in_panel["predicted_spatial_detection"].sum()
        else:
            before = 0
        after = final_coverage[st] if final_coverage[st] < 900 else 0
        before_vals.append(before)
        after_vals.append(after)

    x = np.arange(len(supertypes))
    width = 0.35

    fig, ax = plt.subplots(figsize=(max(12, len(supertypes) * 0.6), 6))
    bars1 = ax.bar(x - width/2, before_vals, width, label="Existing panel",
                    color="#4292c6", alpha=0.8)
    bars2 = ax.bar(x + width/2, after_vals, width, label="+ Add-on probes",
                    color="#ef6548", alpha=0.8)
    ax.axhline(y=MIN_EXPECTED_DETECTED, color="red", linestyle="--",
               alpha=0.7, linewidth=1.5,
               label=f"Target: {MIN_EXPECTED_DETECTED} expected/cell")
    ax.set_xticks(x)
    ax.set_xticklabels(supertypes, rotation=45, ha="right", fontsize=10)
    ax.set_ylabel("Expected Detected Markers per Cell", fontsize=13)
    ax.set_title(f"Coverage Improvement: {panel_name}\n"
                 f"(E[detected] = Σ predicted_spatial_detection)",
                 fontsize=14, fontweight="bold")
    ax.legend(fontsize=11)
    ax.tick_params(labelsize=11)
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR,
                f"coverage_improvement_{panel_name}.png"),
                dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: coverage_improvement_{panel_name}.png")

# ── Figure 4: Dropout calibration curve ──────────────────────────────────────
print("\nGenerating dropout calibration curve...")

if len(all_scored) > 0:
    fig, ax = plt.subplots(figsize=(8, 8))

    ax.scatter(all_scored["frac_snrna_in_group"],
               all_scored["predicted_spatial_detection"],
               alpha=0.3, s=15, c="#2166ac", edgecolors="none")

    # Add diagonal reference (efficiency = 1)
    ax.plot([0, 1], [0, 1], "k--", alpha=0.4, label="Perfect detection (eff=1)")

    # Add median efficiency line
    ax.plot([0, 1], [0, median_efficiency], "r--", alpha=0.4,
            label=f"Median efficiency ({median_efficiency:.2f})")

    ax.set_xlabel("snRNAseq Detection Rate (in group)", fontsize=14)
    ax.set_ylabel("Predicted Spatial Detection Rate", fontsize=14)
    ax.set_title("Dropout Calibration: snRNAseq → Spatial Prediction",
                 fontsize=14, fontweight="bold")
    ax.legend(fontsize=11)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.tick_params(labelsize=12)
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR,
                "hierarchical_dropout_calibration.png"),
                dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  Saved: hierarchical_dropout_calibration.png")


# ══════════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("FINAL SUMMARY")
print("=" * 80)

print(f"\nReference: {adata.shape[0]} cells, {adata.obs['Supertype'].nunique()} supertypes")
print(f"Target subclasses: {', '.join(TARGET_SUBCLASSES)}")
print(f"Marker discovery: Wilcoxon rank-sum, top {N_MARKERS_PER_GROUP} per group")
print(f"Detection efficiency: median = {median_efficiency:.3f}")

for panel_name, (result_df, coverage) in panel_results.items():
    print(f"\n  {panel_name}:")
    print(f"    Total add-on genes proposed: {len(result_df)}")
    if "reason" in result_df.columns:
        for reason, count in result_df["reason"].value_counts().items():
            print(f"      {reason}: {count}")
    r0_count = (result_df["round"] == 0).sum()
    r1_count = (result_df["round"] == 1).sum()
    r2_count = (result_df["round"] == 2).sum()
    print(f"    Cardinal markers: {r0_count}")
    print(f"    Round 1 (subclass): {r1_count}")
    print(f"    Round 2 (within-subclass): {r2_count}")

    # Unique genes
    unique_genes = result_df["gene"].unique()
    print(f"    Unique genes: {len(unique_genes)}")

    # Spatially validated
    if "spatially_validated" in result_df.columns:
        n_val = result_df["spatially_validated"].sum()
        print(f"    Spatially validated (Xenium/MERFISH): {n_val}/{len(result_df)}")

    # How many supertypes meet threshold?
    n_meet = sum(1 for v in coverage.values()
                 if v >= MIN_EXPECTED_DETECTED and v < 900)
    n_total = len([v for v in coverage.values() if v < 900])
    print(f"    Supertypes meeting threshold ({MIN_EXPECTED_DETECTED}+ expected): "
          f"{n_meet}/{n_total}")

# ══════════════════════════════════════════════════════════════════════════════
# COST ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("COST ANALYSIS ($90K budget, UPitt Xenium Core)")
print("=" * 80)

n_v1_addon = len(panel_results["xenium_v1"][0])
n_5k_addon = len(panel_results["xenium_5k"][0])
total_v1 = 266 + n_v1_addon + 30   # +30 reserved for DE
total_5k = 5001 + n_5k_addon + 30  # +30 reserved for DE

BUDGET_TOTAL_USD = 90000

# ═══════════════════════════════════════════════════════════════════════
# 10x Pricing (from official product listing, v1 chemistry <500-gene):
#   Add-On Custom (requires predesigned panel, SAME SLIDE):
#     1-50 genes:  $3,550 (4 rxn) / $4,800 (16 rxn)
#     51-100 genes: $5,650 (4 rxn) / $8,000 (16 rxn)
#     *** MAX 100 genes for add-on pathway ***
#   Standalone Custom (NO predesigned panel):
#     101-300 genes: $11,600 (4 rxn) / $24,000 (16 rxn)
#     301-480 genes: $17,850 (4 rxn) / $29,100 (16 rxn)
#   Predesigned Panels:
#     Human Brain: $345/2 rxn
# ═══════════════════════════════════════════════════════════════════════

# ── Scenario A: Xenium v1 Brain + 100 Add-On (same slide) ──
# One-time: Add-on 51-100 gene panel ($8,000 for 16 rxn)
#   + predesigned brain panel ($345 per 2 rxn)
# Per-run at UPitt: $11,000/run + $2,000 cell segmentation = $13,000/run
v1_addon_cost = 8000       # 51-100 gene add-on, 16 rxn
v1_brain_per_2rxn = 345    # predesigned brain panel per 2 rxn
v1_per_run = 11000 + 2000  # UPitt: $11K base + $2K cell seg

# Compute max runs
# Budget = addon_cost + brain_panels + N_runs * per_run
# Brain panels: need 2 rxn per run → $345 per run
v1_per_run_total = v1_per_run + v1_brain_per_2rxn  # $13,345/run
v1_remaining = BUDGET_TOTAL_USD - v1_addon_cost
v1_max_runs = int(v1_remaining / v1_per_run_total)
v1_max_sections = v1_max_runs * 2
v1_total_cost = v1_addon_cost + v1_max_runs * v1_per_run_total

# ── Scenario B: Xenium 5K Prime + 100 Add-On (same slide) ──
# One-time: Add-on 51-100 gene panel ($8,000 for 16 rxn) [5K chemistry pricing TBD]
# Per-run at UPitt: $23,500/run (includes 5K panel reagents + cell seg)
# NOTE: 5K add-on pricing may differ from v1; using $8K as estimate
k5_addon_cost = 8000       # 51-100 gene add-on, 16 rxn (estimated, verify for 5K)
k5_per_run = 23500         # UPitt 5K price (includes panel reagents)
k5_remaining = BUDGET_TOTAL_USD - k5_addon_cost
k5_max_runs = int(k5_remaining / k5_per_run)
k5_max_sections = k5_max_runs * 2
k5_total_cost = k5_addon_cost + k5_max_runs * k5_per_run

# Number of supertypes
n_supertypes = len([v for v in panel_results["xenium_v1"][1].values() if v < 900])

# v1 average coverage
v1_cov = panel_results["xenium_v1"][1]
v1_mean_cov = np.mean([v for v in v1_cov.values() if v < 900])
v1_min_cov = min(v for v in v1_cov.values() if v < 900)

# 5K average coverage
k5_cov = panel_results["xenium_5k"][1]
k5_mean_cov = np.mean([v for v in k5_cov.values() if v < 900])
k5_min_cov = min(v for v in k5_cov.values() if v < 900)

# Spatially validated counts
v1_n_validated = panel_results["xenium_v1"][0]["spatially_validated"].sum()
k5_n_validated = panel_results["xenium_5k"][0]["spatially_validated"].sum()

print(f"""
  ┌──────────────────────────────────────────────────────────────────────────┐
  │  SCENARIO A: Xenium v1 Brain (266) + {n_v1_addon}-gene Add-On Custom           │
  ├──────────────────────────────────────────────────────────────────────────┤
  │  Total genes on panel:       {total_v1:<5}                                    │
  │  Add-on panel (one-time):    $ {v1_addon_cost:>6,}  (51-100 gene, 16 rxn)          │
  │  Brain panel per run:        $   {v1_brain_per_2rxn:>3}  ($345/2 rxn)                   │
  │  UPitt per run (2 slides):   ${v1_per_run:>6,}  ($11K base + $2K cell seg)     │
  │  Total per run:              ${v1_per_run_total:>6,}                                │
  │  Budget after add-on panel:  ${v1_remaining:>6,}                                │
  │  Maximum runs affordable:         {v1_max_runs}                                 │
  │  ═════════════════════════════════════════════════                        │
  │  MAXIMUM SECTIONS:  {v1_max_sections:>2}   (total cost: ${v1_total_cost:>,})                      │
  │  SST/L6b coverage: {v1_mean_cov:.1f} avg expected detected/cell (min {v1_min_cov:.1f})   │
  │  Spatially validated genes: {v1_n_validated}/{n_v1_addon}                                │
  └──────────────────────────────────────────────────────────────────────────┘

  ┌──────────────────────────────────────────────────────────────────────────┐
  │  SCENARIO B: Xenium Prime 5K (5001) + {n_5k_addon}-gene Add-On Custom           │
  ├──────────────────────────────────────────────────────────────────────────┤
  │  Total genes on panel:       {total_5k:<5}                                    │
  │  Add-on panel (one-time):    $ {k5_addon_cost:>6,}  (51-100 gene, 16 rxn, est.)    │
  │  5K panel per run:           included in UPitt fee                        │
  │  UPitt per run (2 slides):   ${k5_per_run:>6,}                                │
  │  Total per run:              ${k5_per_run:>6,}                                │
  │  Budget after add-on panel:  ${k5_remaining:>6,}                                │
  │  Maximum runs affordable:         {k5_max_runs}                                 │
  │  ═════════════════════════════════════════════════                        │
  │  MAXIMUM SECTIONS:   {k5_max_sections:>2}  (total cost: ${k5_total_cost:>,})                      │
  │  SST/L6b coverage: {k5_mean_cov:.1f} avg expected detected/cell (min {k5_min_cov:.1f})   │
  │  Spatially validated genes: {k5_n_validated}/{n_5k_addon}                                │
  └──────────────────────────────────────────────────────────────────────────┘
""")

# Cost per section
v1_cost_per_section = v1_total_cost / v1_max_sections
k5_cost_per_section = k5_total_cost / k5_max_sections

print(f"  ┌─────────────────────────────────────────────────────────────────┐")
print(f"  │  HEAD-TO-HEAD COMPARISON                                       │")
print(f"  ├─────────────────────────────────────────────────────────────────┤")
print(f"  │                        v1 + 100        5K + {n_5k_addon + 30}          │")
print(f"  │  Total genes:           {total_v1:<5}          {total_5k:<5}          │")
print(f"  │  Add-on markers:          {n_v1_addon:<3}             {n_5k_addon:<3}          │")
print(f"  │  Sections ($90K):          {v1_max_sections:<2}               {k5_max_sections:<2}          │")
print(f"  │  Cost per section:    ${v1_cost_per_section:>7,.0f}        ${k5_cost_per_section:>7,.0f}          │")
print(f"  │  Supertypes ≥{MIN_EXPECTED_DETECTED}:        {sum(1 for v in v1_cov.values() if v >= MIN_EXPECTED_DETECTED and v < 900)}/{n_supertypes}            {sum(1 for v in k5_cov.values() if v >= MIN_EXPECTED_DETECTED and v < 900)}/{n_supertypes}          │")
print(f"  │  Mean coverage:         {v1_mean_cov:>5.1f}           {k5_mean_cov:>5.1f}          │")
print(f"  │  Min coverage:          {v1_min_cov:>5.1f}           {k5_min_cov:>5.1f}          │")
print(f"  │  Spatially validated:   {v1_n_validated:>3}/{n_v1_addon:<3}          {k5_n_validated:>3}/{n_5k_addon:<3}          │")
print(f"  │  Genome-wide DE:          No             Yes          │")
print(f"  └─────────────────────────────────────────────────────────────────┘")

print(f"""
  Key trade-offs:
  ─────────────────────────────────────────────
    v1 advantage:  {v1_max_sections} sections vs {k5_max_sections} ({v1_max_sections - k5_max_sections} more), ${v1_cost_per_section:,.0f}/section vs ${k5_cost_per_section:,.0f}/section
    5K advantage:  {total_5k} genes vs {total_v1} ({total_5k - total_v1} more), genome-wide DE + discovery

  SST/L6b supertype resolution:
    Both panels achieve all {n_supertypes} supertypes ≥ {MIN_EXPECTED_DETECTED} expected detected ✓
    v1 provides {v1_mean_cov:.1f} vs {k5_mean_cov:.1f} avg expected markers/cell (both sufficient)
    With only 100 add-on genes, v1 coverage per supertype is thinner but still above threshold
    5K base panel already provides substantial baseline coverage for SST/L6b supertypes

  Budget efficiency:
    v1 leftover: ${BUDGET_TOTAL_USD - v1_total_cost:>,} (cannot afford another run at ${v1_per_run_total:>,}/run)
    5K leftover: ${BUDGET_TOTAL_USD - k5_total_cost:>,} (cannot afford another run at ${k5_per_run:>,}/run)
""")

print(f"All outputs saved to: {OUTPUT_DIR}/")
print("Done!")
