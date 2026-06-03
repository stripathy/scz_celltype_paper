"""
Supertype marker gene identification and Xenium panel overlap analysis.

Pipeline position: 5 of 7 (probe design)
Upstream: none (independent)
Downstream: none (informational, results in output/marker_analysis/)

Uses scanpy rank_genes_groups (Wilcoxon) to find markers at supertype level,
then checks how many are in the Xenium 5K Prime and v1 panels.

Two analyses:
1. Unrestricted markers (all genes) → check panel coverage
2. Panel-restricted markers → best separating genes within each panel
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import scanpy as sc

# ── Shared module imports ─────────────────────────────────────────────────
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "modules"))
from panel_utils import load_xenium_panels
from reference_utils import load_and_normalize_reference, subsample_by_group

# ── Configuration ──────────────────────────────────────────────────────────
BASE_DIR = os.path.expanduser("~/Github/SCZ_Xenium")
REF_PATH = os.path.join(BASE_DIR, "data", "reference",
                         "nicole_sea_ad_snrnaseq_reference.h5ad")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "marker_analysis")
os.makedirs(OUTPUT_DIR, exist_ok=True)

MAX_CELLS_PER_TYPE = 500
MIN_CELLS_PER_TYPE = 20
N_TOP_MARKERS = 10  # top markers per supertype to examine

# ── Load panel gene lists ─────────────────────────────────────────────────
print("Loading Xenium panel gene lists...")
panels = load_xenium_panels()
panel_5k_genes = panels["xenium_5k"]
panel_v1_genes = panels["xenium_v1"]

# ── Load reference data ──────────────────────────────────────────────────
adata = load_and_normalize_reference(REF_PATH, normalize=True, min_cells=10)
adata = subsample_by_group(adata, "Supertype",
                           max_cells=MAX_CELLS_PER_TYPE,
                           min_cells=MIN_CELLS_PER_TYPE)

# ══════════════════════════════════════════════════════════════════════════
# ANALYSIS 1: Unrestricted marker genes → panel overlap
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("ANALYSIS 1: UNRESTRICTED MARKER GENES")
print("="*80)

print("Running rank_genes_groups (Wilcoxon, 1-vs-rest)...")
t0 = time.time()
sc.tl.rank_genes_groups(adata, groupby="Supertype", method="wilcoxon",
                         n_genes=N_TOP_MARKERS, use_raw=False)
print(f"  Done in {time.time()-t0:.0f}s")

# Extract results
result = adata.uns["rank_genes_groups"]
groups = result["names"].dtype.names
records = []
for group in groups:
    for rank in range(N_TOP_MARKERS):
        gene = result["names"][rank][group]
        score = result["scores"][rank][group]
        pval = result["pvals_adj"][rank][group]
        logfc = result["logfoldchanges"][rank][group]
        records.append({
            "supertype": group,
            "rank": rank + 1,
            "gene": gene,
            "score": score,
            "logfoldchange": logfc,
            "pval_adj": pval,
            "in_5k": gene in panel_5k_genes,
            "in_v1": gene in panel_v1_genes,
        })

markers_df = pd.DataFrame(records)
markers_df.to_csv(os.path.join(OUTPUT_DIR, "unrestricted_markers_top10.csv"), index=False)

# Summary: per-supertype coverage
print(f"\nTop {N_TOP_MARKERS} markers per supertype — panel coverage:")
coverage_records = []
for group in sorted(groups):
    gdf = markers_df[markers_df["supertype"] == group]
    genes = list(gdf["gene"])
    n_in_5k = gdf["in_5k"].sum()
    n_in_v1 = gdf["in_v1"].sum()
    missing_5k = list(gdf[~gdf["in_5k"]]["gene"])
    missing_v1 = list(gdf[~gdf["in_v1"]]["gene"])
    coverage_records.append({
        "supertype": group,
        "markers": ", ".join(genes),
        "n_markers": len(genes),
        "n_in_5k": n_in_5k,
        "pct_in_5k": 100 * n_in_5k / len(genes),
        "missing_from_5k": ", ".join(missing_5k),
        "n_in_v1": n_in_v1,
        "pct_in_v1": 100 * n_in_v1 / len(genes),
        "missing_from_v1": ", ".join(missing_v1),
    })

coverage_df = pd.DataFrame(coverage_records)
coverage_df.to_csv(os.path.join(OUTPUT_DIR, "unrestricted_markers_coverage.csv"), index=False)

# Print aggregate stats
all_marker_genes = set(markers_df["gene"])
in_5k = all_marker_genes & panel_5k_genes
in_v1 = all_marker_genes & panel_v1_genes
print(f"\n  Unique marker genes across all supertypes: {len(all_marker_genes)}")
print(f"  In 5K panel: {len(in_5k)}/{len(all_marker_genes)} ({100*len(in_5k)/len(all_marker_genes):.1f}%)")
print(f"  In v1 panel: {len(in_v1)}/{len(all_marker_genes)} ({100*len(in_v1)/len(all_marker_genes):.1f}%)")

# Per-supertype: how many have ALL top markers covered?
full_5k = (coverage_df["n_in_5k"] == coverage_df["n_markers"]).sum()
full_v1 = (coverage_df["n_in_v1"] == coverage_df["n_markers"]).sum()
print(f"\n  Supertypes with ALL top-{N_TOP_MARKERS} markers in 5K: {full_5k}/{len(coverage_df)}")
print(f"  Supertypes with ALL top-{N_TOP_MARKERS} markers in v1: {full_v1}/{len(coverage_df)}")

# Distribution of coverage
for threshold in [10, 8, 5, 3, 1]:
    n_5k = (coverage_df["n_in_5k"] >= threshold).sum()
    n_v1 = (coverage_df["n_in_v1"] >= threshold).sum()
    print(f"  Supertypes with >= {threshold}/{N_TOP_MARKERS} markers:  5K={n_5k}  v1={n_v1}")

# ── SST supertype details ────────────────────────────────────────────────
print("\n" + "-"*60)
print("SST SUPERTYPE DETAILS (unrestricted markers)")
print("-"*60)
sst_coverage = coverage_df[coverage_df["supertype"].str.startswith("Sst")]
for _, row in sst_coverage.iterrows():
    print(f"\n  {row['supertype']}:")
    print(f"    Top markers: {row['markers']}")
    print(f"    In 5K: {row['n_in_5k']}/{row['n_markers']}  |  Missing: {row['missing_from_5k'] or 'none'}")
    print(f"    In v1: {row['n_in_v1']}/{row['n_markers']}  |  Missing: {row['missing_from_v1'] or 'none'}")

# ══════════════════════════════════════════════════════════════════════════
# ANALYSIS 2: Panel-restricted markers (what's the best you can do?)
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("ANALYSIS 2: PANEL-RESTRICTED MARKERS")
print("="*80)

for panel_name, panel_genes in [("5K_Prime", panel_5k_genes), ("v1_Brain", panel_v1_genes)]:
    print(f"\n--- {panel_name} panel ({len(panel_genes)} genes) ---")

    # Subset to panel genes
    panel_genes_in_ref = [g for g in adata.var_names if g in panel_genes]
    print(f"  Panel genes found in reference: {len(panel_genes_in_ref)}")

    adata_panel = adata[:, panel_genes_in_ref].copy()

    # Run marker analysis
    sc.tl.rank_genes_groups(adata_panel, groupby="Supertype", method="wilcoxon",
                             n_genes=min(N_TOP_MARKERS, len(panel_genes_in_ref)),
                             use_raw=False)

    result_panel = adata_panel.uns["rank_genes_groups"]
    panel_records = []
    for group in groups:
        for rank in range(min(N_TOP_MARKERS, len(panel_genes_in_ref))):
            gene = result_panel["names"][rank][group]
            score = result_panel["scores"][rank][group]
            pval = result_panel["pvals_adj"][rank][group]
            logfc = result_panel["logfoldchanges"][rank][group]
            panel_records.append({
                "supertype": group,
                "rank": rank + 1,
                "gene": gene,
                "score": score,
                "logfoldchange": logfc,
                "pval_adj": pval,
            })

    panel_markers_df = pd.DataFrame(panel_records)
    panel_markers_df.to_csv(
        os.path.join(OUTPUT_DIR, f"panel_restricted_markers_{panel_name}.csv"), index=False)

    # Compare: how similar are panel-restricted top markers to unrestricted top markers?
    overlap_records = []
    for group in sorted(groups):
        unrestricted = set(markers_df[markers_df["supertype"] == group]["gene"])
        restricted = set(panel_markers_df[panel_markers_df["supertype"] == group]["gene"])
        overlap = unrestricted & restricted
        overlap_records.append({
            "supertype": group,
            "unrestricted_markers": ", ".join(sorted(unrestricted)),
            "panel_markers": ", ".join(sorted(restricted)),
            "n_overlap": len(overlap),
            "overlap_genes": ", ".join(sorted(overlap)),
        })

    overlap_df = pd.DataFrame(overlap_records)
    overlap_df.to_csv(
        os.path.join(OUTPUT_DIR, f"marker_overlap_unrestricted_vs_{panel_name}.csv"), index=False)

    n_perfect = (overlap_df["n_overlap"] == N_TOP_MARKERS).sum()
    mean_overlap = overlap_df["n_overlap"].mean()
    print(f"  Overlap with unrestricted top-{N_TOP_MARKERS}: mean {mean_overlap:.1f} genes")
    print(f"  Perfect overlap (all {N_TOP_MARKERS}): {n_perfect}/{len(overlap_df)}")

    # SST details for this panel
    sst_panel = panel_markers_df[panel_markers_df["supertype"].str.startswith("Sst")]
    print(f"\n  SST supertypes — best markers in {panel_name}:")
    for st in sorted(sst_panel["supertype"].unique()):
        st_genes = list(sst_panel[sst_panel["supertype"] == st].head(5)["gene"])
        st_scores = list(sst_panel[sst_panel["supertype"] == st].head(5)["score"])
        gene_str = ", ".join([f"{g} ({s:.0f})" for g, s in zip(st_genes, st_scores)])
        print(f"    {st}: {gene_str}")

# ══════════════════════════════════════════════════════════════════════════
# ANALYSIS 3: Missing genes most needed for supertype resolution
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("MOST-NEEDED MISSING GENES FROM 5K PANEL")
print("="*80)

# Count how many supertypes need each missing gene (in their top N markers)
missing_counts = {}
for _, row in markers_df[~markers_df["in_5k"]].iterrows():
    gene = row["gene"]
    if gene not in missing_counts:
        missing_counts[gene] = {"count": 0, "types": [], "avg_rank": []}
    missing_counts[gene]["count"] += 1
    missing_counts[gene]["types"].append(row["supertype"])
    missing_counts[gene]["avg_rank"].append(row["rank"])

missing_summary = []
for gene, info in missing_counts.items():
    missing_summary.append({
        "gene": gene,
        "n_types_needing": info["count"],
        "avg_rank": np.mean(info["avg_rank"]),
        "types": ", ".join(info["types"][:10]) + ("..." if len(info["types"]) > 10 else ""),
    })

missing_summary_df = pd.DataFrame(missing_summary).sort_values("n_types_needing", ascending=False)
missing_summary_df.to_csv(os.path.join(OUTPUT_DIR, "missing_from_5k_ranked.csv"), index=False)

print(f"\n{'Gene':<20} {'# types':<10} {'Avg rank':<10} {'Types (sample)'}")
print("-" * 80)
for _, row in missing_summary_df.head(30).iterrows():
    print(f"{row['gene']:<20} {row['n_types_needing']:<10} {row['avg_rank']:<10.1f} {row['types'][:50]}")

# Same for v1
print("\n" + "="*80)
print("MOST-NEEDED MISSING GENES FROM v1 PANEL")
print("="*80)
missing_v1_counts = {}
for _, row in markers_df[~markers_df["in_v1"]].iterrows():
    gene = row["gene"]
    if gene not in missing_v1_counts:
        missing_v1_counts[gene] = {"count": 0}
    missing_v1_counts[gene]["count"] += 1

missing_v1_df = pd.DataFrame([
    {"gene": g, "n_types": info["count"], "in_5k": g in panel_5k_genes}
    for g, info in missing_v1_counts.items()
]).sort_values("n_types", ascending=False)
missing_v1_df.to_csv(os.path.join(OUTPUT_DIR, "missing_from_v1_ranked.csv"), index=False)

print(f"\n{'Gene':<20} {'# types':<10} {'In 5K?'}")
print("-" * 40)
for _, row in missing_v1_df.head(30).iterrows():
    print(f"{row['gene']:<20} {row['n_types']:<10} {'Yes' if row['in_5k'] else 'No'}")

# ── Final summary ────────────────────────────────────────────────────────
print("\n" + "="*80)
print("FINAL SUMMARY")
print("="*80)
print(f"\nReference: {adata.shape[0]} cells, {adata.obs['Supertype'].nunique()} supertypes")
print(f"Marker method: Wilcoxon rank-sum (1-vs-rest), top {N_TOP_MARKERS} per supertype")
print(f"\nUnique marker genes: {len(all_marker_genes)}")
print(f"  In 5K Prime ({len(panel_5k_genes)} genes): {len(in_5k)} ({100*len(in_5k)/len(all_marker_genes):.1f}%)")
print(f"  In v1 Brain ({len(panel_v1_genes)} genes): {len(in_v1)} ({100*len(in_v1)/len(all_marker_genes):.1f}%)")
print(f"\nSupertypes fully covered (all top-{N_TOP_MARKERS} in panel):")
print(f"  5K: {full_5k}/{len(coverage_df)} ({100*full_5k/len(coverage_df):.1f}%)")
print(f"  v1: {full_v1}/{len(coverage_df)} ({100*full_v1/len(coverage_df):.1f}%)")

# SST vulnerable types check
print("\nVulnerable SST supertypes (Sst_2, Sst_3, Sst_20, Sst_22, Sst_25):")
vulnerable = ["Sst_2", "Sst_3", "Sst_20", "Sst_22", "Sst_25"]
for st in vulnerable:
    row = coverage_df[coverage_df["supertype"] == st]
    if len(row) == 0:
        print(f"  {st}: NOT FOUND in reference")
        continue
    row = row.iloc[0]
    print(f"  {st}: {row['n_in_5k']}/{row['n_markers']} in 5K, {row['n_in_v1']}/{row['n_markers']} in v1")
    if row['missing_from_5k']:
        print(f"         Missing from 5K: {row['missing_from_5k']}")

print(f"\nAll results saved to: {OUTPUT_DIR}/")
print("Done!")
