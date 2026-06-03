"""
NS-Forest marker gene identification for SEA-AD supertypes.

Derives minimal combinatorial marker gene sets for each supertype using NS-Forest,
then checks overlap with Xenium 5K Prime and v1 panels.

NS-Forest (Aevermann et al., 2021) finds the minimal combination of genes whose
expression patterns are necessary and sufficient to distinguish each cell type.
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
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "nsforest")
os.makedirs(OUTPUT_DIR, exist_ok=True)

MAX_CELLS_PER_TYPE = 200  # Subsample for tractability
MIN_CELLS_PER_TYPE = 20   # Skip very rare types

# ── Load panel gene lists ─────────────────────────────────────────────────
print("Loading Xenium panel gene lists...")
panels = load_xenium_panels()
panel_5k_genes = panels["xenium_5k"]
panel_v1_genes = panels["xenium_v1"]

# ── Load and prepare reference data ───────────────────────────────────────
adata = load_and_normalize_reference(REF_PATH, normalize=True, min_cells=10)
adata = subsample_by_group(adata, "Supertype",
                           max_cells=MAX_CELLS_PER_TYPE,
                           min_cells=MIN_CELLS_PER_TYPE)
print(f"  After gene filtering: {adata.shape[1]} genes")

# ── NS-Forest preprocessing ──────────────────────────────────────────────
print("\nRunning NS-Forest preprocessing...")
import nsforest as nsf

# Compute medians and binary scores (required before NSForest)
nsf.pp.dendrogram(adata, "Supertype")

# ── Run NS-Forest ─────────────────────────────────────────────────────────
print("\nRunning NS-Forest (this will take a while)...")
t0 = time.time()
results = nsf.nsforesting.NSForest(
    adata,
    cluster_header="Supertype",
    n_trees=1000,
    n_top_genes=15,
    n_genes_eval=6,
    n_jobs=-1,
    save=True,
    save_supplementary=True,
    output_folder=OUTPUT_DIR,
    outputfilename_prefix="seaad_supertype_"
)
elapsed = time.time() - t0
print(f"  NS-Forest completed in {elapsed/60:.1f} minutes")

# ── Analyze results ──────────────────────────────────────────────────────
print("\n" + "="*80)
print("NS-FOREST RESULTS SUMMARY")
print("="*80)

print(f"\nResult columns: {list(results.columns)}")
print(f"\nResults shape: {results.shape}")
print(f"\nFirst few rows:")
print(results.head(10).to_string())

# Save results
results.to_csv(os.path.join(OUTPUT_DIR, "nsforest_results_all_supertypes.csv"), index=True)

# Extract all unique marker genes
if 'NSForest_markers' in results.columns:
    marker_col = 'NSForest_markers'
elif 'clusterName' in results.columns:
    # Try to find the markers column
    marker_col = [c for c in results.columns if 'marker' in c.lower()]
    marker_col = marker_col[0] if marker_col else None

print(f"\nMarker column: {marker_col}")

# Collect all markers
all_markers = set()
markers_by_type = {}
for idx, row in results.iterrows():
    cluster = idx if isinstance(idx, str) else row.get('clusterName', idx)
    markers = row[marker_col] if marker_col else []
    if isinstance(markers, str):
        markers = [m.strip() for m in markers.split(",")]
    elif isinstance(markers, (list, tuple)):
        markers = list(markers)
    else:
        markers = []
    markers_by_type[cluster] = markers
    all_markers.update(markers)

print(f"\nTotal unique NS-Forest markers: {len(all_markers)}")
print(f"Average markers per supertype: {np.mean([len(v) for v in markers_by_type.values()]):.1f}")

# ── Panel overlap analysis ────────────────────────────────────────────────
print("\n" + "="*80)
print("PANEL OVERLAP ANALYSIS")
print("="*80)

# Xenium 5K overlap
in_5k = all_markers & panel_5k_genes
not_in_5k = all_markers - panel_5k_genes
print(f"\n--- Xenium 5K Prime Panel ---")
print(f"  NS-Forest markers in 5K panel: {len(in_5k)}/{len(all_markers)} ({100*len(in_5k)/len(all_markers):.1f}%)")
print(f"  Missing from 5K panel: {len(not_in_5k)}")

# Xenium v1 overlap
in_v1 = all_markers & panel_v1_genes
not_in_v1 = all_markers - panel_v1_genes
print(f"\n--- Xenium v1 Brain Panel ---")
print(f"  NS-Forest markers in v1 panel: {len(in_v1)}/{len(all_markers)} ({100*len(in_v1)/len(all_markers):.1f}%)")
print(f"  Missing from v1 panel: {len(not_in_v1)}")

# Per-supertype coverage
print("\n" + "="*80)
print("PER-SUPERTYPE PANEL COVERAGE")
print("="*80)

coverage_records = []
for cluster, markers in sorted(markers_by_type.items()):
    n_markers = len(markers)
    n_in_5k = len(set(markers) & panel_5k_genes)
    n_in_v1 = len(set(markers) & panel_v1_genes)
    missing_5k = sorted(set(markers) - panel_5k_genes)
    missing_v1 = sorted(set(markers) - panel_v1_genes)

    coverage_records.append({
        "supertype": cluster,
        "n_markers": n_markers,
        "markers": ", ".join(markers),
        "n_in_5k": n_in_5k,
        "pct_in_5k": 100 * n_in_5k / n_markers if n_markers > 0 else 0,
        "missing_5k": ", ".join(missing_5k) if missing_5k else "",
        "n_in_v1": n_in_v1,
        "pct_in_v1": 100 * n_in_v1 / n_markers if n_markers > 0 else 0,
        "missing_v1": ", ".join(missing_v1) if missing_v1 else "",
        "fully_covered_5k": n_in_5k == n_markers,
        "fully_covered_v1": n_in_v1 == n_markers,
    })

coverage_df = pd.DataFrame(coverage_records)
coverage_df.to_csv(os.path.join(OUTPUT_DIR, "nsforest_panel_coverage.csv"), index=False)

# Print summary
fully_covered_5k = coverage_df["fully_covered_5k"].sum()
fully_covered_v1 = coverage_df["fully_covered_v1"].sum()
print(f"\nSupertypes with ALL markers in 5K panel: {fully_covered_5k}/{len(coverage_df)} ({100*fully_covered_5k/len(coverage_df):.1f}%)")
print(f"Supertypes with ALL markers in v1 panel: {fully_covered_v1}/{len(coverage_df)} ({100*fully_covered_v1/len(coverage_df):.1f}%)")

# Focus on SST types
print("\n" + "="*80)
print("SST SUPERTYPE DETAILS")
print("="*80)
sst_df = coverage_df[coverage_df["supertype"].str.startswith("Sst")]
for _, row in sst_df.iterrows():
    print(f"\n  {row['supertype']}:")
    print(f"    Markers: {row['markers']}")
    print(f"    In 5K: {row['n_in_5k']}/{row['n_markers']} | Missing: {row['missing_5k'] or 'none'}")
    print(f"    In v1: {row['n_in_v1']}/{row['n_markers']} | Missing: {row['missing_v1'] or 'none'}")

# Print missing genes sorted by frequency (most needed)
print("\n" + "="*80)
print("MOST-NEEDED MISSING GENES (not in 5K panel)")
print("="*80)
missing_gene_counts = {}
for cluster, markers in markers_by_type.items():
    for m in markers:
        if m not in panel_5k_genes:
            missing_gene_counts[m] = missing_gene_counts.get(m, 0) + 1

if missing_gene_counts:
    missing_sorted = sorted(missing_gene_counts.items(), key=lambda x: -x[1])
    print(f"\n{'Gene':<20} {'# types needing it':<20}")
    print("-" * 40)
    for gene, count in missing_sorted[:30]:
        print(f"{gene:<20} {count:<20}")

print("\n" + "="*80)
print("MOST-NEEDED MISSING GENES (not in v1 panel)")
print("="*80)
missing_v1_counts = {}
for cluster, markers in markers_by_type.items():
    for m in markers:
        if m not in panel_v1_genes:
            missing_v1_counts[m] = missing_v1_counts.get(m, 0) + 1

if missing_v1_counts:
    missing_sorted_v1 = sorted(missing_v1_counts.items(), key=lambda x: -x[1])
    print(f"\n{'Gene':<20} {'# types needing it':<20}")
    print("-" * 40)
    for gene, count in missing_sorted_v1[:30]:
        print(f"{gene:<20} {count:<20}")

# Save missing genes list
missing_df = pd.DataFrame([
    {"gene": g, "n_types_needing": c, "in_5k": g in panel_5k_genes, "in_v1": g in panel_v1_genes}
    for g, c in sorted(missing_gene_counts.items(), key=lambda x: -x[1])
])
missing_df.to_csv(os.path.join(OUTPUT_DIR, "nsforest_missing_from_5k.csv"), index=False)

print(f"\n\nAll results saved to {OUTPUT_DIR}")
print("Done!")
