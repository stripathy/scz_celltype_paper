"""
MERSCOPE 4K/250 panel assessment for SST supertype marker coverage.

Loads Fang et al. MERSCOPE data (250 and 4000 gene panels, human MTG),
computes per-gene detection rates, and compares all available panels:
  - MERSCOPE 250-gene brain panel
  - MERSCOPE 4000-gene brain panel
  - Xenium 5K Prime pan-tissue panel
  - Xenium v1 Brain panel (266 genes)

Then calibrates snRNAseq → MERSCOPE detection efficiency using the 4K data
(much better powered than the 140-gene MERFISH calibration).
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt

# ── Shared module imports ─────────────────────────────────────────────────
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "modules"))
from panel_utils import load_xenium_panels

# ── Configuration ──────────────────────────────────────────────────────────
BASE_DIR = os.path.expanduser("~/Github/SCZ_Xenium")
MERSCOPE_DIR = os.path.join(BASE_DIR, "data", "merscope_4k_probe_testing")
REF_PATH = os.path.join(BASE_DIR, "data", "reference",
                         "nicole_sea_ad_snrnaseq_reference.h5ad")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "marker_analysis")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Load Xenium panels ────────────────────────────────────────────────────
print("Loading Xenium panel gene lists...")
panels = load_xenium_panels()
panel_5k_genes = panels["xenium_5k"]
panel_v1_genes = panels["xenium_v1"]


# ══════════════════════════════════════════════════════════════════════════
# PART 1: Load MERSCOPE data and extract panel gene lists
# ══════════════════════════════════════════════════════════════════════════

def load_merscope_sample(prefix):
    """Load a MERSCOPE sample from genes/features/matrix CSVs into AnnData."""
    genes_df = pd.read_csv(f"{prefix}.genes.csv")
    gene_names = list(genes_df["name"])

    features_df = pd.read_csv(f"{prefix}.features.csv", index_col=0)
    n_cells = len(features_df)
    n_genes = len(gene_names)

    matrix_df = pd.read_csv(f"{prefix}.matrix.csv")
    # Sparse triplet format: row (gene index, 1-based), col (cell index, 1-based), val
    rows = matrix_df["row"].values - 1  # to 0-based
    cols = matrix_df["col"].values - 1
    vals = matrix_df["val"].values

    X = sparse.csr_matrix((vals, (cols, rows)), shape=(n_cells, n_genes))

    adata = ad.AnnData(X=X, obs=features_df)
    adata.var_names = gene_names

    return adata


print("Loading MERSCOPE datasets...")

# Discover all available samples
import glob
gene_files = sorted(glob.glob(os.path.join(MERSCOPE_DIR, "*.genes.csv")))
samples = {}
for gf in gene_files:
    prefix = gf.replace(".genes.csv", "")
    basename = os.path.basename(prefix)
    parts = basename.split(".")
    # e.g. H18.06.006.MTG.4000.expand.rep1
    donor = ".".join(parts[:3])
    region = parts[3]
    panel_size = parts[4]
    expand = parts[5]
    rep = parts[6]
    key = f"{donor}.{region}.{panel_size}.{expand}.{rep}"
    samples[key] = {
        "prefix": prefix,
        "donor": donor,
        "region": region,
        "panel_size": int(panel_size),
        "expand": expand,
        "rep": rep,
    }

print(f"  Found {len(samples)} samples:")
for key, info in sorted(samples.items()):
    print(f"    {key}: {info['region']} {info['panel_size']}-gene {info['expand']} {info['rep']}")

# Load all samples
adatas = {}
for key, info in sorted(samples.items()):
    t0 = time.time()
    adata = load_merscope_sample(info["prefix"])
    adata.obs["sample_key"] = key
    adata.obs["donor"] = info["donor"]
    adata.obs["region"] = info["region"]
    adata.obs["panel_size"] = info["panel_size"]
    adatas[key] = adata
    print(f"    {key}: {adata.shape[0]} cells x {adata.shape[1]} genes ({time.time()-t0:.1f}s)")

# Extract gene lists
merscope_250_genes = set(adatas[[k for k in adatas if ".250." in k][0]].var_names)
merscope_4k_genes = set(adatas[[k for k in adatas if ".4000." in k][0]].var_names)

print(f"\n  MERSCOPE 250-gene panel: {len(merscope_250_genes)} genes")
print(f"  MERSCOPE 4000-gene panel: {len(merscope_4k_genes)} genes")

# ══════════════════════════════════════════════════════════════════════════
# PART 2: Panel overlap comparison
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("PANEL OVERLAP MATRIX")
print("="*80)

panels = {
    "MERSCOPE_250": merscope_250_genes,
    "MERSCOPE_4K": merscope_4k_genes,
    "Xenium_5K": panel_5k_genes,
    "Xenium_v1": panel_v1_genes,
}

# Pairwise overlap
print(f"\n{'':>15}", end="")
for name in panels:
    print(f"  {name:>14}", end="")
print()

for name1, genes1 in panels.items():
    print(f"{name1:>15}", end="")
    for name2, genes2 in panels.items():
        overlap = len(genes1 & genes2)
        print(f"  {overlap:>14}", end="")
    print()

# ══════════════════════════════════════════════════════════════════════════
# PART 3: Load supertype markers and check coverage per panel
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("SUPERTYPE MARKER COVERAGE PER PANEL")
print("="*80)

# Load the unrestricted markers from previous analysis
markers_path = os.path.join(OUTPUT_DIR, "unrestricted_markers_top10.csv")
markers_df = pd.read_csv(markers_path)
all_supertypes = sorted(markers_df["supertype"].unique())

coverage_records = []
for st in all_supertypes:
    st_markers = list(markers_df[markers_df["supertype"] == st]["gene"])
    n = len(st_markers)
    record = {"supertype": st, "n_markers": n, "markers": ", ".join(st_markers)}
    for panel_name, panel_genes in panels.items():
        n_in = len(set(st_markers) & panel_genes)
        missing = sorted(set(st_markers) - panel_genes)
        record[f"n_in_{panel_name}"] = n_in
        record[f"pct_{panel_name}"] = 100 * n_in / n if n > 0 else 0
        record[f"missing_{panel_name}"] = ", ".join(missing)
    coverage_records.append(record)

coverage_df = pd.DataFrame(coverage_records)
coverage_df.to_csv(os.path.join(OUTPUT_DIR, "supertype_marker_coverage_all_panels.csv"), index=False)

# Summary stats per panel
print(f"\nTop-10 marker coverage across {len(all_supertypes)} supertypes:")
print(f"\n{'Panel':>15} {'All 10':>8} {'>=8':>8} {'>=5':>8} {'>=3':>8} {'>=1':>8} {'Median':>8}")
print("-" * 63)
for panel_name in panels:
    col = f"n_in_{panel_name}"
    vals = coverage_df[col]
    all10 = (vals == 10).sum()
    ge8 = (vals >= 8).sum()
    ge5 = (vals >= 5).sum()
    ge3 = (vals >= 3).sum()
    ge1 = (vals >= 1).sum()
    med = vals.median()
    print(f"{panel_name:>15} {all10:>8} {ge8:>8} {ge5:>8} {ge3:>8} {ge1:>8} {med:>8.1f}")

# SST-specific
print(f"\nSST SUPERTYPE COVERAGE:")
sst_df = coverage_df[coverage_df["supertype"].str.startswith("Sst")]
print(f"\n{'Supertype':<15}", end="")
for panel_name in panels:
    print(f"  {panel_name:>14}", end="")
print()
for _, row in sst_df.iterrows():
    print(f"{row['supertype']:<15}", end="")
    for panel_name in panels:
        n = row[f"n_in_{panel_name}"]
        print(f"  {n:>11}/10", end="")
    print()

# Vulnerable SST types detail
vulnerable = ["Sst_2", "Sst_3", "Sst_20", "Sst_22", "Sst_25"]
print(f"\nVulnerable SST details:")
for st in vulnerable:
    row = coverage_df[coverage_df["supertype"] == st]
    if len(row) == 0:
        continue
    row = row.iloc[0]
    print(f"\n  {st}: {row['markers']}")
    for panel_name in panels:
        n = row[f"n_in_{panel_name}"]
        missing = row[f"missing_{panel_name}"]
        print(f"    {panel_name:>15}: {n}/10  missing: {missing or 'none'}")

# ══════════════════════════════════════════════════════════════════════════
# PART 4: Detection rate calibration (snRNAseq vs MERSCOPE 4K)
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("DETECTION CALIBRATION: snRNAseq vs MERSCOPE 4K")
print("="*80)

# Concatenate all 4K MTG replicates for better power
mtg_4k_keys = [k for k in adatas if ".MTG.4000." in k]
print(f"\nUsing {len(mtg_4k_keys)} MTG 4K replicates: {mtg_4k_keys}")

merscope_concat = ad.concat([adatas[k] for k in mtg_4k_keys])
print(f"  Combined: {merscope_concat.shape[0]} cells x {merscope_concat.shape[1]} genes")

# Check cell type labels
print(f"\n  Cell type labels available:")
print(f"    cluster_L1: {merscope_concat.obs['cluster_L1'].nunique()} types — {sorted(merscope_concat.obs['cluster_L1'].unique())[:10]}")
print(f"    cluster_L2: {merscope_concat.obs['cluster_L2'].nunique()} types — {sorted(merscope_concat.obs['cluster_L2'].dropna().unique())[:10]}")
print(f"    cluster_L3: {merscope_concat.obs['cluster_L3'].nunique()} types — {sorted(merscope_concat.obs['cluster_L3'].dropna().unique())[:10]}")

# Use L2 labels for calibration (closest to subclass level)
# Show all L2 labels
l2_counts = merscope_concat.obs["cluster_L2"].value_counts()
print(f"\n  L2 cluster distribution:")
for ct, n in l2_counts.items():
    print(f"    {ct}: {n}")

# Load snRNAseq reference
print("\nLoading snRNAseq reference...")
snrna = ad.read_h5ad(REF_PATH)
print(f"  {snrna.shape[0]} cells x {snrna.shape[1]} genes")

# Find shared genes between MERSCOPE 4K and snRNAseq
shared_genes = sorted(set(merscope_concat.var_names) & set(snrna.var_names))
print(f"  Shared genes: {len(shared_genes)}")

# Subsample snRNAseq
np.random.seed(42)
keep_idx = []
for st in snrna.obs["Subclass"].unique():
    idx = np.where(snrna.obs["Subclass"] == st)[0]
    if len(idx) > 500:
        idx = np.random.choice(idx, 500, replace=False)
    keep_idx.extend(idx)
snrna_sub = snrna[sorted(keep_idx)][:, shared_genes].copy()
if sparse.issparse(snrna_sub.X):
    snrna_X = snrna_sub.X.toarray()
else:
    snrna_X = np.array(snrna_sub.X)
snrna_subclass = snrna_sub.obs["Subclass"].values
print(f"  snRNAseq subsampled: {snrna_X.shape[0]} cells")

# MERSCOPE expression matrix
merscope_X = merscope_concat[:, shared_genes].X
if sparse.issparse(merscope_X):
    merscope_X = merscope_X.toarray()

# Map MERSCOPE L2 labels to SEA-AD subclass (rough mapping)
# First, let's see the L2 labels
merscope_l2 = merscope_concat.obs["cluster_L2"].values

# Build per-gene global detection stats (not per cell type — simpler, more robust)
print("\nComputing global per-gene detection rates...")
n_cells_snrna = snrna_X.shape[0]
n_cells_merscope = merscope_X.shape[0]

gene_stats = []
for j, gene in enumerate(shared_genes):
    frac_snrna = np.mean(snrna_X[:, j] > 0)
    mean_snrna = np.mean(snrna_X[:, j])
    frac_merscope = np.mean(merscope_X[:, j] > 0)
    mean_merscope = np.mean(merscope_X[:, j])

    gene_stats.append({
        "gene": gene,
        "frac_snrna": frac_snrna,
        "mean_snrna": mean_snrna,
        "frac_merscope_4k": frac_merscope,
        "mean_merscope_4k": mean_merscope,
        "in_xenium_5k": gene in panel_5k_genes,
        "in_xenium_v1": gene in panel_v1_genes,
        "in_merscope_250": gene in merscope_250_genes,
    })

gene_stats_df = pd.DataFrame(gene_stats).sort_values("frac_merscope_4k", ascending=False)
gene_stats_df["detection_efficiency"] = (
    gene_stats_df["frac_merscope_4k"] / gene_stats_df["frac_snrna"].clip(lower=0.001)
)
gene_stats_df.to_csv(os.path.join(OUTPUT_DIR, "snrnaseq_vs_merscope4k_detection.csv"), index=False)

# Correlation
valid = gene_stats_df[gene_stats_df["frac_snrna"] > 0.01]
r_frac, _ = pearsonr(valid["frac_snrna"], valid["frac_merscope_4k"])
rho_frac, _ = spearmanr(valid["frac_snrna"], valid["frac_merscope_4k"])
r_mean, _ = pearsonr(np.log1p(valid["mean_snrna"]), np.log1p(valid["mean_merscope_4k"]))

print(f"\nCalibration (n={len(valid)} genes with >1% snRNAseq detection):")
print(f"  Detection fraction: Pearson r={r_frac:.3f}, Spearman ρ={rho_frac:.3f}")
print(f"  Log mean counts: Pearson r={r_mean:.3f}")

median_eff = valid["detection_efficiency"].median()
iqr_25 = valid["detection_efficiency"].quantile(0.25)
iqr_75 = valid["detection_efficiency"].quantile(0.75)
print(f"\n  Detection efficiency (MERSCOPE/snRNAseq):")
print(f"    Median: {median_eff:.3f}")
print(f"    IQR: {iqr_25:.3f} - {iqr_75:.3f}")

# ── Predict spatial detectability for SST markers ─────────────────────────
print("\n" + "="*80)
print("PREDICTED SPATIAL DETECTABILITY OF SST MARKERS (MERSCOPE 4K calibration)")
print("="*80)

# Get snRNAseq per-supertype detection for SST markers
sst_markers_all = markers_df[markers_df["supertype"].str.startswith("Sst")]
sst_genes = sorted(sst_markers_all["gene"].unique())
sst_genes_in_ref = [g for g in sst_genes if g in snrna.var_names]

# Use MERSCOPE 4K calibration where available, else use median efficiency
calibration_lookup = dict(zip(gene_stats_df["gene"], gene_stats_df["detection_efficiency"]))

sst_predict_records = []
for st in sorted(sst_markers_all["supertype"].unique()):
    mask = snrna.obs["Supertype"] == st
    if mask.sum() < 10:
        continue
    st_genes_avail = [g for g in sst_genes_in_ref if g in snrna.var_names]
    X_st = snrna[mask][:, st_genes_avail].X
    if sparse.issparse(X_st):
        X_st = X_st.toarray()

    st_markers = list(sst_markers_all[sst_markers_all["supertype"] == st]["gene"])

    for gene in st_markers:
        if gene not in st_genes_avail:
            continue
        j = st_genes_avail.index(gene)
        frac_snrna = np.mean(X_st[:, j] > 0)

        if gene in calibration_lookup and not np.isnan(calibration_lookup[gene]):
            eff = calibration_lookup[gene]
            source = "MERSCOPE_4K"
        else:
            eff = median_eff
            source = "estimated"

        predicted = min(frac_snrna * eff, 1.0)

        sst_predict_records.append({
            "supertype": st,
            "gene": gene,
            "frac_snrna": frac_snrna,
            "detection_efficiency": eff,
            "predicted_spatial": predicted,
            "calibration_source": source,
            "in_merscope_4k": gene in merscope_4k_genes,
            "in_merscope_250": gene in merscope_250_genes,
            "in_xenium_5k": gene in panel_5k_genes,
            "in_xenium_v1": gene in panel_v1_genes,
        })

sst_predict_df = pd.DataFrame(sst_predict_records)
sst_predict_df.to_csv(os.path.join(OUTPUT_DIR, "sst_marker_detectability_merscope4k_calibration.csv"), index=False)

# Print vulnerable types
for st in vulnerable:
    st_data = sst_predict_df[sst_predict_df["supertype"] == st].sort_values("predicted_spatial", ascending=False)
    if len(st_data) == 0:
        continue
    print(f"\n  {st}:")
    for _, row in st_data.iterrows():
        panels_str = []
        if row["in_merscope_4k"]: panels_str.append("M4K")
        if row["in_merscope_250"]: panels_str.append("M250")
        if row["in_xenium_5k"]: panels_str.append("X5K")
        if row["in_xenium_v1"]: panels_str.append("Xv1")
        panel_tag = ",".join(panels_str) if panels_str else "NONE"

        print(f"    {row['gene']:<15} snRNA={row['frac_snrna']:.3f}  "
              f"eff={row['detection_efficiency']:.3f} ({row['calibration_source']:<12})  "
              f"spatial≈{row['predicted_spatial']:.3f}  [{panel_tag}]")

# ══════════════════════════════════════════════════════════════════════════
# PART 5: Plot calibration
# ══════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 3, figsize=(21, 6))

# Panel 1: snRNAseq vs MERSCOPE detection fraction
ax = axes[0]
colors = ['#e41a1c' if row['in_xenium_5k'] else '#377eb8'
          for _, row in valid.iterrows()]
ax.scatter(valid["frac_snrna"], valid["frac_merscope_4k"], s=8, alpha=0.4,
           c=colors, rasterized=True)
ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, lw=1)
ax.set_xlabel("snRNAseq fraction detected", fontsize=16)
ax.set_ylabel("MERSCOPE 4K fraction detected", fontsize=16)
ax.set_title(f"Detection calibration (n={len(valid)} genes)\nr={r_frac:.3f}, ρ={rho_frac:.3f}", fontsize=18)
ax.tick_params(labelsize=14)
# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#e41a1c', markersize=8, label='In Xenium 5K'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#377eb8', markersize=8, label='Not in Xenium 5K'),
]
ax.legend(handles=legend_elements, fontsize=12)

# Panel 2: Detection efficiency histogram
ax = axes[1]
ax.hist(valid["detection_efficiency"].clip(upper=5), bins=50, color='seagreen', edgecolor='white')
ax.axvline(1.0, color='k', linestyle='--', alpha=0.5, label='1:1')
ax.axvline(median_eff, color='red', linestyle='-', label=f'Median={median_eff:.2f}')
ax.set_xlabel("Detection efficiency (MERSCOPE/snRNAseq)", fontsize=16)
ax.set_ylabel("Number of genes", fontsize=16)
ax.set_title("Per-gene detection efficiency\n(4K-gene calibration)", fontsize=18)
ax.tick_params(labelsize=14)
ax.legend(fontsize=14)

# Panel 3: Per-panel marker coverage for all supertypes
ax = axes[2]
panel_names_short = ["MERSCOPE_250", "MERSCOPE_4K", "Xenium_5K", "Xenium_v1"]
panel_colors = ["#4daf4a", "#984ea3", "#e41a1c", "#377eb8"]
x_pos = np.arange(len(panel_names_short))
medians = [coverage_df[f"n_in_{p}"].median() for p in panel_names_short]
means = [coverage_df[f"n_in_{p}"].mean() for p in panel_names_short]
q25 = [coverage_df[f"n_in_{p}"].quantile(0.25) for p in panel_names_short]
q75 = [coverage_df[f"n_in_{p}"].quantile(0.75) for p in panel_names_short]

bars = ax.bar(x_pos, means, color=panel_colors, edgecolor='white', alpha=0.8)
ax.errorbar(x_pos, means, yerr=[np.array(means)-np.array(q25), np.array(q75)-np.array(means)],
            fmt='none', color='black', capsize=5, capthick=2)
ax.set_xticks(x_pos)
ax.set_xticklabels(["MERSCOPE\n250", "MERSCOPE\n4K", "Xenium\n5K", "Xenium\nv1"], fontsize=14)
ax.set_ylabel("Markers covered (of top 10)", fontsize=16)
ax.set_title("Mean marker coverage\nper supertype", fontsize=18)
ax.tick_params(labelsize=14)
ax.set_ylim(0, 10)
ax.axhline(10, color='gray', linestyle=':', alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "panel_comparison_merscope4k.png"), dpi=150, bbox_inches='tight')
plt.close()

# ══════════════════════════════════════════════════════════════════════════
# PART 6: Genes with poor spatial detection efficiency
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("GENES WITH POOR SPATIAL DETECTION (efficiency < 0.5)")
print("="*80)

poor_detect = gene_stats_df[
    (gene_stats_df["detection_efficiency"] < 0.5) &
    (gene_stats_df["frac_snrna"] > 0.1)
].sort_values("detection_efficiency")

print(f"\n{'Gene':<15} {'snRNA det':<10} {'MERSCOPE':<10} {'Efficiency':<12} {'In 5K':<8}")
print("-" * 55)
for _, row in poor_detect.head(30).iterrows():
    print(f"{row['gene']:<15} {row['frac_snrna']:<10.3f} {row['frac_merscope_4k']:<10.3f} "
          f"{row['detection_efficiency']:<12.3f} {'Yes' if row['in_xenium_5k'] else 'No':<8}")

# ══════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("FINAL SUMMARY")
print("="*80)

# Panel sizes and brain relevance
print(f"\nPanel comparison:")
print(f"  {'Panel':<20} {'Total genes':<12} {'In snRNAseq ref':<16} {'Brain-relevant*'}")
print(f"  {'':<20} {'':<12} {'':<16} {'(>5% det in MTG)'}")
print(f"  {'-'*64}")
for name, genes in panels.items():
    in_ref = len(genes & set(snrna.var_names))
    # For brain-relevant, use snRNAseq detection > 5%
    # (we computed this for shared genes)
    if name in ["MERSCOPE_250", "MERSCOPE_4K"]:
        brain_rel = len([g for g in genes if g in gene_stats_df["gene"].values and
                        gene_stats_df[gene_stats_df["gene"]==g]["frac_snrna"].values[0] > 0.05])
    else:
        brain_rel = "—"
    print(f"  {name:<20} {len(genes):<12} {in_ref:<16} {brain_rel}")

print(f"\nMarker coverage (top-10 Wilcoxon markers per supertype):")
print(f"  {'Panel':<20} {'Mean coverage':<15} {'Median':<10} {'All 10':<10} {'>=5':<10}")
print(f"  {'-'*55}")
for panel_name in panels:
    col = f"n_in_{panel_name}"
    vals = coverage_df[col]
    print(f"  {panel_name:<20} {vals.mean():<15.1f} {vals.median():<10.1f} "
          f"{(vals==10).sum():<10} {(vals>=5).sum():<10}")

print(f"\nCalibration: snRNAseq → MERSCOPE 4K (n={len(valid)} shared genes)")
print(f"  Pearson r = {r_frac:.3f}, Spearman ρ = {rho_frac:.3f}")
print(f"  Median detection efficiency = {median_eff:.3f}")

print(f"\nAll results saved to {OUTPUT_DIR}/")
print("Done!")
