#!/usr/bin/env python3
"""
Cross-platform marker adequacy assessment for supertype classification.

Computes within-subclass Wilcoxon DE markers using the snRNAseq reference,
then checks coverage across 5 spatial panels:
  1. Xenium v1 Brain (~266 genes)
  2. SEA-AD MERFISH (~180 genes)
  3. MERSCOPE 250 (~250 genes)
  4. MERSCOPE 4K (~4000 genes)
  5. Xenium 5K Prime (~5000 genes)

Also assesses layer specificity as an independent axis of classification
confidence.

Usage:
    python3 -u cross_platform_marker_adequacy.py
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    BASE_DIR, SNRNASEQ_REF_PATH, PANEL_V1_PATH, PANEL_5K_PATH,
    CRUMBLR_DIR,
)

sys.path.insert(0, os.path.join(BASE_DIR, "code", "modules"))
from reference_utils import load_and_normalize_reference, subsample_by_group

OUTPUT_DIR = os.path.join(BASE_DIR, "output", "marker_analysis")
DEPTH_FILE = os.path.join(BASE_DIR, "output", "presentation",
                           "median_depth_supertype.csv")
DENSITY_RESULTS = os.path.join(BASE_DIR, "output", "density_analysis",
                                "density_results_supertype_cortical.csv")
CRUMBLR_RESULTS = os.path.join(CRUMBLR_DIR, "crumblr_results_supertype.csv")

# MERSCOPE / MERFISH paths
MERFISH_REF = os.path.join(BASE_DIR, "data", "reference",
                            "SEAAD_MTG_MERFISH.2024-12-11.h5ad")
MERSCOPE_DIR = os.path.join(BASE_DIR, "data", "merscope_4k_probe_testing")
MERSCOPE_250_GENES = os.path.join(
    MERSCOPE_DIR, "H18.06.006.MTG.250.expand.rep1.genes.csv")
MERSCOPE_4K_GENES = os.path.join(
    MERSCOPE_DIR, "H18.06.006.MTG.4000.expand.rep1.genes.csv")

N_MARKERS = 50
MAX_CELLS_PER_TYPE = 500
MIN_CELLS_PER_TYPE = 20


def load_all_panels():
    """Load gene sets for all 5 spatial panels."""
    panels = {}

    # 1. Xenium v1
    v1 = pd.read_csv(PANEL_V1_PATH)
    panels["Xenium_v1"] = set(v1["Genes"].str.strip())

    # 2. Xenium 5K
    x5k = pd.read_csv(PANEL_5K_PATH)
    panels["Xenium_5K"] = set(x5k["gene_name"].str.strip())

    # 3. MERFISH
    import anndata as ad
    merfish = ad.read_h5ad(MERFISH_REF, backed="r")
    panels["MERFISH_180"] = set(merfish.var_names)
    merfish.file.close()

    # 4. MERSCOPE 250
    m250 = pd.read_csv(MERSCOPE_250_GENES)
    panels["MERSCOPE_250"] = set(m250.iloc[:, 0].str.strip())

    # 5. MERSCOPE 4K
    m4k = pd.read_csv(MERSCOPE_4K_GENES)
    panels["MERSCOPE_4K"] = set(m4k.iloc[:, 0].str.strip())

    for name, genes in panels.items():
        print(f"  {name}: {len(genes)} genes")

    return panels


def extract_markers(adata_sub, groupby, n_genes):
    """Run Wilcoxon rank-sum and extract markers."""
    sc.tl.rank_genes_groups(adata_sub, groupby=groupby,
                             method="wilcoxon", n_genes=n_genes,
                             use_raw=False)
    result = adata_sub.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    records = []
    for group in groups:
        for rank in range(n_genes):
            try:
                records.append({
                    "group": group,
                    "rank": rank + 1,
                    "gene": result["names"][rank][group],
                    "wilcoxon_score": float(result["scores"][rank][group]),
                    "logfoldchange": float(result["logfoldchanges"][rank][group]),
                    "pval_adj": float(result["pvals_adj"][rank][group]),
                })
            except (IndexError, KeyError):
                break
    return pd.DataFrame(records)


def infer_subclass(supertype):
    """Infer subclass from supertype name."""
    for prefix in ["L2/3 IT", "L4 IT", "L5 ET", "L5 IT", "L5/6 NP",
                    "L6 CT", "L6 IT Car3", "L6 IT", "L6b",
                    "Sst Chodl", "Micro-PVM"]:
        if supertype.startswith(prefix):
            return prefix
    parts = supertype.rsplit("_", 1)
    if len(parts) == 2 and parts[1].replace("-SEAAD", "").isdigit():
        return parts[0]
    return supertype


def main():
    t0 = time.time()
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ── Load panels ──
    print("Loading spatial panels...")
    panels = load_all_panels()
    panel_names = ["Xenium_v1", "MERFISH_180", "MERSCOPE_250",
                    "MERSCOPE_4K", "Xenium_5K"]

    # ── Load reference ──
    print("\nLoading snRNAseq reference...")
    adata = load_and_normalize_reference(SNRNASEQ_REF_PATH,
                                          normalize=True, min_cells=10)
    adata = subsample_by_group(adata, "Supertype",
                                max_cells=MAX_CELLS_PER_TYPE,
                                min_cells=MIN_CELLS_PER_TYPE)
    for col in ["Supertype", "Subclass"]:
        if hasattr(adata.obs[col], "cat"):
            adata.obs[col] = adata.obs[col].cat.remove_unused_categories()
        else:
            adata.obs[col] = pd.Categorical(adata.obs[col])
    print(f"  {adata.n_obs} cells, {adata.obs['Supertype'].nunique()} supertypes")

    # ── Load depth data ──
    depth_df = pd.read_csv(DEPTH_FILE)
    depth_df["subclass"] = depth_df["supertype"].apply(infer_subclass)

    # ── Load significance results ──
    crumblr_df = pd.read_csv(CRUMBLR_RESULTS)
    density_df = pd.read_csv(DENSITY_RESULTS) if os.path.exists(DENSITY_RESULTS) else None

    # ══════════════════════════════════════════════════════════════════
    # SUBCLASS-LEVEL MARKERS (each subclass vs all others)
    # ══════════════════════════════════════════════════════════════════
    print("\n" + "=" * 80)
    print("SUBCLASS-LEVEL MARKERS (subclass vs all)")
    print("=" * 80)

    subclass_markers = extract_markers(adata, groupby="Subclass",
                                        n_genes=N_MARKERS)
    subclass_markers["level"] = "subclass"

    # Coverage per subclass per panel
    print(f"\n{'Subclass':20s}", end="")
    for pn in panel_names:
        print(f" {pn:>12s}", end="")
    print()
    print("-" * (20 + 13 * len(panel_names)))

    subclass_coverage = []
    for sc_name in sorted(subclass_markers["group"].unique()):
        sc_m = subclass_markers[subclass_markers["group"] == sc_name]
        print(f"{sc_name:20s}", end="")
        for pn in panel_names:
            for topn in [10, 20, 50]:
                top = sc_m.head(topn)
                n_in = top["gene"].isin(panels[pn]).sum()
                genes_in = list(top[top["gene"].isin(panels[pn])]["gene"])
                subclass_coverage.append({
                    "celltype": sc_name,
                    "level": "subclass",
                    "subclass": sc_name,
                    "panel": pn,
                    "topN": topn,
                    "n_in_panel": n_in,
                    "genes_in_panel": "; ".join(genes_in),
                })
            # Print top10 count
            n10 = sc_m.head(10)["gene"].isin(panels[pn]).sum()
            print(f" {n10:>12d}", end="")
        print()

    # ══════════════════════════════════════════════════════════════════
    # WITHIN-SUBCLASS SUPERTYPE MARKERS
    # ══════════════════════════════════════════════════════════════════
    print("\n" + "=" * 80)
    print("WITHIN-SUBCLASS SUPERTYPE MARKERS")
    print("=" * 80)

    all_supertype_markers = []
    supertype_coverage = []

    for sc_name in sorted(adata.obs["Subclass"].unique()):
        mask = adata.obs["Subclass"] == sc_name
        adata_sub = adata[mask].copy()
        adata_sub.obs["Supertype"] = (
            adata_sub.obs["Supertype"].cat.remove_unused_categories()
        )
        n_types = adata_sub.obs["Supertype"].nunique()
        if n_types < 2:
            continue

        print(f"\n  {sc_name} ({adata_sub.n_obs} cells, {n_types} supertypes)")
        markers = extract_markers(adata_sub, groupby="Supertype",
                                   n_genes=N_MARKERS)
        markers["subclass"] = sc_name
        markers["level"] = "supertype"
        all_supertype_markers.append(markers)

        for st in sorted(markers["group"].unique()):
            st_m = markers[markers["group"] == st]
            for pn in panel_names:
                for topn in [5, 10, 20, 50]:
                    top = st_m.head(topn)
                    n_in = top["gene"].isin(panels[pn]).sum()
                    genes_in = list(top[top["gene"].isin(panels[pn])]["gene"])
                    supertype_coverage.append({
                        "celltype": st,
                        "level": "supertype",
                        "subclass": sc_name,
                        "panel": pn,
                        "topN": topn,
                        "n_in_panel": n_in,
                        "genes_in_panel": "; ".join(genes_in),
                    })

    supertype_markers_df = pd.concat(all_supertype_markers, ignore_index=True)

    # ── Combine coverage ──
    coverage_df = pd.DataFrame(subclass_coverage + supertype_coverage)
    cov_path = os.path.join(OUTPUT_DIR,
                             "cross_platform_marker_coverage.csv")
    coverage_df.to_csv(cov_path, index=False)
    print(f"\nSaved: {cov_path}")

    # ══════════════════════════════════════════════════════════════════
    # SUMMARY TABLE: Supertype marker coverage across panels
    # ══════════════════════════════════════════════════════════════════
    print("\n" + "=" * 100)
    print("SUPERTYPE MARKER COVERAGE COMPARISON (top-10 within-subclass markers)")
    print("=" * 100)

    # Pivot: supertype x panel → n_in_panel (for topN=10)
    sup_cov = coverage_df[
        (coverage_df["level"] == "supertype") &
        (coverage_df["topN"] == 10)
    ].copy()
    pivot = sup_cov.pivot_table(index=["subclass", "celltype"],
                                 columns="panel",
                                 values="n_in_panel",
                                 aggfunc="first")
    pivot = pivot[panel_names]  # order columns

    # Add depth info
    depth_map = dict(zip(depth_df["supertype"], depth_df["median_depth_merfish"]))
    pivot["depth"] = pivot.index.get_level_values("celltype").map(depth_map)

    # Add layer specificity
    layer_distinct_map = {}
    for subclass in depth_df["subclass"].unique():
        sub_d = depth_df[depth_df["subclass"] == subclass]
        if len(sub_d) < 2:
            for st in sub_d["supertype"]:
                layer_distinct_map[st] = False
            continue
        depths = sub_d["median_depth_merfish"].values
        for _, row in sub_d.iterrows():
            own = row["median_depth_merfish"]
            others = depths[depths != own]
            min_sep = np.min(np.abs(others - own)) if len(others) > 0 else 0
            layer_distinct_map[row["supertype"]] = min_sep > 0.10
    pivot["layer_distinct"] = pivot.index.get_level_values("celltype").map(
        layer_distinct_map)

    # Print
    print(f"\n{'Supertype':20s} {'Sub':12s}", end="")
    for pn in panel_names:
        short = pn.replace("MERSCOPE_", "MS").replace("MERFISH_", "MF").replace("Xenium_", "X")
        print(f" {short:>6s}", end="")
    print(f" {'depth':>6s} {'layer':>6s}")
    print("-" * (20 + 12 + 7 * len(panel_names) + 14))

    for (subclass, st), row in pivot.iterrows():
        print(f"{st:20s} {subclass:12s}", end="")
        for pn in panel_names:
            val = int(row[pn]) if pd.notna(row[pn]) else 0
            print(f" {val:>6d}", end="")
        depth_str = f"{row['depth']:.2f}" if pd.notna(row['depth']) else "n/a"
        layer_str = "✓" if row.get("layer_distinct") else " "
        print(f" {depth_str:>6s} {layer_str:>6s}")

    # ══════════════════════════════════════════════════════════════════
    # AGGREGATE STATISTICS
    # ══════════════════════════════════════════════════════════════════
    print("\n" + "=" * 80)
    print("AGGREGATE STATISTICS")
    print("=" * 80)

    for topn in [10, 20, 50]:
        print(f"\n  Top-{topn} within-subclass markers:")
        sup_topn = coverage_df[
            (coverage_df["level"] == "supertype") &
            (coverage_df["topN"] == topn)
        ]
        for pn in panel_names:
            pn_data = sup_topn[sup_topn["panel"] == pn]
            mean_cov = pn_data["n_in_panel"].mean()
            zero_count = (pn_data["n_in_panel"] == 0).sum()
            n_total = len(pn_data)
            pct_zero = 100 * zero_count / n_total if n_total > 0 else 0
            ge3 = (pn_data["n_in_panel"] >= 3).sum()
            pct_ge3 = 100 * ge3 / n_total if n_total > 0 else 0
            print(f"    {pn:15s}: mean={mean_cov:.1f}/{topn}, "
                  f"zero={zero_count}/{n_total} ({pct_zero:.0f}%), "
                  f"≥3={ge3}/{n_total} ({pct_ge3:.0f}%)")

    # ── Subclass-level coverage ──
    print(f"\n  Subclass-level markers (subclass vs all):")
    for topn in [10, 20]:
        print(f"\n    Top-{topn}:")
        sc_topn = coverage_df[
            (coverage_df["level"] == "subclass") &
            (coverage_df["topN"] == topn)
        ]
        for pn in panel_names:
            pn_data = sc_topn[sc_topn["panel"] == pn]
            mean_cov = pn_data["n_in_panel"].mean()
            zero_count = (pn_data["n_in_panel"] == 0).sum()
            n_total = len(pn_data)
            ge3 = (pn_data["n_in_panel"] >= 3).sum()
            print(f"      {pn:15s}: mean={mean_cov:.1f}/{topn}, "
                  f"zero={zero_count}/{n_total}, ≥3={ge3}/{n_total}")

    # ══════════════════════════════════════════════════════════════════
    # CONFIDENCE RATING WITH ALL PANELS
    # ══════════════════════════════════════════════════════════════════
    print("\n" + "=" * 80)
    print("CLASSIFICATION CONFIDENCE BY PANEL")
    print("=" * 80)

    # Build confidence table for each panel
    for pn in panel_names:
        sup_10 = coverage_df[
            (coverage_df["level"] == "supertype") &
            (coverage_df["topN"] == 10) &
            (coverage_df["panel"] == pn)
        ].copy()

        def rate(row):
            m = row["n_in_panel"]
            ld = layer_distinct_map.get(row["celltype"], False)
            if m >= 3:
                return "HIGH"
            elif m >= 2:
                return "HIGH" if ld else "MEDIUM"
            elif m == 1:
                return "MEDIUM" if ld else "LOW"
            else:
                return "MEDIUM" if ld else "LOW"

        sup_10["confidence"] = sup_10.apply(rate, axis=1)
        counts = sup_10["confidence"].value_counts()
        n_total = len(sup_10)
        print(f"\n  {pn} ({len(panels[pn])} genes):")
        for level in ["HIGH", "MEDIUM", "LOW"]:
            n = counts.get(level, 0)
            print(f"    {level:8s}: {n:3d}/{n_total} supertypes "
                  f"({100*n/n_total:.0f}%)")

    # ══════════════════════════════════════════════════════════════════
    # FIGURE: Heatmap of marker coverage across panels
    # ══════════════════════════════════════════════════════════════════
    print("\nGenerating coverage heatmap...")

    # Build matrix: supertypes x panels (top-10 coverage)
    sup_cov10 = coverage_df[
        (coverage_df["level"] == "supertype") &
        (coverage_df["topN"] == 10)
    ]
    heatmap_pivot = sup_cov10.pivot_table(
        index="celltype", columns="panel", values="n_in_panel", aggfunc="first"
    )[panel_names]

    # Also add top-20 for Xenium v1
    sup_cov20_v1 = coverage_df[
        (coverage_df["level"] == "supertype") &
        (coverage_df["topN"] == 20) &
        (coverage_df["panel"] == "Xenium_v1")
    ].set_index("celltype")["n_in_panel"]
    heatmap_pivot["Xenium_v1_top20"] = sup_cov20_v1

    # Sort by subclass then depth
    subclass_map = dict(zip(
        coverage_df[coverage_df["level"] == "supertype"]["celltype"],
        coverage_df[coverage_df["level"] == "supertype"]["subclass"]))
    heatmap_pivot["subclass"] = heatmap_pivot.index.map(subclass_map)
    heatmap_pivot["depth"] = heatmap_pivot.index.map(depth_map)
    heatmap_pivot = heatmap_pivot.sort_values(["subclass", "depth"])

    # Create figure
    fig, ax = plt.subplots(figsize=(12, max(20, len(heatmap_pivot) * 0.22)))
    plot_cols = panel_names + ["Xenium_v1_top20"]
    data = heatmap_pivot[plot_cols].values.astype(float)

    im = ax.imshow(data, aspect="auto", cmap="YlOrRd", vmin=0, vmax=10,
                    interpolation="nearest")

    ax.set_yticks(range(len(heatmap_pivot)))
    ax.set_yticklabels(heatmap_pivot.index, fontsize=7)
    col_labels = ["Xen v1\n(266)", "MERFISH\n(180)", "MERSCOPE\n(250)",
                   "MERSCOPE\n(4K)", "Xen 5K\n(5001)", "Xen v1\ntop20"]
    ax.set_xticks(range(len(plot_cols)))
    ax.set_xticklabels(col_labels, fontsize=10, ha="center")

    # Add text annotations
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = int(data[i, j]) if not np.isnan(data[i, j]) else 0
            color = "white" if val >= 5 else "black"
            ax.text(j, i, str(val), ha="center", va="center",
                     fontsize=6, color=color)

    # Add subclass separators
    prev_sc = None
    for i, (st, row) in enumerate(heatmap_pivot.iterrows()):
        sc = row["subclass"]
        if prev_sc is not None and sc != prev_sc:
            ax.axhline(i - 0.5, color="white", linewidth=1.5)
        prev_sc = sc

    plt.colorbar(im, ax=ax, label="# top-10 markers in panel", shrink=0.5)
    ax.set_title("Within-subclass marker coverage across spatial panels",
                  fontsize=14, fontweight="bold", pad=10)
    plt.tight_layout()

    fig_path = os.path.join(OUTPUT_DIR,
                             "cross_platform_marker_coverage_heatmap.png")
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {fig_path}")

    # ══════════════════════════════════════════════════════════════════
    # FIGURE 2: Bar chart - mean coverage by subclass across panels
    # ══════════════════════════════════════════════════════════════════
    print("Generating subclass mean coverage bar chart...")

    subclass_means = sup_cov10.groupby(["subclass", "panel"])["n_in_panel"].mean()
    subclass_means = subclass_means.unstack("panel")[panel_names]

    fig, ax = plt.subplots(figsize=(14, 7))
    x = np.arange(len(subclass_means))
    width = 0.15
    colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]

    for i, (pn, color) in enumerate(zip(panel_names, colors)):
        short = pn.replace("MERSCOPE_", "MS ").replace("MERFISH_", "MF ").replace("Xenium_", "X ")
        ax.bar(x + i * width, subclass_means[pn], width, label=short,
               color=color, alpha=0.85)

    ax.set_xticks(x + width * 2)
    ax.set_xticklabels(subclass_means.index, rotation=45, ha="right",
                        fontsize=10)
    ax.set_ylabel("Mean # of top-10 markers in panel", fontsize=14)
    ax.set_title("Within-subclass marker coverage by subclass and panel",
                  fontsize=16, fontweight="bold")
    ax.legend(fontsize=10)
    ax.axhline(3, color="gray", linestyle="--", alpha=0.5, label="≥3 threshold")
    ax.set_ylim(0, 10)
    plt.tight_layout()

    fig2_path = os.path.join(OUTPUT_DIR,
                              "cross_platform_marker_coverage_bars.png")
    plt.savefig(fig2_path, dpi=150)
    plt.close()
    print(f"Saved: {fig2_path}")

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.0f}s")


if __name__ == "__main__":
    main()
