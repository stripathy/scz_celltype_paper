#!/usr/bin/env python3
"""
Visualize quality of selected add-on probe markers using snRNAseq reference data.

Generates grouped dot plots (fraction expressing × mean expression) for:
  1. Round 2 SST supertype markers across all SST supertypes (depth-ordered)
  2. Round 2 L6b supertype markers across all L6b supertypes (depth-ordered)
  3. Round 1 subclass-level SST markers across major subclasses
  4. SST matrix plot (mean expression heatmap, depth-ordered)
  5. SST stacked violin (top 2 markers per supertype)

Pipeline position: 8 of 8 (probe design visualization)
Upstream: hierarchical_probe_selection.py (reads v1_addon_markers_clean.csv)
          plot_median_depth_by_celltype.py (reads median_depth_supertype.csv)
Downstream: none (terminal visualization)

Usage:
    python visualize_marker_quality.py
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=FutureWarning)

# Global font sizes for legibility
plt.rcParams.update({
    "font.size": 14,
    "axes.titlesize": 20,
    "axes.labelsize": 16,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
})

# ── Shared module imports ────────────────────────────────────────────────────
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "modules"))
from reference_utils import load_and_normalize_reference, subsample_by_group

# ── Configuration ────────────────────────────────────────────────────────────
BASE_DIR = os.path.expanduser("~/Github/SCZ_Xenium")
REF_PATH = os.path.join(BASE_DIR, "data", "reference",
                         "nicole_sea_ad_snrnaseq_reference.h5ad")
MARKERS_PATH = os.path.join(BASE_DIR, "output", "marker_analysis",
                             "v1_addon_markers_clean.csv")
DEPTH_PATH = os.path.join(BASE_DIR, "output", "presentation",
                           "median_depth_supertype.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "marker_analysis")
os.makedirs(OUTPUT_DIR, exist_ok=True)

MAX_CELLS_PER_TYPE = 500
MIN_CELLS_PER_TYPE = 20


def natural_sort_key(s):
    """Sort SST/L6b supertypes numerically: Sst_1, Sst_2, ..., Sst_25."""
    import re
    parts = re.split(r'(\d+)', s)
    return [int(p) if p.isdigit() else p for p in parts]


def main():
    # ── Load marker gene list ────────────────────────────────────────────────
    print("Loading marker gene list...")
    markers_df = pd.read_csv(MARKERS_PATH)
    print(f"  {len(markers_df)} total markers loaded")

    # Split by category
    cardinal = markers_df[markers_df["selection_category"] == "1_Cardinal"]
    round1 = markers_df[markers_df["selection_category"] == "2_Subclass_Round1"]
    round2 = markers_df[markers_df["selection_category"] == "3_Supertype_Round2"]

    # Split Round 2 into SST vs L6b markers
    round2_sst = round2[round2["target_group"].str.startswith("Sst_")]
    round2_l6b = round2[round2["target_group"].str.startswith("L6b_")]

    print(f"  Cardinal: {len(cardinal)}, Round1 (subclass): {len(round1)}, "
          f"Round2 SST: {len(round2_sst)}, Round2 L6b: {len(round2_l6b)}")

    # ── Load depth data ─────────────────────────────────────────────────────
    print("Loading supertype depth data...")
    depth_df = pd.read_csv(DEPTH_PATH)
    depth_lookup = dict(zip(depth_df["supertype"], depth_df["median_depth_merfish"]))
    print(f"  Loaded median MERFISH depth for {len(depth_lookup)} supertypes")

    # ── Load reference ───────────────────────────────────────────────────────
    adata = load_and_normalize_reference(REF_PATH, normalize=True, min_cells=10)
    adata = subsample_by_group(adata, "Supertype", max_cells=MAX_CELLS_PER_TYPE,
                                min_cells=MIN_CELLS_PER_TYPE)

    # ══════════════════════════════════════════════════════════════════════════
    # PLOT 1: SST supertype dotplot — Round 2 markers grouped by target type
    # ══════════════════════════════════════════════════════════════════════════
    print("\n── Plot 1: SST supertype dot plot ──")

    # Subset to SST cells
    adata_sst = adata[adata.obs["Subclass"] == "Sst"].copy()
    adata_sst.obs["Supertype"] = adata_sst.obs["Supertype"].cat.remove_unused_categories()
    sst_supertypes = sorted(adata_sst.obs["Supertype"].unique().tolist(),
                             key=lambda x: depth_lookup.get(x, 0.5))
    print(f"  SST cells: {adata_sst.shape[0]}, supertypes: {len(sst_supertypes)}")
    print("  Depth ordering (superficial → deep):")
    for st in sst_supertypes:
        d = depth_lookup.get(st, float("nan"))
        print(f"    {st:>8s}: {d:.3f}")

    # Order genes by target supertype depth, then by composite score
    round2_sst_sorted = round2_sst.copy()
    round2_sst_sorted["_depth"] = round2_sst_sorted["target_group"].map(depth_lookup)
    round2_sst_sorted = round2_sst_sorted.sort_values(
        ["_depth", "dropout_adjusted_score"],
        ascending=[True, False])
    round2_sst_sorted = round2_sst_sorted.drop(columns=["_depth"])

    # Filter to genes present in adata
    sst_marker_genes = [g for g in round2_sst_sorted["gene"].tolist()
                        if g in adata_sst.var_names]
    missing = set(round2_sst_sorted["gene"]) - set(sst_marker_genes)
    if missing:
        print(f"  Warning: {len(missing)} genes not in reference: {missing}")

    # Build var_group_positions for gene grouping by target supertype
    sst_gene_to_target = dict(zip(round2_sst_sorted["gene"],
                                   round2_sst_sorted["target_group"]))
    gene_groups = []
    group_positions = []
    group_labels = []
    current_group = None
    group_start = 0
    for i, gene in enumerate(sst_marker_genes):
        target = sst_gene_to_target.get(gene, "Unknown")
        if target != current_group:
            if current_group is not None:
                group_positions.append((group_start, i - 1))
                group_labels.append(current_group)
            current_group = target
            group_start = i
    if current_group is not None:
        group_positions.append((group_start, len(sst_marker_genes) - 1))
        group_labels.append(current_group)

    print(f"  Plotting {len(sst_marker_genes)} genes across {len(group_labels)} "
          f"target groups")

    # Order supertypes by depth from pia (superficial at top)
    adata_sst.obs["Supertype"] = pd.Categorical(
        adata_sst.obs["Supertype"], categories=sst_supertypes, ordered=True)

    dp = sc.pl.dotplot(
        adata_sst,
        var_names=sst_marker_genes,
        groupby="Supertype",
        var_group_positions=group_positions,
        var_group_labels=group_labels,
        var_group_rotation=45,
        standard_scale="var",
        cmap="Blues",
        show=False,
        return_fig=True,
        figsize=(24, 8),
    )
    dp.style(dot_edge_color="black", dot_edge_lw=0.5)
    dp.savefig(os.path.join(OUTPUT_DIR, "dotplot_sst_supertype_markers.png"),
               dpi=150, bbox_inches="tight")
    plt.close("all")
    print("  Saved: dotplot_sst_supertype_markers.png")

    # ══════════════════════════════════════════════════════════════════════════
    # PLOT 2: L6b supertype dotplot — Round 2 markers grouped by target type
    # ══════════════════════════════════════════════════════════════════════════
    print("\n── Plot 2: L6b supertype dot plot ──")

    adata_l6b = adata[adata.obs["Subclass"] == "L6b"].copy()
    adata_l6b.obs["Supertype"] = adata_l6b.obs["Supertype"].cat.remove_unused_categories()
    l6b_supertypes = sorted(adata_l6b.obs["Supertype"].unique().tolist(),
                             key=lambda x: depth_lookup.get(x, 0.8))
    print(f"  L6b cells: {adata_l6b.shape[0]}, supertypes: {len(l6b_supertypes)}")

    round2_l6b_sorted = round2_l6b.copy()
    round2_l6b_sorted["_depth"] = round2_l6b_sorted["target_group"].map(depth_lookup)
    round2_l6b_sorted = round2_l6b_sorted.sort_values(
        ["_depth", "dropout_adjusted_score"],
        ascending=[True, False])
    round2_l6b_sorted = round2_l6b_sorted.drop(columns=["_depth"])

    l6b_marker_genes = [g for g in round2_l6b_sorted["gene"].tolist()
                        if g in adata_l6b.var_names]
    l6b_gene_to_target = dict(zip(round2_l6b_sorted["gene"],
                                   round2_l6b_sorted["target_group"]))

    # Build var_group_positions for L6b
    l6b_group_positions = []
    l6b_group_labels = []
    current_group = None
    group_start = 0
    for i, gene in enumerate(l6b_marker_genes):
        target = l6b_gene_to_target.get(gene, "Unknown")
        if target != current_group:
            if current_group is not None:
                l6b_group_positions.append((group_start, i - 1))
                l6b_group_labels.append(current_group)
            current_group = target
            group_start = i
    if current_group is not None:
        l6b_group_positions.append((group_start, len(l6b_marker_genes) - 1))
        l6b_group_labels.append(current_group)

    print(f"  Plotting {len(l6b_marker_genes)} genes across {len(l6b_group_labels)} "
          f"target groups")

    # Order supertypes by depth from pia
    adata_l6b.obs["Supertype"] = pd.Categorical(
        adata_l6b.obs["Supertype"], categories=l6b_supertypes, ordered=True)

    dp_l6b = sc.pl.dotplot(
        adata_l6b,
        var_names=l6b_marker_genes,
        groupby="Supertype",
        var_group_positions=l6b_group_positions,
        var_group_labels=l6b_group_labels,
        var_group_rotation=45,
        standard_scale="var",
        cmap="Blues",
        show=False,
        return_fig=True,
    )
    dp_l6b.style(dot_edge_color="black", dot_edge_lw=0.5)
    dp_l6b.savefig(os.path.join(OUTPUT_DIR, "dotplot_l6b_supertype_markers.png"),
                    dpi=150, bbox_inches="tight")
    plt.close("all")
    print("  Saved: dotplot_l6b_supertype_markers.png")

    # ══════════════════════════════════════════════════════════════════════════
    # PLOT 3: Round 1 subclass markers — SST-vs-all across major subclasses
    # ══════════════════════════════════════════════════════════════════════════
    print("\n── Plot 3: Subclass-level marker dot plot ──")

    # Use all cells, grouped by subclass
    adata_sub = subsample_by_group(adata, "Subclass", max_cells=MAX_CELLS_PER_TYPE,
                                    min_cells=MIN_CELLS_PER_TYPE)

    round1_genes = [g for g in round1["gene"].tolist()
                    if g in adata_sub.var_names]

    # Also include cardinal markers for context
    cardinal_genes = [g for g in cardinal["gene"].tolist()
                      if g in adata_sub.var_names]

    subclass_genes = cardinal_genes + round1_genes
    var_group_positions_sub = []
    var_group_labels_sub = []
    if cardinal_genes:
        var_group_positions_sub.append((0, len(cardinal_genes) - 1))
        var_group_labels_sub.append("Cardinal")
    if round1_genes:
        var_group_positions_sub.append(
            (len(cardinal_genes), len(subclass_genes) - 1))
        var_group_labels_sub.append("Sst Round 1")

    # Order subclasses: GABAergic first (SST highlighted), then Glutamatergic, then Non-neuronal
    gaba_order = ["Sst", "Pvalb", "Vip", "Lamp5", "Sncg", "Pax6"]
    glut_order = ["L2/3 IT", "L4 IT", "L5 IT", "L5 ET", "L5/6 NP",
                  "L6 IT", "L6 IT Car3", "L6 CT", "L6b"]
    non_neuronal = ["Astro", "Oligo", "OPC", "Micro-PVM", "Endo", "VLMC"]

    available_subclasses = adata_sub.obs["Subclass"].unique().tolist()
    ordered_subclasses = []
    for s in gaba_order + glut_order + non_neuronal:
        if s in available_subclasses:
            ordered_subclasses.append(s)
    for s in available_subclasses:
        if s not in ordered_subclasses:
            ordered_subclasses.append(s)

    adata_sub.obs["Subclass"] = pd.Categorical(
        adata_sub.obs["Subclass"], categories=ordered_subclasses, ordered=True)

    print(f"  Plotting {len(subclass_genes)} genes across "
          f"{len(ordered_subclasses)} subclasses")

    dp_sub = sc.pl.dotplot(
        adata_sub,
        var_names=subclass_genes,
        groupby="Subclass",
        var_group_positions=var_group_positions_sub,
        var_group_labels=var_group_labels_sub,
        var_group_rotation=45,
        standard_scale="var",
        cmap="Blues",
        show=False,
        return_fig=True,
        figsize=(16, 10),
    )
    dp_sub.style(dot_edge_color="black", dot_edge_lw=0.5)
    dp_sub.savefig(os.path.join(OUTPUT_DIR, "dotplot_subclass_markers.png"),
                    dpi=150, bbox_inches="tight")
    plt.close("all")
    print("  Saved: dotplot_subclass_markers.png")

    # ══════════════════════════════════════════════════════════════════════════
    # PLOT 4: Matrix plot — mean expression heatmap for SST markers
    # ══════════════════════════════════════════════════════════════════════════
    print("\n── Plot 4: SST matrixplot (mean expression) ──")

    mp = sc.pl.matrixplot(
        adata_sst,
        var_names=sst_marker_genes,
        groupby="Supertype",
        var_group_positions=group_positions,
        var_group_labels=group_labels,
        var_group_rotation=45,
        standard_scale="var",
        cmap="viridis",
        show=False,
        return_fig=True,
        figsize=(24, 8),
    )
    mp.savefig(os.path.join(OUTPUT_DIR, "matrixplot_sst_supertype_markers.png"),
               dpi=150, bbox_inches="tight")
    plt.close("all")
    print("  Saved: matrixplot_sst_supertype_markers.png")

    # ══════════════════════════════════════════════════════════════════════════
    # PLOT 5: Stacked violin — top markers per SST supertype
    # ══════════════════════════════════════════════════════════════════════════
    print("\n── Plot 5: SST stacked violin (top 2 per supertype) ──")

    # Pick top 2 markers per target supertype for a cleaner stacked violin
    top2_per_type = (round2_sst_sorted
                     .groupby("target_group")
                     .head(2)["gene"]
                     .tolist())
    top2_genes = [g for g in top2_per_type if g in adata_sst.var_names]

    if len(top2_genes) > 0:
        sc.pl.stacked_violin(
            adata_sst,
            var_names=top2_genes,
            groupby="Supertype",
            standard_scale="var",
            cmap="Blues",
            show=False,
            swap_axes=True,
            figsize=(14, max(4, len(top2_genes) * 0.4)),
        )
        plt.savefig(os.path.join(OUTPUT_DIR,
                                  "stacked_violin_sst_top2_markers.png"),
                     dpi=150, bbox_inches="tight")
        plt.close("all")
        print("  Saved: stacked_violin_sst_top2_markers.png")
    else:
        print("  Skipped: no top2 genes found in reference")

    # ══════════════════════════════════════════════════════════════════════════
    # SUMMARY
    # ══════════════════════════════════════════════════════════════════════════
    print("\n" + "=" * 80)
    print("VISUALIZATION COMPLETE")
    print("=" * 80)
    print(f"Output directory: {OUTPUT_DIR}")
    print("Files generated:")
    for f in ["dotplot_sst_supertype_markers.png",
              "dotplot_l6b_supertype_markers.png",
              "dotplot_subclass_markers.png",
              "matrixplot_sst_supertype_markers.png",
              "stacked_violin_sst_top2_markers.png"]:
        print(f"  - {f}")


if __name__ == "__main__":
    main()
