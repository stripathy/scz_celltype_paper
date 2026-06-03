#!/usr/bin/env python3
"""
Assess supertype classification confidence for all subclasses.

For each subclass, computes within-subclass Wilcoxon DE markers using the
snRNAseq reference, then checks how many of those markers are in the Xenium
v1 brain panel (~300 genes). Also computes layer specificity of supertypes
as a second axis of classification confidence.

Combines both signals into a per-supertype confidence rating.

Usage:
    python3 -u supertype_classification_confidence.py
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import matplotlib.pyplot as plt
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    BASE_DIR, SNRNASEQ_REF_PATH, PANEL_V1_PATH, CRUMBLR_DIR,
)

sys.path.insert(0, os.path.join(BASE_DIR, "code", "modules"))
from reference_utils import load_and_normalize_reference, subsample_by_group

OUTPUT_DIR = os.path.join(BASE_DIR, "output", "marker_analysis")
DEPTH_FILE = os.path.join(BASE_DIR, "output", "presentation",
                           "median_depth_supertype.csv")
DENSITY_RESULTS = os.path.join(BASE_DIR, "output", "density_analysis",
                                "density_results_supertype_cortical.csv")
CRUMBLR_RESULTS = os.path.join(CRUMBLR_DIR, "crumblr_results_supertype.csv")

# Xenium v1 panel genes
V1_PANEL = None  # loaded in main()

N_MARKERS = 50  # markers per supertype
MAX_CELLS_PER_TYPE = 500
MIN_CELLS_PER_TYPE = 20


def extract_markers(adata_sub, groupby, n_genes):
    """Run Wilcoxon rank-sum and extract markers as DataFrame."""
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


def compute_layer_specificity(depth_df):
    """Compute layer specificity metrics for each supertype.

    Returns depth spread (IQR proxy using MERFISH depth) and whether
    supertypes within a subclass occupy distinct depth niches.
    """
    # We only have median depth, but that's enough to compute
    # within-subclass depth separation
    return depth_df


def infer_subclass(supertype):
    """Infer subclass from supertype name."""
    # Handle special cases
    if supertype.startswith("L2/3 IT"):
        return "L2/3 IT"
    if supertype.startswith("L4 IT"):
        return "L4 IT"
    if supertype.startswith("L5 ET"):
        return "L5 ET"
    if supertype.startswith("L5 IT"):
        return "L5 IT"
    if supertype.startswith("L5/6 NP"):
        return "L5/6 NP"
    if supertype.startswith("L6 CT"):
        return "L6 CT"
    if supertype.startswith("L6 IT Car3"):
        return "L6 IT Car3"
    if supertype.startswith("L6 IT"):
        return "L6 IT"
    if supertype.startswith("L6b"):
        return "L6b"
    if supertype.startswith("Sst Chodl"):
        return "Sst Chodl"
    if supertype.startswith("Micro-PVM"):
        return "Micro-PVM"
    # Default: split on underscore, take prefix
    parts = supertype.rsplit("_", 1)
    if len(parts) == 2 and parts[1].replace("-SEAAD", "").isdigit():
        return parts[0]
    return supertype


def main():
    global V1_PANEL
    t0 = time.time()
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ── Load Xenium v1 panel ──
    print("Loading Xenium v1 panel...")
    panel_df = pd.read_csv(PANEL_V1_PATH)
    # Try common column names
    for col in ["Gene", "gene", "Name", "name", "gene_name"]:
        if col in panel_df.columns:
            V1_PANEL = set(panel_df[col].dropna().astype(str))
            break
    if V1_PANEL is None:
        # First column likely has gene names
        V1_PANEL = set(panel_df.iloc[:, 0].dropna().astype(str))
    print(f"  {len(V1_PANEL)} genes in v1 panel")

    # ── Load reference ──
    print("\nLoading snRNAseq reference...")
    adata = load_and_normalize_reference(SNRNASEQ_REF_PATH,
                                          normalize=True, min_cells=10)
    adata = subsample_by_group(adata, "Supertype",
                                max_cells=MAX_CELLS_PER_TYPE,
                                min_cells=MIN_CELLS_PER_TYPE)

    # Clean categories
    for col in ["Supertype", "Subclass"]:
        if hasattr(adata.obs[col], "cat"):
            adata.obs[col] = adata.obs[col].cat.remove_unused_categories()
        else:
            adata.obs[col] = pd.Categorical(adata.obs[col])

    print(f"  {adata.n_obs} cells, {adata.obs['Supertype'].nunique()} supertypes, "
          f"{adata.obs['Subclass'].nunique()} subclasses")

    # ── Load depth data ──
    print("\nLoading depth data...")
    depth_df = pd.read_csv(DEPTH_FILE)
    depth_df["subclass"] = depth_df["supertype"].apply(infer_subclass)
    print(f"  {len(depth_df)} supertypes with depth info")

    # ── Load results files ──
    crumblr_df = pd.read_csv(CRUMBLR_RESULTS)
    density_df = pd.read_csv(DENSITY_RESULTS) if os.path.exists(DENSITY_RESULTS) else None

    # ══════════════════════════════════════════════════════════════════
    # WITHIN-SUBCLASS MARKER DE FOR ALL SUBCLASSES
    # ══════════════════════════════════════════════════════════════════
    print("\n" + "=" * 80)
    print("WITHIN-SUBCLASS WILCOXON DE FOR ALL SUBCLASSES")
    print("=" * 80)

    all_markers = []
    subclasses = sorted(adata.obs["Subclass"].unique())

    for sc_name in subclasses:
        mask = adata.obs["Subclass"] == sc_name
        adata_sub = adata[mask].copy()
        adata_sub.obs["Supertype"] = (
            adata_sub.obs["Supertype"].cat.remove_unused_categories()
        )
        n_types = adata_sub.obs["Supertype"].nunique()

        if n_types < 2:
            print(f"\n  {sc_name}: only {n_types} supertype(s), skipping")
            continue

        n_cells = adata_sub.n_obs
        print(f"\n  {sc_name}: {n_cells} cells, {n_types} supertypes")

        markers = extract_markers(adata_sub, groupby="Supertype",
                                   n_genes=N_MARKERS)
        markers["subclass"] = sc_name

        # Check panel coverage for top-10
        for st in markers["group"].unique():
            top10 = markers[markers["group"] == st].head(10)
            in_panel = top10["gene"].isin(V1_PANEL).sum()
            genes_in = list(top10[top10["gene"].isin(V1_PANEL)]["gene"])
            genes_str = ", ".join(genes_in) if genes_in else "none"
            print(f"    {st}: {in_panel}/10 top markers in v1 panel ({genes_str})")

        all_markers.append(markers)

    markers_df = pd.concat(all_markers, ignore_index=True)

    # Save all within-subclass markers
    markers_path = os.path.join(OUTPUT_DIR,
                                 "within_subclass_markers_all.csv")
    markers_df.to_csv(markers_path, index=False)
    print(f"\nSaved all markers: {markers_path}")

    # ══════════════════════════════════════════════════════════════════
    # PANEL COVERAGE SUMMARY TABLE
    # ══════════════════════════════════════════════════════════════════
    print("\n" + "=" * 80)
    print("PANEL COVERAGE SUMMARY")
    print("=" * 80)

    coverage_records = []
    for st in markers_df["group"].unique():
        st_markers = markers_df[markers_df["group"] == st]
        subclass = st_markers["subclass"].iloc[0]

        for topn in [5, 10, 20]:
            top = st_markers.head(topn)
            in_panel = top["gene"].isin(V1_PANEL).sum()
            genes_in = list(top[top["gene"].isin(V1_PANEL)]["gene"])
            coverage_records.append({
                "supertype": st,
                "subclass": subclass,
                f"top{topn}_in_panel": in_panel,
                f"top{topn}_genes": "; ".join(genes_in),
            })

    # Merge the topN columns
    cov_5 = pd.DataFrame([r for r in coverage_records
                           if "top5_in_panel" in r])[
        ["supertype", "subclass", "top5_in_panel", "top5_genes"]]
    cov_10 = pd.DataFrame([r for r in coverage_records
                            if "top10_in_panel" in r])[
        ["supertype", "top10_in_panel", "top10_genes"]]
    cov_20 = pd.DataFrame([r for r in coverage_records
                            if "top20_in_panel" in r])[
        ["supertype", "top20_in_panel", "top20_genes"]]

    coverage_df = cov_5.merge(cov_10, on="supertype").merge(cov_20, on="supertype")

    # ══════════════════════════════════════════════════════════════════
    # LAYER SPECIFICITY ANALYSIS
    # ══════════════════════════════════════════════════════════════════
    print("\n" + "=" * 80)
    print("LAYER SPECIFICITY ANALYSIS")
    print("=" * 80)

    # For each subclass with multiple supertypes, compute:
    # 1. Depth range of the subclass (max - min median depth among supertypes)
    # 2. Mean pairwise depth separation between supertypes
    # 3. Whether a supertype is at an extreme depth within its subclass

    layer_records = []
    for subclass in depth_df["subclass"].unique():
        sub_depth = depth_df[depth_df["subclass"] == subclass].copy()
        if len(sub_depth) < 2:
            for _, row in sub_depth.iterrows():
                layer_records.append({
                    "supertype": row["supertype"],
                    "subclass": subclass,
                    "median_depth_merfish": row["median_depth_merfish"],
                    "median_depth_xenium": row["median_depth_xenium"],
                    "n_cells_xenium": row["n_cells_xenium"],
                    "depth_range_subclass": 0,
                    "min_separation": 0,
                    "depth_rank_in_subclass": "only",
                    "layer_distinct": False,
                })
            continue

        depths = sub_depth["median_depth_merfish"].values
        depth_range = depths.max() - depths.min()

        for _, row in sub_depth.iterrows():
            own_depth = row["median_depth_merfish"]
            other_depths = depths[depths != own_depth]
            if len(other_depths) == 0:
                min_sep = 0
            else:
                min_sep = np.min(np.abs(other_depths - own_depth))

            # Rank within subclass (shallow→deep)
            rank = int(np.searchsorted(np.sort(depths), own_depth)) + 1

            # "Layer distinct" if min separation > 0.1 normalized depth
            # (roughly one cortical layer worth of separation)
            layer_distinct = min_sep > 0.10

            layer_records.append({
                "supertype": row["supertype"],
                "subclass": subclass,
                "median_depth_merfish": row["median_depth_merfish"],
                "median_depth_xenium": row["median_depth_xenium"],
                "n_cells_xenium": int(row["n_cells_xenium"]),
                "depth_range_subclass": depth_range,
                "min_separation": min_sep,
                "depth_rank_in_subclass": f"{rank}/{len(sub_depth)}",
                "layer_distinct": layer_distinct,
            })

    layer_df = pd.DataFrame(layer_records)

    print("\n  Supertypes with layer-distinct depth positions:")
    distinct = layer_df[layer_df["layer_distinct"]]
    for _, r in distinct.iterrows():
        print(f"    {r['supertype']:20s} (depth={r['median_depth_merfish']:.2f}, "
              f"min_sep={r['min_separation']:.2f}, "
              f"subclass range={r['depth_range_subclass']:.2f})")

    # ══════════════════════════════════════════════════════════════════
    # COMBINED CONFIDENCE TABLE
    # ══════════════════════════════════════════════════════════════════
    print("\n" + "=" * 80)
    print("COMBINED SUPERTYPE CLASSIFICATION CONFIDENCE")
    print("=" * 80)

    # Merge coverage + layer specificity
    conf_df = coverage_df.merge(
        layer_df[["supertype", "median_depth_merfish", "median_depth_xenium",
                  "n_cells_xenium", "depth_range_subclass", "min_separation",
                  "layer_distinct"]],
        on="supertype", how="left"
    )

    # Add crumblr significance
    crumblr_sig = crumblr_df[["celltype", "logFC", "adj.P.Val"]].rename(
        columns={"celltype": "supertype", "adj.P.Val": "crumblr_fdr",
                 "logFC": "crumblr_logFC"})
    conf_df = conf_df.merge(crumblr_sig, on="supertype", how="left")

    # Add density significance if available
    if density_df is not None:
        dens_sig = density_df[["supertype", "logFC", "fdr"]].rename(
            columns={"fdr": "density_fdr", "logFC": "density_logFC"})
        conf_df = conf_df.merge(dens_sig, on="supertype", how="left")

    # Compute confidence rating
    def rate_confidence(row):
        markers = row.get("top10_in_panel", 0)
        layer = row.get("layer_distinct", False)

        if markers >= 3:
            return "HIGH"
        elif markers >= 2:
            return "HIGH" if layer else "MEDIUM"
        elif markers == 1:
            return "MEDIUM" if layer else "LOW"
        else:
            return "MEDIUM" if layer else "LOW"

    conf_df["confidence"] = conf_df.apply(rate_confidence, axis=1)

    # Sort by significance
    conf_df = conf_df.sort_values("crumblr_fdr", na_position="last")

    # Save full table
    conf_path = os.path.join(OUTPUT_DIR,
                              "supertype_classification_confidence.csv")
    conf_df.to_csv(conf_path, index=False)
    print(f"\nSaved: {conf_path}")

    # ── Print results for significant supertypes ──
    # Focus on those with crumblr FDR < 0.1 or density FDR < 0.1
    sig_mask = conf_df["crumblr_fdr"] < 0.1
    if density_df is not None:
        sig_mask = sig_mask | (conf_df["density_fdr"] < 0.1)

    sig_df = conf_df[sig_mask].copy()

    print(f"\n{'='*100}")
    print(f"SIGNIFICANT SUPERTYPES (crumblr or density FDR < 0.1)")
    print(f"{'='*100}")
    print(f"\n{'Supertype':20s} {'Sub':12s} {'Top10':>5s} {'Layer':>6s} "
          f"{'MinSep':>7s} {'Conf':>6s} {'CR_FDR':>8s} {'CR_lFC':>7s} "
          f"{'Top10 genes in panel'}")
    print("-" * 100)

    for _, r in sig_df.iterrows():
        layer_str = "✓" if r.get("layer_distinct", False) else "✗"
        min_sep = r.get("min_separation", 0)
        cr_fdr = f"{r['crumblr_fdr']:.3f}" if pd.notna(r.get("crumblr_fdr")) else "n/a"
        cr_lfc = f"{r['crumblr_logFC']:+.3f}" if pd.notna(r.get("crumblr_logFC")) else "n/a"
        genes = r.get("top10_genes", "")
        print(f"{r['supertype']:20s} {r['subclass']:12s} {r['top10_in_panel']:5d} "
              f"{layer_str:>6s} {min_sep:7.3f} {r['confidence']:>6s} "
              f"{cr_fdr:>8s} {cr_lfc:>7s} {genes}")

    # ── Summary by confidence level ──
    print(f"\n{'='*80}")
    print("SUMMARY: Significant supertypes by confidence level")
    print(f"{'='*80}")

    for conf_level in ["HIGH", "MEDIUM", "LOW"]:
        subset = sig_df[sig_df["confidence"] == conf_level]
        if len(subset) == 0:
            continue
        print(f"\n  {conf_level} confidence ({len(subset)} supertypes):")
        for _, r in subset.iterrows():
            markers_str = r.get("top10_genes", "none") or "none"
            layer_str = f"depth={r['median_depth_merfish']:.2f}, sep={r['min_separation']:.2f}" if pd.notna(r.get("min_separation")) else ""
            cr_fdr = f"FDR={r['crumblr_fdr']:.3f}" if pd.notna(r.get("crumblr_fdr")) else ""
            cr_lfc = f"logFC={r['crumblr_logFC']:+.2f}" if pd.notna(r.get("crumblr_logFC")) else ""
            print(f"    {r['supertype']:20s} markers: {r['top10_in_panel']}/10 [{markers_str}]  "
                  f"layer_distinct={'Y' if r.get('layer_distinct') else 'N'}  "
                  f"{cr_fdr} {cr_lfc}")

    # ── All supertypes summary ──
    print(f"\n{'='*80}")
    print("ALL SUPERTYPES: Panel marker coverage summary by subclass")
    print(f"{'='*80}")

    for subclass in sorted(conf_df["subclass"].unique()):
        sub = conf_df[conf_df["subclass"] == subclass]
        if len(sub) < 2:
            continue
        mean_cov = sub["top10_in_panel"].mean()
        n_zero = (sub["top10_in_panel"] == 0).sum()
        n_types = len(sub)
        depth_range = sub["depth_range_subclass"].iloc[0] if pd.notna(
            sub["depth_range_subclass"].iloc[0]) else 0

        print(f"\n  {subclass} ({n_types} supertypes, "
              f"depth range={depth_range:.2f}):")
        for _, r in sub.sort_values("median_depth_merfish").iterrows():
            layer_str = "✓" if r.get("layer_distinct") else " "
            depth = f"{r['median_depth_merfish']:.2f}" if pd.notna(
                r.get("median_depth_merfish")) else "n/a"
            genes = r.get("top10_genes", "") or ""
            print(f"    {r['supertype']:20s} {r['top10_in_panel']:2d}/10  "
                  f"depth={depth} {layer_str}  [{genes}]")

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.0f}s")


if __name__ == "__main__":
    main()
