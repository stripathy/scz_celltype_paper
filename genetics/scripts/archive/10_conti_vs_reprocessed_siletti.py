#!/usr/bin/env python
"""
10_conti_vs_reprocessed_siletti.py — Apples-to-apples comparison of enrichment
from the original Siletti processing (Conti/Duncan et al.) vs our reprocessed
Siletti data (updated Allen Institute release).

Both analyses now use:
  - All 461 Siletti clusters (not filtered to neocortical)
  - ENTREZ gene IDs directly (no symbol mapping loss)
  - Independently computed specificity (not a subset of combined matrix)

The only difference is the source data version:
  - Conti: original 2022 Siletti loom file
  - Reprocessed: updated Allen Institute data

Also generates a Duncan-style Manhattan plot for both processings.

Outputs:
  - results/tables/conti_siletti_enrichment.csv
  - results/tables/reprocessed_siletti_461_enrichment.csv
  - results/tables/conti_vs_reprocessed_461_comparison.csv
  - results/figures/conti_vs_reprocessed_461_enrichment.png
  - results/figures/duncan_manhattan_replication.png
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats

from scz_celltype_enrichment.config import (
    MAGMA_GENES_OUT, GENE_LOC_FILE, TABLES_DIR, FIGURES_DIR, INTERMEDIATES_DIR,
    FIGURE_DPI, TITLE_FONTSIZE, LABEL_FONTSIZE, TICK_FONTSIZE, LEGEND_FONTSIZE,
    FDR_THRESHOLD, PROJECT_ROOT,
)
from scz_celltype_enrichment.enrichment.gwas import (
    load_magma_genes, build_covariates,
)
from scz_celltype_enrichment.enrichment.celltype import run_enrichment_all_types
from scz_celltype_enrichment.utils import ensure_dirs, Timer


CONTI_SPECIFICITY = PROJECT_ROOT / "data" / "conti_specificity_matrix.txt"
REPRO_SPECIFICITY = INTERMEDIATES_DIR / "reprocessed_siletti_specificity_entrez.csv"


def build_cluster_number_map():
    """
    Build bidirectional mappings between cluster names and Conti-style numbers.
    Conti's Cluster0 → sorted label index 0 (=Bcell_0), etc.
    """
    d = np.load(
        str(INTERMEDIATES_DIR / "siletti_cluster_level_stats.npz"),
        allow_pickle=True,
    )
    names = d["cluster_names"]
    labels = d["cluster_labels"]

    label_nums = np.array([int(l.split("_")[-1]) for l in labels])
    sorted_idx = np.argsort(label_nums)
    sorted_names = names[sorted_idx]

    name_to_num = {sorted_names[i]: i for i in range(len(sorted_names))}
    num_to_name = {i: sorted_names[i] for i in range(len(sorted_names))}

    return name_to_num, num_to_name


def load_conti_specificity(num_to_name):
    """Load Conti specificity and rename columns to cluster names."""
    with Timer("Loading Conti specificity"):
        spec = pd.read_csv(str(CONTI_SPECIFICITY), sep="\t", index_col=0)
    # Rename: Cluster0 -> name
    col_map = {f"Cluster{i}": num_to_name[i] for i in range(len(num_to_name))}
    spec.columns = [col_map.get(c, c) for c in spec.columns]
    print(f"  Shape: {spec.shape[0]:,} genes × {spec.shape[1]} clusters (ENTREZ index)")
    return spec


def load_reprocessed_specificity():
    """Load reprocessed Siletti specificity (ENTREZ-indexed, 461 clusters)."""
    with Timer("Loading reprocessed specificity"):
        spec = pd.read_csv(str(REPRO_SPECIFICITY), index_col=0)
    print(f"  Shape: {spec.shape[0]:,} genes × {spec.shape[1]} clusters (ENTREZ index)")
    return spec


def run_enrichment_entrez(specificity, label):
    """Run enrichment on an ENTREZ-indexed specificity matrix."""
    print(f"\n  Loading MAGMA gene-level results...")
    magma_df = load_magma_genes(MAGMA_GENES_OUT)

    # Intersect on ENTREZ IDs
    spec_genes = set(specificity.index)
    magma_genes = set(magma_df["GENE"].values)
    overlap = sorted(spec_genes & magma_genes)
    print(f"  Overlapping genes (ENTREZ): {len(overlap):,}")

    spec_matched = specificity.loc[overlap]
    magma_matched = magma_df.set_index("GENE").loc[overlap].reset_index()

    covariates = build_covariates(magma_matched)
    zstat = magma_matched["ZSTAT"].values

    print(f"\n  Running enrichment for {spec_matched.shape[1]} clusters ({label})...")
    results = run_enrichment_all_types(spec_matched, zstat, covariates)
    results["source"] = label

    return results


def compare_enrichments(conti_enrich, repro_enrich):
    """Compare enrichment results for the same set of cluster names."""
    conti_types = set(conti_enrich["supertype"])
    repro_types = set(repro_enrich["supertype"])
    shared = sorted(conti_types & repro_types)
    print(f"\n  Conti clusters: {len(conti_types)}")
    print(f"  Reprocessed clusters: {len(repro_types)}")
    print(f"  Shared: {len(shared)}")

    rows = []
    for name in shared:
        c = conti_enrich[conti_enrich["supertype"] == name].iloc[0]
        r = repro_enrich[repro_enrich["supertype"] == name].iloc[0]
        rows.append({
            "cluster": name,
            "conti_beta": c["beta"], "conti_pvalue": c["p_value"],
            "conti_fdr": c["p_fdr"], "conti_sig": c["p_fdr"] < FDR_THRESHOLD,
            "repro_beta": r["beta"], "repro_pvalue": r["p_value"],
            "repro_fdr": r["p_fdr"], "repro_sig": r["p_fdr"] < FDR_THRESHOLD,
            "both_sig": (c["p_fdr"] < FDR_THRESHOLD) and (r["p_fdr"] < FDR_THRESHOLD),
        })

    comp = pd.DataFrame(rows).sort_values("conti_pvalue").reset_index(drop=True)

    n_c = comp["conti_sig"].sum()
    n_r = comp["repro_sig"].sum()
    n_b = comp["both_sig"].sum()
    print(f"\n  FDR-significant:")
    print(f"    Conti: {n_c}")
    print(f"    Reprocessed: {n_r}")
    print(f"    Both: {n_b}")
    print(f"    Conti only: {(comp['conti_sig'] & ~comp['repro_sig']).sum()}")
    print(f"    Reprocessed only: {(~comp['conti_sig'] & comp['repro_sig']).sum()}")

    return comp


def load_supercluster_info():
    """Load supercluster assignments for coloring the Manhattan plot."""
    d = np.load(
        str(INTERMEDIATES_DIR / "siletti_cluster_level_stats.npz"),
        allow_pickle=True,
    )
    names = d["cluster_names"]
    superclusters = d["cluster_superclusters"]
    return dict(zip(names, superclusters))


def plot_duncan_manhattan(conti_enrich, repro_enrich, name_to_num,
                          supercluster_map, output_path):
    """
    Replicate Duncan et al. Manhattan-style plot: cell type number on x-axis,
    -log10(p) on y-axis, colored by supercluster.
    """
    fig, axes = plt.subplots(2, 1, figsize=(24, 16), sharex=True)

    # Assign colors to superclusters
    all_sc = sorted(set(supercluster_map.values()))
    n_sc = len(all_sc)
    cmap = plt.cm.get_cmap("tab20", min(n_sc, 20))
    # Extend with tab20b/tab20c for >20 superclusters
    if n_sc > 20:
        cmap2 = plt.cm.get_cmap("tab20b", 20)
        colors_list = [cmap(i % 20) for i in range(20)] + [cmap2(i % 20) for i in range(n_sc - 20)]
    else:
        colors_list = [cmap(i) for i in range(n_sc)]
    sc_colors = {sc: colors_list[i] for i, sc in enumerate(all_sc)}

    for ax, (enrich, title) in zip(axes, [
        (conti_enrich, "Conti (Duncan et al. 2025) — original Siletti processing"),
        (repro_enrich, "Reprocessed — updated Allen Institute Siletti data"),
    ]):
        enrich = enrich.copy()
        enrich["cluster_num"] = enrich["supertype"].map(name_to_num)
        enrich["logp"] = -np.log10(np.clip(enrich["p_value"], 1e-50, 1))
        enrich["supercluster"] = enrich["supertype"].map(supercluster_map)
        enrich["color"] = enrich["supercluster"].map(sc_colors)
        enrich = enrich.dropna(subset=["cluster_num"])
        enrich = enrich.sort_values("cluster_num")

        # Scatter plot
        for sc in all_sc:
            mask = enrich["supercluster"] == sc
            if mask.sum() == 0:
                continue
            sub = enrich[mask]
            ax.scatter(sub["cluster_num"], sub["logp"],
                       c=[sc_colors[sc]] * len(sub), s=30, alpha=0.7,
                       edgecolors="black", linewidth=0.3)

        # Significance line
        bonf_line = -np.log10(0.05 / len(enrich))
        ax.axhline(bonf_line, color="gray", linestyle="-", linewidth=1, alpha=0.5)

        # Label top hits
        top = enrich.nsmallest(15, "p_value")
        for _, row in top.iterrows():
            num = int(row["cluster_num"])
            ax.annotate(f"{num}", (row["cluster_num"], row["logp"]),
                        fontsize=8, fontweight="bold",
                        xytext=(0, 6), textcoords="offset points",
                        ha="center")

        n_sig = (enrich["p_fdr"] < FDR_THRESHOLD).sum()
        n_bonf = (enrich["p_bonferroni"] < 0.05).sum()
        ax.set_title(f"{title}\n"
                     f"{n_bonf} Bonferroni-sig, {n_sig} FDR-sig / {len(enrich)} total",
                     fontsize=TITLE_FONTSIZE)
        ax.set_ylabel("-log$_{10}$(P)", fontsize=LABEL_FONTSIZE)
        ax.tick_params(labelsize=TICK_FONTSIZE)
        ax.set_xlim(-5, 465)

    axes[1].set_xlabel("Cell type number from Siletti et al.", fontsize=LABEL_FONTSIZE)

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def plot_comparison(comparison, output_path):
    """Scatter comparison of betas and p-values."""
    fig, axes = plt.subplots(1, 2, figsize=(18, 8))

    mask_both = comparison["both_sig"]
    mask_conti = comparison["conti_sig"] & ~comparison["repro_sig"]
    mask_repro = ~comparison["conti_sig"] & comparison["repro_sig"]
    mask_neither = ~comparison["conti_sig"] & ~comparison["repro_sig"]

    for ax, (x_col, y_col, xlabel, ylabel, title_prefix) in zip(axes, [
        ("conti_beta", "repro_beta",
         "Conti (Duncan et al.) beta", "Reprocessed beta", "Enrichment beta"),
        ("conti_pvalue", "repro_pvalue",
         "Conti -log10(p)", "Reprocessed -log10(p)", "Significance"),
    ]):
        if "pvalue" in x_col:
            xv = -np.log10(np.clip(comparison[x_col], 1e-50, 1))
            yv = -np.log10(np.clip(comparison[y_col], 1e-50, 1))
        else:
            xv = comparison[x_col]
            yv = comparison[y_col]

        ax.scatter(xv[mask_neither], yv[mask_neither],
                   c="gray", alpha=0.4, s=40, label=f"Neither ({mask_neither.sum()})")
        ax.scatter(xv[mask_conti], yv[mask_conti],
                   c="tab:blue", alpha=0.6, s=60, label=f"Conti only ({mask_conti.sum()})")
        ax.scatter(xv[mask_repro], yv[mask_repro],
                   c="tab:orange", alpha=0.6, s=60, label=f"Repro only ({mask_repro.sum()})")
        ax.scatter(xv[mask_both], yv[mask_both],
                   c="tab:red", alpha=0.8, s=80, label=f"Both ({mask_both.sum()})")

        lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
                max(ax.get_xlim()[1], ax.get_ylim()[1])]
        ax.plot(lims, lims, "k--", alpha=0.3, zorder=0)

        r, p = stats.spearmanr(xv, yv)
        ax.set_title(f"{title_prefix} comparison\nSpearman r={r:.3f}, p={p:.2e}",
                     fontsize=TITLE_FONTSIZE)
        ax.set_xlabel(xlabel, fontsize=LABEL_FONTSIZE)
        ax.set_ylabel(ylabel, fontsize=LABEL_FONTSIZE)
        ax.tick_params(labelsize=TICK_FONTSIZE)
        ax.legend(fontsize=LEGEND_FONTSIZE - 1, loc="upper left")

        # Label top concordant hits
        if "pvalue" in x_col:
            top = comparison[mask_both].nsmallest(8, "conti_pvalue")
            for _, row in top.iterrows():
                x = -np.log10(max(row["conti_pvalue"], 1e-50))
                y = -np.log10(max(row["repro_pvalue"], 1e-50))
                ax.annotate(row["cluster"], (x, y), fontsize=8,
                            xytext=(5, 5), textcoords="offset points")

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def main():
    ensure_dirs(TABLES_DIR, FIGURES_DIR)

    print("=" * 60)
    print("Step 10: Conti vs Reprocessed Siletti (Apples-to-Apples)")
    print("=" * 60)
    t0 = time.time()

    # --- Cluster mappings ---
    print("\n1. Building cluster mappings...")
    name_to_num, num_to_name = build_cluster_number_map()
    supercluster_map = load_supercluster_info()
    print(f"  {len(name_to_num)} clusters, {len(set(supercluster_map.values()))} superclusters")

    # --- Conti enrichment ---
    conti_output = TABLES_DIR / "conti_siletti_enrichment.csv"
    if conti_output.exists():
        print(f"\n2. Loading cached Conti enrichment...")
        conti_enrich = pd.read_csv(str(conti_output))
    else:
        print("\n2. Loading Conti specificity and running enrichment...")
        conti_spec = load_conti_specificity(num_to_name)
        conti_enrich = run_enrichment_entrez(conti_spec, "Conti")
        conti_enrich.to_csv(str(conti_output), index=False)
    n_sig = (conti_enrich["p_fdr"] < FDR_THRESHOLD).sum()
    print(f"  Conti: {n_sig}/{len(conti_enrich)} FDR-significant")

    # --- Reprocessed enrichment (461 clusters, ENTREZ) ---
    repro_output = TABLES_DIR / "reprocessed_siletti_461_enrichment.csv"
    if repro_output.exists():
        print(f"\n3. Loading cached reprocessed enrichment...")
        repro_enrich = pd.read_csv(str(repro_output))
    else:
        print("\n3. Loading reprocessed specificity and running enrichment...")
        repro_spec = load_reprocessed_specificity()
        repro_enrich = run_enrichment_entrez(repro_spec, "Reprocessed")
        repro_enrich.to_csv(str(repro_output), index=False)
        print(f"  Saved: {repro_output}")
    n_sig2 = (repro_enrich["p_fdr"] < FDR_THRESHOLD).sum()
    print(f"  Reprocessed: {n_sig2}/{len(repro_enrich)} FDR-significant")

    # --- Compare ---
    print("\n4. Comparing enrichments (apples-to-apples)...")
    comparison = compare_enrichments(conti_enrich, repro_enrich)
    comp_output = TABLES_DIR / "conti_vs_reprocessed_461_comparison.csv"
    comparison.to_csv(str(comp_output), index=False)
    print(f"  Saved: {comp_output}")

    # --- Figures ---
    print("\n5. Generating figures...")

    plot_comparison(
        comparison,
        FIGURES_DIR / "conti_vs_reprocessed_461_enrichment.png",
    )

    plot_duncan_manhattan(
        conti_enrich, repro_enrich, name_to_num, supercluster_map,
        FIGURES_DIR / "duncan_manhattan_replication.png",
    )

    print(f"\nStep 10 complete. Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
