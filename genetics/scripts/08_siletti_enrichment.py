#!/usr/bin/env python
"""
08_siletti_enrichment.py — Compute and compare SCZ enrichment for Siletti cell types.

Runs enrichment on two versions of the Siletti et al. (2023) atlas:
  1. Conti/Duncan original (precomputed specificity, ENTREZ IDs)
  2. Reprocessed from updated Allen Institute data (computed here, ENTREZ IDs)

Both use all 461 clusters and ENTREZ gene IDs directly (no symbol mapping loss).
Compares enrichment between the two versions.

Inputs:
  - data/conti_specificity_matrix.txt (original Siletti specificity)
  - results/intermediates/siletti_cluster_level_stats.npz (reprocessed cluster means)
  - data/adult_human_20221007.loom (gene ID mapping)
  - linking_.../NCBI37.3.gene.loc.extendedMHCexcluded (ENTREZ mapping)
  - linking_.../PGC3_SCZ...genes.out (MAGMA results)

Outputs:
  - results/tables/conti_siletti_enrichment.csv
  - results/tables/reprocessed_siletti_461_enrichment.csv
  - results/intermediates/reprocessed_siletti_specificity_entrez.csv
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
from scipy import stats

from scz_celltype_enrichment.config import (
    MAGMA_GENES_OUT, GENE_LOC_FILE, TABLES_DIR, FIGURES_DIR, INTERMEDIATES_DIR,
    FIGURE_DPI, TITLE_FONTSIZE, LABEL_FONTSIZE, TICK_FONTSIZE, LEGEND_FONTSIZE,
    FDR_THRESHOLD, PROJECT_ROOT,
)
from scz_celltype_enrichment.enrichment.gwas import load_magma_genes, build_covariates
from scz_celltype_enrichment.enrichment.celltype import run_enrichment_all_types
from scz_celltype_enrichment.siletti.specificity import (
    build_ensembl_to_entrez_mapping,
    compute_siletti_specificity_entrez,
    load_conti_specificity,
)
from scz_celltype_enrichment.siletti.clusters import (
    build_cluster_number_map,
    load_supercluster_info,
)
from scz_celltype_enrichment.utils import ensure_dirs, Timer


CONTI_SPECIFICITY = PROJECT_ROOT / "data" / "conti_specificity_matrix.txt"
LOOM_PATH = PROJECT_ROOT / "data" / "adult_human_20221007.loom"
CLUSTER_STATS = INTERMEDIATES_DIR / "siletti_cluster_level_stats.npz"


def run_enrichment_entrez(specificity, label):
    """Run MAGMA-style enrichment on an ENTREZ-indexed specificity matrix."""
    magma_df = load_magma_genes(MAGMA_GENES_OUT)
    overlap = sorted(set(specificity.index) & set(magma_df["GENE"]))
    print(f"  Overlapping genes (ENTREZ): {len(overlap):,}")

    spec_matched = specificity.loc[overlap]
    magma_matched = magma_df.set_index("GENE").loc[overlap].reset_index()
    covariates = build_covariates(magma_matched)
    zstat = magma_matched["ZSTAT"].values

    results = run_enrichment_all_types(spec_matched, zstat, covariates)
    results["source"] = label
    return results


def compare_enrichments(conti_enrich, repro_enrich):
    """Compare enrichment between Conti and reprocessed Siletti."""
    shared = sorted(set(conti_enrich["supertype"]) & set(repro_enrich["supertype"]))
    print(f"  Shared clusters: {len(shared)}")

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
    n_b = comp["both_sig"].sum()
    print(f"  Both FDR-sig: {n_b}, Conti-only: {(comp['conti_sig'] & ~comp['repro_sig']).sum()}, "
          f"Repro-only: {(~comp['conti_sig'] & comp['repro_sig']).sum()}")
    return comp


def plot_comparison(comparison, output_path):
    """Scatter comparison of betas and p-values."""
    fig, axes = plt.subplots(1, 2, figsize=(18, 8))

    mask_both = comparison["both_sig"]
    mask_conti = comparison["conti_sig"] & ~comparison["repro_sig"]
    mask_repro = ~comparison["conti_sig"] & comparison["repro_sig"]
    mask_neither = ~comparison["conti_sig"] & ~comparison["repro_sig"]

    for ax, (x_col, y_col, xlabel, ylabel) in zip(axes, [
        ("conti_beta", "repro_beta", "Conti (Duncan et al.) beta", "Reprocessed beta"),
        ("conti_pvalue", "repro_pvalue", "Conti -log10(p)", "Reprocessed -log10(p)"),
    ]):
        if "pvalue" in x_col:
            xv = -np.log10(np.clip(comparison[x_col], 1e-50, 1))
            yv = -np.log10(np.clip(comparison[y_col], 1e-50, 1))
        else:
            xv, yv = comparison[x_col], comparison[y_col]

        for mask, color, sz, label in [
            (mask_neither, "gray", 40, f"Neither ({mask_neither.sum()})"),
            (mask_conti, "tab:blue", 60, f"Conti only ({mask_conti.sum()})"),
            (mask_repro, "tab:orange", 60, f"Repro only ({mask_repro.sum()})"),
            (mask_both, "tab:red", 80, f"Both ({mask_both.sum()})"),
        ]:
            ax.scatter(xv[mask], yv[mask], c=color, s=sz, alpha=0.6, label=label)

        lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
                max(ax.get_xlim()[1], ax.get_ylim()[1])]
        ax.plot(lims, lims, "k--", alpha=0.3, zorder=0)

        r, p = stats.spearmanr(xv, yv)
        prefix = "Enrichment beta" if "beta" in x_col else "Significance"
        ax.set_title(f"{prefix} comparison\nSpearman r={r:.3f}, p={p:.2e}",
                     fontsize=TITLE_FONTSIZE)
        ax.set_xlabel(xlabel, fontsize=LABEL_FONTSIZE)
        ax.set_ylabel(ylabel, fontsize=LABEL_FONTSIZE)
        ax.tick_params(labelsize=TICK_FONTSIZE)
        ax.legend(fontsize=LEGEND_FONTSIZE - 1, loc="upper left")

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def plot_duncan_manhattan(conti_enrich, repro_enrich, name_to_num,
                          supercluster_map, output_path):
    """Replicate Duncan et al. Manhattan-style plot for both data versions."""
    fig, axes = plt.subplots(2, 1, figsize=(24, 16), sharex=True)

    all_sc = sorted(set(supercluster_map.values()))
    cmap = plt.colormaps.get_cmap("tab20")
    cmap2 = plt.colormaps.get_cmap("tab20b")
    colors_list = [cmap(i % 20) for i in range(20)]
    colors_list += [cmap2(i % 20) for i in range(max(0, len(all_sc) - 20))]
    sc_colors = {sc: colors_list[i] for i, sc in enumerate(all_sc)}

    for ax, (enrich, title) in zip(axes, [
        (conti_enrich, "Conti (Duncan et al. 2025) \u2014 original Siletti processing"),
        (repro_enrich, "Reprocessed \u2014 updated Allen Institute Siletti data"),
    ]):
        df = enrich.copy()
        df["cluster_num"] = df["supertype"].map(name_to_num)
        df["logp"] = -np.log10(np.clip(df["p_value"], 1e-300, 1)).clip(upper=50)
        df["sc"] = df["supertype"].map(supercluster_map)
        df = df.dropna(subset=["cluster_num"]).sort_values("cluster_num")

        for sc in all_sc:
            sub = df[df["sc"] == sc]
            if len(sub) == 0:
                continue
            ax.scatter(sub["cluster_num"], sub["logp"],
                       c=[sc_colors[sc]] * len(sub), s=30, alpha=0.7,
                       edgecolors="black", linewidth=0.3)

        bonf = -np.log10(0.05 / len(df))
        ax.axhline(bonf, color="gray", linestyle="-", linewidth=1, alpha=0.5)

        for _, row in df.nsmallest(15, "p_value").iterrows():
            ax.annotate(f"{int(row['cluster_num'])}", (row["cluster_num"], row["logp"]),
                        fontsize=8, fontweight="bold",
                        xytext=(0, 6), textcoords="offset points", ha="center")

        n_fdr = (df["p_fdr"] < FDR_THRESHOLD).sum()
        n_bonf = (df["p_bonferroni"] < 0.05).sum()
        ax.set_title(f"{title}\n{n_bonf} Bonferroni-sig, {n_fdr} FDR-sig / {len(df)} total",
                     fontsize=TITLE_FONTSIZE)
        ax.set_ylabel("-log$_{10}$(P)", fontsize=LABEL_FONTSIZE)
        ax.tick_params(labelsize=TICK_FONTSIZE)
        ax.set_xlim(-5, 465)

    axes[1].set_xlabel("Cell type number from Siletti et al.", fontsize=LABEL_FONTSIZE)
    plt.tight_layout()
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def main():
    ensure_dirs(TABLES_DIR, FIGURES_DIR, INTERMEDIATES_DIR)

    print("=" * 60)
    print("Step 08: Siletti Enrichment (Conti vs Reprocessed)")
    print("=" * 60)
    t0 = time.time()

    # --- Cluster mappings ---
    print("\n1. Building cluster mappings...")
    name_to_num, num_to_name = build_cluster_number_map(CLUSTER_STATS)
    supercluster_map = load_supercluster_info(CLUSTER_STATS)

    # --- Conti enrichment ---
    conti_output = TABLES_DIR / "conti_siletti_enrichment.csv"
    if conti_output.exists():
        print(f"\n2. Loading cached Conti enrichment...")
        conti_enrich = pd.read_csv(str(conti_output))
    else:
        print("\n2. Running Conti enrichment...")
        conti_spec = load_conti_specificity(CONTI_SPECIFICITY, num_to_name)
        conti_enrich = run_enrichment_entrez(conti_spec, "Conti")
        conti_enrich.to_csv(str(conti_output), index=False)
        print(f"  Saved: {conti_output}")
    print(f"  Conti: {(conti_enrich['p_fdr'] < FDR_THRESHOLD).sum()}/{len(conti_enrich)} FDR-sig")

    # --- Reprocessed enrichment ---
    repro_output = TABLES_DIR / "reprocessed_siletti_461_enrichment.csv"
    if repro_output.exists():
        print(f"\n3. Loading cached reprocessed enrichment...")
        repro_enrich = pd.read_csv(str(repro_output))
    else:
        print("\n3. Computing reprocessed Siletti specificity and enrichment...")
        ens_to_entrez = build_ensembl_to_entrez_mapping(LOOM_PATH, GENE_LOC_FILE)
        repro_spec = compute_siletti_specificity_entrez(CLUSTER_STATS, ens_to_entrez)
        repro_spec.to_csv(str(INTERMEDIATES_DIR / "reprocessed_siletti_specificity_entrez.csv"))
        repro_enrich = run_enrichment_entrez(repro_spec, "Reprocessed")
        repro_enrich.to_csv(str(repro_output), index=False)
        print(f"  Saved: {repro_output}")
    print(f"  Reprocessed: {(repro_enrich['p_fdr'] < FDR_THRESHOLD).sum()}/{len(repro_enrich)} FDR-sig")

    # --- Compare ---
    print("\n4. Comparing enrichments...")
    comparison = compare_enrichments(conti_enrich, repro_enrich)
    comparison.to_csv(str(TABLES_DIR / "conti_vs_reprocessed_461_comparison.csv"), index=False)

    # --- Figures ---
    print("\n5. Generating figures...")
    plot_comparison(comparison, FIGURES_DIR / "conti_vs_reprocessed_461_enrichment.png")
    plot_duncan_manhattan(conti_enrich, repro_enrich, name_to_num, supercluster_map,
                          FIGURES_DIR / "duncan_manhattan_replication.png")

    print(f"\nStep 08 complete. Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
