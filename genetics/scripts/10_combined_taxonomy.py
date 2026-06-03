#!/usr/bin/env python
"""
10_combined_taxonomy.py — Build and evaluate the RBH combined taxonomy.

Uses MetaNeighbor reciprocal best hits to merge SEA-AD and Siletti:
  - Keep all 137 SEA-AD types
  - Exclude Siletti clusters that are RBH matches (already represented)
  - Include remaining Siletti clusters as novel types

Runs enrichment on the combined taxonomy, compares with the old 20%-cortical
approach, and generates publication figures + webapp data.

Inputs:
  - results/tables/metaneighbor_full461_all_matches.csv
  - results/intermediates/full_seaad_siletti_461_mean_expression.csv
  - linking_.../PGC3_SCZ...genes.out (MAGMA)

Outputs:
  - results/tables/rbh_combined_enrichment.csv
  - results/tables/rbh_combined_enrichment_app.csv (webapp format)
  - results/tables/rbh_vs_old_combined_comparison.csv
  - results/intermediates/rbh_combined_type_order_info.csv
  - results/figures/rbh_combined_taxonomy_manhattan.png
  - results/figures/rbh_combined_taxonomy_top_enrichments.png
  - results/figures/rbh_vs_old_combined_taxonomy.png
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
    FDR_THRESHOLD,
)
from scz_celltype_enrichment.enrichment.gwas import (
    load_magma_genes, load_entrez_to_symbol, map_and_merge, build_covariates,
)
from scz_celltype_enrichment.enrichment.celltype import run_enrichment_all_types
from scz_celltype_enrichment.siletti.metaneighbor import build_rbh_combined_taxonomy
from scz_celltype_enrichment.siletti.clusters import (
    load_supercluster_info, classify_siletti_cluster,
)
from scz_celltype_enrichment.utils import (
    ensure_dirs, Timer, classify_supertypes_by_class, extract_subclass,
)


def compute_combined_specificity(mean_expr, seaad_cols, novel_siletti):
    """Compute specificity for the combined taxonomy from mean expression."""
    cols = list(seaad_cols) + [f"Siletti_{c}" for c in novel_siletti]
    sub = mean_expr[cols]

    # Column normalize (sum to 1000), then row normalize (sum to 1.0)
    col_sums = sub.sum(axis=0)
    col_norm = sub * (1000.0 / col_sums)
    row_sums = col_norm.sum(axis=1)
    nonzero = row_sums > 0
    spec = col_norm.loc[nonzero].div(row_sums[nonzero], axis=0)

    # Clean column names
    spec.columns = [c.replace("Siletti_", "") if c.startswith("Siletti_") else c
                    for c in spec.columns]
    return spec


def build_type_order(seaad_types, novel_siletti, supercluster_map):
    """Build type ordering metadata for the webapp.

    SEA-AD types are ordered by subclass (GABAergic → Glutamatergic → Non-neuronal).
    Siletti types are ordered by supercluster, then alphabetically within each.
    """
    rows = []
    for t in seaad_types:
        cls = classify_supertypes_by_class([t])
        rows.append({
            "cell_type": t, "source": "SEA-AD",
            "class": cls.get(t, "Other"),
            "subclass": extract_subclass(t),
            "supercluster": "",
        })

    # Sort Siletti by supercluster, then by name within each supercluster
    siletti_with_sc = [(t, supercluster_map.get(t, "ZZZ_Unknown")) for t in novel_siletti]
    siletti_with_sc.sort(key=lambda x: (x[1], x[0]))

    for t, sc in siletti_with_sc:
        rows.append({
            "cell_type": t, "source": "Siletti",
            "class": "Siletti",
            "subclass": "",
            "supercluster": sc,
        })
    return pd.DataFrame(rows)


def plot_manhattan(enrich, novel_siletti, supercluster_map, output_path):
    """Manhattan-style plot of the RBH combined taxonomy."""
    from scz_celltype_enrichment.utils import extract_subclass

    seaad_colors = {"GABAergic": "#e74c3c", "Glutamatergic": "#3498db", "Non-neuronal": "#2ecc71"}
    gaba_prefixes = ["Chandelier", "Lamp5", "Lamp5 Lhx6", "Pax6", "Pvalb",
                     "Sncg", "Sst", "Sst Chodl", "Vip"]
    glut_prefixes = ["L2/3 IT", "L4 IT", "L5 ET", "L5 IT", "L5/6 NP",
                     "L6 CT", "L6 IT", "L6 IT Car3", "L6b"]

    def seaad_sort_key(name):
        sc = extract_subclass(name)
        for i, p in enumerate(gaba_prefixes):
            if sc == p:
                return (0, i, name)
        for i, p in enumerate(glut_prefixes):
            if sc == p:
                return (1, i, name)
        return (2, 0, name)

    df = enrich.copy()
    df["logp"] = -np.log10(np.clip(df["p_value"], 1e-300, 1)).clip(upper=50)

    seaad_df = df[df["source"] == "SEA-AD"].copy()
    seaad_df = seaad_df.sort_values("supertype", key=lambda s: s.map(seaad_sort_key))
    siletti_df = df[df["source"] == "Siletti"].copy()
    siletti_df["_sc"] = siletti_df["supertype"].map(lambda x: supercluster_map.get(x, "ZZZ"))
    siletti_df = siletti_df.sort_values(["_sc", "supertype"]).drop(columns=["_sc"])
    ordered = pd.concat([seaad_df, siletti_df]).reset_index(drop=True)
    ordered["x"] = range(len(ordered))

    def get_color(row):
        if row["source"] == "SEA-AD":
            cls = classify_supertypes_by_class([row["supertype"]])
            return seaad_colors.get(cls.get(row["supertype"], ""), "gray")
        return "#8b5cf6"

    ordered["color"] = ordered.apply(get_color, axis=1)

    fig, ax = plt.subplots(figsize=(24, 10))
    for _, row in ordered.iterrows():
        marker = "o" if row["source"] == "SEA-AD" else "s"
        size = 50 if row["p_fdr"] < FDR_THRESHOLD else 20
        alpha = 0.8 if row["p_fdr"] < FDR_THRESHOLD else 0.4
        ax.scatter(row["x"], row["logp"], c=[row["color"]], s=size, alpha=alpha,
                   marker=marker, edgecolors="black", linewidth=0.3)

    bonf = -np.log10(0.05 / len(ordered))
    ax.axhline(bonf, color="gray", linestyle="-", linewidth=1, alpha=0.5)
    sep = len(seaad_df) - 0.5
    ax.axvline(sep, color="black", linestyle="--", linewidth=1.5, alpha=0.5)
    ax.text(sep / 2, 49, "SEA-AD (137)", ha="center", fontsize=LABEL_FONTSIZE, fontweight="bold")
    ax.text(sep + (len(ordered) - sep) / 2, 49, f"Novel Siletti ({len(siletti_df)})",
            ha="center", fontsize=LABEL_FONTSIZE, fontweight="bold")

    for _, row in ordered.nsmallest(20, "p_value").iterrows():
        y_off = 8 if row["logp"] < 45 else -12
        ax.annotate(row["supertype"], (row["x"], row["logp"]), fontsize=7,
                    fontweight="bold", rotation=45, ha="left", va="bottom",
                    xytext=(3, y_off), textcoords="offset points")

    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#e74c3c", markersize=10, label="SEA-AD GABAergic"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#3498db", markersize=10, label="SEA-AD Glutamatergic"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#2ecc71", markersize=10, label="SEA-AD Non-neuronal"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#8b5cf6", markersize=10, label="Novel Siletti"),
    ]
    ax.legend(handles=legend_elements, fontsize=LEGEND_FONTSIZE, loc="upper right")

    n_sig = (ordered["p_fdr"] < FDR_THRESHOLD).sum()
    n_bonf = (ordered["p_bonferroni"] < 0.05).sum()
    ax.set_title(f"RBH Combined Taxonomy: SCZ Cell-Type Enrichment\n"
                 f"{n_bonf} Bonferroni-sig, {n_sig} FDR-sig / {len(ordered)} total",
                 fontsize=TITLE_FONTSIZE)
    ax.set_xlabel("Cell type (SEA-AD \u2192 Novel Siletti)", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("-log$_{10}$(P)", fontsize=LABEL_FONTSIZE)
    ax.tick_params(labelsize=TICK_FONTSIZE)
    ax.set_xticks([])
    ax.set_xlim(-3, len(ordered) + 3)
    ax.set_ylim(-1, 53)

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def main():
    ensure_dirs(TABLES_DIR, FIGURES_DIR, INTERMEDIATES_DIR)

    print("=" * 60)
    print("Step 10: RBH Combined Taxonomy")
    print("=" * 60)
    t0 = time.time()

    # --- Load MetaNeighbor results ---
    print("\n1. Loading MetaNeighbor matches...")
    all_matches = pd.read_csv(str(TABLES_DIR / "metaneighbor_full461_all_matches.csv"))

    # --- Load mean expression ---
    print("\n2. Loading combined mean expression...")
    with Timer("Loading"):
        mean_expr = pd.read_csv(
            str(INTERMEDIATES_DIR / "full_seaad_siletti_461_mean_expression.csv"), index_col=0
        )
    seaad_cols = [c for c in mean_expr.columns if not c.startswith("Siletti_")]
    siletti_cols = [c.replace("Siletti_", "") for c in mean_expr.columns if c.startswith("Siletti_")]

    supercluster_map = load_supercluster_info(INTERMEDIATES_DIR / "siletti_cluster_level_stats.npz")

    # --- Build RBH taxonomy ---
    print("\n3. Building RBH combined taxonomy...")
    novel_siletti, claimed = build_rbh_combined_taxonomy(all_matches, seaad_cols, siletti_cols)

    # --- Compute specificity and run enrichment ---
    rbh_output = TABLES_DIR / "rbh_combined_enrichment.csv"
    if rbh_output.exists():
        print(f"\n4. Loading cached RBH enrichment...")
        rbh_enrich = pd.read_csv(str(rbh_output))
    else:
        print("\n4. Computing specificity and running enrichment...")
        spec = compute_combined_specificity(mean_expr, seaad_cols, novel_siletti)
        magma_df = load_magma_genes(MAGMA_GENES_OUT)
        entrez_to_symbol = load_entrez_to_symbol(GENE_LOC_FILE)
        spec_matched, gwas_matched = map_and_merge(magma_df, spec, entrez_to_symbol)
        covariates = build_covariates(gwas_matched)
        zstat = gwas_matched["ZSTAT"].values
        rbh_enrich = run_enrichment_all_types(spec_matched, zstat, covariates)
        rbh_enrich["source"] = rbh_enrich["supertype"].apply(
            lambda x: "Siletti" if x in novel_siletti else "SEA-AD"
        )
        rbh_enrich.to_csv(str(rbh_output), index=False)
        print(f"  Saved: {rbh_output}")

    n_sig = (rbh_enrich["p_fdr"] < FDR_THRESHOLD).sum()
    print(f"  RBH combined: {n_sig} FDR-sig / {len(rbh_enrich)} total")

    # --- Prepare webapp data ---
    print("\n5. Preparing webapp data...")
    app_enrich = rbh_enrich.copy().rename(columns={"supertype": "cell_type"})
    app_enrich["supercluster"] = app_enrich.apply(
        lambda r: supercluster_map.get(r["cell_type"], "") if r["source"] == "Siletti" else "",
        axis=1,
    )
    app_enrich.to_csv(str(TABLES_DIR / "rbh_combined_enrichment_app.csv"), index=False)

    type_order = build_type_order(seaad_cols, novel_siletti, supercluster_map)
    type_order.to_csv(str(INTERMEDIATES_DIR / "rbh_combined_type_order_info.csv"), index=False)

    # --- Compare with old taxonomy ---
    print("\n6. Comparing with old combined taxonomy...")
    old_path = TABLES_DIR / "combined_cluster_enrichment.csv"
    if old_path.exists():
        old_enrich = pd.read_csv(str(old_path))
        # Compare SEA-AD betas
        comp_rows = []
        for t in seaad_cols:
            old_row = old_enrich[old_enrich["cell_type"] == t]
            rbh_row = rbh_enrich[rbh_enrich["supertype"] == t]
            if len(old_row) > 0 and len(rbh_row) > 0:
                comp_rows.append({
                    "type": t, "source": "SEA-AD",
                    "old_beta": old_row.iloc[0]["beta"],
                    "old_fdr": old_row.iloc[0]["p_fdr"],
                    "rbh_beta": rbh_row.iloc[0]["beta"],
                    "rbh_fdr": rbh_row.iloc[0]["p_fdr"],
                })
        comp_df = pd.DataFrame(comp_rows)
        r, p = stats.spearmanr(comp_df["old_beta"], comp_df["rbh_beta"])
        print(f"  SEA-AD beta correlation (old vs RBH): r={r:.3f}")
        comp_df.to_csv(str(TABLES_DIR / "rbh_vs_old_combined_comparison.csv"), index=False)

    # --- Figures ---
    print("\n7. Generating figures...")
    plot_manhattan(rbh_enrich, novel_siletti, supercluster_map,
                   FIGURES_DIR / "rbh_combined_taxonomy_manhattan.png")

    print(f"\nStep 10 complete. Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
