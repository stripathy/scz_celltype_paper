#!/usr/bin/env python
"""
08_three_taxonomy_enrichment.py — Compare SCZ enrichment across 3 taxonomies.

Runs enrichment analysis using three independent cell-type taxonomies:
  1. SEA-AD only (137 supertypes) — loaded from existing results
  2. Siletti only (348 clusters) — computed here from combined specificity
  3. Combined SEA-AD + Siletti (485 types) — loaded from existing results

Then compares enrichment results across taxonomies using matched cell types
from MetaNeighbor.

Outputs:
  - results/tables/siletti_only_enrichment.csv
  - results/tables/three_taxonomy_enrichment_comparison.csv
  - results/figures/three_taxonomy_enrichment_comparison.png
  - results/figures/taxonomy_enrichment_scatter.png
  - results/figures/taxonomy_significant_types_barplot.png
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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
from scz_celltype_enrichment.utils import ensure_dirs, Timer


def compute_siletti_only_enrichment():
    """Extract Siletti-only columns, re-normalize, and run enrichment."""
    print("=" * 60)
    print("Computing Siletti-only enrichment")
    print("=" * 60)

    # Load combined specificity
    print("\n1. Loading combined specificity matrix...")
    combined_spec = pd.read_csv(
        str(INTERMEDIATES_DIR / "combined_cluster_specificity.csv"), index_col=0
    )
    print(f"  Combined shape: {combined_spec.shape}")

    # Extract Siletti-only columns
    siletti_cols = [c for c in combined_spec.columns if c.startswith("Siletti_")]
    print(f"  Siletti clusters: {len(siletti_cols)}")

    siletti_spec = combined_spec[siletti_cols].copy()
    # Strip Siletti_ prefix for clean names
    siletti_spec.columns = [c.replace("Siletti_", "") for c in siletti_spec.columns]

    # Re-normalize rows to sum to 1.0
    print("\n2. Re-normalizing specificity (rows sum to 1.0)...")
    row_sums = siletti_spec.sum(axis=1)
    nonzero = row_sums > 0
    n_dropped = (~nonzero).sum()
    if n_dropped > 0:
        print(f"  Dropping {n_dropped} genes with zero expression across Siletti clusters")
    siletti_spec = siletti_spec.loc[nonzero].div(row_sums[nonzero], axis=0)

    # Validate
    row_check = siletti_spec.sum(axis=1)
    assert np.allclose(row_check, 1.0, atol=1e-6), "Row normalization failed"
    print(f"  Siletti specificity: {siletti_spec.shape[0]:,} genes × {siletti_spec.shape[1]} clusters")

    # Load GWAS data and run enrichment
    print("\n3. Loading GWAS data...")
    magma_df = load_magma_genes(MAGMA_GENES_OUT)
    entrez_to_symbol = load_entrez_to_symbol(GENE_LOC_FILE)

    print("\n4. Mapping and intersecting...")
    spec_matched, gwas_matched = map_and_merge(magma_df, siletti_spec, entrez_to_symbol)

    print("\n5. Building covariates...")
    covariates = build_covariates(gwas_matched)
    zstat = gwas_matched["ZSTAT"].values

    print("\n6. Running enrichment for all Siletti clusters...")
    results = run_enrichment_all_types(spec_matched, zstat, covariates)

    # Add source column
    results["source"] = "Siletti"

    return results


def load_seaad_enrichment():
    """Load existing SEA-AD enrichment results."""
    path = TABLES_DIR / "seaad_magma_scz_enrichment_all_supertypes.csv"
    df = pd.read_csv(str(path))
    df["source"] = "SEA-AD"
    print(f"  SEA-AD enrichment: {len(df)} types, "
          f"{(df['p_fdr'] < FDR_THRESHOLD).sum()} FDR-significant")
    return df


def load_combined_enrichment():
    """Load existing combined enrichment results."""
    path = TABLES_DIR / "combined_cluster_enrichment.csv"
    df = pd.read_csv(str(path))
    print(f"  Combined enrichment: {len(df)} types, "
          f"{(df['p_fdr'] < FDR_THRESHOLD).sum()} FDR-significant")
    return df


def load_metaneighbor_matches():
    """Load MetaNeighbor best matches for cross-taxonomy comparison."""
    path = INTERMEDIATES_DIR / "seaad_to_siletti_metaneighbor.csv"
    df = pd.read_csv(str(path))
    return df


def build_comparison_table(seaad_enrich, siletti_enrich, combined_enrich, mn_matches):
    """
    Build a comparison table linking SEA-AD types to their best Siletti match.
    """
    rows = []
    for _, match in mn_matches.iterrows():
        seaad_name = match["seaad_supertype"]
        siletti_name = match["best_siletti_match"]
        auroc = match["best_siletti_auroc"]

        # Get SEA-AD enrichment stats
        seaad_row = seaad_enrich[seaad_enrich["supertype"] == seaad_name]
        if len(seaad_row) == 0:
            continue
        seaad_row = seaad_row.iloc[0]

        # Get Siletti enrichment stats
        siletti_row = siletti_enrich[siletti_enrich["supertype"] == siletti_name]
        if len(siletti_row) == 0:
            continue
        siletti_row = siletti_row.iloc[0]

        # Get combined enrichment stats for both
        comb_seaad = combined_enrich[combined_enrich["cell_type"] == seaad_name]
        comb_siletti = combined_enrich[
            combined_enrich["cell_type"] == f"Siletti_{siletti_name}"
        ]

        rows.append({
            "seaad_type": seaad_name,
            "siletti_match": siletti_name,
            "metaneighbor_auroc": auroc,
            "seaad_beta": seaad_row["beta"],
            "seaad_pvalue": seaad_row["p_value"],
            "seaad_fdr": seaad_row["p_fdr"],
            "siletti_beta": siletti_row["beta"],
            "siletti_pvalue": siletti_row["p_value"],
            "siletti_fdr": siletti_row["p_fdr"],
            "combined_seaad_beta": comb_seaad.iloc[0]["beta"] if len(comb_seaad) > 0 else np.nan,
            "combined_seaad_fdr": comb_seaad.iloc[0]["p_fdr"] if len(comb_seaad) > 0 else np.nan,
            "combined_siletti_beta": comb_siletti.iloc[0]["beta"] if len(comb_siletti) > 0 else np.nan,
            "combined_siletti_fdr": comb_siletti.iloc[0]["p_fdr"] if len(comb_siletti) > 0 else np.nan,
            "seaad_sig": seaad_row["p_fdr"] < FDR_THRESHOLD,
            "siletti_sig": siletti_row["p_fdr"] < FDR_THRESHOLD,
            "both_sig": (seaad_row["p_fdr"] < FDR_THRESHOLD) and (siletti_row["p_fdr"] < FDR_THRESHOLD),
        })

    comparison = pd.DataFrame(rows)
    comparison = comparison.sort_values("seaad_pvalue").reset_index(drop=True)

    # Summary stats
    n_total = len(comparison)
    n_seaad_sig = comparison["seaad_sig"].sum()
    n_siletti_sig = comparison["siletti_sig"].sum()
    n_both_sig = comparison["both_sig"].sum()
    print(f"\n  Comparison table: {n_total} matched pairs")
    print(f"  SEA-AD significant: {n_seaad_sig}")
    print(f"  Siletti significant: {n_siletti_sig}")
    print(f"  Both significant: {n_both_sig}")

    return comparison


def plot_enrichment_scatter(comparison, output_path):
    """Scatter plot of SEA-AD vs Siletti enrichment betas for matched types."""
    fig, axes = plt.subplots(1, 2, figsize=(18, 8))

    # --- Panel A: Beta comparison ---
    ax = axes[0]
    mask_both = comparison["both_sig"]
    mask_seaad_only = comparison["seaad_sig"] & ~comparison["siletti_sig"]
    mask_siletti_only = ~comparison["seaad_sig"] & comparison["siletti_sig"]
    mask_neither = ~comparison["seaad_sig"] & ~comparison["siletti_sig"]

    ax.scatter(comparison.loc[mask_neither, "seaad_beta"],
               comparison.loc[mask_neither, "siletti_beta"],
               c="gray", alpha=0.4, s=40, label=f"Neither (n={mask_neither.sum()})")
    ax.scatter(comparison.loc[mask_seaad_only, "seaad_beta"],
               comparison.loc[mask_seaad_only, "siletti_beta"],
               c="tab:blue", alpha=0.6, s=60, label=f"SEA-AD only (n={mask_seaad_only.sum()})")
    ax.scatter(comparison.loc[mask_siletti_only, "seaad_beta"],
               comparison.loc[mask_siletti_only, "siletti_beta"],
               c="tab:orange", alpha=0.6, s=60, label=f"Siletti only (n={mask_siletti_only.sum()})")
    ax.scatter(comparison.loc[mask_both, "seaad_beta"],
               comparison.loc[mask_both, "siletti_beta"],
               c="tab:red", alpha=0.8, s=80, label=f"Both (n={mask_both.sum()})")

    # Add identity line
    lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, "k--", alpha=0.3, zorder=0)

    # Correlation
    r, p = stats.spearmanr(comparison["seaad_beta"], comparison["siletti_beta"])
    ax.set_title(f"Enrichment beta comparison\nSpearman r={r:.3f}, p={p:.2e}",
                 fontsize=TITLE_FONTSIZE)
    ax.set_xlabel("SEA-AD enrichment beta", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("Siletti enrichment beta", fontsize=LABEL_FONTSIZE)
    ax.tick_params(labelsize=TICK_FONTSIZE)
    ax.legend(fontsize=LEGEND_FONTSIZE - 1, loc="upper left")

    # --- Panel B: -log10(p) comparison ---
    ax = axes[1]
    seaad_logp = -np.log10(np.clip(comparison["seaad_pvalue"], 1e-50, 1))
    siletti_logp = -np.log10(np.clip(comparison["siletti_pvalue"], 1e-50, 1))

    ax.scatter(seaad_logp[mask_neither], siletti_logp[mask_neither],
               c="gray", alpha=0.4, s=40, label=f"Neither (n={mask_neither.sum()})")
    ax.scatter(seaad_logp[mask_seaad_only], siletti_logp[mask_seaad_only],
               c="tab:blue", alpha=0.6, s=60, label=f"SEA-AD only (n={mask_seaad_only.sum()})")
    ax.scatter(seaad_logp[mask_siletti_only], siletti_logp[mask_siletti_only],
               c="tab:orange", alpha=0.6, s=60, label=f"Siletti only (n={mask_siletti_only.sum()})")
    ax.scatter(seaad_logp[mask_both], siletti_logp[mask_both],
               c="tab:red", alpha=0.8, s=80, label=f"Both (n={mask_both.sum()})")

    # Label top hits
    top_both = comparison[mask_both].nsmallest(5, "seaad_pvalue")
    for _, row in top_both.iterrows():
        x = -np.log10(max(row["seaad_pvalue"], 1e-50))
        y = -np.log10(max(row["siletti_pvalue"], 1e-50))
        ax.annotate(row["seaad_type"], (x, y), fontsize=9,
                    xytext=(5, 5), textcoords="offset points")

    r2, p2 = stats.spearmanr(seaad_logp, siletti_logp)
    ax.set_title(f"Significance comparison\nSpearman r={r2:.3f}, p={p2:.2e}",
                 fontsize=TITLE_FONTSIZE)
    ax.set_xlabel("SEA-AD -log10(p)", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("Siletti -log10(p)", fontsize=LABEL_FONTSIZE)
    ax.tick_params(labelsize=TICK_FONTSIZE)
    ax.legend(fontsize=LEGEND_FONTSIZE - 1, loc="upper left")

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def plot_significant_types_by_class(seaad_enrich, siletti_enrich, combined_enrich, output_path):
    """Bar plot showing number of FDR-significant types per class per taxonomy."""
    from scz_celltype_enrichment.utils import extract_subclass

    fig, axes = plt.subplots(1, 3, figsize=(24, 8))
    taxonomies = [
        ("SEA-AD only", seaad_enrich, "supertype"),
        ("Siletti only", siletti_enrich, "supertype"),
        ("Combined", combined_enrich, "cell_type"),
    ]

    for ax, (name, df, col) in zip(axes, taxonomies):
        sig = df[df["p_fdr"] < FDR_THRESHOLD].copy()

        if name == "Siletti only":
            # Map Siletti cluster names to broad classes
            siletti_class_map = {
                "MGE": "GABAergic (MGE)", "LLC": "GABAergic (MGE)",
                "Thex": "GABAergic (MGE)",
                "CGE": "GABAergic (CGE)", "Splat": "GABAergic (mixed/other)",
                "Midi": "GABAergic (mixed/other)",
                "Amex": "Glutamatergic", "Misc": "Glutamatergic",
                "DLIT": "Glutamatergic", "DLCT6b": "Glutamatergic",
                "L5ET": "Glutamatergic", "DG": "Glutamatergic",
                "CA13": "Glutamatergic", "CA4": "Glutamatergic",
                "MSN": "Glutamatergic", "EMSN": "Glutamatergic",
                "URL": "Glutamatergic", "CBI": "GABAergic (mixed/other)",
                "LRL": "GABAergic (mixed/other)", "Mmb": "GABAergic (mixed/other)",
                "Astro": "Non-neuronal", "Oligo": "Non-neuronal",
                "OPC": "Non-neuronal", "COP": "Non-neuronal",
                "Mgl": "Non-neuronal", "Epen": "Non-neuronal",
                "Fbl": "Non-neuronal", "Vend": "Non-neuronal",
                "Vsmc": "Non-neuronal", "Per": "Non-neuronal",
                "Bgl": "Non-neuronal", "Chrp": "Non-neuronal",
            }
            classes = []
            for t in sig[col]:
                prefix = t.split("_")[0]
                classes.append(siletti_class_map.get(prefix, "Other"))
            sig["class"] = classes
        elif name == "Combined":
            # Handle combined names (Siletti_ prefix or SEA-AD names)
            classes = []
            siletti_class_map = {
                "MGE": "GABAergic (MGE)", "LLC": "GABAergic (MGE)",
                "Thex": "GABAergic (MGE)",
                "CGE": "GABAergic (CGE)", "Splat": "GABAergic (mixed/other)",
                "Midi": "GABAergic (mixed/other)",
                "Amex": "Glutamatergic", "Misc": "Glutamatergic",
                "DLIT": "Glutamatergic", "DLCT6b": "Glutamatergic",
                "L5ET": "Glutamatergic", "DG": "Glutamatergic",
                "CA13": "Glutamatergic", "CA4": "Glutamatergic",
                "MSN": "Glutamatergic", "EMSN": "Glutamatergic",
                "URL": "Glutamatergic", "CBI": "GABAergic (mixed/other)",
                "LRL": "GABAergic (mixed/other)", "Mmb": "GABAergic (mixed/other)",
                "Astro": "Non-neuronal", "Oligo": "Non-neuronal",
                "OPC": "Non-neuronal", "COP": "Non-neuronal",
                "Mgl": "Non-neuronal", "Epen": "Non-neuronal",
                "Fbl": "Non-neuronal", "Vend": "Non-neuronal",
                "Vsmc": "Non-neuronal", "Per": "Non-neuronal",
                "Bgl": "Non-neuronal", "Chrp": "Non-neuronal",
            }
            from scz_celltype_enrichment.utils import classify_supertypes_by_class
            for t in sig[col]:
                if t.startswith("Siletti_"):
                    prefix = t.replace("Siletti_", "").split("_")[0]
                    classes.append(siletti_class_map.get(prefix, "Other"))
                else:
                    cls = classify_supertypes_by_class([t])
                    classes.append(cls.get(t, "Other"))
            sig["class"] = classes
        else:
            from scz_celltype_enrichment.utils import classify_supertypes_by_class
            cls_map = classify_supertypes_by_class(sig[col].tolist())
            sig["class"] = [cls_map.get(t, "Other") for t in sig[col]]

        class_counts = sig["class"].value_counts().sort_values(ascending=True)
        colors = {
            "GABAergic": "#e74c3c", "GABAergic (MGE)": "#e74c3c",
            "GABAergic (CGE)": "#e67e22", "GABAergic (mixed/other)": "#f39c12",
            "Glutamatergic": "#3498db", "Non-neuronal": "#2ecc71", "Other": "#95a5a6",
        }
        bar_colors = [colors.get(c, "#95a5a6") for c in class_counts.index]

        class_counts.plot(kind="barh", ax=ax, color=bar_colors)
        ax.set_title(f"{name}\n({len(sig)} FDR-significant)", fontsize=TITLE_FONTSIZE)
        ax.set_xlabel("Number of significant types", fontsize=LABEL_FONTSIZE)
        ax.tick_params(labelsize=TICK_FONTSIZE)

        # Add count labels
        for i, (count, label) in enumerate(zip(class_counts.values, class_counts.index)):
            ax.text(count + 0.3, i, str(count), va="center", fontsize=TICK_FONTSIZE)

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def plot_top_enrichment_comparison(seaad_enrich, siletti_enrich, output_path):
    """Side-by-side barplots of top enriched types from each taxonomy."""
    fig, axes = plt.subplots(1, 2, figsize=(20, 10))

    n_show = 25

    # SEA-AD top types
    ax = axes[0]
    top_seaad = seaad_enrich.nsmallest(n_show, "p_value")
    logp = -np.log10(np.clip(top_seaad["p_value"].values, 1e-50, 1))
    colors = ["tab:red" if fdr < FDR_THRESHOLD else "gray"
              for fdr in top_seaad["p_fdr"]]
    ax.barh(range(n_show), logp[::-1], color=colors[::-1])
    ax.set_yticks(range(n_show))
    ax.set_yticklabels(top_seaad["supertype"].values[::-1], fontsize=TICK_FONTSIZE)
    ax.set_xlabel("-log10(p-value)", fontsize=LABEL_FONTSIZE)
    ax.set_title(f"Top {n_show} SEA-AD enrichments", fontsize=TITLE_FONTSIZE)
    ax.tick_params(labelsize=TICK_FONTSIZE)

    # Siletti top types
    ax = axes[1]
    top_siletti = siletti_enrich.nsmallest(n_show, "p_value")
    logp = -np.log10(np.clip(top_siletti["p_value"].values, 1e-50, 1))
    colors = ["tab:red" if fdr < FDR_THRESHOLD else "gray"
              for fdr in top_siletti["p_fdr"]]
    ax.barh(range(n_show), logp[::-1], color=colors[::-1])
    ax.set_yticks(range(n_show))
    ax.set_yticklabels(top_siletti["supertype"].values[::-1], fontsize=TICK_FONTSIZE)
    ax.set_xlabel("-log10(p-value)", fontsize=LABEL_FONTSIZE)
    ax.set_title(f"Top {n_show} Siletti enrichments", fontsize=TITLE_FONTSIZE)
    ax.tick_params(labelsize=TICK_FONTSIZE)

    plt.tight_layout()
    fig.savefig(str(output_path), dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def main():
    ensure_dirs(TABLES_DIR, FIGURES_DIR)

    print("=" * 60)
    print("Step 08: Three-Taxonomy Enrichment Comparison")
    print("=" * 60)
    t0 = time.time()

    # --- Load existing enrichments ---
    print("\n1. Loading SEA-AD enrichment...")
    seaad_enrich = load_seaad_enrichment()

    print("\n2. Loading combined enrichment...")
    combined_enrich = load_combined_enrichment()

    # --- Compute Siletti-only enrichment ---
    siletti_output = TABLES_DIR / "siletti_only_enrichment.csv"
    if siletti_output.exists():
        print(f"\n3. Loading cached Siletti-only enrichment from {siletti_output}")
        siletti_enrich = pd.read_csv(str(siletti_output))
        siletti_enrich["source"] = "Siletti"
        print(f"  Siletti enrichment: {len(siletti_enrich)} types, "
              f"{(siletti_enrich['p_fdr'] < FDR_THRESHOLD).sum()} FDR-significant")
    else:
        print("\n3. Computing Siletti-only enrichment...")
        siletti_enrich = compute_siletti_only_enrichment()
        siletti_enrich.to_csv(str(siletti_output), index=False)
        print(f"  Saved: {siletti_output}")

    # --- Summary ---
    print("\n" + "=" * 60)
    print("Summary of FDR-significant types across taxonomies:")
    print(f"  SEA-AD only:  {(seaad_enrich['p_fdr'] < FDR_THRESHOLD).sum()}/{len(seaad_enrich)} types")
    print(f"  Siletti only: {(siletti_enrich['p_fdr'] < FDR_THRESHOLD).sum()}/{len(siletti_enrich)} types")
    print(f"  Combined:     {(combined_enrich['p_fdr'] < FDR_THRESHOLD).sum()}/{len(combined_enrich)} types")

    # --- Build comparison table ---
    print("\n4. Building cross-taxonomy comparison table...")
    mn_matches = load_metaneighbor_matches()
    comparison = build_comparison_table(seaad_enrich, siletti_enrich, combined_enrich, mn_matches)
    comparison_output = TABLES_DIR / "three_taxonomy_enrichment_comparison.csv"
    comparison.to_csv(str(comparison_output), index=False)
    print(f"  Saved: {comparison_output}")

    # --- Figures ---
    print("\n5. Generating comparison figures...")

    plot_enrichment_scatter(
        comparison,
        FIGURES_DIR / "taxonomy_enrichment_scatter.png"
    )

    plot_significant_types_by_class(
        seaad_enrich, siletti_enrich, combined_enrich,
        FIGURES_DIR / "taxonomy_significant_types_barplot.png"
    )

    plot_top_enrichment_comparison(
        seaad_enrich, siletti_enrich,
        FIGURES_DIR / "three_taxonomy_top_enrichments.png"
    )

    print(f"\nStep 08 complete. Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
