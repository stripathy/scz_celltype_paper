#!/usr/bin/env python
"""
sst_depth_volcano.py — Volcano plot of gene-depth correlations in SST supertypes.

For each gene, computes Spearman correlation between its expression across 18 SST
supertypes and the continuous cortical depth of those supertypes (from spatial
transcriptomics: MERFISH + Xenium weighted average).

Genes on the left (negative rho) are enriched in upper layers (near pial surface).
Genes on the right (positive rho) are enriched in deeper layers (near white matter).

Points are colored by risk direction from S-PrediXcan (GTEx Brain_Cortex):
  Red = genetically predicted UPregulation increases SCZ risk (z > 0)
  Blue = genetically predicted UPregulation is protective (z < 0)
For genes without S-PrediXcan models, falls back to lead-SNP BETA direction.

Two versions are generated:
  1. Using SEA-AD specificity scores (row-normalized: fraction of total expression)
  2. Using SEA-AD mean expression (log-normalized, not row-normalized)

Output:
  results/tables/sst_gene_depth_correlations_specificity.csv
  results/tables/sst_gene_depth_correlations_expression.csv
  results/figures/sst_depth_volcano_specificity.png
  results/figures/sst_depth_volcano_expression.png
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
import numpy as np

from scz_celltype_enrichment.config import (
    PROJECT_ROOT, SPECIFICITY_CSV, MEAN_EXPRESSION_CSV,
    SCZ_GWAS_GENE_SET_CSV, TABLES_DIR, FIGURES_DIR,
    MAGMA_GENES_OUT, GENE_LOC_FILE, PGC3_SNP_SUMSTATS,
    GENE_DIRECTION_CSV, INTERMEDIATES_DIR,
)
from scz_celltype_enrichment.enrichment.gwas import (
    load_scz_gwas_gene_set,
    load_entrez_to_symbol,
    compute_gene_effect_direction,
)
from scz_celltype_enrichment.integration.gene_ephys import compute_supertype_correlations
from scz_celltype_enrichment.integration.spatial import load_depth_data, get_sst_supertypes
from scz_celltype_enrichment.plotting import plot_depth_volcano
from scz_celltype_enrichment.utils import ensure_dirs


# Key genes to always label on the volcano plot
HIGHLIGHT_GENES = [
    "HCN1", "CACNA1I", "TPH2", "BCL11B", "ZNF536",
    "CALB1", "SST", "ERBB4", "TCF4", "RERE",
]

# S-PrediXcan results
SPREDIXCAN_DIR = PROJECT_ROOT / "results" / "spredixcan"
SPREDIXCAN_TISSUE = "Brain_Cortex"  # primary tissue for direction annotation


def load_spredixcan_results(tissue=SPREDIXCAN_TISSUE):
    """Load S-PrediXcan results for a brain tissue."""
    path = SPREDIXCAN_DIR / f"SCZ_PGC3_{tissue}.csv"
    if not path.exists():
        print(f"  WARNING: S-PrediXcan results not found: {path}")
        return None
    df = pd.read_csv(str(path))
    print(f"  Loaded S-PrediXcan ({tissue}): {len(df)} genes")
    n_sig = (df["pvalue_fdr"] < 0.05).sum() if "pvalue_fdr" in df.columns else 0
    print(f"    FDR < 0.05: {n_sig}")
    return df


def annotate_depth_correlations(corr_df, scz_genes, gwas_matched,
                                 gene_direction=None, spredixcan=None):
    """
    Add SCZ GWAS, MAGMA, effect direction, and S-PrediXcan annotations.

    Parameters
    ----------
    corr_df : pd.DataFrame
        From compute_supertype_correlations: gene, spearman_rho, pval, pval_fdr.
    scz_genes : set
        SCZ GWAS gene set.
    gwas_matched : pd.DataFrame
        MAGMA gene-level results with SYMBOL, ZSTAT columns.
    gene_direction : pd.DataFrame, optional
        Gene-level effect directions with gene, lead_snp_beta columns.
    spredixcan : pd.DataFrame, optional
        S-PrediXcan results with gene_name, zscore, pvalue, pvalue_fdr.

    Returns
    -------
    pd.DataFrame
        Annotated with is_scz_gwas, magma_zstat, spredixcan_z, direction_source,
        signed_magma_z (composite direction metric).
    """
    df = corr_df.copy()

    # SCZ GWAS gene set
    df["is_scz_gwas"] = df["gene"].isin(scz_genes)

    # MAGMA gene-level Z-stat (unsigned)
    zstat_map = gwas_matched.set_index("SYMBOL")["ZSTAT"].to_dict()
    df["magma_zstat"] = df["gene"].map(zstat_map)

    # S-PrediXcan z-score (SIGNED: positive = upregulation is risk)
    if spredixcan is not None:
        spx_map = spredixcan.set_index("gene_name")["zscore"].to_dict()
        spx_fdr_map = spredixcan.set_index("gene_name")["pvalue_fdr"].to_dict() if "pvalue_fdr" in spredixcan.columns else {}
        df["spredixcan_z"] = df["gene"].map(spx_map)
        df["spredixcan_fdr"] = df["gene"].map(spx_fdr_map)
        n_spx = df["spredixcan_z"].notna().sum()
        print(f"    S-PrediXcan matched: {n_spx}/{len(df)}")

    # Lead-SNP BETA direction (fallback)
    if gene_direction is not None:
        beta_map = gene_direction.set_index("gene")["lead_snp_beta"].to_dict()
        df["lead_snp_beta"] = df["gene"].map(beta_map)

    # Build composite signed MAGMA Z:
    # Priority 1: S-PrediXcan z-score (already signed, well-calibrated)
    # Priority 2: sign(lead_snp_beta) × |MAGMA_ZSTAT| (heuristic)
    df["signed_magma_z"] = np.nan
    df["direction_source"] = "none"

    # First: fill with lead-SNP approach where available
    if gene_direction is not None:
        has_both = df["lead_snp_beta"].notna() & df["magma_zstat"].notna()
        df.loc[has_both, "signed_magma_z"] = (
            np.sign(df.loc[has_both, "lead_snp_beta"]) * df.loc[has_both, "magma_zstat"].abs()
        )
        df.loc[has_both, "direction_source"] = "lead_snp"

    # Then: override with S-PrediXcan where available (better signal)
    if spredixcan is not None:
        has_spx = df["spredixcan_z"].notna()
        # Use S-PrediXcan z directly as the signed value
        # (it's already a z-score, comparable to MAGMA Z in scale for FDR-sig genes)
        df.loc[has_spx, "signed_magma_z"] = df.loc[has_spx, "spredixcan_z"]
        df.loc[has_spx, "direction_source"] = "spredixcan"

    # Summary
    n_scz = df["is_scz_gwas"].sum()
    n_magma = (df["magma_zstat"].abs() > 3).sum()
    print(f"    SCZ GWAS genes: {n_scz}/{len(df)}")
    print(f"    MAGMA |Z|>3: {n_magma}/{len(df)}")

    n_spx_dir = (df["direction_source"] == "spredixcan").sum()
    n_snp_dir = (df["direction_source"] == "lead_snp").sum()
    n_none = (df["direction_source"] == "none").sum()
    print(f"    Direction source: S-PrediXcan={n_spx_dir}, lead_SNP={n_snp_dir}, none={n_none}")

    if df["signed_magma_z"].notna().any():
        has_dir = df["signed_magma_z"].notna() & (df["magma_zstat"].abs() > 3)
        n_risk = (df.loc[has_dir, "signed_magma_z"] > 0).sum()
        n_prot = (df.loc[has_dir, "signed_magma_z"] < 0).sum()
        print(f"    High-signal (|MAGMA Z|>3) risk: {n_risk}, protective: {n_prot}")

    return df


def main():
    ensure_dirs(TABLES_DIR, FIGURES_DIR, INTERMEDIATES_DIR)

    print("=" * 60)
    print("SST Depth Volcano: Gene-Depth Correlations")
    print("=" * 60)

    t0 = time.time()

    # ----------------------------------------------------------------
    # 1. Load depth data
    # ----------------------------------------------------------------
    print("\n1. Loading SST depth data...")
    depth_csv = TABLES_DIR / "sst_supertype_layer_info.csv"
    if not depth_csv.exists():
        from scz_celltype_enrichment.config import DEPTH_SUPERTYPE_CSV
        depth_df = load_depth_data(DEPTH_SUPERTYPE_CSV)
        sst_depth = get_sst_supertypes(depth_df)
    else:
        sst_depth = pd.read_csv(str(depth_csv))

    depth_by_supertype = sst_depth.set_index("supertype")["avg_depth"]
    sst_types = depth_by_supertype.index.tolist()
    print(f"  {len(sst_types)} SST supertypes, depth range: "
          f"[{depth_by_supertype.min():.3f}, {depth_by_supertype.max():.3f}]")

    # ----------------------------------------------------------------
    # 2. Load GWAS annotation data + effect direction + S-PrediXcan
    # ----------------------------------------------------------------
    print("\n2. Loading GWAS data...")
    scz_genes = load_scz_gwas_gene_set(SCZ_GWAS_GENE_SET_CSV)
    gwas_matched = pd.read_csv(str(PROJECT_ROOT / "results" / "intermediates" / "gwas_matched.csv"))
    print(f"  SCZ GWAS gene set: {len(scz_genes)} genes")
    print(f"  MAGMA results: {len(gwas_matched)} genes")

    # Gene-level effect direction from SNP-level data (fallback)
    print("\n  Loading gene-level effect direction from PGC3 SNP-level data...")
    gene_direction_path = Path(GENE_DIRECTION_CSV)
    if gene_direction_path.exists():
        gene_direction = pd.read_csv(str(gene_direction_path))
        print(f"  Loaded cached direction for {len(gene_direction)} genes")
    elif Path(PGC3_SNP_SUMSTATS).exists():
        entrez_map = load_entrez_to_symbol(GENE_LOC_FILE)
        gene_direction = compute_gene_effect_direction(
            PGC3_SNP_SUMSTATS, MAGMA_GENES_OUT, entrez_map,
            output_path=gene_direction_path,
        )
    else:
        print("  WARNING: PGC3 SNP-level summary stats not found.")
        gene_direction = None

    # S-PrediXcan results (primary direction source)
    print("\n  Loading S-PrediXcan results...")
    spredixcan = load_spredixcan_results(SPREDIXCAN_TISSUE)

    # Print key genes direction from both sources
    if spredixcan is not None:
        print("\n  Effect direction for key genes (S-PrediXcan):")
        for g in HIGHLIGHT_GENES:
            spx_row = spredixcan[spredixcan["gene_name"] == g]
            snp_row = gene_direction[gene_direction["gene"] == g] if gene_direction is not None else pd.DataFrame()
            if len(spx_row) > 0:
                r = spx_row.iloc[0]
                direction = "UP->risk" if r["zscore"] > 0 else "DOWN->risk"
                fdr_str = f"FDR={r.get('pvalue_fdr', float('nan')):.1e}"
                print(f"    {g:12s}  SPX z={r['zscore']:+7.3f}  {fdr_str}  ({direction})")
            elif len(snp_row) > 0:
                r = snp_row.iloc[0]
                direction = "RISK" if r["lead_snp_beta"] > 0 else "PROTECTIVE"
                print(f"    {g:12s}  [no SPX model] lead_beta={r['lead_snp_beta']:+.4f} ({direction})")
            else:
                print(f"    {g:12s}  [no data]")

    # ----------------------------------------------------------------
    # 3. Process SPECIFICITY version
    # ----------------------------------------------------------------
    print("\n3. Computing depth correlations (SPECIFICITY)...")
    specificity = pd.read_csv(str(SPECIFICITY_CSV), index_col=0)
    sst_spec_types = [t for t in sst_types if t in specificity.columns]
    sst_spec = specificity[sst_spec_types]
    print(f"  Genes: {sst_spec.shape[0]:,}, SST types: {len(sst_spec_types)}")

    spec_corr = compute_supertype_correlations(
        sst_spec, depth_by_supertype.loc[sst_spec_types],
        feature_name="depth (specificity)"
    )

    spec_corr = annotate_depth_correlations(
        spec_corr, scz_genes, gwas_matched, gene_direction, spredixcan
    )
    spec_corr.to_csv(TABLES_DIR / "sst_gene_depth_correlations_specificity.csv", index=False)
    print(f"  Saved: sst_gene_depth_correlations_specificity.csv ({len(spec_corr):,} genes)")

    # ----------------------------------------------------------------
    # 4. Process MEAN EXPRESSION version
    # ----------------------------------------------------------------
    mean_expr_path = Path(MEAN_EXPRESSION_CSV)
    if mean_expr_path.exists():
        print("\n4. Computing depth correlations (MEAN EXPRESSION)...")
        mean_expr = pd.read_csv(str(mean_expr_path), index_col=0)
        sst_expr_types = [t for t in sst_types if t in mean_expr.columns]
        sst_expr = mean_expr[sst_expr_types]
        print(f"  Genes: {sst_expr.shape[0]:,}, SST types: {len(sst_expr_types)}")

        expr_corr = compute_supertype_correlations(
            sst_expr, depth_by_supertype.loc[sst_expr_types],
            feature_name="depth (expression)"
        )

        expr_corr = annotate_depth_correlations(
            expr_corr, scz_genes, gwas_matched, gene_direction, spredixcan
        )
        expr_corr.to_csv(TABLES_DIR / "sst_gene_depth_correlations_expression.csv", index=False)
        print(f"  Saved: sst_gene_depth_correlations_expression.csv ({len(expr_corr):,} genes)")
    else:
        print("\n4. SKIPPING mean expression version (not yet computed).")
        expr_corr = None

    # ----------------------------------------------------------------
    # 5. Generate volcano plots
    # ----------------------------------------------------------------
    print("\n5. Generating volcano plots...")

    plot_depth_volcano(
        spec_corr,
        FIGURES_DIR / "sst_depth_volcano_specificity.png",
        title="Gene-Depth Correlations (Specificity) Across SST Supertypes",
        n_label=12,
        highlight_genes=HIGHLIGHT_GENES,
    )

    if expr_corr is not None:
        plot_depth_volcano(
            expr_corr,
            FIGURES_DIR / "sst_depth_volcano_expression.png",
            title="Gene-Depth Correlations (Mean Expression) Across SST Supertypes",
            n_label=12,
            highlight_genes=HIGHLIGHT_GENES,
        )

    # ----------------------------------------------------------------
    # 6. Summary
    # ----------------------------------------------------------------
    if expr_corr is not None:
        print("\n6. Comparing specificity vs expression correlations...")
        merged = spec_corr[["gene", "spearman_rho", "pval_fdr"]].merge(
            expr_corr[["gene", "spearman_rho", "pval_fdr"]],
            on="gene", suffixes=("_spec", "_expr")
        )
        from scipy.stats import spearmanr
        rho_corr, _ = spearmanr(merged["spearman_rho_spec"], merged["spearman_rho_expr"])
        print(f"  Spearman correlation of rho values: {rho_corr:.3f}")

    # Direction summary for FDR-significant depth genes
    print("\n7. Directional summary for FDR-significant depth genes...")
    sig = spec_corr[spec_corr["pval_fdr"] < 0.05].copy()
    upper_sig = sig[sig["spearman_rho"] < 0]
    deep_sig = sig[sig["spearman_rho"] > 0]

    for label, subset in [("Upper-layer enriched", upper_sig), ("Deep-layer enriched", deep_sig)]:
        has_dir = subset["signed_magma_z"].notna()
        n_risk = (subset.loc[has_dir, "signed_magma_z"] > 0).sum()
        n_prot = (subset.loc[has_dir, "signed_magma_z"] < 0).sum()
        print(f"  {label} (FDR<0.05): {len(subset)} genes")
        print(f"    With direction data: {has_dir.sum()}")
        print(f"    Risk (signed Z > 0): {n_risk}")
        print(f"    Protective (signed Z < 0): {n_prot}")

        # By direction source
        n_spx = (subset["direction_source"] == "spredixcan").sum()
        n_snp = (subset["direction_source"] == "lead_snp").sum()
        print(f"    Source: S-PrediXcan={n_spx}, lead_SNP={n_snp}")

    print(f"\nComplete. Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
