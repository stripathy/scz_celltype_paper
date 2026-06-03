#!/usr/bin/env python3
"""
Custom S-PrediXcan implementation for PGC3 SCZ × GTEx brain tissues.

For each gene with a prediction model, computes:
    z_gene = (w' · z) / sqrt(w' · Σ · w)

where w = eQTL prediction weights, z = GWAS z-scores, Σ = LD covariance.

Positive z_gene → genetically predicted upregulation associated with SCZ risk.
Negative z_gene → genetically predicted upregulation is protective.

Vectorized implementation: merges all weights with GWAS once, then groups by gene.
"""
import pandas as pd
import numpy as np
import sqlite3
import gzip
import time
import sys
from pathlib import Path
from scipy import stats
from collections import defaultdict

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# ============================================================================
# Configuration
# ============================================================================

GWAS_FILE = PROJECT_ROOT / "results" / "intermediates" / "pgc3_clean.tsv.gz"
MODEL_DIR = PROJECT_ROOT / "data" / "predixcan_models" / "eqtl" / "mashr"
OUTPUT_DIR = PROJECT_ROOT / "results" / "spredixcan"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

BRAIN_TISSUES = [
    "Brain_Cortex",
    "Brain_Frontal_Cortex_BA9",
]


# ============================================================================
# Core functions
# ============================================================================

def load_gwas():
    """Load PGC3 SCZ GWAS and compute z-scores."""
    print("Loading GWAS summary statistics...")
    gwas = pd.read_csv(
        str(GWAS_FILE), sep="\t", compression="gzip",
        usecols=["ID", "A1", "A2", "BETA", "SE", "PVAL"],
    )
    gwas["zscore"] = gwas["BETA"] / gwas["SE"]
    print(f"  {len(gwas):,} SNPs loaded")
    return gwas


def load_model_and_merge(tissue, gwas):
    """
    Load prediction model + covariance for a tissue.
    Merge weights with GWAS in one vectorized operation.
    Returns merged DataFrame and covariance dict.
    """
    db_path = MODEL_DIR / f"mashr_{tissue}.db"
    cov_path = MODEL_DIR / f"mashr_{tissue}.txt.gz"

    # Load weights + gene names
    conn = sqlite3.connect(str(db_path))
    weights = pd.read_sql(
        "SELECT gene, rsid, varID, ref_allele, eff_allele, weight FROM weights",
        conn,
    )
    extra = pd.read_sql("SELECT gene, genename FROM extra", conn)
    conn.close()
    gene_names = dict(zip(extra["gene"], extra["genename"]))

    n_model_genes = weights["gene"].nunique()
    n_model_snps = weights["rsid"].nunique()

    # --- Vectorized merge: weights × GWAS ---
    merged = weights.merge(gwas, left_on="rsid", right_on="ID", how="inner")
    print(f"  Merged: {len(merged):,} SNP-gene pairs "
          f"({merged['gene'].nunique():,} genes, "
          f"{merged['rsid'].nunique():,}/{n_model_snps} model SNPs matched)")

    # Allele alignment: flip z-score if effect alleles are swapped
    # Case 1: model eff_allele == GWAS A1 → aligned (no flip)
    # Case 2: model eff_allele == GWAS A2 → flip z
    # Case 3: neither → allele mismatch, drop
    aligned = merged["eff_allele"] == merged["A1"]
    flipped = merged["eff_allele"] == merged["A2"]
    keep = aligned | flipped

    merged = merged[keep].copy()
    merged["z_aligned"] = merged["zscore"].copy()
    flip_mask = merged["eff_allele"] == merged["A2"]
    merged.loc[flip_mask, "z_aligned"] = -merged.loc[flip_mask, "zscore"]

    print(f"  After allele alignment: {len(merged):,} SNP-gene pairs "
          f"({flip_mask.sum():,} flipped)")

    # Load covariance, indexed by gene
    print(f"  Loading covariance...")
    t0 = time.time()
    cov_by_gene = defaultdict(list)
    with gzip.open(str(cov_path), "rt") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4 or parts[0] == "GENE":
                continue
            cov_by_gene[parts[0]].append((parts[1], parts[2], float(parts[3])))
    print(f"    {len(cov_by_gene):,} genes in {time.time()-t0:.1f}s")

    return merged, gene_names, cov_by_gene, n_model_genes


def compute_spredixcan(merged, gene_names, cov_by_gene):
    """Compute S-PrediXcan z-score for all genes using pre-merged data.

    Note: covariance files use varID (chr_pos_ref_alt_b38) not rsID,
    so we index on varID for covariance matching.
    """
    results = []
    grouped = merged.groupby("gene")
    n_genes = len(grouped)
    n_no_cov = 0
    n_neg_var = 0

    for i, (gene, grp) in enumerate(grouped):
        if (i + 1) % 5000 == 0:
            print(f"    {i+1}/{n_genes} genes...")

        # Use varID for covariance matching
        varids = grp["varID"].values
        w_arr = grp["weight"].values
        z_arr = grp["z_aligned"].values
        n_snps = len(varids)

        # Numerator: w' · z
        numerator = np.dot(w_arr, z_arr)

        # Denominator: sqrt(w' · Σ · w)
        # Covariance file uses varID format
        cov_entries = cov_by_gene.get(gene, [])

        if len(cov_entries) == 0:
            # No covariance → assume independent (diagonal = 1)
            var_g = np.sum(w_arr ** 2)
            n_no_cov += 1
        else:
            snp_idx = {s: j for j, s in enumerate(varids)}
            sigma = np.zeros((n_snps, n_snps))

            for vid1, vid2, val in cov_entries:
                ii = snp_idx.get(vid1)
                jj = snp_idx.get(vid2)
                if ii is not None and jj is not None:
                    sigma[ii, jj] = val
                    sigma[jj, ii] = val

            var_g = w_arr @ sigma @ w_arr

        if var_g <= 0:
            n_neg_var += 1
            continue

        z_gene = numerator / np.sqrt(var_g)
        p_gene = 2 * stats.norm.sf(abs(z_gene))

        results.append({
            "gene_id": gene,
            "gene_name": gene_names.get(gene, gene),
            "zscore": z_gene,
            "pvalue": p_gene,
            "n_snps_used": n_snps,
        })

    print(f"    Genes without covariance data: {n_no_cov}")
    print(f"    Genes with var_g <= 0 (skipped): {n_neg_var}")

    if len(results) == 0:
        return pd.DataFrame(columns=["gene_id", "gene_name", "zscore", "pvalue", "n_snps_used"])
    return pd.DataFrame(results).sort_values("pvalue")


# ============================================================================
# Main
# ============================================================================

def main():
    t_start = time.time()
    print("=" * 60)
    print("S-PrediXcan: PGC3 SCZ x GTEx Brain Tissues")
    print("=" * 60)

    gwas = load_gwas()

    for tissue in BRAIN_TISSUES:
        print(f"\n{'=' * 50}")
        print(f"Tissue: {tissue}")
        print(f"{'=' * 50}")

        t0 = time.time()
        merged, gene_names, cov_by_gene, n_model = load_model_and_merge(tissue, gwas)

        print(f"\n  Computing gene-level z-scores...")
        results = compute_spredixcan(merged, gene_names, cov_by_gene)
        elapsed = time.time() - t0
        print(f"  Done: {len(results):,} genes with results in {elapsed:.1f}s")

        # FDR correction
        n_bonf = (results["pvalue"] < 0.05 / len(results)).sum()
        n_fdr = 0
        if len(results) > 0:
            from statsmodels.stats.multitest import multipletests
            _, pval_fdr, _, _ = multipletests(results["pvalue"], method="fdr_bh")
            results["pvalue_fdr"] = pval_fdr
            n_fdr = (pval_fdr < 0.05).sum()

        print(f"\n  Bonferroni significant (p < {0.05/len(results):.1e}): {n_bonf}")
        print(f"  FDR < 0.05: {n_fdr}")

        # Directionality summary
        sig = results[results["pvalue_fdr"] < 0.05]
        if len(sig) > 0:
            n_up = (sig["zscore"] > 0).sum()
            n_down = (sig["zscore"] < 0).sum()
            print(f"  Among FDR-sig: {n_up} upregulated-in-risk, {n_down} downregulated-in-risk")

        # Save
        outpath = OUTPUT_DIR / f"SCZ_PGC3_{tissue}.csv"
        results.to_csv(str(outpath), index=False)
        print(f"  Saved: {outpath}")

        # Top hits
        print(f"\n  Top 20 hits:")
        top = results.head(20)
        for _, row in top.iterrows():
            direction = "UP" if row["zscore"] > 0 else "DOWN"
            fdr = row.get("pvalue_fdr", float("nan"))
            print(f"    {row['gene_name']:15s}  z={row['zscore']:+7.3f}  "
                  f"p={row['pvalue']:.2e}  FDR={fdr:.2e}  "
                  f"({direction}, {row['n_snps_used']:.0f} SNPs)")

    print(f"\n{'=' * 60}")
    print(f"Total time: {time.time()-t_start:.1f}s")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
