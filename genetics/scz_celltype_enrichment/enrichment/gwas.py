"""
gwas.py — Load and process MAGMA gene-level association results.

MAGMA aggregates all SNPs within/near a gene (35kb upstream, 10kb downstream)
into a single gene-level Z-statistic, correcting for linkage disequilibrium.
This captures the full polygenic signal, not just genome-wide significant hits.

Input files:
  - MAGMA .genes.out: space-delimited, columns GENE CHR START STOP NSNPS NPARAM N ZSTAT P
    (GENE column is ENTREZ ID, 18,449 genes for PGC3 SCZ)
  - NCBI gene location file: tab-delimited, no header,
    columns ENTREZ CHR START STOP STRAND SYMBOL (19,175 genes, MHC excluded)
"""
import numpy as np
import pandas as pd

from ..utils import Timer


def load_magma_genes(genes_out_path):
    """
    Load MAGMA step2 .genes.out file with gene-level GWAS statistics.

    Parameters
    ----------
    genes_out_path : str or Path
        Path to .genes.out file.

    Returns
    -------
    pd.DataFrame
        Columns: GENE (int, ENTREZ), CHR, START, STOP, NSNPS, NPARAM, N, ZSTAT, P
    """
    with Timer("Loading MAGMA gene-level results"):
        df = pd.read_csv(str(genes_out_path), sep=r"\s+")

    # GENE column is ENTREZ ID
    df["GENE"] = df["GENE"].astype(int)
    print(f"    Loaded {len(df):,} genes")
    print(f"    ZSTAT range: [{df['ZSTAT'].min():.2f}, {df['ZSTAT'].max():.2f}]")
    print(f"    Top 5 by ZSTAT: {df.nlargest(5, 'ZSTAT')['GENE'].tolist()}")
    return df


def load_entrez_to_symbol(gene_loc_path):
    """
    Load NCBI gene location file and build ENTREZ → SYMBOL mapping.

    File format: tab-delimited, no header
    Columns: ENTREZ_ID  CHR  START  STOP  STRAND  SYMBOL

    Parameters
    ----------
    gene_loc_path : str or Path
        Path to NCBI37.3.gene.loc.extendedMHCexcluded.

    Returns
    -------
    dict
        Mapping from ENTREZ ID (int) to gene SYMBOL (str).
    """
    with Timer("Loading ENTREZ→SYMBOL mapping"):
        gene_loc = pd.read_csv(
            str(gene_loc_path),
            sep="\t",
            header=None,
            names=["ENTREZ", "CHR", "START", "STOP", "STRAND", "SYMBOL"],
        )
        mapping = dict(zip(gene_loc["ENTREZ"].astype(int), gene_loc["SYMBOL"]))

    print(f"    Loaded {len(mapping):,} gene mappings")
    return mapping


def load_scz_gwas_gene_set(gene_set_path):
    """
    Load the SCZ GWAS gene set (genes within ±100kb of significant loci).

    Parameters
    ----------
    gene_set_path : str or Path
        Path to scz_gwas_gene_set.csv.

    Returns
    -------
    set
        Set of gene symbols in the SCZ GWAS gene set.
    """
    df = pd.read_csv(str(gene_set_path))
    genes = set(df["gene"].dropna().unique())
    print(f"  Loaded SCZ GWAS gene set: {len(genes)} genes")
    return genes


def map_and_merge(magma_df, specificity_df, entrez_to_symbol):
    """
    Map ENTREZ IDs to gene symbols and intersect with specificity matrix.

    Parameters
    ----------
    magma_df : pd.DataFrame
        MAGMA gene-level results with GENE (ENTREZ) column.
    specificity_df : pd.DataFrame
        Specificity matrix with gene symbols as index.
    entrez_to_symbol : dict
        ENTREZ → SYMBOL mapping.

    Returns
    -------
    spec_matched : pd.DataFrame
        Specificity matrix filtered to overlapping genes, gene symbols as index.
    gwas_matched : pd.DataFrame
        GWAS results for overlapping genes, with SYMBOL column added.
    """
    with Timer("Mapping ENTREZ to SYMBOL and intersecting"):
        # Map ENTREZ to SYMBOL
        magma_df = magma_df.copy()
        magma_df["SYMBOL"] = magma_df["GENE"].map(entrez_to_symbol)

        # Drop genes without mapping
        n_before = len(magma_df)
        magma_df = magma_df.dropna(subset=["SYMBOL"])
        n_mapped = len(magma_df)
        print(f"    ENTREZ→SYMBOL mapping: {n_mapped:,}/{n_before:,} genes mapped "
              f"({n_mapped/n_before*100:.1f}%)")

        # Drop duplicate symbols (keep first occurrence)
        magma_df = magma_df.drop_duplicates(subset="SYMBOL", keep="first")

        # Intersect with specificity matrix
        overlap = sorted(set(magma_df["SYMBOL"]) & set(specificity_df.index))
        n_overlap = len(overlap)
        print(f"    Overlapping genes: {n_overlap:,}")

        spec_matched = specificity_df.loc[overlap]
        gwas_matched = magma_df.set_index("SYMBOL").loc[overlap].reset_index()

    return spec_matched, gwas_matched


def compute_gene_effect_direction(
    snp_sumstats_path,
    genes_out_path,
    entrez_to_symbol_map,
    output_path=None,
):
    """
    Compute per-gene effect direction from SNP-level GWAS summary statistics.

    For each gene, identifies the lead SNP (most significant) within the gene
    boundaries and extracts its signed BETA (log-OR). Positive BETA = risk
    allele increases SCZ liability; negative BETA = protective.

    Also computes a significance-weighted mean BETA across all gene SNPs as
    a more robust directional estimate.

    Parameters
    ----------
    snp_sumstats_path : str or Path
        Path to PGC3 SCZ summary statistics (gzipped TSV with ## comment lines).
        Must have columns: CHROM, POS, BETA, PVAL.
    genes_out_path : str or Path
        Path to MAGMA .genes.out file (for gene boundaries).
    entrez_to_symbol_map : dict
        ENTREZ ID → gene SYMBOL mapping.
    output_path : str or Path, optional
        If provided, save results to CSV.

    Returns
    -------
    pd.DataFrame
        Columns: gene, lead_snp_beta, lead_snp_pval, n_snps_in_gene,
                 weighted_mean_beta
    """
    with Timer("Loading SNP-level GWAS summary statistics"):
        snps = pd.read_csv(
            str(snp_sumstats_path),
            sep="\t", compression="gzip", comment="#",
            usecols=["CHROM", "POS", "BETA", "PVAL"],
            dtype={"CHROM": str, "POS": int, "BETA": float, "PVAL": float},
        )
    print(f"    {len(snps):,} SNPs, BETA range: [{snps['BETA'].min():.4f}, {snps['BETA'].max():.4f}]")

    # Load gene boundaries from MAGMA
    genes_out = pd.read_csv(str(genes_out_path), sep=r"\s+")
    genes_out["SYMBOL"] = genes_out["GENE"].map(entrez_to_symbol_map)
    genes_out = genes_out.dropna(subset=["SYMBOL"])
    genes_out["CHR"] = genes_out["CHR"].astype(str)

    # For each gene, find lead SNP and compute directional summary
    print(f"    Mapping SNPs to {len(genes_out)} genes...")
    results = []
    for chrom in sorted(snps["CHROM"].unique(), key=lambda x: int(x) if x.isdigit() else 99):
        chr_snps = snps[snps["CHROM"] == chrom].sort_values("POS")
        chr_genes = genes_out[genes_out["CHR"] == chrom]

        if len(chr_genes) == 0:
            continue

        snp_pos = chr_snps["POS"].values
        snp_beta = chr_snps["BETA"].values
        snp_pval = chr_snps["PVAL"].values

        for _, gene in chr_genes.iterrows():
            mask = (snp_pos >= gene["START"]) & (snp_pos <= gene["STOP"])
            if not mask.any():
                continue

            gene_pvals = snp_pval[mask]
            gene_betas = snp_beta[mask]

            # Lead SNP: most significant
            lead_idx = np.argmin(gene_pvals)
            lead_beta = gene_betas[lead_idx]
            lead_pval = gene_pvals[lead_idx]

            # Weighted mean beta (weighted by -log10(p))
            weights = -np.log10(np.clip(gene_pvals, 1e-300, 1))
            weighted_beta = np.average(gene_betas, weights=weights)

            results.append({
                "gene": gene["SYMBOL"],
                "lead_snp_beta": float(lead_beta),
                "lead_snp_pval": float(lead_pval),
                "n_snps_in_gene": int(mask.sum()),
                "weighted_mean_beta": float(weighted_beta),
            })

    result_df = pd.DataFrame(results)

    n_risk = (result_df["lead_snp_beta"] > 0).sum()
    n_prot = (result_df["lead_snp_beta"] < 0).sum()
    print(f"    Computed direction for {len(result_df)} genes: "
          f"{n_risk} risk (beta>0), {n_prot} protective (beta<0)")

    if output_path is not None:
        from pathlib import Path
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        result_df.to_csv(str(output_path), index=False)
        print(f"    Saved: {output_path}")

    return result_df


def build_covariates(gwas_matched):
    """
    Build the covariate matrix for MAGMA-style regression.

    Covariates: intercept, log(gene_size), log(nsnps), log(N)
    where gene_size = STOP - START.

    These control for known confounds: larger genes have more SNPs,
    genes with more SNPs have inflated test statistics, and sample
    size affects statistical power.

    Parameters
    ----------
    gwas_matched : pd.DataFrame
        GWAS results with START, STOP, NSNPS, N columns.

    Returns
    -------
    np.ndarray, shape (n_genes, 4)
        Covariate matrix: [intercept, log_gene_size, log_nsnps, log_N]
    """
    gene_size = (gwas_matched["STOP"] - gwas_matched["START"]).values.astype(float)
    gene_size = np.maximum(gene_size, 1)  # Avoid log(0)

    nsnps = gwas_matched["NSNPS"].values.astype(float)
    nsnps = np.maximum(nsnps, 1)

    sample_n = gwas_matched["N"].values.astype(float)

    intercept = np.ones(len(gwas_matched))
    log_gene_size = np.log(gene_size)
    log_nsnps = np.log(nsnps)
    log_n = np.log(sample_n)

    covariates = np.column_stack([intercept, log_gene_size, log_nsnps, log_n])

    print(f"  Built covariate matrix: {covariates.shape}")
    return covariates
