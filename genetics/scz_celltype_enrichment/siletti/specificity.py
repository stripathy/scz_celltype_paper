"""
specificity.py — Compute cell-type specificity from Siletti cluster-level data.

Handles two specificity sources:
  1. Conti/Duncan precomputed (ENTREZ-indexed, numeric cluster IDs)
  2. Reprocessed from Allen Institute data (Ensembl → ENTREZ mapping)
"""
import numpy as np
import pandas as pd
import loompy

from ..utils import Timer


def build_ensembl_to_entrez_mapping(loom_path, gene_loc_path):
    """
    Build Ensembl → ENTREZ gene ID mapping by chaining:
      Ensembl → Symbol (from loom Gene attribute)
      Symbol → ENTREZ (from NCBI gene loc file)

    Parameters
    ----------
    loom_path : str or Path
        Path to Siletti loom file (for Accession → Gene mapping).
    gene_loc_path : str or Path
        Path to NCBI37.3.gene.loc.extendedMHCexcluded.

    Returns
    -------
    dict
        Mapping from Ensembl base ID (no version) to ENTREZ ID (int).
    """
    with Timer("Building Ensembl→ENTREZ mapping"):
        # Ensembl → Symbol from loom
        with loompy.connect(str(loom_path), "r") as ds:
            accessions = ds.ra["Accession"]
            gene_names = ds.ra["Gene"]

        ens_base = [a.split(".")[0] for a in accessions]
        ens_to_symbol = dict(zip(ens_base, gene_names))

        # Symbol → ENTREZ from NCBI gene loc
        gene_loc = pd.read_csv(
            str(gene_loc_path), sep="\t", header=None,
            names=["ENTREZ", "CHR", "START", "STOP", "STRAND", "SYMBOL"],
        )
        symbol_to_entrez = dict(zip(gene_loc["SYMBOL"], gene_loc["ENTREZ"].astype(int)))

        # Chain: Ensembl → Symbol → ENTREZ
        mapping = {}
        for ens, sym in ens_to_symbol.items():
            if sym in symbol_to_entrez:
                mapping[ens] = symbol_to_entrez[sym]

    print(f"    Ensembl→Symbol: {len(ens_to_symbol):,}")
    print(f"    Symbol→ENTREZ: {len(symbol_to_entrez):,}")
    print(f"    Chained Ensembl→ENTREZ: {len(mapping):,}")
    return mapping


def compute_siletti_specificity_entrez(cluster_stats_path, ens_to_entrez):
    """
    Compute specificity matrix for all 461 Siletti clusters using ENTREZ IDs.

    Steps:
      1. Load cluster-level mean expression from NPZ
      2. Map Ensembl gene IDs → ENTREZ (using pre-built mapping)
      3. Deduplicate genes (keep highest-expressed per ENTREZ ID)
      4. Column-normalize (each cluster sums to 1000)
      5. Row-normalize (each gene sums to 1.0)

    Parameters
    ----------
    cluster_stats_path : str or Path
        Path to siletti_cluster_level_stats.npz.
    ens_to_entrez : dict
        Ensembl → ENTREZ mapping from build_ensembl_to_entrez_mapping().

    Returns
    -------
    pd.DataFrame
        Specificity matrix with ENTREZ IDs as index, cluster names as columns.
        Shape: (n_genes, 461). Rows sum to 1.0.
    """
    with Timer("Computing reprocessed Siletti specificity"):
        d = np.load(str(cluster_stats_path), allow_pickle=True)
        gene_ids = d["gene_ids"]
        mean_raw = d["mean_raw"]  # (461, n_genes)
        cluster_names = d["cluster_names"]

        print(f"    Raw: {mean_raw.shape[0]} clusters × {mean_raw.shape[1]} genes")

        # Map to ENTREZ
        entrez_ids = np.array([ens_to_entrez.get(g, -1) for g in gene_ids])
        has_entrez = entrez_ids > 0
        mean_mapped = mean_raw[:, has_entrez]
        entrez_mapped = entrez_ids[has_entrez]
        print(f"    ENTREZ-mapped: {mean_mapped.shape[1]} genes")

        # Deduplicate: keep gene with highest total expression per ENTREZ
        entrez_to_best = {}
        for i, eid in enumerate(entrez_mapped):
            total = mean_mapped[:, i].sum()
            if eid not in entrez_to_best or total > entrez_to_best[eid][1]:
                entrez_to_best[eid] = (i, total)

        keep_idx = [v[0] for v in entrez_to_best.values()]
        keep_entrez = list(entrez_to_best.keys())
        mean_dedup = mean_mapped[:, keep_idx]
        print(f"    After dedup: {mean_dedup.shape[1]} unique ENTREZ genes")

        # Column normalize (each cluster sums to 1000)
        col_sums = mean_dedup.sum(axis=1, keepdims=True)
        col_sums[col_sums == 0] = 1
        col_norm = mean_dedup * (1000.0 / col_sums)

        # Row normalize (each gene sums to 1.0)
        row_sums = col_norm.sum(axis=0, keepdims=True)
        nonzero = row_sums.flatten() > 0
        specificity = col_norm[:, nonzero] / row_sums[:, nonzero]
        keep_entrez_filt = [keep_entrez[i] for i in range(len(keep_entrez)) if nonzero[i]]

        # Build DataFrame (genes × clusters)
        spec_df = pd.DataFrame(
            specificity.T,
            index=keep_entrez_filt,
            columns=cluster_names,
        )
        spec_df.index.name = "ENTREZ"

    print(f"    Final: {spec_df.shape[0]:,} genes × {spec_df.shape[1]} clusters")

    # Validate
    gene_sums = spec_df.sum(axis=1)
    assert np.allclose(gene_sums, 1.0, atol=1e-6), "Row normalization failed"
    print("    Validated: all rows sum to 1.0")

    return spec_df


def load_conti_specificity(conti_path, num_to_name):
    """
    Load the precomputed Conti/Duncan specificity matrix and rename columns.

    Parameters
    ----------
    conti_path : str or Path
        Path to conti_specificity_matrix.txt.
    num_to_name : dict
        Mapping from integer cluster index → cluster name
        (from clusters.build_cluster_number_map).

    Returns
    -------
    pd.DataFrame
        Specificity matrix with ENTREZ IDs as index, named clusters as columns.
    """
    with Timer("Loading Conti specificity"):
        spec = pd.read_csv(str(conti_path), sep="\t", index_col=0)

    # Rename: Cluster0 → name, Cluster1 → name, ...
    col_map = {f"Cluster{i}": num_to_name[i] for i in range(len(num_to_name))}
    spec.columns = [col_map.get(c, c) for c in spec.columns]

    print(f"    Shape: {spec.shape[0]:,} genes × {spec.shape[1]} clusters (ENTREZ index)")
    return spec


def build_combined_mean_expression(seaad_mean_expr_path, cluster_stats_path, loom_path):
    """
    Build a combined mean expression matrix (gene symbols) for MetaNeighbor.

    Combines SEA-AD log-normalized mean expression with Siletti cluster means
    (normalized from raw counts: target_sum=1e4 + log1p).

    Parameters
    ----------
    seaad_mean_expr_path : str or Path
        Path to seaad_supertype_mean_expression.csv.
    cluster_stats_path : str or Path
        Path to siletti_cluster_level_stats.npz.
    loom_path : str or Path
        Path to loom file (for Ensembl → Symbol mapping).

    Returns
    -------
    pd.DataFrame
        Combined mean expression (genes × cell types), gene symbols as index.
        SEA-AD columns are plain names, Siletti columns have "Siletti_" prefix.
    """
    with Timer("Building combined mean expression"):
        # Ensembl → Symbol from loom
        with loompy.connect(str(loom_path), "r") as ds:
            accessions = ds.ra["Accession"]
            gene_names = ds.ra["Gene"]

        ens_to_symbol = dict(zip(
            [a.split(".")[0] for a in accessions],
            gene_names,
        ))

        # Load Siletti raw means
        d = np.load(str(cluster_stats_path), allow_pickle=True)
        gene_ids = d["gene_ids"]
        mean_raw = d["mean_raw"]
        cluster_names = d["cluster_names"]

        # Map to symbols
        symbols = np.array([ens_to_symbol.get(g, "") for g in gene_ids])
        has_symbol = symbols != ""
        mean_mapped = mean_raw[:, has_symbol]
        symbols_mapped = symbols[has_symbol]

        # Deduplicate symbols (keep highest expressed)
        sym_to_best = {}
        for i, sym in enumerate(symbols_mapped):
            total = mean_mapped[:, i].sum()
            if sym not in sym_to_best or total > sym_to_best[sym][1]:
                sym_to_best[sym] = (i, total)

        keep_idx = [v[0] for v in sym_to_best.values()]
        keep_syms = list(sym_to_best.keys())
        mean_dedup = mean_mapped[:, keep_idx]

        # Normalize: target_sum=1e4, then log1p
        cluster_sums = mean_dedup.sum(axis=1, keepdims=True)
        cluster_sums[cluster_sums == 0] = 1
        log_norm = np.log1p(mean_dedup * (1e4 / cluster_sums))

        siletti_df = pd.DataFrame(
            log_norm.T,
            index=keep_syms,
            columns=[f"Siletti_{n}" for n in cluster_names],
        )

        # Load SEA-AD mean expression
        seaad_df = pd.read_csv(str(seaad_mean_expr_path), index_col=0)

        # Intersect genes and combine
        shared_genes = sorted(set(siletti_df.index) & set(seaad_df.index))
        combined = pd.concat([seaad_df.loc[shared_genes], siletti_df.loc[shared_genes]], axis=1)

    seaad_n = len([c for c in combined.columns if not c.startswith("Siletti_")])
    siletti_n = len([c for c in combined.columns if c.startswith("Siletti_")])
    print(f"    Combined: {combined.shape[0]:,} genes × {combined.shape[1]} types "
          f"({seaad_n} SEA-AD + {siletti_n} Siletti)")

    return combined
