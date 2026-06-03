#!/usr/bin/env python
"""
12_atac_peak_overlap.py — ATAC-seq peak overlap with FINEMAP credible sets.

For each SCZ GWAS fine-mapped variant, measures chromatin accessibility in each
cell type from SEA-AD snATAC-seq. Computes PIP-weighted accessibility scores
per cell type and compares with expression-based MAGMA enrichment.

Inputs:
  - data/fine_mapping/pgc3_finemap_credible_sets.csv (FINEMAP results)
  - data/atac/SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad (snATAC-seq)
  - results/tables/rbh_combined_enrichment.csv (MAGMA enrichment for comparison)

Outputs:
  - results/tables/atac_pip_weighted_scores.csv
  - results/intermediates/atac_finemap_peak_overlaps.csv
  - results/intermediates/atac_mean_accessibility_finemap_peaks.csv
  - results/figures/atac_vs_gwas_enrichment.png
  - results/figures/atac_vs_gwas_manhattan.png
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats

from scz_celltype_enrichment.config import (
    TABLES_DIR, FIGURES_DIR, INTERMEDIATES_DIR, PROJECT_ROOT,
    FIGURE_DPI, TITLE_FONTSIZE, LABEL_FONTSIZE, TICK_FONTSIZE, LEGEND_FONTSIZE,
)
from scz_celltype_enrichment.utils import (
    ensure_dirs, classify_supertypes_by_class, extract_subclass,
)


ATAC_H5AD = PROJECT_ROOT / "data" / "atac" / "SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad"
FINEMAP_CSV = PROJECT_ROOT / "data" / "fine_mapping" / "pgc3_finemap_credible_sets.csv"


def find_peak_variant_overlaps(finemap_df, h5ad_file):
    """Find which ATAC peaks contain FINEMAP credible set variants."""
    f = h5py.File(str(h5ad_file), "r")

    # Parse peak coordinates
    peaks_raw = f["var"]["_index"][:]
    peaks_str = [p.decode() if isinstance(p, bytes) else str(p) for p in peaks_raw]

    peak_chr = []
    peak_start = []
    peak_end = []
    for p in peaks_str:
        parts = p.split(":")
        coords = parts[1].split("-")
        peak_chr.append(parts[0])
        peak_start.append(int(coords[0]))
        peak_end.append(int(coords[1]))

    peak_chr = np.array(peak_chr)
    peak_start = np.array(peak_start)
    peak_end = np.array(peak_end)

    # Build chromosome index
    chr_peak_idx = {}
    for chrom in np.unique(peak_chr):
        chr_peak_idx[chrom] = np.where(peak_chr == chrom)[0]

    # Find overlaps
    overlaps = []
    for _, var in finemap_df.iterrows():
        chrom = "chr" + str(var["chr"])
        pos = var["pos"]
        if chrom not in chr_peak_idx:
            continue
        idx = chr_peak_idx[chrom]
        hits = idx[(peak_start[idx] <= pos) & (peak_end[idx] >= pos)]
        for peak_i in hits:
            overlaps.append({
                "variant_rsid": var["rsid"],
                "pip": var["pip"],
                "locus": var["locus"],
                "gwas_p": var["gwas_p"],
                "peak_idx": int(peak_i),
                "peak": peaks_str[peak_i],
            })

    f.close()
    return pd.DataFrame(overlaps), peaks_str


def compute_per_supertype_accessibility(h5ad_file, overlap_peak_idx, peaks_str):
    """Compute mean accessibility per supertype at overlapping peaks."""
    f = h5py.File(str(h5ad_file), "r")

    n_peaks = len(peaks_str)
    n_overlap = len(overlap_peak_idx)

    # Cell type annotations
    supertype_cats = f["obs"]["Supertype"]["categories"][:]
    supertype_codes = f["obs"]["Supertype"]["codes"][:]
    supertype_names = [c.decode() if isinstance(c, bytes) else str(c) for c in supertype_cats]
    st_codes = supertype_codes[:]

    unique_st = sorted(set(supertype_names))
    n_st = len(unique_st)
    st_to_row = {st: i for i, st in enumerate(unique_st)}
    code_to_row = np.array([st_to_row[supertype_names[c]] for c in range(len(supertype_names))])

    indptr = f["X"]["indptr"][:]
    n_cells = len(indptr) - 1

    # Cell counts per supertype
    st_counts = np.zeros(n_st, dtype=np.int64)
    for c in st_codes:
        st_counts[code_to_row[c]] += 1

    # Accumulate accessibility in chunks
    st_sums = np.zeros((n_st, n_overlap), dtype=np.float64)
    CHUNK = 50000
    n_chunks = (n_cells + CHUNK - 1) // CHUNK

    print(f"  Processing {n_cells:,} cells in {n_chunks} chunks...")
    for chunk_i in range(n_chunks):
        start = chunk_i * CHUNK
        end = min(start + CHUNK, n_cells)

        ptr_start = int(indptr[start])
        ptr_end = int(indptr[end])

        chunk_data = f["X"]["data"][ptr_start:ptr_end]
        chunk_indices = f["X"]["indices"][ptr_start:ptr_end]
        chunk_indptr = indptr[start:end + 1] - ptr_start

        chunk_csr = sparse.csr_matrix(
            (chunk_data, chunk_indices, chunk_indptr),
            shape=(end - start, n_peaks),
        )

        # Extract overlapping peak columns
        chunk_overlap = chunk_csr[:, overlap_peak_idx].toarray()

        # Aggregate by supertype
        chunk_st = code_to_row[st_codes[start:end]]
        for st_row in range(n_st):
            mask = chunk_st == st_row
            if mask.any():
                st_sums[st_row] += chunk_overlap[mask].sum(axis=0)

        if (chunk_i + 1) % 5 == 0 or chunk_i == n_chunks - 1:
            print(f"    Chunk {chunk_i + 1}/{n_chunks}")

    f.close()

    st_means = st_sums / np.maximum(st_counts[:, np.newaxis], 1)
    mean_acc = pd.DataFrame(
        st_means, index=unique_st,
        columns=[peaks_str[i] for i in overlap_peak_idx],
    )
    return mean_acc, st_counts, unique_st


def compute_pip_weighted_scores(mean_acc, overlaps_df, peaks_str, overlap_peak_idx,
                                 unique_st, st_counts):
    """Compute PIP-weighted cell-type accessibility scores."""
    # Sum PIP per peak
    peak_pip = {}
    for _, row in overlaps_df.iterrows():
        p = peaks_str[row["peak_idx"]]
        peak_pip[p] = peak_pip.get(p, 0) + row["pip"]

    pip_vec = np.array([peak_pip.get(peaks_str[i], 0) for i in overlap_peak_idx])
    scores = mean_acc.values @ pip_vec

    return pd.DataFrame({
        "supertype": unique_st,
        "pip_weighted_atac_score": scores,
        "n_cells": st_counts.astype(int),
    }).sort_values("pip_weighted_atac_score", ascending=False).reset_index(drop=True)


def plot_comparison(scores_df, enrich_df, output_dir):
    """Generate comparison figures: ATAC vs GWAS enrichment."""
    mpl.rcParams.update({
        "figure.facecolor": "white", "axes.facecolor": "white",
        "axes.grid": False, "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 1.2, "font.size": 14,
    })

    seaad_enrich = enrich_df[enrich_df["source"] == "SEA-AD"].copy()
    merged = scores_df.merge(seaad_enrich, on="supertype", how="inner")
    merged["gwas_logp"] = -np.log10(np.clip(merged["p_value"], 1e-300, 1))

    cls_map = classify_supertypes_by_class(merged["supertype"].tolist())
    merged["class"] = [cls_map.get(ct, "Other") for ct in merged["supertype"]]

    class_colors = {"GABAergic": "#D55E00", "Glutamatergic": "#0072B2", "Non-neuronal": "#009E73"}
    LS = 18; TS = 22; TKS = 15; LGS = 13; ANNOT = 10

    r_all, p_all = stats.spearmanr(merged["pip_weighted_atac_score"], merged["gwas_logp"])

    # Manhattan comparison
    fig, axes = plt.subplots(2, 1, figsize=(20, 14), gridspec_kw={"height_ratios": [1, 1]})

    gaba_pref = ["Chandelier", "Lamp5", "Pax6", "Pvalb", "Sncg", "Sst", "Sst Chodl", "Vip"]
    def sort_key(name):
        sc = extract_subclass(name)
        cls = cls_map.get(name, "Other")
        if cls == "GABAergic":
            idx = gaba_pref.index(sc) if sc in gaba_pref else 99
            return (0, idx, name)
        elif cls == "Glutamatergic":
            return (1, 0, name)
        return (2, 0, name)

    ordered = merged.sort_values("supertype", key=lambda s: s.map(sort_key)).reset_index(drop=True)
    ordered["x"] = range(len(ordered))

    for panel, (ax, ycol, ylabel, title) in enumerate(zip(
        axes,
        ["gwas_logp", "pip_weighted_atac_score"],
        [r"GWAS $-\log_{10}(p)$", "PIP-weighted ATAC score"],
        ["Expression-based GWAS enrichment", "Chromatin accessibility at fine-mapped variants"],
    )):
        for cls in ["Non-neuronal", "Glutamatergic", "GABAergic"]:
            mask = ordered["class"] == cls
            s = 50 if panel == 1 else np.where(ordered.loc[mask, "p_fdr"] < 0.05, 50, 20)
            ax.scatter(ordered.loc[mask, "x"], ordered.loc[mask, ycol],
                       c=class_colors[cls], s=s, alpha=0.7, label=cls,
                       edgecolors="black", linewidth=0.3)

        top = ordered.nlargest(10, ycol)
        for _, row in top.iterrows():
            ax.annotate(row["supertype"], (row["x"], row[ycol]),
                        fontsize=8, rotation=30, ha="left", va="bottom",
                        xytext=(3, 5), textcoords="offset points")

        ax.set_ylabel(ylabel, fontsize=LS)
        ax.set_title(title, fontsize=TS, fontweight="bold")
        ax.set_xticks([])
        ax.tick_params(labelsize=TKS)
        if panel == 1:
            ax.legend(fontsize=LGS, frameon=False)
            ax.set_xlabel("Cell type (GABAergic → Glutamatergic → Non-neuronal)", fontsize=LS)

    plt.tight_layout()
    fig.savefig(str(output_dir / "atac_vs_gwas_manhattan.png"), dpi=150,
                bbox_inches="tight", facecolor="white")
    plt.close()

    # Scatter comparison
    fig, ax = plt.subplots(figsize=(10, 9))
    for cls in ["Non-neuronal", "Glutamatergic", "GABAergic"]:
        mask = merged["class"] == cls
        ax.scatter(merged.loc[mask, "pip_weighted_atac_score"],
                   merged.loc[mask, "gwas_logp"],
                   c=class_colors[cls], s=60, alpha=0.7, label=cls,
                   edgecolors="black", linewidth=0.3)

    sig = merged[(merged["p_fdr"] < 0.05) | (merged["pip_weighted_atac_score"] > 50)]
    for _, row in sig.iterrows():
        ax.annotate(row["supertype"], (row["pip_weighted_atac_score"], row["gwas_logp"]),
                    fontsize=ANNOT, xytext=(4, 4), textcoords="offset points")

    ax.text(0.03, 0.97, f"r = {r_all:.2f}\np = {p_all:.1e}",
            transform=ax.transAxes, fontsize=LGS, va="top", fontstyle="italic", color="#555")
    ax.set_xlabel("PIP-weighted ATAC score", fontsize=LS)
    ax.set_ylabel(r"GWAS enrichment $-\log_{10}(p)$", fontsize=LS)
    ax.set_title("ATAC vs expression-based enrichment", fontsize=TS, fontweight="bold")
    ax.tick_params(labelsize=TKS)
    ax.legend(fontsize=LGS, frameon=False)

    plt.tight_layout()
    fig.savefig(str(output_dir / "atac_vs_gwas_enrichment.png"), dpi=150,
                bbox_inches="tight", facecolor="white")
    plt.close()
    mpl.rcParams.update(mpl.rcParamsDefault)


def main():
    ensure_dirs(TABLES_DIR, FIGURES_DIR, INTERMEDIATES_DIR)

    print("=" * 60)
    print("Step 12: ATAC Peak Overlap with FINEMAP Credible Sets")
    print("=" * 60)
    t0 = time.time()

    # Check inputs
    if not ATAC_H5AD.exists():
        print(f"ERROR: ATAC h5ad not found: {ATAC_H5AD}")
        return
    if not FINEMAP_CSV.exists():
        print(f"ERROR: FINEMAP CSV not found: {FINEMAP_CSV}")
        return

    # Load FINEMAP
    print("\n1. Loading FINEMAP credible sets...")
    finemap = pd.read_csv(str(FINEMAP_CSV))
    print(f"  {len(finemap)} variants, {finemap['locus'].nunique()} loci")

    # Find overlaps (fast)
    print("\n2. Finding peak-variant overlaps...")
    overlaps_df, peaks_str = find_peak_variant_overlaps(finemap, ATAC_H5AD)
    overlap_peak_idx = sorted(overlaps_df["peak_idx"].unique())
    print(f"  {len(overlaps_df)} overlaps, {len(overlap_peak_idx)} peaks, "
          f"{overlaps_df['locus'].nunique()} loci")

    overlaps_df.to_csv(str(INTERMEDIATES_DIR / "atac_finemap_peak_overlaps.csv"), index=False)

    # Compute per-supertype accessibility (slow — reads 18GB h5ad)
    cached = INTERMEDIATES_DIR / "atac_mean_accessibility_finemap_peaks.csv"
    if cached.exists():
        print(f"\n3. Loading cached mean accessibility...")
        mean_acc = pd.read_csv(str(cached), index_col=0)
        unique_st = mean_acc.index.tolist()
        # Reconstruct st_counts from h5ad
        f = h5py.File(str(ATAC_H5AD), "r")
        cats = f["obs"]["Supertype"]["categories"][:]
        codes = f["obs"]["Supertype"]["codes"][:]
        names = [c.decode() if isinstance(c, bytes) else str(c) for c in cats]
        st_to_row = {st: i for i, st in enumerate(unique_st)}
        st_counts = np.zeros(len(unique_st), dtype=np.int64)
        for c in codes[:]:
            n = names[c]
            if n in st_to_row:
                st_counts[st_to_row[n]] += 1
        f.close()
    else:
        print(f"\n3. Computing per-supertype accessibility...")
        mean_acc, st_counts, unique_st = compute_per_supertype_accessibility(
            ATAC_H5AD, overlap_peak_idx, peaks_str
        )
        mean_acc.to_csv(str(cached))

    print(f"  Mean accessibility: {mean_acc.shape}")

    # PIP-weighted scores
    print("\n4. Computing PIP-weighted scores...")
    scores_df = compute_pip_weighted_scores(
        mean_acc, overlaps_df, peaks_str, overlap_peak_idx, unique_st, st_counts
    )
    scores_df.to_csv(str(TABLES_DIR / "atac_pip_weighted_scores.csv"), index=False)

    print(f"\n  Top 15 by PIP-weighted ATAC score:")
    for i, (_, row) in enumerate(scores_df.head(15).iterrows()):
        print(f"    {i + 1:>2d}. {row['supertype']:25s}  score={row['pip_weighted_atac_score']:.2f}")

    # Compare with MAGMA enrichment
    print("\n5. Generating comparison figures...")
    enrich_path = TABLES_DIR / "rbh_combined_enrichment.csv"
    if enrich_path.exists():
        enrich = pd.read_csv(str(enrich_path))
        plot_comparison(scores_df, enrich, FIGURES_DIR)
        print("  Saved figures")
    else:
        print("  WARNING: rbh_combined_enrichment.csv not found, skipping figures")

    print(f"\nStep 12 complete. Total time: {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
