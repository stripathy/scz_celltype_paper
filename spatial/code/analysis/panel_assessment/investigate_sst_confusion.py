#!/usr/bin/env python3
"""
Investigate Sst supertype classification confusion in Xenium.

For each Sst supertype, computes:
1. Correlation margin (how confidently cells are classified)
2. Second-best supertype assignment (what would cells be called otherwise)
3. Whether margins differ between SCZ and Control (suggesting disease-driven
   misclassification)

Also computes per-gene expression differences between SCZ and Control within
each Sst supertype for all panel genes, to check whether SCZ downregulates
markers that distinguish subtypes.

Usage:
    python3 -u investigate_sst_confusion.py
"""

import os
import sys
import glob
import time
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import BASE_DIR, H5AD_DIR, EXCLUDE_SAMPLES, load_sample_adata

sys.path.insert(0, os.path.join(BASE_DIR, "code"))
from modules.metadata import get_subject_info
from modules.correlation_classifier import (
    build_supertype_centroids, correlate, build_flat_centroids
)

METADATA_PATH = os.path.join(BASE_DIR, "data", "sample_metadata.xlsx")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "density_analysis")

SST_TYPES = ["Sst_2", "Sst_3", "Sst_20", "Sst_22", "Sst_25"]


def main():
    t0 = time.time()
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    meta = get_subject_info(METADATA_PATH).set_index("sample_id")

    # Load all Sst cells with expression data
    h5ad_files = sorted(glob.glob(os.path.join(H5AD_DIR, "*_annotated.h5ad")))

    all_adatas = []
    for fpath in h5ad_files:
        sid = os.path.basename(fpath).replace("_annotated.h5ad", "")
        if sid in EXCLUDE_SAMPLES:
            continue
        adata = load_sample_adata(sid, cortical_only=True, qc_mode="corr")
        # Keep only Sst cells
        sst_mask = adata.obs["subclass_label"] == "Sst"
        if sst_mask.sum() > 0:
            adata_sst = adata[sst_mask].copy()
            adata_sst.obs["diagnosis"] = meta.loc[sid, "diagnosis"]
            all_adatas.append(adata_sst)
            print(f"  {sid}: {sst_mask.sum():,} Sst cells "
                  f"({meta.loc[sid, 'diagnosis']})")

    import anndata as ad
    adata_all = ad.concat(all_adatas, merge="same")
    print(f"\nTotal: {adata_all.n_obs:,} Sst cells")
    print(f"Supertypes: {adata_all.obs['supertype_label'].value_counts().to_dict()}")

    # ── Part 1: Build Sst supertype centroids and re-correlate ──
    print("\n" + "=" * 70)
    print("Part 1: Classification confidence by diagnosis")
    print("=" * 70)

    # Build centroids from all cells (same approach as the classifier)
    # Normalize
    from scipy.sparse import issparse
    X = adata_all.X
    if issparse(X):
        X = X.toarray()
    X = X.astype(np.float64)

    # log CP10K
    total_counts = X.sum(axis=1, keepdims=True)
    total_counts[total_counts == 0] = 1
    X_norm = np.log1p(X * 10000 / total_counts)

    gene_names = list(adata_all.var_names)
    sup_labels = adata_all.obs["supertype_label"].values.astype(str)
    diagnoses = adata_all.obs["diagnosis"].values.astype(str)

    # Build centroids per supertype (using ALL cells as reference)
    unique_sups = sorted(set(sup_labels))
    centroids_dict = {}
    for sup in unique_sups:
        mask = sup_labels == sup
        centroids_dict[sup] = X_norm[mask].mean(axis=0)
    centroids_df = pd.DataFrame(centroids_dict, index=gene_names).T
    print(f"Built {len(centroids_df)} Sst supertype centroids")

    # Correlate all cells against all Sst centroids
    corr_matrix, type_names = correlate(X_norm, centroids_df)

    # For each cell: best corr, margin, second-best type
    sorted_corr = np.sort(corr_matrix, axis=1)
    best_idx = np.argmax(corr_matrix, axis=1)
    best_corr = sorted_corr[:, -1]
    second_corr = sorted_corr[:, -2]
    margin = best_corr - second_corr

    # Second-best type
    n_cells = corr_matrix.shape[0]
    second_best = []
    for i in range(n_cells):
        corrs = corr_matrix[i].copy()
        corrs[best_idx[i]] = -np.inf
        second_best.append(type_names[np.argmax(corrs)])
    second_best = np.array(second_best)

    # Report margins by diagnosis for each vulnerable Sst type
    print("\n── Classification margins by diagnosis ──")
    margin_results = []
    for st in SST_TYPES:
        st_mask = sup_labels == st
        if st_mask.sum() < 10:
            continue

        ctrl_margins = margin[st_mask & (diagnoses == "Control")]
        scz_margins = margin[st_mask & (diagnoses == "SCZ")]

        if len(ctrl_margins) < 5 or len(scz_margins) < 5:
            continue

        t_stat, t_pval = ttest_ind(ctrl_margins, scz_margins)

        # Second-best confusion partners
        second_counts = pd.Series(second_best[st_mask]).value_counts()
        top_confusions = second_counts.head(3)

        print(f"\n  {st} (n_ctrl={len(ctrl_margins)}, n_scz={len(scz_margins)}):")
        print(f"    Ctrl margin: {ctrl_margins.mean():.4f} ± {ctrl_margins.std():.4f}")
        print(f"    SCZ  margin: {scz_margins.mean():.4f} ± {scz_margins.std():.4f}")
        print(f"    t-test p = {t_pval:.4g} (SCZ {'lower' if scz_margins.mean() < ctrl_margins.mean() else 'higher'})")
        print(f"    Top confusion partners: {dict(top_confusions)}")

        margin_results.append({
            "supertype": st,
            "n_ctrl": len(ctrl_margins),
            "n_scz": len(scz_margins),
            "ctrl_margin_mean": ctrl_margins.mean(),
            "scz_margin_mean": scz_margins.mean(),
            "margin_diff": scz_margins.mean() - ctrl_margins.mean(),
            "tstat": t_stat,
            "pval": t_pval,
            "top_confusion_1": top_confusions.index[0] if len(top_confusions) > 0 else "",
            "top_confusion_2": top_confusions.index[1] if len(top_confusions) > 1 else "",
        })

    margin_df = pd.DataFrame(margin_results)
    margin_df.to_csv(os.path.join(OUTPUT_DIR, "sst_margin_by_diagnosis.csv"), index=False)
    print(f"\nSaved: sst_margin_by_diagnosis.csv")

    # ── Part 2: Per-gene SCZ vs Control expression within each Sst type ──
    print("\n" + "=" * 70)
    print("Part 2: Per-gene SCZ effects within Sst supertypes")
    print("=" * 70)

    gene_results = []
    for st in SST_TYPES:
        st_mask = sup_labels == st
        ctrl_mask = st_mask & (diagnoses == "Control")
        scz_mask = st_mask & (diagnoses == "SCZ")

        if ctrl_mask.sum() < 20 or scz_mask.sum() < 20:
            continue

        X_ctrl = X_norm[ctrl_mask]
        X_scz = X_norm[scz_mask]

        for gi, gene in enumerate(gene_names):
            ctrl_expr = X_ctrl[:, gi]
            scz_expr = X_scz[:, gi]

            # Only test genes with some expression
            if ctrl_expr.mean() < 0.01 and scz_expr.mean() < 0.01:
                continue

            try:
                u_stat, u_pval = mannwhitneyu(ctrl_expr, scz_expr,
                                               alternative="two-sided")
            except Exception:
                continue

            logfc = np.log(scz_expr.mean() / ctrl_expr.mean()) if ctrl_expr.mean() > 0 else np.nan

            gene_results.append({
                "supertype": st,
                "gene": gene,
                "ctrl_mean": ctrl_expr.mean(),
                "scz_mean": scz_expr.mean(),
                "logFC": logfc,
                "pval": u_pval,
                "n_ctrl": int(ctrl_mask.sum()),
                "n_scz": int(scz_mask.sum()),
            })

    gene_df = pd.DataFrame(gene_results)

    # FDR per supertype
    for st in SST_TYPES:
        st_mask = gene_df["supertype"] == st
        if st_mask.sum() > 0:
            _, fdr, _, _ = multipletests(gene_df.loc[st_mask, "pval"], method="fdr_bh")
            gene_df.loc[st_mask, "fdr"] = fdr

    gene_df.to_csv(os.path.join(OUTPUT_DIR, "sst_pergene_scz_effects.csv"), index=False)

    # Report top DE genes for each type (focusing on known markers)
    print("\n── Top SCZ-affected genes per Sst supertype ──")
    for st in SST_TYPES:
        st_genes = gene_df[gene_df["supertype"] == st].sort_values("pval")
        n_sig = (st_genes["fdr"] < 0.05).sum()
        print(f"\n  {st}: {n_sig} genes FDR < 0.05 (of {len(st_genes)} tested)")
        for _, row in st_genes.head(10).iterrows():
            sig = "***" if row["fdr"] < 0.05 else ("*" if row["pval"] < 0.05 else "")
            print(f"    {row['gene']:15s} logFC={row['logFC']:+.3f}  "
                  f"ctrl={row['ctrl_mean']:.3f}  scz={row['scz_mean']:.3f}  "
                  f"p={row['pval']:.4g}  FDR={row['fdr']:.4g} {sig}")

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.0f}s")


if __name__ == "__main__":
    main()
