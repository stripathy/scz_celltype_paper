#!/usr/bin/env python3
"""
Does our reprocessing reproduce Reiner et al.'s own published DE result?

This is the pipeline-validation step, independent of anything about our paper: if our
cell calling, QC, label transfer and limma-voom model recover the DE genes the original
authors reported from the same nuclei, the reprocessing is sound and any disagreement
with our meta-analysis is about biology rather than about our handling of the data.

Their Supplementary Table 4 reports 4,766 DE events in 2,994 genes across 16 of 20
clusters, from a MAST hurdle model on individual nuclei -- a different estimator and a
different (coarser, ad hoc) taxonomy than ours, so the expectation is directional
agreement and enrichment, not identical p-values.
"""
import numpy as np
import pandas as pd
from scipy import stats

XLSX = "/Users/shreejoy/Downloads/media-2 (6).xlsx"
OUTD = "/Users/shreejoy/Github/scz_celltype_paper/external_validation/gse158516/output"

# Reiner cluster ID -> the SEA-AD subclass(es) it most plausibly corresponds to.
# Their taxonomy is coarser and partly laminar; only the unambiguous ones are mapped.
CLUSTER_TO_SUBCLASS = {
    "SST+ interneurons": ["Sst"],
    "PVALB": ["Pvalb"],
    "VIP interneurons": ["Vip"],
    "Layer2/3 CUX2+ excitatory": ["L2/3 IT"],
    "Layer 4/5 excitatory": ["L4 IT", "L5 IT"],
    "Layer 5 excitatory": ["L5 IT"],
    "Layer 5 - HTR2C+ excitatory": ["L5 IT"],
    "Layer 5/6 excitatory": ["L5/6 NP", "L6 IT"],
    "Layer 6 NR4A2+ excitatory": ["L6b"],
    "Oligdendrocytes": ["Oligodendrocyte"],
    "OPC": ["OPC"],
    "Astrocytes": ["Astrocyte"],
    "Microglia": ["Microglia-PVM"],
    "Endothelial cells": ["Endothelial"],
}

reiner = pd.read_excel(XLSX, sheet_name="Supplementary Table 4", header=1)
reiner = reiner.rename(columns={"primerid": "gene", "coef (logFC)": "reiner_logFC",
                                "Pr(>Chisq)": "reiner_p", "fdr": "reiner_fdr",
                                "Cluster ID": "cluster"})
print(f"Reiner Supp. Table 4: {len(reiner)} DE events, {reiner.gene.nunique()} unique genes, "
      f"{reiner.cluster.nunique()} clusters")

ours = pd.read_csv(f"{OUTD}/de_subclass.csv.gz")

rows = []
for cluster, subclasses in CLUSTER_TO_SUBCLASS.items():
    r = reiner[reiner["cluster"] == cluster]
    if not len(r):
        continue
    for sc in subclasses:
        o = ours[ours["cell_type"] == sc]
        m = r.merge(o, on="gene", how="inner")
        if len(m) < 5:
            continue
        same = int(np.sum(np.sign(m["reiner_logFC"]) == np.sign(m["logFC"])))
        sp = stats.binomtest(same, len(m), 0.5).pvalue
        rho, rp = stats.spearmanr(m["reiner_logFC"], m["logFC"])
        # is our p-value distribution for their DE genes shifted vs all our genes?
        bg = o["p"].to_numpy()
        u, up = stats.mannwhitneyu(m["p"], bg, alternative="less")
        rows.append(dict(reiner_cluster=cluster, our_subclass=sc, n_genes=len(m),
                         sign_concordant=same, frac=same / len(m), sign_P=sp,
                         spearman=rho, spearman_P=rp,
                         median_our_p=float(np.median(m["p"])),
                         median_bg_p=float(np.median(bg)), enrich_P=up))

t = pd.DataFrame(rows).sort_values("n_genes", ascending=False)
print("\n=== per-cluster reproduction of Reiner et al.'s own DE genes ===")
print(t.round(4).to_string(index=False))
t.to_csv(f"{OUTD}/reproduce_reiner.csv", index=False)

tot_n = t["n_genes"].sum()
tot_c = t["sign_concordant"].sum()
print(f"\npooled across mapped clusters: {tot_c}/{tot_n} = {tot_c/tot_n:.1%} sign-concordant "
      f"(binomial P = {stats.binomtest(int(tot_c), int(tot_n), 0.5).pvalue:.3g})")

# The Sst cluster specifically -- the one this paper cares about
hdr = reiner[reiner["cluster"] == "SST+ interneurons"].merge(
    ours[ours["cell_type"] == "Sst"], on="gene", how="left")
print("\n=== Reiner's 29 SST+ interneuron DE genes, in our reprocessing ===")
print(hdr[["gene", "reiner_logFC", "reiner_fdr", "logFC", "p", "FDR"]]
      .sort_values("p").round(4).to_string(index=False))
