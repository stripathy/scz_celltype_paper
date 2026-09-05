"""
Rebuild the four Figure 4 panels that still read the MTG reference -- d, e
(HCN1 expression) and h, i (vulnerable-vs-not-depleted volcano and CALB1
violin) -- from the SEA-AD A9/DLPFC neurotypical reference instead.

Statistics and normalisation follow the retired export_fig4_new_panels.py:
CP10K + log1p, cell-level Wilcoxon for the volcano, Seurat-style avg_log2FC on
the un-logged CP10K scale, BH FDR, and a donor-level paired t-test for the
violin annotation.

NOTE the donor-level test drops from n = 5 donors (MTG) to n = 3 (DLPFC),
because A9 has only three neurotypical reference brains. That test is
correspondingly weaker and its p-values should be read with that in mind.
"""
import h5py
import numpy as np
import scipy.sparse as sp
import pandas as pd
from scipy import stats
import glob

import os as _os

# Repo-relative, so the chain runs from any clone. External data that is not in
# the repo is still resolved by absolute path or an environment override below.
GEN = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
PAPER = _os.path.dirname(GEN)

W = f"{GEN}/results/intermediates"
OUT = f"{GEN}/results/figures/r_panels"


VULNERABLE = ["Sst_25", "Sst_22", "Sst_2", "Sst_20", "Sst_3"]
NOT_DEPLETED = ["Sst_19", "Sst_9", "Sst_23", "Sst_11", "Sst_13",
                "Sst_1", "Sst_4", "Sst_5", "Sst_7", "Sst_10", "Sst_12"]
KEEP = set(VULNERABLE + NOT_DEPLETED)

mats, sts, donors, genes = [], [], [], None
for path in sorted(glob.glob(_os.environ.get("SEAAD_A9_GLOB",
        _os.path.expanduser("~/Downloads/*A9_RNAseq_final-nuclei.2024-02-13.h5ad")))):
    with h5py.File(path, "r") as f:
        d = f["obs/Supertype"]
        cats = np.array([x.decode() if isinstance(x, bytes) else x
                         for x in f[d.attrs["categories"]][:]])
        codes = d[:]
        st = cats[np.clip(codes, 0, None)]
        st[codes < 0] = ""
        sel = np.isin(st, list(KEEP))
        if genes is None:
            genes = np.array([g.decode() for g in f["var/_index"][:]])
        n_cells, n_genes = f["X"].attrs["shape"]
        indptr = f["layers/UMIs/indptr"][:]
        rows = np.flatnonzero(sel)
        blocks = []
        for i in rows:
            p0, p1 = int(indptr[i]), int(indptr[i + 1])
            blocks.append(sp.csr_matrix(
                (f["layers/UMIs/data"][p0:p1], f["layers/UMIs/indices"][p0:p1], [0, p1 - p0]),
                shape=(1, n_genes)))
        mats.append(sp.vstack(blocks).tocsr())
        sts.append(st[sel])
        dn = np.array([x.decode() if isinstance(x, bytes) else x
                       for x in f[f["obs/Donor ID"].attrs["categories"]][:]])[0]
        donors.append(np.full(sel.sum(), dn))
        print(f"  {dn}: {sel.sum():,} Sst nuclei", flush=True)

X = sp.vstack(mats).tocsr().astype(np.float64)
supertype = np.concatenate(sts)
donor = np.concatenate(donors)
group = np.where(np.isin(supertype, VULNERABLE), "Vulnerable", "Not depleted")
print(f"total {X.shape[0]:,} Sst nuclei | Vulnerable {(group=='Vulnerable').sum():,} | "
      f"Not depleted {(group=='Not depleted').sum():,} | donors {len(set(donor))}")

# CP10K + log1p, exactly as scanpy normalize_total(1e4) + log1p
tot = np.asarray(X.sum(axis=1)).ravel()
tot[tot == 0] = 1
X = sp.diags(1e4 / tot) @ X
cp10k = X.copy()                       # un-logged CP10K, for Seurat-style FC
X.data = np.log1p(X.data)

dense = np.asarray(X.todense())
c10 = np.asarray(cp10k.todense())
is_v = group == "Vulnerable"

stat, pval = stats.mannwhitneyu(dense[is_v], dense[~is_v], axis=0, alternative="two-sided")
mv, mn = c10[is_v].mean(axis=0), c10[~is_v].mean(axis=0)
log2fc = np.log2((mv + 1.0) / (mn + 1.0))


def bh(p):
    n = len(p)
    order = np.argsort(p)
    ranked = np.asarray(p)[order]
    f = np.minimum.accumulate((ranked * n / (np.arange(n) + 1))[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(f, 0, 1)
    return out


uds = sorted(set(donor))
dmean_v = np.vstack([dense[(donor == d) & is_v].mean(axis=0) for d in uds])
dmean_n = np.vstack([dense[(donor == d) & ~is_v].mean(axis=0) for d in uds])
t_d, p_d = stats.ttest_rel(dmean_v, dmean_n, axis=0)
p_d = np.nan_to_num(p_d, nan=1.0)

res = pd.DataFrame({
    "gene": genes, "log2FC": log2fc,
    "pct_vulnerable": (c10[is_v] > 0).mean(axis=0),
    "pct_not_depleted": (c10[~is_v] > 0).mean(axis=0),
    "mean_vulnerable": dense[is_v].mean(axis=0),
    "mean_not_depleted": dense[~is_v].mean(axis=0),
    "p_wilcoxon": pval, "fdr_wilcoxon": bh(pval),
    "p_donor_paired": p_d, "fdr_donor": bh(p_d),
}).sort_values("p_wilcoxon")
res["neg_log10_p"] = -np.log10(np.clip(res.p_wilcoxon, 1e-300, None))
res.to_csv(f"{OUT}/panel_volcano_vulnerable_vs_notdepleted.csv", index=False)
print(f"volcano: {len(res):,} genes")
for g in ["CALB1", "HCN1", "SST", "NOS1", "CALB2", "NPY", "RELN"]:
    r = res[res.gene == g]
    if len(r):
        r = r.iloc[0]
        print(f"  {g:6s} log2FC={r.log2FC:+.3f}  wilcox p={r.p_wilcoxon:.2e}  "
              f"FDR={r.fdr_wilcoxon:.2e}  donor-paired p={r.p_donor_paired:.4f}")

viol = []
for g in ["CALB1", "SST"]:
    gi = np.where(genes == g)[0]
    if not len(gi):
        continue
    viol.append(pd.DataFrame({"gene": g, "expression": dense[:, int(gi[0])],
                              "group": group, "donor": donor, "supertype": supertype}))
viol = pd.concat(viol, ignore_index=True)
viol.to_csv(f"{OUT}/panel_violin_calb1_percell.csv", index=False)

dsum = (viol.groupby(["gene", "donor", "group"], observed=True).expression.mean()
        .reset_index().pivot_table(index=["gene", "donor"], columns="group",
                                   values="expression").reset_index())
fdr_lookup = dict(zip(res.gene, res.fdr_donor))
rows = []
for g, gg in dsum.groupby("gene"):
    gg = gg.dropna(subset=["Vulnerable", "Not depleted"])
    tt = stats.ttest_rel(gg["Vulnerable"], gg["Not depleted"])
    rows.append(dict(gene=g, n_donors=len(gg), mean_vulnerable=gg["Vulnerable"].mean(),
                     mean_not_depleted=gg["Not depleted"].mean(), t=tt.statistic,
                     p_donor_paired=tt.pvalue, fdr_donor=fdr_lookup.get(g, np.nan)))
    print(f"  violin {g}: n={len(gg)} donors, paired p={tt.pvalue:.4f}, "
          f"FDR={fdr_lookup.get(g, float('nan')):.4f}")
dsum.to_csv(f"{OUT}/panel_violin_donor_means.csv", index=False)

# The violin is annotated with the cell-level BH FDR, not the donor-paired one.
# Both are computed here, but the volcano beside it defines its gene set on the
# cell-level test; annotating the violin with the donor test made the two panels
# report different statistics for the same comparison. The donor-paired p is
# retained in the CSV and quoted in the text, where the three-donor caveat can
# be stated.
cell_fdr = dict(zip(res.gene, res.fdr_wilcoxon))
stats_df = pd.DataFrame(rows)
stats_df["fdr_cell"] = stats_df.gene.map(cell_fdr)
stats_df.to_csv(f"{OUT}/panel_violin_stats.csv", index=False)
for r in stats_df.itertuples():
    print(f"  violin {r.gene}: cell-level FDR={r.fdr_cell:.3g} "
          f"(donor-paired p={r.p_donor_paired:.4f})")

# ---- panels d/e: per-supertype mean expression (same CP10K+log1p scale) ----
means = pd.DataFrame(
    {s: dense[supertype == s].mean(axis=0) for s in sorted(set(supertype))}, index=genes)
means.to_csv(f"{W}/a9_sst_supertype_mean_expression_cp10k.csv")
ephys = pd.read_csv(f"{GEN}/results/tables/sst_supertype_ephys_summary.csv")
ephys["HCN1_expr"] = ephys.supertype.map(means.loc["HCN1"].to_dict())
dfE = ephys.dropna(subset=["HCN1_expr", "mean_sag"]).copy()
old = pd.read_csv(f"{GEN}/results/figures/r_panels/panel_E_hcn1_vs_sag.csv")
dfE = dfE[[c for c in old.columns if c in dfE.columns]]
for c in old.columns:
    if c not in dfE.columns:
        dfE[c] = dfE.supertype.map(dict(zip(old.supertype, old[c])))
dfE[old.columns].to_csv(f"{OUT}/panel_E_hcn1_vs_sag.csv", index=False)
rho, p = stats.spearmanr(dfE.HCN1_expr, dfE.mean_sag)
print(f"panel E (DLPFC HCN1 vs sag, n={len(dfE)}): rho={rho:+.3f}, p={p:.5f}")
print(f"  HCN1 in Sst_25 = {means.loc['HCN1','Sst_25']:.4f} (MTG value was 2.4404)")
