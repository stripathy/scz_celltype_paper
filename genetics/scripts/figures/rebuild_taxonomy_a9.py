"""
Rebuild the SEA-AD/Siletti merge on the Duncan-normalized expression
(per-cell ln(1+raw) averaged within type, ENSEMBL space).

Reproduces the existing matching rule exactly -- CV-based HVG selection,
Spearman correlation between type profiles, the row-wise rank statistic the
codebase calls "AUROC", and reciprocal best hits on it -- but additionally
records the underlying correlation for each pair.

That matters because the "AUROC" is 1.0 by construction for every best hit
(it is the rank of the top match among the other candidates), so the documented
"mean AUROC >= 0.9" bar cannot fail and carries no quality information. The
correlation is the quantity that actually distinguishes a redundant pair from a
forced one.
"""
import pandas as pd
import numpy as np
from scipy.stats import rankdata
import h5py

W = "/Users/shreejoy/Github/scz_celltype_paper/genetics/results/intermediates"
N_HVG = 5000

ens_map = pd.read_csv(f"{W}/map_entrez_ensembl.csv", dtype=str)
sym_map = pd.read_csv(f"{W}/map_entrez_symbol.csv", dtype=str)

# A9/DLPFC reference: the figure's expression data is A9 throughout, so the
# taxonomy match is made in the same data rather than in MTG.
sea = pd.read_csv(f"{W}/a9_supertype_log1p_mean.csv", index_col=0)
s2e = sym_map.merge(ens_map, on="ENTREZ")[["SYMBOL", "ENSEMBL"]]
s2e = s2e[~s2e.SYMBOL.duplicated(keep=False) & ~s2e.ENSEMBL.duplicated(keep=False)]
sea = sea[~sea.index.duplicated(keep=False)].join(
    s2e.set_index("SYMBOL"), how="inner").set_index("ENSEMBL")

with h5py.File(f"{W}/Siletti_L2-cluster-log1p_matrix.h5", "r") as f:
    mat, acc = f["matrix"][:], np.array([a.decode().split(".")[0] for a in f["Accession"][:]])
names = pd.read_csv(f"{W}/siletti_cluster_names.csv")
n2n = dict(zip(names.cluster_num, names.cluster_name))
sil = pd.DataFrame(mat, index=acc, columns=[n2n.get(i, f"Cluster{i}") for i in range(mat.shape[1])])
sil = sil[~sil.index.duplicated(keep=False)]
print(f"SEA-AD {sea.shape} | Siletti {sil.shape}")

# --- HVGs: shared genes, mean > 0.1, top-N by CV (existing rule) ---
shared = sorted(set(sea.index) & set(sil.index))
comb = pd.concat([sea.loc[shared], sil.loc[shared]], axis=1)
gm, gs = comb.mean(axis=1), comb.std(axis=1)
expressed = gm > 0.1
cv = (gs[expressed] / gm[expressed]).replace([np.inf, -np.inf], np.nan).dropna()
hvgs = cv.nlargest(min(N_HVG, len(cv))).index.tolist()
print(f"HVGs: {len(shared):,} shared -> {int(expressed.sum()):,} expressed -> {len(hvgs)} HVG")

# --- Spearman correlation between type profiles ---
A, B = sea.loc[hvgs], sil.loc[hvgs]
corr = np.corrcoef(A.rank(axis=0).values.T, B.rank(axis=0).values.T)[:A.shape[1], A.shape[1]:]
corr_df = pd.DataFrame(corr, index=A.columns, columns=B.columns)
print(f"correlation range [{corr.min():.3f}, {corr.max():.3f}]")

# --- the codebase's row-wise rank statistic ("AUROC") ---
nb = corr.shape[1]
auroc = np.apply_along_axis(lambda r: (rankdata(r) - 1) / (nb - 1), 1, corr)
auroc_df = pd.DataFrame(auroc, index=A.columns, columns=B.columns)

a_best, b_best = auroc_df.idxmax(axis=1), auroc_df.idxmax(axis=0)
rows = []
for a in auroc_df.index:
    b = a_best[a]
    if b_best[b] == a:
        rows.append({"seaad_type": a, "siletti_cluster": b,
                     "mean_auroc": (auroc_df.loc[a, b] + auroc_df[b].max()) / 2,
                     "spearman_r": corr_df.loc[a, b],
                     "second_best_r": corr_df.loc[a].nlargest(2).iloc[1]})
rbh = pd.DataFrame(rows).sort_values("spearman_r", ascending=False).reset_index(drop=True)
rbh["r_margin"] = rbh.spearman_r - rbh.second_best_r
rbh.to_csv(f"{W}/franken_rbh_A9.csv", index=False)

old = pd.read_csv("/Users/shreejoy/Github/scz_celltype_paper/genetics/franken_taxonomy/franken_rbh.csv")
print(f"\nRBH pairs: {len(rbh)} (previously {len(old)})")
print(f"  mean_auroc: min {rbh.mean_auroc.min():.4f}  <- saturated by construction")
print(f"\nSpearman r across the {len(rbh)} RBH pairs -- the real quality measure:")
print(rbh.spearman_r.describe()[["min", "25%", "50%", "75%", "max"]].to_string())
for t in (0.5, 0.6, 0.7, 0.8, 0.9):
    print(f"  r >= {t}: {(rbh.spearman_r >= t).sum():3d} of {len(rbh)}")
print(f"\nmargin over 2nd-best correlation: median {rbh.r_margin.median():.4f}, "
      f"< 0.01 for {(rbh.r_margin < 0.01).sum()} of {len(rbh)}")

same = set(zip(old.seaad_type, old.siletti_cluster)) & set(zip(rbh.seaad_type, rbh.siletti_cluster))
print(f"\npairs identical to the old taxonomy: {len(same)} of {len(rbh)}")
print(f"Siletti clusters removed then vs now: {len(set(old.siletti_cluster) & set(rbh.siletti_cluster))} shared")
print("\nweakest 10 merges by correlation:")
print(rbh.nsmallest(10, "spearman_r")[
    ["seaad_type", "siletti_cluster", "spearman_r", "second_best_r", "r_margin"]
].to_string(index=False, float_format=lambda v: f"{v:.4f}"))
