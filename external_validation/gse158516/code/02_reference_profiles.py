#!/usr/bin/env python3
"""
Build the SEA-AD MTG reference artifacts used to annotate GSE158516.

  * ref_profiles_{subclass,supertype}.parquet -- mean log-normalized expression per
    label, accumulated by streaming so the full 137k x 36.6k matrix is never held.
  * ref_pca.npz + ref_embed.npz -- a 50-PC space fit on HVGs of a supertype-stratified
    subsample of the reference, plus those cells' coordinates. Query nuclei get
    projected in and labelled by kNN vote -- the light-weight stand-in for the Seurat
    label transfer the main paper's pipeline uses.

Subsampling for the PCA step (not for the profiles) is what keeps this inside RAM:
a kNN vote does not need all 137k reference cells, it needs every *label* densely
represented, which CAP_PER_SUPERTYPE guarantees.
"""
import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

REF = "/Users/shreejoy/Github/shared_data/nicole_sea_ad_snrnaseq_reference.h5ad"
OUT = "/Users/shreejoy/Github/shared_data/GSE158516/reference"
QUERY_VAR = "/Users/shreejoy/Github/shared_data/GSE158516/raw/GSM4802078_CT_01_features.tsv.gz"
BLOCK = 10_000
CAP_PER_SUPERTYPE = 400
N_HVG = 3000
N_PCS = 50
SEED = 0

os.makedirs(OUT, exist_ok=True)

# ---- gene space: intersect reference (2020-A, 36,601) with query (3.0.0, 33,538) ----
qfeat = pd.read_csv(QUERY_VAR, sep="\t", header=None,
                    names=["gene_id", "gene_symbol", "feature_type"])
qsym = pd.Index(qfeat["gene_symbol"].astype(str))
r = ad.read_h5ad(REF, backed="r")
shared = r.var_names[r.var_names.isin(qsym[~qsym.duplicated()])]
gene_pos = pd.Series(np.arange(r.n_vars), index=r.var_names).loc[shared].to_numpy()
print(f"reference genes {r.n_vars}, query genes {qsym.nunique()}, shared {len(shared)}", flush=True)

obs = r.obs[["Class", "Subclass", "Supertype", "donor_id"]].astype(str).reset_index(drop=True)
n = r.n_obs


def norm_block(start, stop):
    """CP10K + log1p on the shared gene space for reference rows [start, stop)."""
    X = r.X[start:stop][:, gene_pos].astype(np.float32)
    tot = np.asarray(X.sum(1)).ravel()
    tot[tot == 0] = 1.0
    X = (sp.diags(1e4 / tot) @ X).tocsr()
    X.data = np.log1p(X.data)
    return X


# ---- pass 1: streaming mean profiles per label ------------------------------
levels = ["Subclass", "Supertype"]
cats = {lv: pd.Index(sorted(obs[lv].unique())) for lv in levels}
codes = {lv: cats[lv].get_indexer(obs[lv]) for lv in levels}
sums = {lv: np.zeros((len(cats[lv]), len(shared)), dtype=np.float64) for lv in levels}
cnts = {lv: np.zeros(len(cats[lv]), dtype=np.int64) for lv in levels}

for start in range(0, n, BLOCK):
    stop = min(start + BLOCK, n)
    X = norm_block(start, stop)
    for lv in levels:
        c = codes[lv][start:stop]
        # one-hot indicator @ X gives per-label sums in a single sparse product
        ind = sp.csr_matrix((np.ones(stop - start, dtype=np.float32),
                             (c, np.arange(stop - start))),
                            shape=(len(cats[lv]), stop - start))
        sums[lv] += np.asarray((ind @ X).todense(), dtype=np.float64)
        cnts[lv] += np.bincount(c, minlength=len(cats[lv]))
    del X
    if start % 50_000 == 0:
        print(f"  profiles block {start}/{n}", flush=True)

for lv in levels:
    M = (sums[lv] / np.maximum(cnts[lv], 1)[:, None]).astype(np.float32)
    df = pd.DataFrame(M, index=pd.Index(cats[lv], name=lv), columns=shared)
    df.to_parquet(f"{OUT}/ref_profiles_{lv.lower()}.parquet")
    pd.Series(cnts[lv], index=cats[lv], name="n_cells").to_csv(f"{OUT}/ref_ncells_{lv.lower()}.csv")
    print(f"wrote ref_profiles_{lv.lower()}.parquet {df.shape}", flush=True)
del sums

# ---- pass 2: stratified subsample -> HVG -> PCA ------------------------------
rng = np.random.default_rng(SEED)
take = []
for st, idx in obs.groupby("Supertype", observed=True).groups.items():
    idx = np.asarray(idx)
    take.append(idx if idx.size <= CAP_PER_SUPERTYPE
                else rng.choice(idx, CAP_PER_SUPERTYPE, replace=False))
take = np.sort(np.concatenate(take))
print(f"PCA subsample: {take.size} of {n} reference cells "
      f"({obs.loc[take, 'Supertype'].nunique()} supertypes)", flush=True)

blocks, kept = [], []
for start in range(0, n, BLOCK):
    stop = min(start + BLOCK, n)
    sel = take[(take >= start) & (take < stop)]
    if sel.size == 0:
        continue
    blocks.append(norm_block(start, stop)[sel - start])
    kept.append(sel)
Xs = sp.vstack(blocks).tocsr()
kept = np.concatenate(kept)
del blocks

a = ad.AnnData(X=Xs, obs=obs.loc[kept].reset_index(drop=True))
a.var_names = shared
sc.pp.highly_variable_genes(a, n_top_genes=N_HVG, flavor="seurat", batch_key="donor_id")
hvg = a.var_names[a.var["highly_variable"]].to_numpy()
a = a[:, hvg].copy()
sc.pp.scale(a, max_value=10)
sc.tl.pca(a, n_comps=N_PCS, svd_solver="arpack")

np.savez_compressed(f"{OUT}/ref_pca.npz",
                    hvg=hvg.astype(object),
                    mean=a.var["mean"].to_numpy().astype(np.float32),
                    std=a.var["std"].to_numpy().astype(np.float32),
                    components=a.varm["PCs"].astype(np.float32),
                    variance_ratio=a.uns["pca"]["variance_ratio"].astype(np.float32))
np.savez_compressed(f"{OUT}/ref_embed.npz",
                    X_pca=a.obsm["X_pca"].astype(np.float32),
                    Class=a.obs["Class"].to_numpy().astype(object),
                    Subclass=a.obs["Subclass"].to_numpy().astype(object),
                    Supertype=a.obs["Supertype"].to_numpy().astype(object),
                    donor_id=a.obs["donor_id"].to_numpy().astype(object))
print("wrote ref_pca.npz / ref_embed.npz; HVGs:", len(hvg), flush=True)
