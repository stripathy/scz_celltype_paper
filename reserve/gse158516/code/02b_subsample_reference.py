#!/usr/bin/env python3
"""
Write a subsample of the SEA-AD MTG reference for kNN label transfer.

The full reference is 137,303 x 36,601 with int64 counts -- 0.78B non-zeros at 12
bytes each is ~9.4 GB in memory, which alone blows the budget before the query is even
loaded. Counts are downcast to float32 (they are small integers) and the gene space is
cut to what the query shares.

MODE matters, and getting it wrong corrupts the composition analysis:

  "proportional" (default) keeps each supertype's share of the reference intact, taking
      FRACTION of every supertype with a FLOOR so rare types survive. kNN voting is
      sensitive to the reference's class priors, so this is what a compositional
      comparison requires.
  "capped" takes up to CAP cells per supertype. It represents rare labels more densely,
      but it *flattens the class balance*: at CAP=400 it inflates Vip 2.1x and
      Sncg/Pax6/Chandelier 3.6x while shrinking L4 IT to 0.34x, pushing GABAergic from
      32.7% to 60.6% of reference neurons. Retained only as a sensitivity check.
"""
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

REF = "/Users/shreejoy/Github/shared_data/nicole_sea_ad_snrnaseq_reference.h5ad"
QUERY_VAR = "/Users/shreejoy/Github/shared_data/GSE158516/raw/GSM4802078_CT_01_features.tsv.gz"
OUTDIR = "/Users/shreejoy/Github/shared_data/GSE158516/reference"

MODE = sys.argv[1] if len(sys.argv) > 1 else "proportional"
FRACTION = 0.45          # proportional mode: ~62k cells
FLOOR = 25               # ...but never fewer than this per supertype
CAP = 400                # capped mode
OUT = f"{OUTDIR}/seaad_mtg_{MODE}.h5ad"
BLOCK = 10_000
SEED = 0

r = ad.read_h5ad(REF, backed="r")
qsym = pd.Index(pd.read_csv(QUERY_VAR, sep="\t", header=None,
                            names=["id", "sym", "type"])["sym"].astype(str))
shared = r.var_names[r.var_names.isin(qsym[~qsym.duplicated()])]
gene_pos = pd.Series(np.arange(r.n_vars), index=r.var_names).loc[shared].to_numpy()

obs = r.obs[["Class", "Subclass", "Supertype", "donor_id"]].astype(str).reset_index(drop=True)
rng = np.random.default_rng(SEED)
take = []
for _, idx in obs.groupby("Supertype", observed=True).groups.items():
    idx = np.asarray(idx)
    n = (min(idx.size, CAP) if MODE == "capped"
         else min(idx.size, max(FLOOR, int(round(FRACTION * idx.size)))))
    take.append(idx if n >= idx.size else rng.choice(idx, n, replace=False))
take = np.sort(np.concatenate(take))
print(f"[{MODE}] subsampling {take.size:,} of {r.n_obs:,} cells, "
      f"{len(shared):,} shared genes", flush=True)

blocks = []
for start in range(0, r.n_obs, BLOCK):
    stop = min(start + BLOCK, r.n_obs)
    sel = take[(take >= start) & (take < stop)]
    if sel.size:
        blocks.append(r.X[start:stop][:, gene_pos][sel - start].astype(np.float32))
X = sp.vstack(blocks).tocsr()
del blocks

a = ad.AnnData(X=X, obs=obs.loc[take].reset_index(drop=True))
a.obs_names = r.obs_names[take]
a.var_names = shared
a.write_h5ad(OUT, compression="gzip")
print(f"wrote {OUT}: {a.shape}, {X.nnz/1e6:.0f}M nnz", flush=True)
print(a.obs["Subclass"].value_counts().to_string())
