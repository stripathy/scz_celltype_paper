#!/usr/bin/env python
"""Independent-cohort check: stratum pseudobulks from the Jens snRNA-seq dataset.

The Jens dataset (83 donors, 42 CTRL / 41 SCZ; NOT one of the 7 meta-analysis
cohorts, and not the Parse cohort) was brisc-mapped onto the SEA-AD DFC
supertype taxonomy in scz_parse_taxonomy. Counts live in ext_full_jens.h5ad
(cells x genes, obs_names == the label parquet's index); labels + donor
metadata in full_mapped_jens.parquet.

Builds per-donor pseudobulk counts for the three Sst depletion strata, in two
variants: all Sst-labeled cells (primary) and supertype_confidence >= 0.5
(sensitivity). Also audits stratum-level label ambiguity (does the runner-up
supertype fall in the same stratum?) because 4/5 depleted-stratum supertypes
are flagged CONFIDENTLY_WRONG in this cohort's mapping assessment.

Output -> transcriptomic/results/sst_strata_gsea/jens/
"""
import os, time
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

t0 = time.time()
H5   = "/Users/shreejoy/Github/scz_parse_taxonomy/data/ext_full_jens.h5ad"
PARQ = "/Users/shreejoy/Github/scz_parse_taxonomy/output/full_mapped_jens.parquet"
OUT  = "transcriptomic/results/sst_strata_gsea/jens"
os.makedirs(OUT, exist_ok=True)

STRATA = {
    "depleted":     ["Sst_2", "Sst_3", "Sst_20", "Sst_22", "Sst_25"],
    "intermediate": ["Sst_9", "Sst_11", "Sst_12", "Sst_13", "Sst_19", "Sst_23"],
    "non_depleted": ["Sst_1", "Sst_4", "Sst_5", "Sst_7", "Sst_10"],
}
st_of = {s: k for k, v in STRATA.items() for s in v}
say = lambda *a: print(f"[{time.time()-t0:6.1f}s]", *a, flush=True)

lab = pd.read_parquet(PARQ)
sst = lab[lab["supertype"].isin(st_of)].copy()
sst["stratum"] = sst["supertype"].map(st_of)
say(f"Sst cells: {len(sst):,} | donors {sst.donor_id.nunique()} | "
    f"dx {sst.diagnosis.value_counts().to_dict()}")

# ---- stratum-level ambiguity audit ------------------------------------------
nxt = sst["supertype_next"].map(st_of)          # NaN if runner-up not Sst
tab = (pd.DataFrame({"stratum": sst["stratum"], "next": nxt.fillna("non-Sst")})
       .groupby(["stratum", "next"], observed=True).size().unstack(fill_value=0))
say("runner-up stratum by assigned stratum (row-normalized):")
print((tab.T / tab.sum(axis=1)).T.round(3).to_string())
say("mean supertype confidence by stratum:")
print(sst.groupby("stratum", observed=True)["supertype_confidence"].mean().round(3).to_string())

# ---- counts ------------------------------------------------------------------
a = ad.read_h5ad(H5, backed="r")
say(f"h5ad: {a.shape} (backed)")
common = sst.index.intersection(a.obs_names)
assert len(common) == len(sst), f"index mismatch: {len(common)} of {len(sst)}"
# verify metadata agreement on a sample
pos = a.obs_names.get_indexer(sst.index[:2000])
agree = (a.obs["donor_id"].values[pos].astype(str) ==
         sst["donor_id"].astype(str).values[:2000]).mean()
say(f"alignment: donor agreement on 2000 cells = {agree:.4f}")
assert agree > 0.999

sub = a[sst.index].to_memory()
X = sparse.csr_matrix(sub.X)
say(f"Sst count matrix: {X.shape}, nnz {X.nnz:,}")

def pseudobulk(mask, tag):
    s = sst[mask]
    grp = s["donor_id"].astype(str) + "|" + s["stratum"]
    keys = sorted(grp.unique())
    kidx = grp.map({k: i for i, k in enumerate(keys)}).values
    onehot = sparse.csr_matrix(
        (np.ones(len(s)), (np.arange(len(s)), kidx)),
        shape=(len(s), len(keys)))
    pb = (X[np.flatnonzero(mask.values)].T @ onehot).toarray()
    dfpb = pd.DataFrame(pb, index=sub.var_names, columns=keys)
    dfpb = dfpb.groupby(level=0).sum()
    dfpb.to_parquet(f"{OUT}/jens_stratum_pseudobulk{tag}.parquet")
    meta = (s.groupby(["donor_id", "stratum"], observed=True)
            .agg(n_cells=("supertype", "size"), diagnosis=("diagnosis", "first"),
                 sex=("sex", "first"), age=("age", "first"), pmi=("pmi", "first"))
            .reset_index())
    meta["key"] = meta["donor_id"].astype(str) + "|" + meta["stratum"]
    meta = meta.set_index("key").loc[keys].reset_index()
    meta.to_csv(f"{OUT}/jens_stratum_pseudobulk_meta{tag}.csv", index=False)
    say(f"{tag or 'primary'}: {dfpb.shape[0]:,} genes x {len(keys)} groups; "
        f"median cells/group {int(meta.n_cells.median())}")
    print(meta.groupby(["stratum", "diagnosis"], observed=True)["n_cells"]
          .sum().unstack().to_string())

pseudobulk(pd.Series(True, index=sst.index), "")
pseudobulk(sst["supertype_confidence"] >= 0.5, "_conf05")
say("done")
