#!/usr/bin/env python
"""Sst-strata pseudobulks from SEA-AD DFC (A9), for DE along CPS.

Applies OUR SCZ depletion strata to the SEA-AD DFC snRNA-seq release
(~/Downloads/SEAAD_DFC_RNAseq_final-nuclei.2026-06-22.h5ad, 1.36M nuclei) so
the AD state axis can be compared to the SCZ one with the same pipeline.
Counts from layers/UMIs (raw). Donors joined to the committed CPS lookup
(crossdisorder/data/seaad_donor_cps.csv, 80 donors); the severely-affected
flag is carried in the metadata so the DE step can exclude those donors, as
Gabitto et al. did for gene-expression tests.

Output -> transcriptomic/results/sst_strata_gsea/seaad_a9/
"""
import os, time
import h5py
import numpy as np
import pandas as pd

t0 = time.time()
H5  = "/Users/shreejoy/Downloads/SEAAD_DFC_RNAseq_final-nuclei.2026-06-22.h5ad"
CPS = "crossdisorder/data/seaad_donor_cps.csv"
OUT = "transcriptomic/results/sst_strata_gsea/seaad_a9"
os.makedirs(OUT, exist_ok=True)

STRATA = {
    "depleted":     ["Sst_2", "Sst_3", "Sst_20", "Sst_22", "Sst_25"],
    "intermediate": ["Sst_9", "Sst_11", "Sst_12", "Sst_13", "Sst_19", "Sst_23"],
    "non_depleted": ["Sst_1", "Sst_4", "Sst_5", "Sst_7", "Sst_10"],
}
st_of = {s: k for k, v in STRATA.items() for s in v}
say = lambda *a: print(f"[{time.time()-t0:7.1f}s]", *a, flush=True)

def cat_col(obs, name):
    g = obs[name]
    if isinstance(g, h5py.Group):
        cats = np.array([c.decode() if isinstance(c, bytes) else c
                         for c in g["categories"][:]])
        return cats[g["codes"][:]]
    v = g[:]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in v]) \
        if v.dtype.kind in "OS" else v

with h5py.File(H5, "r") as f:
    obs = f["obs"]
    sup = cat_col(obs, "Supertype")
    donor = cat_col(obs, "Donor ID")
    say("obs loaded")

    cps = pd.read_csv(CPS)
    cps_map = dict(zip(cps.donor, cps.CPS))
    is_sst = np.isin(sup, list(st_of))
    has_cps = np.isin(donor, list(cps_map))
    sel = np.flatnonzero(is_sst & has_cps)
    say(f"Sst cells with CPS donor: {len(sel):,} "
        f"({len(np.unique(donor[sel]))} donors)")

    strat = np.array([st_of[s] for s in sup[sel]])
    grp = np.char.add(np.char.add(donor[sel].astype(str), "|"), strat)
    keys, kinv = np.unique(grp, return_inverse=True)
    say(f"{len(keys)} donor-stratum groups")

    var = f["var"]
    vidx = var.attrs.get("_index", "_index")
    vidx = vidx.decode() if isinstance(vidx, bytes) else vidx
    genes = np.array([g.decode() if isinstance(g, bytes) else g
                      for g in var[vidx][:]])
    n_genes = len(genes)

    U = f["layers"]["UMIs"]
    indptr = U["indptr"][:]
    data, indices = U["data"], U["indices"]
    pb = np.zeros((n_genes, len(keys)), dtype=np.float64)
    order = np.argsort(sel)
    for j, (r, g) in enumerate(zip(sel[order], kinv[order])):
        s, e = indptr[r], indptr[r + 1]
        np.add.at(pb[:, g], indices[s:e], data[s:e])
        if j % 10000 == 0:
            say(f"  {j:,}/{len(sel):,} cells")
    say("pseudobulk accumulated")

meta = pd.DataFrame({"key": grp}).groupby("key").size().rename("n_cells").reset_index()
info = pd.DataFrame({
    "key": grp,
    "donor": donor[sel].astype(str),
    "stratum": strat,
})
with h5py.File(H5, "r") as f:
    obs = f["obs"]
    aux = pd.DataFrame({
        "donor": donor, "sex": cat_col(obs, "Sex"),
        "age": cat_col(obs, "Age at Death"),
        "pmi": cat_col(obs, "PMI"),
        "severe": cat_col(obs, "Severely Affected Donor"),
    }).drop_duplicates("donor")
meta = (meta.merge(info.drop_duplicates("key"), on="key")
        .merge(aux, on="donor", how="left"))
meta["CPS"] = meta.donor.map(cps_map)
meta = meta.set_index("key").loc[keys].reset_index()
meta.to_csv(f"{OUT}/a9_stratum_pseudobulk_meta.csv", index=False)

out = pd.DataFrame(pb, index=genes, columns=keys).groupby(level=0).sum()
out.to_parquet(f"{OUT}/a9_stratum_pseudobulk.parquet")
say(f"wrote {out.shape[0]:,} genes x {out.shape[1]} groups")
print(meta.groupby(["stratum"], observed=True)["n_cells"].sum().to_string())
print("severely affected donors:", meta.loc[meta.severe.astype(str).isin(['True','TRUE','1']), 'donor'].nunique())
