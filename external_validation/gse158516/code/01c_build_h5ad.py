#!/usr/bin/env python3
"""
Pass 2 of the GSE158516 build: subset the raw matrices to called cells -> per-sample h5ad.

Cell-calling rule (chosen in 01b from the barcode-rank curves, see knee_plots.png):

  * 26 samples in the authors' supplement -> top-N barcodes by total UMI, with N the
    authors' own reported nuclei count. This reproduces their cell set almost exactly:
    the median UMI per nucleus we recover matches their published Supplementary Table 2
    value to <1% in 23/26 samples (max 3.2%). Using their N means our replication tests
    the biology rather than our cell-calling.
  * 6 samples the authors dropped (CT_03, CT_08, SZ_03... see `authors_qc_pass`) have no
    reported N and no age/PMI, so they get the assumption-free inflection-point rule and
    are flagged. They are built for completeness but cannot enter the covariate models.

CellRanger v2 "ordmag" was rejected: it keeps only 0.56x the authors' nuclei.
"""
import os
import sys
import time
from multiprocessing import Pool

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

RAW = "/Users/shreejoy/Github/shared_data/GSE158516/raw"
TOT = "/Users/shreejoy/Github/shared_data/GSE158516/barcode_totals"
OUT = "/Users/shreejoy/Github/shared_data/GSE158516/per_sample"
META = "/Users/shreejoy/Github/scz_celltype_paper/external_validation/gse158516/data/sample_metadata.csv"
CALLS = "/Users/shreejoy/Github/scz_celltype_paper/external_validation/gse158516/output/cell_calling_comparison.csv"
CHUNK = 20_000_000


def mtx_header(path):
    import gzip
    with gzip.open(path, "rt") as fh:
        n = 0
        for line in fh:
            n += 1
            if line.startswith("%"):
                continue
            g, b, nnz = (int(v) for v in line.split())
            return g, b, nnz, n
    raise ValueError(path)


def build_one(args):
    gsm, sid, n_target, rule = args
    out = f"{OUT}/{sid}.h5ad"
    if os.path.exists(out):
        return f"{sid}: exists, skipped"
    t0 = time.time()

    z = np.load(f"{TOT}/{sid}.npz", allow_pickle=True)
    bc_idx, umi, ngene = z["bc_idx"], z["umi"], z["ngene"]
    order = np.argsort(umi)[::-1]
    sel = order[:n_target]
    keep_idx = bc_idx[sel]                    # 1-based barcode indices in the raw mtx
    keep_umi, keep_ngene = umi[sel], ngene[sel]
    srt = np.argsort(keep_idx)                # store in barcode order for reproducibility
    keep_idx, keep_umi, keep_ngene = keep_idx[srt], keep_umi[srt], keep_ngene[srt]
    n_cells = keep_idx.size

    stem = f"{RAW}/{gsm}_{sid}"
    mtx = f"{stem}_matrix.mtx.gz"
    n_genes, n_bc, nnz, n_header = mtx_header(mtx)

    remap = np.full(n_bc + 1, -1, dtype=np.int32)
    remap[keep_idx] = np.arange(n_cells, dtype=np.int32)

    gi, ci, vi = [], [], []
    for ch in pd.read_csv(mtx, sep=" ", header=None, skiprows=n_header,
                          usecols=[0, 1, 2], names=["g", "b", "c"],
                          dtype={0: np.int32, 1: np.int64, 2: np.int32},
                          chunksize=CHUNK, engine="c", compression="gzip"):
        newcol = remap[ch["b"].to_numpy()]
        m = newcol >= 0
        if m.any():
            gi.append(ch["g"].to_numpy()[m] - 1)
            ci.append(newcol[m])
            vi.append(ch["c"].to_numpy()[m])
    gi, ci, vi = np.concatenate(gi), np.concatenate(ci), np.concatenate(vi)
    X = sp.coo_matrix((vi.astype(np.float32), (ci, gi)), shape=(n_cells, n_genes)).tocsr()

    barcodes = pd.read_csv(f"{stem}_barcodes.tsv.gz", header=None)[0].to_numpy()[keep_idx - 1]
    feat = pd.read_csv(f"{stem}_features.tsv.gz", sep="\t", header=None,
                       names=["gene_id", "gene_symbol", "feature_type"])

    adata = ad.AnnData(X=X)
    adata.obs_names = pd.Index([f"{sid}_{b}" for b in barcodes])
    adata.obs["sample_id"] = sid
    adata.obs["gsm"] = gsm
    adata.obs["barcode"] = barcodes
    adata.obs["n_umi"] = keep_umi
    adata.obs["n_genes"] = keep_ngene
    # .to_numpy() so the Index does not inherit the Series name -- an index named
    # "gene_symbol" collides with the (pre-uniquify, hence different) var column below.
    adata.var_names = pd.Index(feat["gene_symbol"].astype(str).to_numpy())
    adata.var_names_make_unique()
    adata.var["gene_id"] = feat["gene_id"].to_numpy()
    adata.var["gene_symbol_raw"] = feat["gene_symbol"].to_numpy()
    adata.uns["cell_call"] = {"rule": rule, "n_target": int(n_target),
                              "min_umi_kept": int(keep_umi.min()),
                              "n_barcodes_total": int(n_bc)}
    adata.write_h5ad(out, compression="gzip")
    return (f"{sid}: {n_cells:6d} cells  rule={rule:8s}  minUMI={int(keep_umi.min()):6d}  "
            f"medUMI={int(np.median(keep_umi)):6d}  {time.time()-t0:4.0f}s")


def main():
    os.makedirs(OUT, exist_ok=True)
    meta = pd.read_csv(META).set_index("Deidentified ID")
    calls = pd.read_csv(CALLS).set_index("sample")

    jobs = []
    for f in sorted(os.listdir(RAW)):
        if not f.endswith("_matrix.mtx.gz"):
            continue
        gsm, sid = f.split("_")[0], "_".join(f.split("_")[1:3])
        if sid in meta.index:
            jobs.append((gsm, sid, int(meta.loc[sid, "Number of Nuclei"]), "topN"))
        else:
            jobs.append((gsm, sid, int(calls.loc[sid, "n_inflect"]), "inflect"))
    print(f"{len(jobs)} samples to build", flush=True)

    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    with Pool(nproc) as pool:
        for msg in pool.imap_unordered(build_one, jobs):
            print(msg, flush=True)


if __name__ == "__main__":
    main()
