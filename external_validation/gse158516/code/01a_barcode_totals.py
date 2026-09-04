#!/usr/bin/env python3
"""
Pass 1 of the GSE158516 build: per-barcode UMI totals and gene counts.

GEO ships CellRanger *raw* matrices (33,538 features x 6,794,880 whitelist
barcodes), so cell calling is on us. Caching the barcode-rank curve as its own
artifact lets us look at the knee and choose a calling rule once, instead of
re-streaming 10 GB every time we want to change a threshold.

Only barcodes with >= 100 UMI are kept (the rest are unambiguously ambient and
would make the cache 100x bigger for no information).
"""
import gzip
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

RAW = "/Users/shreejoy/Github/shared_data/GSE158516/raw"
OUT = "/Users/shreejoy/Github/shared_data/GSE158516/barcode_totals"
CHUNK = 20_000_000
MIN_KEEP = 100


def mtx_header(path):
    with gzip.open(path, "rt") as fh:
        n = 0
        for line in fh:
            n += 1
            if line.startswith("%"):
                continue
            g, b, nnz = (int(v) for v in line.split())
            return g, b, nnz, n
    raise ValueError(path)


def totals(args):
    gsm, sid = args
    out = f"{OUT}/{sid}.npz"
    if os.path.exists(out):
        return f"{sid}: cached"
    t0 = time.time()
    mtx = f"{RAW}/{gsm}_{sid}_matrix.mtx.gz"
    n_genes, n_bc, nnz, n_header = mtx_header(mtx)

    umi = np.zeros(n_bc + 1, dtype=np.int64)
    ngene = np.zeros(n_bc + 1, dtype=np.int32)
    for ch in pd.read_csv(mtx, sep=" ", header=None, skiprows=n_header,
                          usecols=[1, 2], names=["b", "c"],
                          dtype={1: np.int64, 2: np.int64},
                          chunksize=CHUNK, engine="c", compression="gzip"):
        b = ch["b"].to_numpy()
        umi += np.bincount(b, weights=ch["c"].to_numpy(), minlength=n_bc + 1).astype(np.int64)
        ngene += np.bincount(b, minlength=n_bc + 1).astype(np.int32)

    idx = np.where(umi >= MIN_KEEP)[0]
    idx = idx[idx > 0]
    np.savez_compressed(out, bc_idx=idx.astype(np.int64),
                        umi=umi[idx].astype(np.int64), ngene=ngene[idx].astype(np.int32),
                        n_barcodes=n_bc, n_genes=n_genes, nnz=nnz)
    return f"{sid}: {idx.size:8d} barcodes >= {MIN_KEEP} UMI  ({time.time()-t0:4.0f}s)"


def main():
    os.makedirs(OUT, exist_ok=True)
    jobs = []
    for f in sorted(os.listdir(RAW)):
        if f.endswith("_matrix.mtx.gz"):
            jobs.append((f.split("_")[0], "_".join(f.split("_")[1:3])))
    print(f"{len(jobs)} samples", flush=True)
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 6) as pool:
        for msg in pool.imap_unordered(totals, jobs):
            print(msg, flush=True)


if __name__ == "__main__":
    main()
