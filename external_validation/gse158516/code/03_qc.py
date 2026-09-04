#!/usr/bin/env python3
"""
Per-sample QC for GSE158516 with brisc.

Run per sample rather than on the concatenated object for two reasons: the full query
is ~11.7 GB of non-zeros on a 24 GB machine, and doublets cannot span samples anyway,
so per-sample is the correct scope for doublet detection regardless.

brisc qc() defaults are the snRNA-seq conventions we want (max_mito_fraction=0.05,
nonzero_MALAT1=True -- MALAT1 is ubiquitous in nuclei, so its absence flags an empty
droplet or a cytoplasmic fragment). We raise min_genes from brisc's default 100 to 500,
which is standard for nuclei and consistent across the cohort's depth range.

Writes filtered per-sample h5ads via anndata rather than brisc's own save(), because
brisc 0.1.0 cannot read back what it writes (BRISC_BUGS.md Bug 8).
"""
import os
import sys
import time

import brisc_setup as B
import numpy as np
import pandas as pd
import polars as pl

IN = "/Users/shreejoy/Github/shared_data/GSE158516/per_sample"
OUT = "/Users/shreejoy/Github/shared_data/GSE158516/per_sample_qc"
OUTD = "/Users/shreejoy/Github/scz_celltype_paper/external_validation/gse158516/output"
META = "/Users/shreejoy/Github/scz_celltype_paper/external_validation/gse158516/data/sample_metadata.csv"

MIN_GENES = 500
MAX_MITO = 0.05


def main():
    os.makedirs(OUT, exist_ok=True)
    meta = pd.read_csv(META).set_index("Deidentified ID")
    samples = sorted(f[:-5] for f in os.listdir(IN) if f.endswith(".h5ad"))

    rows = []
    for sid in samples:
        t0 = time.time()
        sc = B.SingleCell(f"{IN}/{sid}.h5ad", num_threads=B.NUM_THREADS)
        B.verify_load(sc, sid)
        n0 = sc.shape[0]

        # qc_metrics first so we can report what each filter cost, then qc() to apply.
        sc = sc.qc_metrics(allow_float=True)
        mito = sc.obs["mito_fraction"].to_numpy()
        ngene = sc.obs["num_genes"].to_numpy()

        sc = sc.qc(remove_doublets=True, min_genes=MIN_GENES,
                   max_mito_fraction=MAX_MITO, subset=False, verbose=False,
                   allow_float=True)
        passed = sc.obs["passed_QC"].to_numpy()

        sub = sc.filter_obs("passed_QC")
        adata = sub.copy().to_scanpy()
        adata.obs["sample_id"] = sid
        adata.obs["diagnosis"] = (meta.loc[sid, "Diagnosis"] if sid in meta.index
                                  else "not_in_supplement")
        adata.write_h5ad(f"{OUT}/{sid}.h5ad", compression="gzip")

        rows.append(dict(
            sample=sid, in_paper=sid in meta.index, n_called=n0,
            n_pass=int(passed.sum()), frac_removed=1 - passed.mean(),
            frac_fail_genes=float((ngene < MIN_GENES).mean()),
            frac_fail_mito=float((mito >= MAX_MITO).mean()),
            med_genes=float(np.median(ngene[passed])),
            med_mito=float(np.median(mito[passed])),
            secs=round(time.time() - t0, 1)))
        print(f"{sid}: {n0:6d} -> {int(passed.sum()):6d} "
              f"({rows[-1]['frac_removed']:5.1%} removed; genes<{MIN_GENES} "
              f"{rows[-1]['frac_fail_genes']:5.1%}, mito {rows[-1]['frac_fail_mito']:5.1%})"
              f"  {rows[-1]['secs']}s", flush=True)
        del sc, sub, adata

    df = pd.DataFrame(rows)
    df.to_csv(f"{OUTD}/qc_per_sample.csv", index=False)
    pd.set_option("display.width", 200)
    print("\n" + df.round(4).to_string(index=False))
    print(f"\ntotal nuclei after QC: {df.n_pass.sum():,} "
          f"({df[df.in_paper].n_pass.sum():,} in the authors' 26 samples)")


if __name__ == "__main__":
    main()
