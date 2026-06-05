#!/usr/bin/env python3
"""
Build pseudobulk DE input for the Xenium SCZ-vs-Control edgeR analysis.

Mirrors build_crumblr_input.py and uses the SAME cell inclusion as the
compositional (crumblr) analysis: cortical cells passing corr_qc_pass, typed by
corr_subclass / corr_supertype (via load_sample_adata, qc_mode='corr'), across
all 24 donors. For each (sample x cell type) with >= MIN_CELLS cells, raw counts
are summed across cells into a genes x donors pseudobulk matrix.

Built at BOTH subclass and supertype level.

Outputs (consumed by run_de.R):
  output/de/pseudobulk_{level}.csv          long: celltype, donor, gene, count
  output/de/pseudobulk_{level}_samples.csv  celltype, donor, n_cells, total,
                                            diagnosis, sex, age, pmi, cell_class
"""
import os
import sys
import glob
import time
import argparse
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import BASE_DIR, H5AD_DIR, EXCLUDE_SAMPLES, load_sample_adata
sys.path.insert(0, os.path.join(BASE_DIR, "code"))
from modules.metadata import get_subject_info

METADATA_PATH = os.path.join(BASE_DIR, "data", "sample_metadata.xlsx")
OUT_DIR = os.path.join(BASE_DIR, "output", "de")
MIN_CELLS = 10                       # min cells per (sample x cell type) for a pseudobulk
LEVELS = [("subclass", "subclass_label"), ("supertype", "supertype_label")]


def main():
    ap = argparse.ArgumentParser(description="Build pseudobulk DE input CSVs")
    ap.add_argument("--qc-mode", default="corr", choices=["corr", "hybrid"],
                    help="QC mode (default 'corr' = corr_qc_pass, matching crumblr)")
    args = ap.parse_args()
    t0 = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"QC mode: {args.qc_mode} | min cells per pseudobulk: {MIN_CELLS}")

    meta = get_subject_info(METADATA_PATH).set_index("sample_id")

    h5ads = sorted(glob.glob(os.path.join(H5AD_DIR, "*_annotated.h5ad")))
    # acc: level -> celltype -> donor -> dict(pb=int64[genes], n_cells, total, cell_class)
    acc = {lvl: {} for lvl, _ in LEVELS}
    genes = None

    for fp in h5ads:
        sid = os.path.basename(fp).replace("_annotated.h5ad", "")
        if sid in EXCLUDE_SAMPLES:
            print(f"  skip {sid} (excluded)")
            continue
        a = load_sample_adata(sid, cortical_only=True, qc_mode=args.qc_mode)
        if genes is None:
            genes = list(a.var_names)
        X = a.X.tocsr() if sp.issparse(a.X) else np.asarray(a.X)
        cls = (a.obs["corr_class"].astype(str).values if "corr_class" in a.obs
               else np.array(["NA"] * a.n_obs))
        for lvl, col in LEVELS:
            labels = a.obs[col].astype(str).values
            for ct in np.unique(labels):
                m = labels == ct
                n = int(m.sum())
                if n < MIN_CELLS:
                    continue
                pb = np.asarray(X[m].sum(axis=0)).ravel().astype(np.int64)
                ccls = pd.Series(cls[m]).mode()
                acc[lvl].setdefault(ct, {})[sid] = dict(
                    pb=pb, n=n, total=int(pb.sum()),
                    cell_class=(ccls.iat[0] if len(ccls) else "NA"))
        print(f"  {sid}: {a.n_obs:,} cortical corr-QC cells")
        del a, X

    G = len(genes)
    for lvl, _ in LEVELS:
        cts, donors, mats, samp_rows = [], [], [], []
        for ct, dct in acc[lvl].items():
            for donor, rec in dct.items():
                cts.append(ct); donors.append(donor); mats.append(rec["pb"])
                samp_rows.append(dict(
                    celltype=ct, donor=donor, n_cells=rec["n"], total=rec["total"],
                    diagnosis=meta.loc[donor, "diagnosis"], sex=meta.loc[donor, "sex"],
                    age=float(meta.loc[donor, "age"]),
                    pmi=(float(meta.loc[donor, "pmi"]) if "pmi" in meta.columns else np.nan),
                    cell_class=rec["cell_class"]))
        M = np.vstack(mats)                                  # n_pb x genes
        long = pd.DataFrame({
            "celltype": np.repeat(cts, G),
            "donor": np.repeat(donors, G),
            "gene": np.tile(genes, len(cts)),
            "count": M.ravel()})
        long.to_csv(f"{OUT_DIR}/pseudobulk_{lvl}.csv", index=False)
        pd.DataFrame(samp_rows).to_csv(f"{OUT_DIR}/pseudobulk_{lvl}_samples.csv", index=False)
        n_donor = pd.DataFrame(samp_rows).groupby("celltype")["donor"].nunique()
        print(f"\n{lvl}: {len(acc[lvl])} cell types, {len(cts)} pseudobulk samples "
              f"(donors/type: median {int(n_donor.median())}, min {int(n_donor.min())})")
        print(f"  -> pseudobulk_{lvl}.csv ({len(long):,} rows) + pseudobulk_{lvl}_samples.csv")

    print(f"\nDone in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
