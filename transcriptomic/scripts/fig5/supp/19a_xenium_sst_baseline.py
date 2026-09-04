#!/usr/bin/env python
"""Baseline SST expression per Sst supertype from Xenium CONTROL cells.

Gate: qc_pass & corr_qc_pass & spatial_domain == 'Cortical' (the paper's
analysis set), diagnosis == Control, corr_supertype in the 16 Sst supertypes.
Per cell: SST CP10K = SST UMIs / total UMIs * 1e4. Aggregation: mean per donor
per supertype, then mean across the 12 control donors (guards against donor
imbalance). Output: transcriptomic/results/sst_strata_gsea/xenium_sst_baseline_controls.csv
"""
import time
import anndata as ad
import numpy as np
import pandas as pd

t0 = time.time()
H5 = "spatial/output/all_samples_annotated.h5ad"
OUT = "transcriptomic/results/sst_strata_gsea/xenium_sst_baseline_controls.csv"

print(f"[{time.time()-t0:6.1f}s] loading {H5} ...", flush=True)
a = ad.read_h5ad(H5)
print(f"[{time.time()-t0:6.1f}s] {a.shape[0]:,} cells x {a.shape[1]} genes", flush=True)

obs = a.obs
gate = (obs["qc_pass"] & obs["corr_qc_pass"]
        & (obs["spatial_domain"] == "Cortical")
        & (obs["diagnosis"] == "Control")
        & obs["corr_supertype"].astype(str).str.startswith("Sst_"))
sub = a[gate.values]
print(f"[{time.time()-t0:6.1f}s] control cortical Sst cells: {sub.shape[0]:,} "
      f"({sub.obs['sample_id'].nunique()} donors)", flush=True)

sst_idx = list(sub.var_names).index("SST")
sst = np.asarray(sub.X[:, sst_idx].todense()).ravel()
total = np.asarray(sub.X.sum(axis=1)).ravel()
cp10k = sst / total * 1e4

df = pd.DataFrame({
    "sample_id": sub.obs["sample_id"].values,
    "supertype": sub.obs["corr_supertype"].astype(str).values,
    "sst_cp10k": cp10k,
})
per_donor = (df.groupby(["supertype", "sample_id"], observed=True)
               .agg(mean_cp10k=("sst_cp10k", "mean"), n_cells=("sst_cp10k", "size"))
               .reset_index())
out = (per_donor.groupby("supertype", observed=True)
       .agg(xen_sst_cp10k=("mean_cp10k", "mean"),
            n_donors=("sample_id", "nunique"),
            n_cells=("n_cells", "sum"))
       .reset_index())
out.to_csv(OUT, index=False)
print(out.to_string(index=False))
print(f"[{time.time()-t0:6.1f}s] wrote {OUT}", flush=True)
