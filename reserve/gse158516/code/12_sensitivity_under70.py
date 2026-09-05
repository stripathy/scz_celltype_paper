#!/usr/bin/env python3
"""
Sensitivity analysis: restrict GSE158516 to donors aged <= 70, matching the paper's
own inclusion rule ("to limit survivorship bias, only cases and controls younger than
70 years at death were included").

This is not just a compliance check. In this cohort the rule drops 5 donors, all of
them controls (ages 71, 73, 76, 80, 87), and in doing so removes the cohort's main
confound: controls are 12.8 y older than cases across all 26 donors (Welch p = 0.064),
versus a 0.0 y difference among the 21 who are <= 70 (p = 1.000). The cost is power --
9 controls vs 12 cases.

Reruns the DE from the cached pseudobulk and writes age-restricted counts for crumblr.
"""
import os

import brisc_setup as B
import numpy as np
import pandas as pd
import polars as pl
from brisc import Pseudobulk

DATA = "/Users/shreejoy/Github/shared_data/GSE158516"
BASE = "/Users/shreejoy/Github/scz_celltype_paper/external_validation/gse158516"
OUTD = f"{BASE}/output"
FORMULA = "~ diagnosis + age + PMI + log2(num_cells) + log2(library_size)"
COEF = "diagnosisSchizophrenia"
MAX_AGE = 70

meta = pd.read_csv(f"{BASE}/data/sample_metadata.csv")
keep = set(meta.loc[meta["Age"] <= MAX_AGE, "Deidentified ID"])
print(f"[12] donors <= {MAX_AGE}: {len(keep)} of {len(meta)} "
      f"({(meta['Age'] <= MAX_AGE).groupby(meta['Diagnosis']).sum().to_dict()})", flush=True)

# ---- age-restricted counts for crumblr --------------------------------------
for level in ["subclass", "supertype", "class"]:
    c = pd.read_csv(f"{OUTD}/counts_{level}.csv", index_col=0)
    c.loc[[i for i in c.index if i in keep]].to_csv(f"{OUTD}/under70_counts_{level}.csv")
print("[12] wrote under70_counts_*.csv", flush=True)

# ---- DE on the restricted donor set -----------------------------------------
pb = Pseudobulk(f"{DATA}/pseudobulk_subclass")
pb = pb.filter_obs(pl.col("sample_id").is_in(list(keep)))
n_per = {ct: obs.height for ct, (X, obs, var) in pb.items()}
print(f"[12] samples per cell type after restriction: "
      f"min={min(n_per.values())}, max={max(n_per.values())}", flush=True)

# re-run pseudobulk QC on the restricted set: the >=80%-of-samples gene filter and the
# outlier check both depend on which donors are present
pb = pb.qc("diagnosis", min_cells=10, allow_float=True, verbose=False)
pb = pb.library_size(allow_float=True, overwrite=True)
print(f"[12] cell types surviving: {len(pb.keys())}", flush=True)

de = B.de(pb, FORMULA, coefficient=COEF, allow_float=True, verbose=False, strict=True)
tab = de.table.to_pandas()
tab.to_csv(f"{OUTD}/under70_de_subclass.csv.gz", index=False, compression="gzip")
print(f"[12] wrote under70_de_subclass.csv.gz {tab.shape}")

full = pd.read_csv(f"{OUTD}/de_subclass.csv.gz")
print("\n=== canonical markers: all 26 donors vs <= 70 only ===")
for gene, ct in [("SST", "Sst"), ("PVALB", "Pvalb"), ("VIP", "Vip"),
                 ("GAD1", "Sst"), ("GAD2", "Sst")]:
    a = full[(full.gene == gene) & (full.cell_type == ct)]
    b = tab[(tab.gene == gene) & (tab.cell_type == ct)]
    fa = f"{a.logFC.iloc[0]:+.3f} (p={a.p.iloc[0]:.3f})" if len(a) else "not tested"
    fb = f"{b.logFC.iloc[0]:+.3f} (p={b.p.iloc[0]:.3f})" if len(b) else "not tested"
    print(f"{gene:6s} in {ct:6s}:  all 26 {fa:24s}   <=70 {fb}")

print("\n=== DE hits at FDR < 0.05 ===")
print(de.get_num_hits(significance_column="FDR", threshold=0.05).to_pandas().to_string(index=False))
