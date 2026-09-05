#!/usr/bin/env python3
"""
Run the DE step from the cached pseudobulk written by 05.

Same model, no re-load of the 9 GB query -- use this to iterate on the design without
paying for the pseudobulk rebuild.
"""
import os

import brisc_setup as B
import numpy as np
import pandas as pd
import polars as pl
from brisc import Pseudobulk

DATA = "/Users/shreejoy/Github/shared_data/GSE158516"
OUTD = "/Users/shreejoy/Github/scz_celltype_paper/reserve/gse158516/output"
FORMULA = "~ diagnosis + age + PMI + log2(num_cells) + log2(library_size)"
COEF = "diagnosisSchizophrenia"

pb = Pseudobulk(f"{DATA}/pseudobulk_subclass")
print(f"[05b] brisc {'.'.join(map(str, B.BRISC_VERSION))}; "
      f"{len(pb.keys())} cell types; DE: {FORMULA}", flush=True)

de = B.de(pb, FORMULA, coefficient=COEF, allow_float=True, verbose=False, strict=True)
os.makedirs(f"{DATA}/de", exist_ok=True)
de.save(f"{DATA}/de/gse158516_subclass_de", overwrite=True)
tab = de.table.to_pandas()
tab.to_csv(f"{OUTD}/de_subclass.csv.gz", index=False, compression="gzip")
print(f"[05b] wrote de_subclass.csv.gz {tab.shape}")

print("\n=== hits per cell type (FDR < 0.05) ===")
print(de.get_num_hits(significance_column="FDR", threshold=0.05).to_pandas().to_string(index=False))

print("\n=== canonical markers ===")
for gene, ct in [("SST", "Sst"), ("PVALB", "Pvalb"), ("VIP", "Vip"),
                 ("GAD1", "Sst"), ("GAD2", "Sst"), ("SST", "Pvalb")]:
    r = tab[(tab["gene"] == gene) & (tab["cell_type"] == ct)]
    if len(r):
        r = r.iloc[0]
        print(f"{gene:6s} in {ct:6s}: logFC={r['logFC']:+.3f}  P={r['p']:.4g}  FDR={r['FDR']:.4g}")
    else:
        print(f"{gene:6s} in {ct:6s}: not tested")
