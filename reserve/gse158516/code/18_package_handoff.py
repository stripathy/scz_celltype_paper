#!/usr/bin/env python3
"""
Package GSE158516 as a label-transfer *input*, not an annotated dataset.

If this cohort is ever used in the paper it must be annotated by the same Seurat
reference-based label transfer as the other seven, so the provisional labels produced
here must not travel inside the object anyone will load. Leaving a column called
`supertype` in `obs` is an invitation to use it by accident.

So: strip every provisional annotation out of the h5ad into a clearly-named sidecar, and
rename the object to say what it actually is. What remains is counts + donor metadata +
QC flags -- exactly what the label-transfer step consumes.

  GSE158516_annotated.h5ad          ->  GSE158516_counts_qc.h5ad
  (provisional labels, confidences, cluster, UMAP, Harmony)
                                    ->  GSE158516_provisional_labels.parquet
                                        GSE158516_provisional_harmony.npy

The h5ad is edited in place with h5py rather than rebuilt: a rebuild means re-staging
~8 GB through concat_on_disk for what is a metadata change. HDF5 does not reclaim the
freed space, which is irrelevant here since the file is dominated by X.
"""
import os
import shutil

import anndata as ad
import h5py
import numpy as np
import pandas as pd

DATA = "/Users/shreejoy/Github/shared_data/GSE158516"
SRC = f"{DATA}/GSE158516_annotated.h5ad"
DST = f"{DATA}/GSE158516_counts_qc.h5ad"
LABELS = f"{DATA}/GSE158516_provisional_labels.parquet"
HARMONY = f"{DATA}/GSE158516_provisional_harmony.npy"

PROVISIONAL = ["class", "subclass", "supertype",
               "class_confidence", "subclass_confidence", "supertype_confidence",
               "subclass_next", "subclass_next_confidence",
               "supertype_next", "supertype_next_confidence",
               "cluster", "UMAP1", "UMAP2"]

if not os.path.exists(SRC):
    raise SystemExit(f"{SRC} not found (already packaged?)")

# ---- 1. lift the provisional annotation out into a sidecar --------------------
a = ad.read_h5ad(SRC, backed="r")
present = [c for c in PROVISIONAL if c in a.obs.columns]
side = a.obs[present].copy()
side.index.name = "cell_id"
side.to_parquet(LABELS)
print(f"wrote {LABELS}  {side.shape}")

harmony = a.obsm["X_harmony"][:] if "X_harmony" in a.obsm else None
cell_ids = a.obs_names.to_numpy()
a.file.close()
if harmony is not None:
    np.save(HARMONY, harmony)
    print(f"wrote {HARMONY}  {harmony.shape}")
np.save(f"{DATA}/GSE158516_provisional_cell_ids.npy", cell_ids.astype(object),
        allow_pickle=True)

# ---- 2. strip them from the object -------------------------------------------
shutil.move(SRC, DST)
with h5py.File(DST, "a") as f:
    order = [c.decode() if isinstance(c, bytes) else c
             for c in f["obs"].attrs["column-order"]]
    for col in present:
        if col in f["obs"]:
            del f["obs"][col]
    f["obs"].attrs["column-order"] = np.array(
        [c for c in order if c not in present], dtype=h5py.special_dtype(vlen=str))
    if "obsm" in f and "X_harmony" in f["obsm"]:
        del f["obsm/X_harmony"]

# ---- 3. verify ----------------------------------------------------------------
b = ad.read_h5ad(DST, backed="r")
leaked = [c for c in PROVISIONAL if c in b.obs.columns]
assert not leaked, f"provisional columns still present: {leaked}"
assert b.n_obs == len(side), f"cell count changed: {b.n_obs} vs {len(side)}"
assert list(b.obs_names[:5]) == list(side.index[:5]), "cell order changed"
print(f"\n{DST}")
print(b)
print(f"\nsize on disk: {os.path.getsize(DST) / 1e9:.2f} GB")
print("\nobs now carries only counts-level QC and donor metadata:")
print("  " + ", ".join(b.obs.columns))
print("\nprovisional annotation lives in:")
print(f"  {os.path.basename(LABELS)}  (per-cell labels + confidences + UMAP)")
print(f"  {os.path.basename(HARMONY)}  (Harmony embedding, cell order matches the h5ad)")
