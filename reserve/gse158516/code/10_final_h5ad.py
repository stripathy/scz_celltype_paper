#!/usr/bin/env python3
"""
Assemble the concatenated object.

NOTE: run `18_package_handoff.py` after this. It lifts the provisional labels out into a
sidecar and renames the result `GSE158516_counts_qc.h5ad`, so that the object anyone
loads cannot be mistaken for a properly annotated dataset. This script deliberately still
writes the labels in, so that 18 has a single place to strip them from.

One object, 32 samples, raw counts, with donor metadata, QC metrics, SEA-AD
Class/Subclass/Supertype labels + transfer confidences, and the UMAP/Harmony
embeddings attached.

Built with concat_on_disk because the concatenated matrix is ~8 GB of non-zeros --
more than we want resident on a 24 GB machine alongside anything else.
"""
import os
import shutil

import anndata as ad
import numpy as np
import pandas as pd

DATA = "/Users/shreejoy/Github/shared_data/GSE158516"
QC = f"{DATA}/per_sample_qc"
STAGE = f"{DATA}/per_sample_final"
OUT = f"{DATA}/GSE158516_annotated.h5ad"
META = "/Users/shreejoy/Github/scz_celltype_paper/external_validation/gse158516/data/sample_metadata.csv"

LABEL_COLS = ["class", "subclass", "supertype", "class_confidence",
              "subclass_confidence", "supertype_confidence",
              "subclass_next", "subclass_next_confidence",
              "supertype_next", "supertype_next_confidence"]

obs_all = pd.read_parquet(f"{DATA}/annotated/query_obs.parquet").set_index("_index")
if "cluster" in obs_all.columns:
    LABEL_COLS.append("cluster")
umap = np.load(f"{DATA}/annotated/query_umap.npy")
harmony = np.load(f"{DATA}/annotated/query_harmony.npy")
obs_all["UMAP1"], obs_all["UMAP2"] = umap[:, 0], umap[:, 1]

meta = pd.read_csv(META).rename(columns={
    "Deidentified ID": "sample_id", "Diagnosis": "diagnosis", "Age": "age",
    "Sex": "sex", "Sequencing Batch": "seq_batch", "PFC_pH": "pH",
    "Number of Nuclei": "authors_n_nuclei"}).set_index("sample_id")

os.makedirs(STAGE, exist_ok=True)
paths = []
for f in sorted(os.listdir(QC)):
    if not f.endswith(".h5ad"):
        continue
    sid = f[:-5]
    a = ad.read_h5ad(f"{QC}/{f}")
    a.obs = a.obs.join(obs_all[LABEL_COLS + ["UMAP1", "UMAP2"]], how="left")
    a.obs["in_authors_analysis"] = sid in meta.index
    # string columns need a string sentinel, not NaN: a mixed str/float object column
    # is not writable to h5ad
    for col, blank in [("diagnosis", "not_in_supplement"), ("sex", "unknown"),
                       ("age", np.nan), ("PMI", np.nan), ("seq_batch", np.nan),
                       ("pH", np.nan), ("authors_n_nuclei", np.nan)]:
        a.obs[col] = meta.loc[sid, col] if sid in meta.index else blank
    a.obs["sample_id"] = sid
    a.uns = {}                    # per-sample cell_call records differ; would block concat
    out = f"{STAGE}/{f}"
    a.write_h5ad(out, compression="gzip")
    paths.append(out)
    print(f"staged {sid}: {a.n_obs} cells", flush=True)
    del a

print("concatenating on disk...", flush=True)
if os.path.exists(OUT):
    os.remove(OUT)
ad.experimental.concat_on_disk(paths, OUT, label=None, index_unique=None)

# harmony embedding is per-cell across the whole cohort; attach after concat
final = ad.read_h5ad(OUT, backed="r")
order = pd.Index(final.obs_names)
h = harmony[obs_all.index.get_indexer(order)]
final.file.close()
import h5py
with h5py.File(OUT, "a") as f:
    if "obsm" not in f:
        f.create_group("obsm")
        f["obsm"].attrs["encoding-type"] = "dict"
        f["obsm"].attrs["encoding-version"] = "0.1.0"
    if "X_harmony" in f["obsm"]:
        del f["obsm/X_harmony"]
    d = f["obsm"].create_dataset("X_harmony", data=h.astype(np.float32))
    d.attrs["encoding-type"] = "array"
    d.attrs["encoding-version"] = "0.2.0"

a = ad.read_h5ad(OUT, backed="r")
print(f"\nwrote {OUT}")
print(a)
print(f"\nsize on disk: {os.path.getsize(OUT)/1e9:.2f} GB")
shutil.rmtree(STAGE)
print("cleaned staging dir")
