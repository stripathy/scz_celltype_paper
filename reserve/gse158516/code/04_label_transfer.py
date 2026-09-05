#!/usr/bin/env python3
"""
Annotate GSE158516 nuclei by label transfer from the SEA-AD MTG reference (brisc).

This is the light-weight stand-in for the Seurat reference-based label transfer the
main paper's pipeline uses, and it puts this cohort in the *same* taxonomy (SEA-AD
Class/Subclass/Supertype) as our seven datasets, which is the whole point -- otherwise
there is nothing to compare against.

Memory discipline (24 GB box, query is ~7-9 GB of non-zeros post-QC):
  * the reference is subsampled by 02b (proportionally, so class priors survive);
  * both datasets are cut to the shared gene space before anything else;
  * normalize() runs on the full gene set (library sizes must be honest) and only then
    do we drop to HVGs, which cuts the working set ~4x before PCA/Harmony.

Also runs de novo clustering so the transferred labels can be checked against
unsupervised structure rather than trusted on faith.
"""
import gc
import os
import sys

import brisc_setup as B
import numpy as np
import pandas as pd
import polars as pl

QC = "/Users/shreejoy/Github/shared_data/GSE158516/per_sample_qc"
# Pre-subsampled by 02b: the full reference is 9.4 GB of int64 non-zeros in memory,
# which blows the budget before the query is even loaded. "proportional" preserves the
# reference's class balance, which kNN voting leaks into the calls -- the "capped"
# variant inflates GABAergic from 33% to 61% of reference neurons and is a sensitivity
# check only. Pass the mode as argv[1].
MODE = sys.argv[1] if len(sys.argv) > 1 else "proportional"
REF = f"/Users/shreejoy/Github/shared_data/GSE158516/reference/seaad_mtg_{MODE}.h5ad"
SUFFIX = "" if MODE == "proportional" else f"_{MODE}"
OUTD = "/Users/shreejoy/Github/scz_celltype_paper/external_validation/gse158516/output"
DATA = "/Users/shreejoy/Github/shared_data/GSE158516"

N_HVG = 3000
SEED = 0


def log(msg):
    print(f"[04] {msg}", flush=True)


# ---- reference ---------------------------------------------------------------
log(f"loading SEA-AD reference ({MODE}, pre-subsampled by 02b)")
ref = B.SingleCell(REF, num_threads=B.NUM_THREADS)
B.verify_load(ref, "reference")
ref = ref.skip_qc()
log(f"reference {ref.shape}, {ref.obs['Supertype'].n_unique()} supertypes, "
    f"{ref.obs['Subclass'].n_unique()} subclasses")

# ---- query -------------------------------------------------------------------
log("loading QC'd query samples")
parts = []
for f in sorted(os.listdir(QC)):
    if not f.endswith(".h5ad"):
        continue
    s = B.SingleCell(f"{QC}/{f}", num_threads=B.NUM_THREADS)
    B.verify_load(s, f)
    parts.append(s)
# flexible=True because each per-sample h5ad carries its own uns['cell_call'] record;
# brisc otherwise refuses to concat datasets whose uns differ. The instance-method and
# single-positional-arg forms of concat_obs both fail outright (see BRISC_BUGS.md 5).
query = B.concat_obs(parts, flexible=True)
del parts
gc.collect()
query = query.skip_qc()          # QC already applied in 03; flag the state for brisc
query.obs = query.obs.with_columns(pl.col("sample_id").cast(pl.String).alias("donor_id"))
log(f"query {query.shape}, {query.obs['sample_id'].n_unique()} samples")

# ---- shared gene space -------------------------------------------------------
shared = sorted(set(ref.var_names) & set(query.var_names))
log(f"shared genes: {len(shared)}")
ref = ref.filter_var(pl.col._index.is_in(shared))
query = query.filter_var(pl.col._index.is_in(shared))
gc.collect()

# ---- joint HVG / normalize / PCA / Harmony -----------------------------------
log("hvg (raw counts, batched by donor)")
ref, query = ref.hvg(query, batch_column="donor_id", num_genes=N_HVG,
                     verbose=False)

log("normalize (full gene set, so library sizes are honest)")
ref = ref.normalize()
query = query.normalize()

log("dropping to HVGs before PCA")
ref = ref.filter_var("highly_variable")
query = query.filter_var("highly_variable")
gc.collect()

log("pca")
ref, query = ref.pca(query, num_PCs=50, seed=SEED, verbose=False)
log("harmonize")
ref, query = ref.harmonize(query, batch_column="donor_id", seed=SEED, verbose=False)

# ---- label transfer ----------------------------------------------------------
for level in ["Class", "Subclass", "Supertype"]:
    log(f"label transfer: {level}")
    query = query.label_transfer_from(
        ref, level,
        cell_type_column=level.lower(),
        confidence_column=f"{level.lower()}_confidence",
        next_best=True,
        next_best_cell_type_column=f"{level.lower()}_next",
        next_best_confidence_column=f"{level.lower()}_next_confidence",
        num_neighbors=20, seed=SEED, verbose=False)

# ---- save the labels before anything optional can fail ------------------------
os.makedirs(f"{DATA}/annotated", exist_ok=True)
query.obs.to_pandas().to_parquet(f"{DATA}/annotated/query_obs{SUFFIX}.parquet")
np.save(f"{DATA}/annotated/query_harmony{SUFFIX}.npy", query.obsm["harmony"])
log("saved transferred labels + harmony embedding")

# ---- de novo clustering, to check the transferred labels against ---------------
log("de novo clustering (neighbors -> SNN -> leiden)")
query = query.neighbors(seed=SEED, PC_key="harmony")
query = query.shared_neighbors()
query = query.cluster(resolution=1.0, seed=SEED)
log("umap")
query = query.umap(seed=SEED, PC_key="harmony")

obs = query.obs.to_pandas()
obs.to_parquet(f"{DATA}/annotated/query_obs{SUFFIX}.parquet")
np.save(f"{DATA}/annotated/query_umap{SUFFIX}.npy", query.obsm["umap"])
log(f"wrote query_obs.parquet {obs.shape}")

print("\n=== transferred subclass composition ===")
print(obs["subclass"].value_counts().to_string())
print("\n=== confidence by level ===")
for level in ["class", "subclass", "supertype"]:
    c = obs[f"{level}_confidence"]
    print(f"{level:10s} median={c.median():.3f}  frac>0.5={np.mean(c > 0.5):.3f}  "
          f"frac>0.8={np.mean(c > 0.8):.3f}")
