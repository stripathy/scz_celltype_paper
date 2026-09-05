#!/usr/bin/env python3
"""
Pseudobulk differential expression in GSE158516, mirroring the paper's snRNA-seq model.

Our meta-analysis fits, per dataset, per SEA-AD subclass:
    pseudobulk counts ~ diagnosis + age + sex + PMI      (limma-voom / edgeR)
Sex drops out here because every one of the 26 donors is male. brisc's DE is
limma-voom over ryp, so this is the same estimator, not an approximation.

Two thresholds differ from the meta-analysis and are deliberate:
  * cells per donor per cell type: 10, not 500. A 26-donor cohort at ~10k nuclei/donor
    puts Sst at a few hundred cells per donor, so the 500 rule would delete the very
    cell type this replication exists to test. 10 is the threshold our own Xenium
    replication used for the same reason.
  * genes: brisc's default (>=1 count in >=80% of samples per group) matches the
    meta-analysis's ">=80% of retained donors" rule.

Also exports per-donor x cell-type counts for the crumblr compositional analysis (06).
"""
import gc
import os

import brisc_setup as B
import numpy as np
import pandas as pd
import polars as pl

QC = "/Users/shreejoy/Github/shared_data/GSE158516/per_sample_qc"
DATA = "/Users/shreejoy/Github/shared_data/GSE158516"
OUTD = "/Users/shreejoy/Github/scz_celltype_paper/external_validation/gse158516/output"
META = "/Users/shreejoy/Github/scz_celltype_paper/external_validation/gse158516/data/sample_metadata.csv"

MIN_CELLS = 10
FORMULA = "~ diagnosis + age + PMI + log2(num_cells) + log2(library_size)"
COEF = "diagnosisSchizophrenia"


def log(m):
    print(f"[05] {m}", flush=True)


meta = pd.read_csv(META).rename(columns={
    "Deidentified ID": "sample_id", "Diagnosis": "diagnosis", "Age": "age",
    "PMI": "PMI", "Sequencing Batch": "batch", "PFC_pH": "pH"})
keep_samples = set(meta["sample_id"])
log(f"{len(keep_samples)} donors with metadata (the authors' QC-passing set)")

ann = pd.read_parquet(f"{DATA}/annotated/query_obs.parquet")
LABEL_COLS = ["_index", "class", "subclass", "supertype",
              "subclass_confidence", "supertype_confidence"]
lab = pl.from_pandas(ann[LABEL_COLS].astype({c: str for c in LABEL_COLS[:4]}))
log(f"annotations: {len(ann):,} cells")

# ---- load query, attach labels ----------------------------------------------
parts = []
for f in sorted(os.listdir(QC)):
    sid = f[:-5]
    if not f.endswith(".h5ad") or sid not in keep_samples:
        continue
    s = B.SingleCell(f"{QC}/{f}", num_threads=B.NUM_THREADS)
    B.verify_load(s, sid)
    parts.append(s)
query = B.concat_obs(parts, flexible=True)
del parts
gc.collect()
query = query.skip_qc()
log(f"query {query.shape}")

md = pl.from_pandas(meta[["sample_id", "diagnosis", "age", "PMI", "batch"]]
                    .astype({"sample_id": str, "diagnosis": str}))
query.obs = (query.obs
             .drop("diagnosis", strict=False)      # re-joined from meta below
             .with_columns(pl.col("sample_id").cast(pl.String))
             .join(lab, on="_index", how="left")
             .join(md, on="sample_id", how="left"))
n_missing = query.obs["subclass"].null_count()
if n_missing:
    log(f"WARNING: {n_missing} cells without a transferred label; dropping")
    query = query.filter_obs(pl.col("subclass").is_not_null())

# ---- per-donor x cell-type counts for crumblr (06) ---------------------------
for level in ["subclass", "supertype", "class"]:
    counts = (query.obs.group_by(["sample_id", level]).len()
              .to_pandas().pivot(index="sample_id", columns=level, values="len")
              .fillna(0).astype(int))
    counts.to_csv(f"{OUTD}/counts_{level}.csv")
    log(f"wrote counts_{level}.csv {counts.shape}")

# Compartment map for the stratified crumblr model. SEA-AD Class is one of
# "Neuronal: GABAergic" / "Neuronal: Glutamatergic" / "Non-neuronal and Non-neural".
comp = []
for level in ["subclass", "supertype"]:
    m = (query.obs.select(["class", level]).unique().to_pandas()
         .rename(columns={level: "cell_type"}))
    m["compartment"] = np.where(m["class"].astype(str).str.startswith("Neuronal"),
                                "neuronal", "non-neuronal")
    comp.append(m[["cell_type", "compartment"]])
pd.concat(comp).drop_duplicates().to_csv(f"{OUTD}/celltype_compartment.csv", index=False)
log("wrote celltype_compartment.csv")

# ---- pseudobulk + DE ---------------------------------------------------------
log("pseudobulk by (sample_id, subclass)")
pb = query.pseudobulk("sample_id", "subclass")
del query
gc.collect()

sizes = pd.DataFrame({ct: dict(n_samples=obs.height,
                               median_cells=float(np.median(obs["num_cells"])),
                               min_cells=int(obs["num_cells"].min()))
                      for ct, (X, obs, var) in pb.items()}).T
sizes.to_csv(f"{OUTD}/pseudobulk_sizes.csv")
print("\npseudobulk sizes per subclass:\n" + sizes.sort_values("median_cells",
                                                              ascending=False).to_string())

log("pseudobulk qc + TMM library sizes")
pb = pb.qc("diagnosis", min_cells=MIN_CELLS, allow_float=True, verbose=False)
pb = pb.library_size(allow_float=True)
log(f"cell types surviving pseudobulk QC: {len(pb.keys())}")

# Persist before DE: rebuilding the pseudobulk means re-loading ~9 GB of query, so a
# failure in the design matrix should not cost that.
pb.save(f"{DATA}/pseudobulk_subclass", overwrite=True)
log(f"saved pseudobulk to {DATA}/pseudobulk_subclass")

log(f"DE: {FORMULA}  coefficient={COEF}")
# B.de() handles both the 0.1.3 DE->de rename and the integer-covariate cast that
# brisc still needs as of 0.1.3 (see brisc_setup.de docstring).
de = B.de(pb, FORMULA, coefficient=COEF, allow_float=True, verbose=False, strict=True)
os.makedirs(f"{DATA}/de", exist_ok=True)
de.save(f"{DATA}/de/gse158516_subclass_de", overwrite=True)
tab = de.table.to_pandas()
tab.to_csv(f"{OUTD}/de_subclass.csv.gz", index=False, compression="gzip")
log(f"wrote de_subclass.csv.gz {tab.shape}")

print("\n=== hits per cell type (FDR < 0.05) ===")
print(de.get_num_hits(significance_column="FDR", threshold=0.05).to_pandas().to_string(index=False))
print("\n=== the two canonical markers ===")
for gene, ct in [("SST", "Sst"), ("PVALB", "Pvalb"), ("VIP", "Vip"), ("GAD1", "Sst"), ("GAD2", "Sst")]:
    r = tab[(tab["gene"] == gene) & (tab["cell_type"] == ct)]
    if len(r):
        r = r.iloc[0]
        print(f"{gene:6s} in {ct:6s}: logFC={r['logFC']:+.3f}  P={r['p']:.3g}  FDR={r['FDR']:.3g}")
    else:
        print(f"{gene:6s} in {ct:6s}: not tested")
