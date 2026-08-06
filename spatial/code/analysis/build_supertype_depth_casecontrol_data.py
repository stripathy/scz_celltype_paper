#!/usr/bin/env python3
"""
Per-cell cortical depth for each NEURONAL supertype, labelled by diagnosis and donor,
for a Control-vs-SCZ depth comparison (violins + per-supertype mixed-effects model in R).

Restricts to neuronal (Glutamatergic + GABAergic) cortical QC-pass cells with a predicted
depth. Exports a tidy per-cell CSV (donor, dx, class, subclass, supertype, depth) plus a
per-supertype power summary so the R step can filter to supertypes estimable with a
donor-random-effect model.

Output: output/depth_casecontrol/supertype_depth_casecontrol.csv.gz (+ _counts.csv)
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import H5AD_DIR, EXCLUDE_SAMPLES, SAMPLE_TO_DX, SUBCLASS_TO_CLASS, load_cells

_SP = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # spatial/
OUT = os.path.join(_SP, "output", "depth_casecontrol"); os.makedirs(OUT, exist_ok=True)

rows = []
for f in sorted(os.listdir(H5AD_DIR)):
    if not f.endswith("_annotated.h5ad"): continue
    sid = f.replace("_annotated.h5ad", "")
    if sid in EXCLUDE_SAMPLES: continue
    dx = SAMPLE_TO_DX.get(sid)
    if dx not in ("Control", "SCZ"): continue
    obs = load_cells(sid, cortical_only=True, extra_obs_columns=["predicted_norm_depth"])
    obs = obs[["sample_id", "subclass_label", "supertype_label", "predicted_norm_depth"]].copy()
    obs["dx"] = dx
    rows.append(obs)
    print(f"  {sid} ({dx}): {len(obs):,} cortical cells", flush=True)

df = pd.concat(rows, ignore_index=True)
df["class"] = df["subclass_label"].map(SUBCLASS_TO_CLASS)
df["depth"] = pd.to_numeric(df["predicted_norm_depth"], errors="coerce")
df = df[df["class"].isin(["Glutamatergic", "GABAergic"]) & df["depth"].notna()]
df = df.rename(columns={"sample_id": "donor", "subclass_label": "subclass", "supertype_label": "supertype"})
df = df[["donor", "dx", "class", "subclass", "supertype", "depth"]]
# gzipped: ~20 MB of per-cell depth compresses to ~4 MB, small enough to commit so
# the figure regenerates from a clean clone. readr/pandas both handle .gz natively.
df.to_csv(os.path.join(OUT, "supertype_depth_casecontrol.csv.gz"), index=False)
print(f"\nsaved supertype_depth_casecontrol.csv.gz: {len(df):,} neuronal cells, "
      f"{df.supertype.nunique()} supertypes, {df.donor.nunique()} donors "
      f"({(df.dx=='Control').sum():,} Control / {(df.dx=='SCZ').sum():,} SCZ cells)")

# per-supertype power summary
g = df.groupby(["class", "subclass", "supertype", "dx"]).agg(n=("depth", "size"),
        ndonor=("donor", "nunique")).reset_index()
piv = g.pivot_table(index=["class", "subclass", "supertype"], columns="dx",
                    values=["n", "ndonor"], fill_value=0)
piv.columns = [f"{a}_{b}" for a, b in piv.columns]; piv = piv.reset_index()
piv.to_csv(os.path.join(OUT, "supertype_depth_casecontrol_counts.csv"), index=False)
for nc, nd in [(50, 3), (100, 3), (30, 2)]:
    ok = piv[(piv.n_Control >= nc) & (piv.n_SCZ >= nc) & (piv.ndonor_Control >= nd) & (piv.ndonor_SCZ >= nd)]
    print(f"  supertypes passing >={nc} cells & >={nd} donors per group: {len(ok)} "
          f"(Glut {sum(ok['class']=='Glutamatergic')}, GABA {sum(ok['class']=='GABAergic')})")
