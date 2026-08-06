#!/usr/bin/env python3
"""
Per-cell cortical depth for each neuronal supertype in BOTH spatial platforms, for a
Xenium-vs-MERFISH depth comparison (violins + per-supertype summary in R).

This is the distributional counterpart to the median-depth scatter in Supplementary
Fig. S2: instead of one point per supertype, it shows the whole depth distribution the
two platforms assign, so a reader can see where they agree in position but differ in
spread.

The two depth values are NOT the same kind of measurement, and the figure exists partly
to make that legible:
  - MERFISH  "Normalized depth from pia" — annotated by hand by SEA-AD.
  - Xenium   predicted_norm_depth        — the neighbourhood-composition model's output,
                                           which was trained on the MERFISH annotations.
So MERFISH is the reference the Xenium model was fit against; agreement is expected, and
the informative signal is where it breaks down (and the spread, which the model shrinks).

Restricts to neuronal (Glutamatergic + GABAergic) cortical cells present in both
platforms. Supertypes too rare to analyse in Xenium are dropped (see SM1).

Outputs to output/depth_platform/:
  supertype_depth_platform.csv.gz       per-cell (subsampled for plotting)
  supertype_depth_platform_summary.csv  per-supertype medians/IQR/n on ALL cells
"""
import os, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (H5AD_DIR, EXCLUDE_SAMPLES, SUBCLASS_TO_CLASS,
                    load_cells, load_merfish_cortical)

_SP = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(_SP, "output", "depth_platform"); os.makedirs(OUT, exist_ok=True)

TOO_RARE = {"Lamp5 Lhx6"}          # SM1: 16 cortical cells in Xenium
NEURONAL = ("Glutamatergic", "GABAergic")
MAX_PER_CELL_GROUP = 4000          # violin plotting only; summaries use every cell
RNG = np.random.default_rng(0)     # fixed so the figure is reproducible

# ---- MERFISH: manual depth annotation ----
print("Loading MERFISH cortical cells (manual depth)...")
mer = load_merfish_cortical().rename(columns={"depth": "depth"})
mer["platform"] = "MERFISH"
mer = mer[["donor", "subclass", "supertype", "depth", "platform"]]
print(f"  {len(mer):,} cells, {mer.donor.nunique()} donors")

# ---- Xenium: model-predicted depth ----
print("Loading Xenium cortical cells (predicted depth)...")
rows = []
for f in sorted(os.listdir(H5AD_DIR)):
    if not f.endswith("_annotated.h5ad"):
        continue
    sid = f.replace("_annotated.h5ad", "")
    if sid in EXCLUDE_SAMPLES:
        continue
    o = load_cells(sid, cortical_only=True, extra_obs_columns=["predicted_norm_depth"])
    rows.append(o[["sample_id", "subclass_label", "supertype_label", "predicted_norm_depth"]])
xen = pd.concat(rows, ignore_index=True).rename(columns={
    "sample_id": "donor", "subclass_label": "subclass",
    "supertype_label": "supertype", "predicted_norm_depth": "depth"})
xen["platform"] = "Xenium"
print(f"  {len(xen):,} cells, {xen.donor.nunique()} donors")

df = pd.concat([mer, xen], ignore_index=True)
df["depth"] = pd.to_numeric(df["depth"], errors="coerce")
df["class"] = df["subclass"].map(SUBCLASS_TO_CLASS)
df = df[df["class"].isin(NEURONAL) & df["depth"].notna() & ~df["subclass"].isin(TOO_RARE)]

# keep only supertypes both platforms actually measure, so every violin pair is complete
both = (df.groupby(["supertype", "platform"]).size().unstack(fill_value=0)
          .pipe(lambda t: t[(t["MERFISH"] >= 50) & (t["Xenium"] >= 50)]).index)
dropped = sorted(set(df.supertype) - set(both))
df = df[df.supertype.isin(both)]
print(f"\n{len(both)} neuronal supertypes with >=50 cells on BOTH platforms")
if dropped:
    print(f"  not shown ({len(dropped)}, <50 cells on one platform): {', '.join(dropped)}")

# ---- per-supertype summary on ALL cells (before any subsampling) ----
s = (df.groupby(["class", "subclass", "supertype", "platform"])["depth"]
       .agg(n="size", median="median",
            q25=lambda x: x.quantile(.25), q75=lambda x: x.quantile(.75))
       .reset_index())
wide = s.pivot_table(index=["class", "subclass", "supertype"], columns="platform",
                     values=["n", "median", "q25", "q75"])
wide.columns = [f"{a}_{b}" for a, b in wide.columns]
wide = wide.reset_index()
wide["delta_median"] = wide["median_Xenium"] - wide["median_MERFISH"]
wide["iqr_MERFISH"] = wide["q75_MERFISH"] - wide["q25_MERFISH"]
wide["iqr_Xenium"] = wide["q75_Xenium"] - wide["q25_Xenium"]
wide.to_csv(os.path.join(OUT, "supertype_depth_platform_summary.csv"), index=False)

r = np.corrcoef(wide["median_MERFISH"], wide["median_Xenium"])[0, 1]
print(f"\nmedian depth across {len(wide)} supertypes: Pearson r = {r:.3f}")
print(f"  |delta median| : median {wide.delta_median.abs().median():.3f}, "
      f"90th pct {wide.delta_median.abs().quantile(.9):.3f}, "
      f"max {wide.delta_median.abs().max():.3f} ({wide.loc[wide.delta_median.abs().idxmax(),'supertype']})")
print(f"  IQR            : MERFISH median {wide.iqr_MERFISH.median():.3f} vs "
      f"Xenium {wide.iqr_Xenium.median():.3f}")

# ---- subsample for plotting ----
def cap(g):
    return g if len(g) <= MAX_PER_CELL_GROUP else g.iloc[
        RNG.choice(len(g), MAX_PER_CELL_GROUP, replace=False)]
plot = (df.groupby(["supertype", "platform"], group_keys=False).apply(cap)
          [["donor", "platform", "class", "subclass", "supertype", "depth"]])
# gzipped: ~20 MB of per-cell depth compresses to ~3 MB, small enough to commit so
# the figure regenerates from a clean clone. readr/pandas both handle .gz natively.
plot.to_csv(os.path.join(OUT, "supertype_depth_platform.csv.gz"), index=False)
print(f"\nsaved supertype_depth_platform.csv.gz: {len(plot):,} cells "
      f"(capped at {MAX_PER_CELL_GROUP:,} per supertype x platform)")
print(f"saved supertype_depth_platform_summary.csv: {len(wide)} supertypes")
