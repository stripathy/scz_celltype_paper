#!/usr/bin/env python3
"""
Compare our predicted cortical layers against the original authors' (Lieber/Kwon)
spatial-domain layer annotations from spatialDLPFC_SCZ_XENIUM.

Their labels: processed-data/04_label_transfer/label_transfer_N24_k50_smoothed_labels.csv
  cellid = "{sample}_{rawpos}" (1-indexed raw Xenium cell position), spd01..spd07.
Both pipelines read the same raw Xenium in raw order and never reorder, so we
match by (sample, raw position): our 0-indexed cell p  <->  their "{sample}_{p+1}".

Outputs (output/depth_validation/lieber_layers/):
  joined_cells.csv         per-cell our_layer/our_depth/their_spd (matched, cortical)
  confusion_spd_layer.csv  crosstab their spd x our layer
  depth_per_spd.csv        our predicted-depth distribution per their spd
"""
import os, sys
import numpy as np
import pandas as pd
import anndata as ad

_HERE = os.path.dirname(os.path.abspath(__file__))
_SPATIAL = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_SPATIAL, "code", "analysis"))
from config import H5AD_DIR, EXCLUDE_SAMPLES

LAB = os.path.join(_SPATIAL, "output/depth_validation/lieber_layers")
CSV = os.path.join(LAB, "label_transfer_N24_k50_smoothed_labels.csv")
OUR_LAYERS = ["L1", "L2/3", "L4", "L5", "L6", "WM"]


def main():
    their = pd.read_csv(CSV)
    their.columns = ["cellid", "spd"]
    their = their.set_index("cellid")["spd"]
    print(f"their labels: {len(their):,} cells, domains {sorted(their.unique())}")

    parts = []
    for f in sorted(os.listdir(H5AD_DIR)):
        if not f.endswith("_annotated.h5ad"):
            continue
        sid = f.replace("_annotated.h5ad", "")
        if sid in EXCLUDE_SAMPLES:
            continue
        o = ad.read_h5ad(os.path.join(H5AD_DIR, f), backed="r").obs
        n = len(o)
        df = pd.DataFrame({
            "cellid": [f"{sid}_{i+1}" for i in range(n)],
            "sample": sid,
            "our_layer": o["layer"].values.astype(str),
            "spatial_domain": o["spatial_domain"].values.astype(str),
            "our_depth": o["predicted_norm_depth"].values.astype(float),
        })
        df["their_spd"] = df["cellid"].map(their)
        parts.append(df)
    ours = pd.concat(parts, ignore_index=True)

    matched = ours["their_spd"].notna()
    print(f"\nmatch rate: {matched.mean()*100:.1f}% of our cells got a Lieber label "
          f"({matched.sum():,}/{len(ours):,}); their cells matched: "
          f"{matched.sum():,}/{len(their):,} = {100*matched.sum()/len(their):.1f}%")

    m = ours[matched].copy()
    m.to_csv(os.path.join(LAB, "joined_cells.csv"), index=False)

    # depth per their spd (reveals laminar ordering of spd's)
    dps = m.groupby("their_spd")["our_depth"].agg(["count", "median", "mean", "std"]).round(4)
    dps.to_csv(os.path.join(LAB, "depth_per_spd.csv"))
    print("\nOur predicted depth per their spatial domain (median, pia=0 -> WM=1):")
    print(dps.sort_values("median").to_string())

    # confusion: their spd x our layer (cortical + WM cells only)
    cm = m[m["our_layer"].isin(OUR_LAYERS)]
    ct = pd.crosstab(cm["their_spd"], cm["our_layer"])
    ct = ct.reindex(columns=[c for c in OUR_LAYERS if c in ct.columns])
    ct.to_csv(os.path.join(LAB, "confusion_spd_layer.csv"))
    # spd -> our-layer correspondence (argmax) + agreement under that mapping
    row_norm = ct.div(ct.sum(axis=1), axis=0)
    spd2layer = row_norm.idxmax(axis=1).to_dict()
    cm = cm.copy()
    cm["spd_layer"] = cm["their_spd"].map(spd2layer)
    agree = (cm["spd_layer"] == cm["our_layer"]).mean()
    print("\nspd -> dominant our-layer correspondence:")
    for spd in sorted(spd2layer): print(f"  {spd} -> {spd2layer[spd]}  ({row_norm.loc[spd, spd2layer[spd]]*100:.0f}% of {spd})")
    print(f"\nOverall agreement (their spd mapped to our layer): {agree*100:.1f}%  (n={len(cm):,} cortical+WM cells)")
    print(f"\nConfusion matrix (rows=their spd, cols=our layer):\n{ct.to_string()}")


if __name__ == "__main__":
    main()
