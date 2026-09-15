#!/usr/bin/env python3
"""Build the two small Xenium inputs Figure 3 needs, from the full cell table.

Figure 3 panels d/e draw two sections, and panel h needs one number per Sst
supertype. Neither needs the 1.34M-row, 332 MB per-cell metadata dump that used
to be committed here, and a repository is the wrong place for that dump anyway
(it is git-ignored now). This script derives the two small inputs from it, so
the reduction is reproducible rather than a hand-made file of unknown origin.

Input (git-ignored; regenerate on the cluster, or export from the Seurat object
with `write.csv(seu@meta.data, ...)` as step 1 of Figure_3.r once did):
    Data/xenium_metadata.csv

Outputs (committed, ~5 MB total):
    Data/xenium_sections_fig3de.csv      all cells of Br5931 and Br1139 (qc_pass kept
                                         as a column: panels d/e draw every cell,
                                         as they always have, so it is not applied)
    Data/xenium_sst_supertype_depth.csv  mean predicted depth per Sst supertype

Run from snrnaseq/Final_figures/.
"""
import os
import sys

import pandas as pd

SRC = "Data/xenium_metadata.csv"
SECTIONS = ["Br5931", "Br1139"]          # the control / SCZ pair drawn in d and e
KEEP = ["sample_id", "x", "y", "subclass", "supertype", "layer", "qc_pass"]

if not os.path.exists(SRC):
    sys.exit(f"{SRC} not found. It is git-ignored; see this file's docstring.")

cols = set(KEEP) | {"qc_pass", "predicted_norm_depth"}
d = pd.read_csv(SRC, usecols=lambda c: c in cols, low_memory=False)
qc = d["qc_pass"].astype(str).isin(["True", "TRUE"])

# No QC filter here: the original panels drew every cell in the window, and
# silently dropping ~8% of them would change the published figure.
sec = d[d["sample_id"].isin(SECTIONS)][KEEP]
sec.to_csv("Data/xenium_sections_fig3de.csv", index=False)
print(f"wrote Data/xenium_sections_fig3de.csv  {len(sec):,} cells, "
      f"{sec.sample_id.nunique()} sections")

dep = (d[qc & (d["subclass"] == "Sst") & d["predicted_norm_depth"].notna()]
       .groupby("supertype", as_index=False)["predicted_norm_depth"].mean()
       .rename(columns={"supertype": "CellType", "predicted_norm_depth": "mean_depth"}))
dep.to_csv("Data/xenium_sst_supertype_depth.csv", index=False)
print(f"wrote Data/xenium_sst_supertype_depth.csv  {len(dep)} Sst supertypes")
