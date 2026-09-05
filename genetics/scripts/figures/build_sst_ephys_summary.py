"""
Rebuild results/tables/sst_supertype_ephys_summary.csv, the per-Sst-supertype
patch-seq layer and intrinsic-property table behind Figure 4d and the depth
gradient.

The table had been a carried-over artifact: scripts/06_gene_ephys_correlations.py
reads a `sst_supertype_layer_info.csv` that does not exist, and its write sits
inside `if sst_layer_info is not None`, so it silently never fired. This rebuilds
the table from the upstream patch-seq release and reproduces the committed values
exactly.

Cell selection, and why:
  subclass_scANVI == "Sst"          the scANVI call is the only label source used
                                    by the figure; the kNN transfer in the same
                                    release is a cross-check and is not read here.
                                    An exact string match also excludes the seven
                                    "Sst Chodl" cells, which are a distinct
                                    subclass rather than a low-confidence tail.
  cortical_layer not null           needed for the depth gradient. It costs five
                                    cells that have sag but no layer call, and
                                    moves the sag correlation by 0.003, so the
                                    same filter is left in place for both panels
                                    rather than letting them diverge.

`n_cells` counts layer-annotated cells. Sag and tau are means over whichever of
those cells carry an extracted feature, which is fewer -- 129 of 150 for sag.
The two counts are reported separately so the distinction is visible.
"""
import pandas as pd

import os as _os

# Repo-relative, so the chain runs from any clone. External data that is not in
# the repo is still resolved by absolute path or an environment override below.
GEN = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
PAPER = _os.path.dirname(GEN)

UP = _os.environ.get("PATCHSEQ_REPO",
                     _os.path.expanduser("~/Github/human_int_patch_seq"))

SRC = f"{UP}/data/patchseq/patchseq_combined.csv"
OUT = f"{GEN}/results/tables/sst_supertype_ephys_summary.csv"

# cortical_layer is an ordinal string, not a depth in microns (soma_depth_um is
# empty for every cell in this release). "5_6" denotes a cell at the layer 5/6
# border and is coded as 5.5.
LAYER_CODE = {"2": 2.0, "3": 3.0, "4": 4.0, "5": 5.0, "5_6": 5.5, "6": 6.0}

cells = pd.read_csv(SRC, low_memory=False)
sst = cells[(cells.subclass_scANVI == "Sst") & cells.cortical_layer.notna()].copy()
sst["layer_num"] = sst.cortical_layer.astype(str).map(LAYER_CODE)
assert sst.layer_num.notna().all(), sorted(set(sst.cortical_layer.astype(str)))

g = sst.groupby("supertype_scANVI")
out = pd.DataFrame({
    "mean_layer": g.layer_num.mean(),
    "median_layer": g.layer_num.median(),
    "n_cells": g.size(),
    "mean_sag": g.sag.mean(),
    # upstream tau is in seconds; the table has always carried milliseconds
    "mean_tau": g.tau.mean() * 1000.0,
}).reset_index().rename(columns={"supertype_scANVI": "supertype"})

# Upper/deep split at the layer 3/4 boundary, matching the committed table.
out["layer_group"] = out.mean_layer.map(
    lambda m: "Upper (L2-3)" if m < 3.5 else "Deep (L4-6)")
out = out.sort_values("mean_layer").reset_index(drop=True)
out.to_csv(OUT, index=False)

n_sag = int(sst.sag.notna().sum())
print(f"{len(out)} Sst supertypes | {int(out.n_cells.sum())} layer-annotated cells "
      f"({sst.donor.nunique()} donors) | {n_sag} with sag")
