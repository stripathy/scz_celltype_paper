"""
Per-cell SST/PVALB grain-density input for composite-figure panel K (the
Dienel-style grains/cell-area metric) and the panel-L exemplar selection.

For each SST-in-Sst and PVALB-in-Pvalb cell on the CANONICAL definition
(corr_subclass + cortical + qc_pass & corr_qc_pass -- identical to the crumblr
compositional analysis and the edgeR DE), across all 24 Xenium donors, records
the raw marker count, the cell area (from the deploy boundary polygon), and the
library size, with donor/diagnosis/sex/age.

Grain density = marker count / cell area (grains/100 um^2), the Dienel
(2023, Am J Psychiatry) per-neuron mRNA measure.

DATA PROVENANCE (SCZ_Xenium repo): output/h5ad/<s>_annotated.h5ad (.X raw counts,
corr_subclass, spatial_domain, qc_pass, corr_qc_pass, total_counts) and
output/deploy/boundaries/<s>.json (cell polygons, obs order).

Output: results/tables/percell_grain_density.csv
        columns: sample, dx, sex, age, gene, subclass, count, total_counts, cell_area_um2
Consumed by: scripts/09_composite_figure.R (panel K boxplots) and
             scripts/10_xenium_exemplar_cells.py (exemplar grain-density target).
"""
import os
import sys
import json
import time
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

XEN  = os.path.expanduser("~/Github/SCZ_Xenium")
H5AD = os.path.join(XEN, "output/h5ad/{s}_annotated.h5ad")
BND  = os.path.join(XEN, "output/deploy/boundaries/{s}.json")
OUT  = "results/tables"
sys.path.insert(0, os.path.join(XEN, "code/analysis"))
sys.path.insert(0, os.path.join(XEN, "code"))
from config import BASE_DIR, SAMPLE_TO_DX             # noqa: E402
from modules.metadata import get_subject_info          # noqa: E402

PAIRS = [("SST", "Sst"), ("PVALB", "Pvalb")]


def canonical_mask(o, subclass, n):
    sc = (o["corr_subclass"] if "corr_subclass" in o else o["subclass_label"]).astype(str).values
    cort = (o["spatial_domain"].astype(str).values == "Cortical")
    qc = o["qc_pass"].values.astype(bool) if "qc_pass" in o else np.ones(n, bool)
    cq = o["corr_qc_pass"].values.astype(bool) if "corr_qc_pass" in o else np.ones(n, bool)
    return (sc == subclass) & cort & qc & cq


def cell_areas(sample, n_obs):
    """Vectorised cell-polygon area (um^2) per cell; index = h5ad obs order."""
    b = json.load(open(BND.format(s=sample)))
    V = b["verts_per_cell"]
    X = np.asarray(b["bx"]).reshape(-1, V) * b["x_scale"] + b["x_offset"]
    Y = np.asarray(b["by"]).reshape(-1, V) * b["y_scale"] + b["y_offset"]
    assert X.shape[0] == n_obs, f"boundary/obs mismatch {sample}: {X.shape[0]} vs {n_obs}"
    return 0.5 * np.abs((X * np.roll(Y, -1, axis=1) - np.roll(X, -1, axis=1) * Y).sum(axis=1))


def xcount(a, gene):
    col = a.X[:, a.var_names.get_loc(gene)]
    return np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()


def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    meta = get_subject_info(os.path.join(BASE_DIR, "data", "sample_metadata.xlsx")).set_index("sample_id")
    rows = []
    for i, s in enumerate(sorted(SAMPLE_TO_DX), 1):
        a = ad.read_h5ad(H5AD.format(s=s))
        o = a.obs
        area = cell_areas(s, a.n_obs)
        tot = o["total_counts"].values.astype(int)
        dx, sex, age = SAMPLE_TO_DX[s], meta.loc[s, "sex"], float(meta.loc[s, "age"])
        for gene, sub in PAIRS:
            cnt = xcount(a, gene).astype(int)
            sel = canonical_mask(o, sub, a.n_obs)
            for j in np.where(sel)[0]:
                rows.append((s, dx, sex, age, gene, sub, int(cnt[j]), int(tot[j]), float(area[j])))
        print(f"[{i:2d}/24] {s} ({dx})")
        del a
    df = pd.DataFrame(rows, columns=["sample", "dx", "sex", "age", "gene", "subclass",
                                     "count", "total_counts", "cell_area_um2"])
    df.to_csv(f"{OUT}/percell_grain_density.csv", index=False)
    print(f"\nSaved {len(df):,} cells -> {OUT}/percell_grain_density.csv  ({time.time()-t0:.0f}s)")
    g = df.assign(gd=df["count"] / df["cell_area_um2"] * 100)
    print(g.groupby(["gene", "dx"]).agg(n=("count", "size"),
          median_grain_density=("gd", "median")).round(2).to_string())


if __name__ == "__main__":
    main()
