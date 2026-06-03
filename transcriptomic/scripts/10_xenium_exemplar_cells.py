"""
Extract exemplar Xenium cells illustrating marker downregulation in SCZ:
  - SST  in Sst        cells: one Control vs one SCZ
  - BDNF in L2/3 IT    cells: one Control vs one SCZ

For each chosen cell we save (in micron coordinates, recentred on the cell):
  - the cell boundary polygon
  - the marker-gene transcript molecules that fall INSIDE that polygon
plus a metadata row (sample, marker count, total counts).

These small CSVs are read by scripts/09_composite_figure.R to draw the
exemplar-cell panels.

DATA PROVENANCE (SCZ_Xenium repo, Kwon 2026):
  ~/Github/SCZ_Xenium/output/h5ad/<sample>_annotated.h5ad
        cells x 300 genes; .X = raw integer counts; obs.subclass_label,
        obs.qc_pass, obs.total_counts; obsm['spatial'] = centroids.
  ~/Github/SCZ_Xenium/output/deploy/boundaries/<sample>.json
        cell polygons in obs order (25 verts/cell); decode:
        micron = quant * x_scale + x_offset.
  ~/Github/SCZ_Xenium/output/deploy/transcripts/<sample>/<GENE>.json (+ gene_index.json)
        per-gene molecule coords; decode with gene_index offsets.
  Diagnosis per sample: SCZ_Xenium code/analysis/config.py SAMPLE_TO_DX.

Output: results/tables/exemplar_<gene>_<dx>_{boundary,dots}.csv
        results/tables/exemplar_cells_meta.csv
"""

import json
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
from matplotlib.path import Path
from scipy.spatial import ConvexHull

XEN  = os.path.expanduser("~/Github/SCZ_Xenium")
H5AD = os.path.join(XEN, "output/h5ad/{s}_annotated.h5ad")
BND  = os.path.join(XEN, "output/deploy/boundaries/{s}.json")
TXG  = os.path.join(XEN, "output/deploy/transcripts/{s}/{g}.json")
TXI  = os.path.join(XEN, "output/deploy/transcripts/{s}/gene_index.json")
OUT  = "results/tables"

CONTROL_SAMPLE = "Br6432"   # control sample with transcript export
SCZ_SAMPLE     = "Br2039"   # SCZ sample with transcript export
# (gene, Xenium subclass_label). Chosen to be (1) strongly DE-down in the
# snRNA-seq meta-analysis, (2) highly enough expressed per cell in Xenium to
# show as molecule dots, and (3) depth-matchable between Control/SCZ cells:
#   SST in Sst interneurons  (meta padj 0.049; 42 vs 8 dots)
#   RASGRF2 in Pvalb interneurons (meta padj 0.022; 14 vs 8 dots)
# Genes like SMAD1/BDNF (strong meta-DE but ~1 dot/cell) and CX3CR1/TF
# (depth-confounded) were rejected — see notes.
PAIRS = [("SST", "Sst"), ("RASGRF2", "Pvalb")]

_adata_cache = {}
def load(sample):
    if sample not in _adata_cache:
        _adata_cache[sample] = ad.read_h5ad(H5AD.format(s=sample))
    return _adata_cache[sample]

def gene_counts(a, gene):
    col = a.X[:, a.var_names.get_loc(gene)]
    return np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()

_bnd_cache = {}
def all_polys(sample):
    """Decode every cell's boundary polygon (micron coords), vectorised.
    Returns (X, Y) arrays of shape (n_cells, verts_per_cell)."""
    if sample not in _bnd_cache:
        b = json.load(open(BND.format(s=sample)))
        V = b["verts_per_cell"]
        X = np.asarray(b["bx"]).reshape(-1, V) * b["x_scale"] + b["x_offset"]
        Y = np.asarray(b["by"]).reshape(-1, V) * b["y_scale"] + b["y_offset"]
        _bnd_cache[sample] = (X, Y)
    return _bnd_cache[sample]

def cell_polygon(sample, idx):
    X, Y = all_polys(sample)
    return np.column_stack([X[idx], Y[idx]])

def poly_metrics(poly):
    """Return (circularity, solidity, area) for a polygon.
    circularity = 4*pi*A/P^2  (1 = circle; lower = elongated/jagged)
    solidity    = A / convex_hull_A  (1 = convex; lower = concave 'bites'
                  from neighbouring-cell segmentation)."""
    # drop consecutive duplicate (padding) vertices
    keep = np.insert(np.any(np.diff(poly, axis=0) != 0, axis=1), 0, True)
    u = poly[keep]
    if len(u) < 3:
        return 0.0, 0.0, 0.0
    x, y = u[:, 0], u[:, 1]
    area = 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
    per  = np.sum(np.hypot(np.diff(np.append(x, x[0])), np.diff(np.append(y, y[0]))))
    circ = 4 * np.pi * area / per**2 if per > 0 else 0.0
    try:
        sol = area / ConvexHull(u).volume        # 2D hull .volume == area
    except Exception:
        sol = 0.0
    return circ, sol, area

def marker_dots_in_poly(sample, gene, poly):
    gi = json.load(open(TXI.format(s=sample)))
    g  = json.load(open(TXG.format(s=sample, g=gene)))
    tx = np.asarray(g["x"]) * gi["x_scale"] + gi["x_offset"]
    ty = np.asarray(g["y"]) * gi["y_scale"] + gi["y_offset"]
    xmin, ymin = poly.min(0); xmax, ymax = poly.max(0)
    win = (tx >= xmin - 1) & (tx <= xmax + 1) & (ty >= ymin - 1) & (ty <= ymax + 1)
    pts = np.column_stack([tx[win], ty[win]])
    if len(pts) == 0:
        return pts
    inside = Path(poly).contains_points(pts)
    return pts[inside]

def group_target(sample, gene, subclass):
    """Median marker count for this subclass in this sample. If the group
    median is 0 (lowly-expressed gene, e.g. BDNF), fall back to the median
    among expressing cells so the exemplar still shows some molecules."""
    a = load(sample)
    sc = a.obs["subclass_label"].astype(str).values
    qc = a.obs["qc_pass"].values if "qc_pass" in a.obs else np.ones(a.n_obs, bool)
    vals = gene_counts(a, gene)[(sc == subclass) & qc]
    med = float(np.median(vals))
    if med < 1:                               # group median is 0
        pos = vals[vals > 0]
        med = float(np.median(pos)) if len(pos) else 1.0
    return med

def pick_cell(sample, gene, subclass, target):
    """Return obs index of a REPRESENTATIVE, well-segmented exemplar cell:
    marker count near the group target (>0), typical-sized, and as close to a
    circle/oval as possible (high circularity + convex) to avoid segmentation
    artifacts from neighbouring cells."""
    a = load(sample)
    sc = a.obs["subclass_label"].astype(str).values
    qc = a.obs["qc_pass"].values if "qc_pass" in a.obs else np.ones(a.n_obs, bool)
    tot = a.obs["total_counts"].values.astype(float)
    mk  = gene_counts(a, gene)
    mask = (sc == subclass) & qc & (mk > 0)            # require some expression
    lo, hi = np.percentile(tot[mask], [20, 90])         # typical-sized cells
    band = mask & (tot >= lo) & (tot <= hi)
    # pool of cells with marker count near the group target (representative)
    near = band & (mk >= 0.5 * target) & (mk <= 1.8 * target)
    pool = np.where(near)[0]
    if len(pool) < 5:
        pool = np.where(band)[0]
    # among the pool, pick the roundest convex cell
    best, best_circ = None, -1.0
    for i in pool:
        circ, sol, area = poly_metrics(cell_polygon(sample, i))
        if area < 40 or sol < 0.93:                     # require convex, sensible size
            continue
        if circ > best_circ:
            best, best_circ = i, circ
    if best is None:                                    # relax convexity if needed
        for i in pool:
            circ, sol, area = poly_metrics(cell_polygon(sample, i))
            if area >= 40 and circ > best_circ:
                best, best_circ = i, circ
    return best if best is not None else pool[0]

os.makedirs(OUT, exist_ok=True)
meta_rows = []
for gene, subclass in PAIRS:
    ct_tgt = group_target(CONTROL_SAMPLE, gene, subclass)
    sz_tgt = group_target(SCZ_SAMPLE, gene, subclass)
    print(f"[{gene}/{subclass}] group medians -> Control~{ct_tgt:.0f}, SCZ~{sz_tgt:.0f}")
    ci = pick_cell(CONTROL_SAMPLE, gene, subclass, ct_tgt)
    si = pick_cell(SCZ_SAMPLE, gene, subclass, sz_tgt)

    for sample, idx, dx in [(CONTROL_SAMPLE, ci, "Control"), (SCZ_SAMPLE, si, "SCZ")]:
        a = load(sample)
        poly = cell_polygon(sample, idx)
        ctr = poly.mean(0)                       # recentre on cell centroid
        polyc = poly - ctr
        dots = marker_dots_in_poly(sample, gene, poly)
        dotsc = (dots - ctr) if len(dots) else dots.reshape(0, 2)
        mk_count = int(gene_counts(a, gene)[idx])

        circ, sol, area = poly_metrics(poly)
        gslug = gene
        pd.DataFrame(polyc, columns=["x", "y"]).to_csv(
            f"{OUT}/exemplar_{gslug}_{dx}_boundary.csv", index=False)
        pd.DataFrame(dotsc, columns=["x", "y"]).to_csv(
            f"{OUT}/exemplar_{gslug}_{dx}_dots.csv", index=False)
        meta_rows.append(dict(gene=gene, subclass=subclass, dx=dx, sample=sample,
                              cell_index=int(idx), marker_count=mk_count,
                              n_dots_in_poly=len(dotsc),
                              total_counts=int(a.obs["total_counts"].values[idx]),
                              circularity=round(circ, 3), solidity=round(sol, 3),
                              area_um2=round(area, 1)))
        print(f"{gene:7} {subclass:8} {dx:7} {sample}: cell {idx}  "
              f"{gene}={mk_count}  dots={len(dotsc)}  total={meta_rows[-1]['total_counts']}  "
              f"circ={circ:.2f} sol={sol:.2f}")

pd.DataFrame(meta_rows).to_csv(f"{OUT}/exemplar_cells_meta.csv", index=False)
print(f"\nSaved exemplar CSVs + {OUT}/exemplar_cells_meta.csv")
