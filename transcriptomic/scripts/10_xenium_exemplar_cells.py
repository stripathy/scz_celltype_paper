"""
Extract exemplar Xenium cells illustrating marker downregulation in SCZ, for the
exemplar panel of the composite (scripts/09_composite_figure.R):
  - SST   in Sst   cells: one Control vs one SCZ
  - FGFR3 in Astrocytes : one Control vs one SCZ

Each exemplar is chosen to be REPRESENTATIVE — its raw marker count is as close
as possible to the POOLED group median for that diagnosis (median across all
qc-pass cells of the subclass in all 24 Xenium donors, not one section). Among
cells that hit the group median we then require a typical, capture-matched
library size and pick the roundest, most convex cell (clean segmentation for a
representative image). Cortical depth/layer is reported but NOT used for
selection (circular morphology is prioritised over laminar matching).

For each chosen cell we save (micron coords, recentred on the cell centroid):
  - the cell boundary polygon            exemplar_<gene>_<dx>_boundary.csv
  - the nucleus boundary polygon         exemplar_<gene>_<dx>_nucleus.csv
  - the marker-gene transcript molecules INSIDE the cell polygon
                                         exemplar_<gene>_<dx>_dots.csv
plus a metadata row (exemplar_cells_meta.csv).

DATA PROVENANCE (SCZ_Xenium repo, Kwon 2026):
  output/h5ad/<sample>_annotated.h5ad      cells x 300 genes; .X = raw counts;
        obs.subclass_label, obs.qc_pass, obs.total_counts, obs.predicted_norm_depth,
        obs.layer.
  output/deploy/boundaries/<sample>.json          cell polygons (25 verts), obs order
  output/deploy/boundaries/<sample>_nucleus.json  nucleus polygons, same indexing
  output/deploy/transcripts/<sample>/<GENE>.json (+ gene_index.json)  molecules
  Diagnosis per sample: SCZ_Xenium code/analysis/config.py SAMPLE_TO_DX.

Drawing sections are chosen PER PAIR (see PAIRS) as the donor whose canonical
grain-density median is closest to the diagnosis group median: SST uses
Br6432 / Br5973; FGFR3 uses Br5400 (control) / Br5973 (SCZ) — Br6432 is
atypically LOW in FGFR3, so a better-matched control section is used.

Output: results/tables/exemplar_*.csv
"""

import json
import os
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
from matplotlib.path import Path
from scipy.spatial import ConvexHull

XEN  = os.path.expanduser("~/Github/SCZ_Xenium")
H5AD = os.path.join(XEN, "output/h5ad/{s}_annotated.h5ad")
BND  = os.path.join(XEN, "output/deploy/boundaries/{s}.json")
NUC  = os.path.join(XEN, "output/deploy/boundaries/{s}_nucleus.json")  # same indexing as BND
TXG  = os.path.join(XEN, "output/deploy/transcripts/{s}/{g}.json")
TXI  = os.path.join(XEN, "output/deploy/transcripts/{s}/gene_index.json")
OUT  = "results/tables"

sys.path.insert(0, os.path.join(XEN, "code/analysis"))
from config import SAMPLE_TO_DX                       # noqa: E402

# (gene, subclass, {dx: drawing section}). Section per diagnosis = the donor whose
# canonical-cell grain-density median is closest to that diagnosis's pooled group
# median, chosen per pair. FGFR3 control = Br5400 (closest); Br6432 is atypically
# low in FGFR3. SCZ = Br5973 for both (closest exported SCZ section).
PAIRS = [
    ("SST",   "Sst",       {"Control": "Br6432", "SCZ": "Br5973"}),
    ("FGFR3", "Astrocyte",  {"Control": "Br5400", "SCZ": "Br5973"}),
    ("PVALB", "Pvalb",      {"Control": "Br6432", "SCZ": "Br5973"}),  # composite panel h (PVALB in Pvalb)
]
ALL_SAMPLES = sorted(SAMPLE_TO_DX)

# ── caches ──
_adata, _bnd, _nuc, _tx = {}, {}, {}, {}
def load(s):
    if s not in _adata:
        _adata[s] = ad.read_h5ad(H5AD.format(s=s))
    return _adata[s]

def xcount(a, gene):
    col = a.X[:, a.var_names.get_loc(gene)]
    return np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()

def _decode(path):
    b = json.load(open(path))
    V = b["verts_per_cell"]
    X = np.asarray(b["bx"]).reshape(-1, V) * b["x_scale"] + b["x_offset"]
    Y = np.asarray(b["by"]).reshape(-1, V) * b["y_scale"] + b["y_offset"]
    return X, Y

def all_polys(s):
    if s not in _bnd:
        _bnd[s] = _decode(BND.format(s=s))
    return _bnd[s]

def all_nuc_polys(s):
    if s not in _nuc:
        _nuc[s] = _decode(NUC.format(s=s))
        assert _nuc[s][0].shape[0] == all_polys(s)[0].shape[0], "nucleus/cell index mismatch"
    return _nuc[s]

def cell_polygon(s, idx):
    X, Y = all_polys(s);     return np.column_stack([X[idx], Y[idx]])
def nucleus_polygon(s, idx):
    X, Y = all_nuc_polys(s); return np.column_stack([X[idx], Y[idx]])

def poly_metrics(poly):
    """(circularity 4*pi*A/P^2, solidity A/hullA, area) after dropping padding verts."""
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

def poly_eccentricity(poly):
    """Eccentricity of the cell outline from its boundary-vertex covariance
    (0 = circle, ->1 = elongated). sqrt(1 - lambda_min/lambda_max)."""
    keep = np.insert(np.any(np.diff(poly, axis=0) != 0, axis=1), 0, True)
    u = poly[keep]
    if len(u) < 3:
        return 1.0
    c = u - u.mean(0)
    ev = np.sort(np.linalg.eigvalsh((c.T @ c) / len(c)))   # [lambda_min, lambda_max]
    return float(np.sqrt(max(0.0, 1.0 - ev[0] / ev[1]))) if ev[1] > 0 else 1.0

def gene_tx(s, gene):
    """Cached, decoded transcript-molecule coords (micron) for one gene/section."""
    if (s, gene) not in _tx:
        gi = json.load(open(TXI.format(s=s)))
        g  = json.load(open(TXG.format(s=s, g=gene)))
        tx = np.asarray(g["x"]) * gi["x_scale"] + gi["x_offset"]
        ty = np.asarray(g["y"]) * gi["y_scale"] + gi["y_offset"]
        _tx[(s, gene)] = (tx, ty)
    return _tx[(s, gene)]

def marker_dots_in_poly(s, gene, poly):
    """Marker molecules geometrically inside the cell polygon (what panel K draws)."""
    tx, ty = gene_tx(s, gene)
    xmn, ymn = poly.min(0); xmx, ymx = poly.max(0)
    w = (tx >= xmn - 1) & (tx <= xmx + 1) & (ty >= ymn - 1) & (ty <= ymx + 1)
    pts = np.column_stack([tx[w], ty[w]])
    return pts[Path(poly).contains_points(pts)] if len(pts) else pts.reshape(0, 2)


def canonical_mask(a, subclass):
    """Canonical cell set = corr_subclass + cortical + (qc_pass & corr_qc_pass),
    matching the crumblr compositional analysis and the edgeR DE (load_cells /
    load_sample_adata, qc_mode='corr'). corr_qc_pass is a subset of qc_pass."""
    o = a.obs
    sc = (o["corr_subclass"] if "corr_subclass" in o else o["subclass_label"]).astype(str).values
    cort = (o["spatial_domain"].astype(str).values == "Cortical")
    qc = o["qc_pass"].values.astype(bool) if "qc_pass" in o else np.ones(a.n_obs, bool)
    cq = o["corr_qc_pass"].values.astype(bool) if "corr_qc_pass" in o else np.ones(a.n_obs, bool)
    return (sc == subclass) & cort & qc & cq


def cell_areas(s):
    """Vectorised cell-polygon area (um^2) per cell; index = obs order (all_polys)."""
    X, Y = all_polys(s)
    return 0.5 * np.abs((X * np.roll(Y, -1, axis=1) - np.roll(X, -1, axis=1) * Y).sum(axis=1))

def pooled_grain_median(gene, subclass, dx):
    """Median GRAIN DENSITY (marker dots per cell area, grains/100 um^2) across ALL
    `dx` donors on the canonical cell set -> the Dienel-style group target the
    exemplar aims at (cell area from the deploy boundary polygons)."""
    vals = []
    for s in ALL_SAMPLES:
        if SAMPLE_TO_DX[s] != dx:
            continue
        a = load(s)
        if gene not in a.var_names:
            continue
        m = canonical_mask(a, subclass)
        gd = xcount(a, gene)[m] / cell_areas(s)[m] * 100.0
        vals.append(gd[np.isfinite(gd)])
    v = np.concatenate(vals)
    return float(np.median(v)), len(v)

def shape_anchors(subclass, sections):
    """Shared typical cell SHAPE over the pair's two drawing sections' canonical
    cells: median area (size-match, so displayed dots reflect density not size)
    and median eccentricity (so exemplars have a representative outline rather
    than being cherry-picked as the roundest cell)."""
    areas, eccs = [], []
    for s in sections.values():
        a = load(s)
        m = canonical_mask(a, subclass)
        areas.append(cell_areas(s)[m])
        for i in np.where(m)[0]:
            eccs.append(poly_eccentricity(cell_polygon(s, i)))
    return float(np.median(np.concatenate(areas))), float(np.median(eccs))

def pick_cell(section, gene, subclass, gd_target, area_anc, ecc_target):
    """Obs index of a REPRESENTATIVE exemplar in `section`: among canonical,
    size-matched, convex (non-jagged) cells, the cell whose GRAIN DENSITY
    (grains/100 um^2) is closest to the group target AND whose ECCENTRICITY is
    closest to the median outline (a typical-shaped cell, not the roundest).
    Stage 1 ranks by grain-density band then eccentricity-near-median; stage 2
    picks, among the best, the cell whose DISPLAYED in-polygon grain density best
    matches the target (so drawn dots are representative and Control>SCZ holds)."""
    a = load(section)
    tot = a.obs["total_counts"].values.astype(float)
    dep = a.obs["predicted_norm_depth"].values.astype(float)
    lay = a.obs["layer"].astype(str).values
    mk  = xcount(a, gene)
    cand = np.where(canonical_mask(a, subclass))[0]
    a_lo, a_hi = 0.8 * area_anc, 1.2 * area_anc            # size-matched to the shared anchor

    def collect(sol_min):
        out = []
        for i in cand:
            circ, sol, area = poly_metrics(cell_polygon(section, i))
            if area < a_lo or area > a_hi or sol < sol_min:   # size-matched + convex (no jagged bites)
                continue
            ecc = poly_eccentricity(cell_polygon(section, i))
            gd  = mk[i] / area * 100.0                          # assigned-count grain density
            gd_band = int(abs(gd - gd_target) // max(1e-6, 0.05 * gd_target))  # within 5% bands
            out.append((gd_band, abs(ecc - ecc_target), i, circ, sol, area, ecc))
        return out
    scored = collect(0.93) or collect(0.90) or collect(0.0)
    scored.sort()
    # stage 2: among the best (grain-density band, then eccentricity-near-median),
    # pick the cell whose DISPLAYED in-polygon grain density best matches the
    # target, tie-broken again by eccentricity-near-median.
    best = None
    for _, _, i, circ, sol, area, ecc in scored[:50]:
        nd = len(marker_dots_in_poly(section, gene, cell_polygon(section, i)))
        key = (abs(nd / area * 100.0 - gd_target), abs(ecc - ecc_target))
        if best is None or key < best[0]:
            best = (key, i, circ, sol, area, ecc, nd)
    _, i, circ, sol, area, ecc, nd = best
    return i, dict(x_count=int(mk[i]), n_dots=int(nd), total_counts=int(tot[i]),
                   grain_density=round(mk[i] / area * 100.0, 2),
                   disp_grain_density=round(nd / area * 100.0, 2),
                   eccentricity=round(ecc, 3), norm_depth=round(float(dep[i]), 3), layer=lay[i],
                   circularity=round(circ, 3), solidity=round(sol, 3), area_um2=round(area, 1))


os.makedirs(OUT, exist_ok=True)
meta_rows = []
for gene, subclass, sections in PAIRS:
    area_anc, ecc_tgt = shape_anchors(subclass, sections)
    for dx in ["Control", "SCZ"]:
        gd_tgt, n = pooled_grain_median(gene, subclass, dx)
        section = sections[dx]
        idx, diag = pick_cell(section, gene, subclass, gd_tgt, area_anc, ecc_tgt)
        print(f"[{gene}/{subclass} {dx}] gd target={gd_tgt:.2f}/100um2, ecc target={ecc_tgt:.2f} (n={n}); "
              f"{section} cell {idx}: dots={diag['n_dots']} gd={diag['disp_grain_density']} "
              f"ecc={diag['eccentricity']} area={diag['area_um2']} depth={diag['norm_depth']} "
              f"layer={diag['layer']} circ={diag['circularity']}")

        poly = cell_polygon(section, idx)
        ctr  = poly.mean(0)                              # recentre everything on cell centroid
        nuc  = nucleus_polygon(section, idx)
        dots = marker_dots_in_poly(section, gene, poly)
        pd.DataFrame(poly - ctr, columns=["x", "y"]).to_csv(
            f"{OUT}/exemplar_{gene}_{dx}_boundary.csv", index=False)
        pd.DataFrame(nuc - ctr, columns=["x", "y"]).to_csv(
            f"{OUT}/exemplar_{gene}_{dx}_nucleus.csv", index=False)
        pd.DataFrame((dots - ctr) if len(dots) else dots, columns=["x", "y"]).to_csv(
            f"{OUT}/exemplar_{gene}_{dx}_dots.csv", index=False)
        meta_rows.append(dict(gene=gene, subclass=subclass, dx=dx, sample=section,
                              cell_index=int(idx), grain_density_target=round(gd_tgt, 2),
                              ecc_target=round(ecc_tgt, 3),
                              marker_count=diag["x_count"], n_dots_in_poly=len(dots),
                              grain_density=diag["grain_density"], disp_grain_density=diag["disp_grain_density"],
                              eccentricity=diag["eccentricity"], total_counts=diag["total_counts"],
                              norm_depth=diag["norm_depth"], layer=diag["layer"],
                              circularity=diag["circularity"], solidity=diag["solidity"],
                              area_um2=diag["area_um2"]))

pd.DataFrame(meta_rows).to_csv(f"{OUT}/exemplar_cells_meta.csv", index=False)
print(f"\nSaved exemplar CSVs (boundary/nucleus/dots) + {OUT}/exemplar_cells_meta.csv")
