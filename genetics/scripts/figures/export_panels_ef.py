"""
Rebuild Figure 4 panels f and g with an Sst_3 exemplar inserted between the
existing Sst_25 (high sag, upper layer) and Sst_5 (no sag, deep) cells, so the
pair becomes a three-cell depth/sag series.

Sst_3 cell 1037461069 was chosen over the only other reconstructed Sst_3 cell
(1079573757) because it has the higher sag (0.357 vs 0.294) and the higher
assignment confidence, and its soma depth (1109 um) still sits between the
other two exemplars.

Target currents are set so the three cells receive a comparable voltage
deflection rather than a comparable current: the existing pair uses -100 pA at
117 MOhm and -30 pA at 400 MOhm, both ~-12 mV, so Sst_3 at 238 MOhm gets -50 pA.

NWB for the Sst_3 cell is not in the upstream cache; it was pulled from
DANDI dandiset 000636 (sub-1036289507_ses-1037460966_icephys.nwb).
"""
import sys
import numpy as np
import json

import pandas as pd

import os as _os

# Repo-relative, so the chain runs from any clone. External data that is not in
# the repo is still resolved by absolute path or an environment override below.
GEN = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
PAPER = _os.path.dirname(GEN)

W = f"{GEN}/results/intermediates"

UP = _os.environ.get("PATCHSEQ_REPO",
                     _os.path.expanduser("~/Github/human_int_patch_seq"))
OUT = f"{GEN}/results/figures/r_panels"

sys.path.insert(0, UP)
from patchseq_builder.morphology.download import parse_swc            # noqa: E402
from patchseq_builder.morphology.orientation import (                  # noqa: E402
    INVERTED_SPECIMEN_IDS, flip_swc_y)

# SEA-AD supertype palette (was read from the retired 501-type enrichment export)
colors = json.load(open(f"{GEN}/data/seaad_supertype_colors.json"))

CELLS = [
    # Sst_20 replaced the Sst_22 exemplar (907585117) on 2026-08-31 so the
    # series shows the supertype that is both most SCZ-enriched and depleted.
    # 1079568285 is the only reconstructed Sst_20 cell (16 Sst_25, 2 Sst_22
    # reconstructions exist); assignment confidence 0.86, the lowest in the
    # series, declared in the legend. Target current from its 214 MOhm input
    # resistance under the same ~-12 mV rule as the others.
    dict(specimen_id=1079568285, supertype="Sst_20", layer="L2",
         swc=f"{GEN}/data/patchseq/swc/1079568285_upright.swc",
         nwb=f"{GEN}/data/patchseq/nwb/1079568285.nwb",
         pia_dist_um=304.031759, sag=0.281601, target_pa=-55),
    dict(specimen_id=819770858, supertype="Sst_25", layer="L2",
         swc=f"{GEN}/data/patchseq/swc/819770858_upright.swc",
         nwb=f"{GEN}/data/patchseq/nwb/819770858.nwb",
         pia_dist_um=405.610339, sag=0.564403, target_pa=-100),
    dict(specimen_id=1037461069, supertype="Sst_3", layer="L3",
         swc=f"{GEN}/data/patchseq/swc/1037461069_upright.swc",
         nwb=f"{GEN}/data/patchseq/nwb/1037461069.nwb",
         pia_dist_um=1108.721055, sag=0.357344, target_pa=-50),
    dict(specimen_id=758996755, supertype="Sst_5", layer="L4",
         swc=f"{GEN}/data/patchseq/swc/758996755_upright.swc",
         nwb=f"{GEN}/data/patchseq/nwb/758996755.nwb",
         # kept at the value the committed panel uses (true value is 1460.08)
         pia_dist_um=1500.0, sag=0.001955, target_pa=-30),
    # second not-depleted exemplar: deepest cell in the series, near-zero sag,
    # and the only Sst_1 reconstruction with perfect assignment confidence.
    dict(specimen_id=797048104, supertype="Sst_1", layer="L5",
         swc=f"{GEN}/data/patchseq/swc/797048104_upright.swc",
         nwb=f"{GEN}/data/patchseq/nwb/797048104.nwb",
         pia_dist_um=1733.511976, sag=0.035902, target_pa=-30),
]
# SEA-AD assigns Sst_3 #f2ad49 and Sst_5 #e6a343 -- visually identical, which
# is unusable in a three-cell panel. Sst_3 is given a mid tone here so the
# series reads dark -> mid -> light with depth. Panels f/g only; every other
# panel keeps the SEA-AD palette.
# SEA-AD gives Sst_3 #f2ad49 / Sst_5 #e6a343 (nearly the same light orange),
# so Sst_3 and Sst_1 take intermediate tones and the series reads as an even
# dark -> light ramp with depth. Sst_20 keeps its SEA-AD colour (#885616),
# which already separates from Sst_25's #693d07. Panels e/f only.
PANEL_FG_COLORS = {"Sst_3": "#c8862e", "Sst_1": "#f5c77e"}
for c in CELLS:
    c["color"] = PANEL_FG_COLORS.get(c["supertype"],
                                     colors.get(c["supertype"], "#888888"))

# ---------------------------------------------------------------- panel F
segs, meta_F = [], []
for ci, cell in enumerate(CELLS):
    nodes = parse_swc(cell["swc"])
    if cell["specimen_id"] in INVERTED_SPECIMEN_IDS:
        nodes = flip_swc_y(nodes)
    soma = next((n for n in nodes.values() if n["type"] == 1), None)
    sx, sy = ((soma["x"], soma["y"]) if soma else
              (np.mean([n["x"] for n in nodes.values()]),
               np.mean([n["y"] for n in nodes.values()])))

    def to_um(n):
        return (n["x"] - sx, cell["pia_dist_um"] + (sy - n["y"]))

    for n in nodes.values():
        if n["parent"] < 0 or n["parent"] not in nodes:
            continue
        px, py = to_um(nodes[n["parent"]])
        cx, cy = to_um(n)
        segs.append(dict(cell=cell["supertype"], comp_type=n["type"],
                         x0=px, y0=py, x1=cx, y1=cy))
    meta_F.append(dict(specimen_id=cell["specimen_id"], supertype=cell["supertype"],
                       layer=cell["layer"], pia_dist_um=cell["pia_dist_um"],
                       color=cell["color"], soma_x_um=0.0,
                       soma_y_um=cell["pia_dist_um"], cell_order=ci))
pd.DataFrame(segs).to_csv(f"{OUT}/panel_F_morphology.csv", index=False)
pd.DataFrame(meta_F).to_csv(f"{OUT}/panel_F_meta.csv", index=False)
print(f"panel F: {len(segs):,} segments across {len(CELLS)} cells")

# ---------------------------------------------------------------- panel G
import pynwb                                                           # noqa: E402
from pynwb import NWBHDF5IO                                            # noqa: E402


def best_hyperpol_sweep(nwb_path, target_pa):
    cands = []
    with NWBHDF5IO(nwb_path, "r", load_namespaces=True) as nio:
        nwb = nio.read()
        for k, v in nwb.acquisition.items():
            if not isinstance(v, pynwb.icephys.CurrentClampSeries):
                continue
            sd = getattr(v, "stimulus_description", "") or ""
            if not any(t in sd for t in ("X1PS_SubThresh", "X3LP_Rheo", "X4PS_SupraThresh")):
                continue
            s = nwb.stimulus.get(k.replace("_AD0", "_DA0"))
            if s is None:
                continue
            v_mv = np.asarray(v.data[:]) * v.conversion * 1000.0
            i_pa = np.asarray(s.data[:]) * s.conversion * 1e12
            rate = float(v.rate)
            t = np.arange(len(v_mv)) / rate
            mask = np.abs(i_pa) > 5
            if not mask.any():
                continue
            edges = np.diff(mask.astype(int))
            starts = np.where(edges == 1)[0] + 1
            ends = np.where(edges == -1)[0] + 1
            if mask[0]:
                starts = np.r_[0, starts]
            if mask[-1]:
                ends = np.r_[ends, len(mask)]
            b = int(np.argmax(ends - starts))
            step_start, step_end = t[starts[b]], t[ends[b] - 1]
            mm = np.zeros_like(mask)
            mm[starts[b]:ends[b]] = True
            amp = float(np.round(i_pa[mm].mean()))
            if amp >= 0:
                continue
            t0 = max(0, int((step_start - 0.15) * rate))
            t1 = min(len(t), int((step_end + 0.30) * rate))
            v_seg, i_seg = v_mv[t0:t1], i_pa[t0:t1]
            req_end = int((step_end + 0.10 - t[t0]) * rate)
            zm = np.abs(v_seg) < 1e-9
            fz = int(np.argmax(zm)) if zm.any() else -1
            if zm.any() and fz < req_end:
                continue
            if zm.any():
                v_seg, i_seg = v_seg[:fz], i_seg[:fz]
            cands.append(dict(step_amp_pa=amp, t=t[t0:t0 + len(v_seg)] - step_start,
                              v=v_seg, i=i_seg))
    return min(cands, key=lambda c: abs(c["step_amp_pa"] - target_pa)) if cands else None


traces, meta_G = [], []
for cell in CELLS:
    sw = best_hyperpol_sweep(cell["nwb"], cell["target_pa"])
    if sw is None:
        print(f"  !! no sweep for {cell['supertype']}")
        continue
    t_ms = sw["t"] * 1000.0
    traces.append(pd.DataFrame({"cell": cell["supertype"], "t_ms": t_ms,
                                "v_mV": sw["v"], "i_pA": sw["i"]}))
    base = sw["v"][t_ms < 0].mean()
    instep = (t_ms >= 0) & (t_ms <= 1000)
    pk_i = int(np.argmin(sw["v"][instep]))
    v_peak = sw["v"][instep][pk_i]
    late = (t_ms > 800) & (t_ms <= 1000)
    v_steady = sw["v"][late].mean()
    meta_G.append(dict(specimen_id=cell["specimen_id"], supertype=cell["supertype"],
                       layer=cell["layer"], sag=cell["sag"], color=cell["color"],
                       step_amp_pa=sw["step_amp_pa"], V_baseline=base,
                       V_peak=v_peak, V_peak_t_ms=t_ms[instep][pk_i],
                       V_steady=v_steady))
    print(f"  {cell['supertype']:6s} step {sw['step_amp_pa']:+.0f} pA  "
          f"base {base:.1f}  peak {v_peak:.1f}  steady {v_steady:.1f} mV")
pd.concat(traces, ignore_index=True).to_csv(f"{OUT}/panel_G_traces.csv", index=False)
pd.DataFrame(meta_G).to_csv(f"{OUT}/panel_G_meta.csv", index=False)
print("panel G written")
