#!/usr/bin/env python3
"""
Zoomed cell-resolution view at an L6/WM boundary, comparing depth and layer
annotations across the three model versions:
  OLD       = predicted_norm_depth_prev   (K=50, 24-donor)
  NEW raw   = predicted_norm_depth_raw     (K=100, 6 low-CPS, unsmoothed)
  NEW smooth= predicted_norm_depth          (K=100 + post-hoc k=30 smoothing)

At full-section scale the smoothing is washed out; zoomed in, the raw
salt-and-pepper vs the smoothed coherent gradient is obvious.

Output: output/depth_validation/layer_compare/{sample}_zoom_boundary.png
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import anndata as ad

_HERE = os.path.dirname(os.path.abspath(__file__))
_SPATIAL = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_SPATIAL, "code", "analysis"))
sys.path.insert(0, os.path.join(_SPATIAL, "code"))
from config import H5AD_DIR
from modules.depth_model import (assign_discrete_layers, smooth_layers_spatial, LAYER_COLORS)

OUT = os.path.join(_SPATIAL, "output/depth_validation/layer_compare")
SAMPLE = "Br5588"
BOX = 600.0                # crop window side length (microns)
LAYER_ORDER = ["L1", "L2/3", "L4", "L5", "L6", "WM", "Vascular"]


def recon_old_layers(adata):
    o = adata.obs
    qc = (o["corr_qc_pass"].values.astype(bool) if "corr_qc_pass" in o else o["qc_pass"].values.astype(bool))
    dom = o["banksy_domain"].values.astype(str)[qc]; is_l1 = o["banksy_is_l1"].values.astype(bool)[qc]
    dold = o["predicted_norm_depth_prev"].values[qc]; coords = adata.obsm["spatial"][qc][:, :2]
    lay = assign_discrete_layers(dold); lay[dom == "Vascular"] = "Vascular"
    sm = smooth_layers_spatial(coords=coords, layers=lay, domains=dom, is_l1_banksy=is_l1, depths=dold, verbose=False)
    full = np.full(adata.n_obs, "Unassigned", dtype=object); full[np.where(qc)[0]] = sm
    return full


def main():
    a = ad.read_h5ad(os.path.join(H5AD_DIR, f"{SAMPLE}_annotated.h5ad"))
    xy = a.obsm["spatial"][:, :2]
    dprev = a.obs["predicted_norm_depth_prev"].values.astype(float)
    draw = a.obs["predicted_norm_depth_raw"].values.astype(float)
    dnew = a.obs["predicted_norm_depth"].values.astype(float)
    new_layer = a.obs["layer"].values.astype(str)
    old_layer = recon_old_layers(a)

    # pick a crop window centred on the L6/WM transition (depth ~0.8), densest area
    sel = np.isfinite(dnew) & (dnew > 0.75) & (dnew < 0.95)
    cx, cy = np.median(xy[sel, 0]), np.median(xy[sel, 1])
    win = (np.abs(xy[:, 0] - cx) < BOX / 2) & (np.abs(xy[:, 1] - cy) < BOX / 2)
    print(f"{SAMPLE}: crop centre ({cx:.0f},{cy:.0f}), {win.sum():,} cells in {BOX:.0f}µm box")

    fig, ax = plt.subplots(2, 3, figsize=(16, 11), facecolor="black")
    s = 26
    depths = [(dprev, "Depth OLD (K=50)"), (draw, "Depth NEW raw (K=100)"),
              (dnew, "Depth NEW smoothed (K=100+k30)")]
    for j, (dd, ttl) in enumerate(depths):
        m = win & np.isfinite(dd)
        sc = ax[0, j].scatter(xy[m, 0], xy[m, 1], c=dd[m], s=s, cmap="turbo", vmin=0.4, vmax=1.0,
                              edgecolors="none")
        ax[0, j].set_title(ttl, color="white", fontsize=14)
        cb = fig.colorbar(sc, ax=ax[0, j], fraction=0.045); cb.set_label("depth", color="white")
        plt.setp(cb.ax.get_yticklabels(), color="white"); cb.ax.yaxis.set_tick_params(color="white")
    # layers: old | new ; leave 3rd blank
    for j, (ll, ttl) in enumerate([(old_layer, "Layer OLD"), (new_layer, "Layer NEW (smoothed)")]):
        for lay in LAYER_ORDER:
            mm = win & (ll == lay)
            if mm.sum():
                ax[1, j].scatter(xy[mm, 0], xy[mm, 1], s=s, color=LAYER_COLORS.get(lay, (.5, .5, .5)),
                                 edgecolors="none", label=lay)
        ax[1, j].set_title(ttl, color="white", fontsize=14)
        ax[1, j].legend(markerscale=2, fontsize=9, loc="upper right", facecolor="black",
                        labelcolor="white", framealpha=0.4)
    ax[1, 2].axis("off")
    for a_ in ax.ravel()[:5]:
        a_.set_facecolor("black"); a_.set_aspect("equal"); a_.set_xticks([]); a_.set_yticks([])
    fig.suptitle(f"{SAMPLE}: zoom at L6/WM boundary ({BOX:.0f}µm) — raw salt-and-pepper vs smoothed",
                 color="white", fontsize=16)
    plt.tight_layout()
    op = os.path.join(OUT, f"{SAMPLE}_zoom_boundary.png")
    fig.savefig(op, dpi=150, facecolor="black", bbox_inches="tight"); plt.close(fig)
    print(f"  saved {op}")


if __name__ == "__main__":
    main()
