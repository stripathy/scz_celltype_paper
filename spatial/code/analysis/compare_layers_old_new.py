#!/usr/bin/env python3
"""
Visual + quantitative comparison of layer annotations BEFORE vs AFTER the
depth-model re-deployment (K=50/24-donor -> K=100/6-bottom-third-donor).

Old depth is preserved per-cell as `predicted_norm_depth_prev`; old layers are
reconstructed with the SAME step-05 smoothing (assign_discrete_layers +
smooth_layers_spatial) using the existing BANKSY domains, so the comparison is
fair (smoothed-vs-smoothed). New layers come straight from the h5ad.

Smoothness metric (from DECISIONS.md): mean over cells of the local std of
depth among the 10 nearest spatial neighbours (lower = smoother).

Outputs: output/depth_validation/layer_compare/{sample}_layer_compare.png
         + prints per-sample roughness old vs new across all samples.
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import anndata as ad
from scipy.spatial import cKDTree

_HERE = os.path.dirname(os.path.abspath(__file__))
_SPATIAL = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_SPATIAL, "code", "analysis"))
sys.path.insert(0, os.path.join(_SPATIAL, "code"))
from config import H5AD_DIR
from modules.depth_model import assign_discrete_layers, smooth_layers_spatial, LAYER_COLORS

OUT = os.path.join(_SPATIAL, "output/depth_validation/layer_compare")
os.makedirs(OUT, exist_ok=True)
PLOT_SAMPLES = ["Br8667", "Br6032", "Br5588"]   # full-cortical-span sections
LAYER_ORDER = ["L1", "L2/3", "L4", "L5", "L6", "WM", "Vascular"]


def roughness(coords, depth):
    """Mean local std of depth among 10 nearest spatial neighbours."""
    ok = np.isfinite(depth)
    c, d = coords[ok], depth[ok]
    if len(d) < 20:
        return np.nan
    _, idx = cKDTree(c).query(c, k=11)         # self + 10
    return float(np.nanmean(np.nanstd(d[idx[:, 1:]], axis=1)))


def recon_old_layers(adata):
    """Reconstruct old smoothed layers from predicted_norm_depth_prev (step-05 logic)."""
    o = adata.obs
    qc = (o["corr_qc_pass"].values.astype(bool) if "corr_qc_pass" in o
          else o["qc_pass"].values.astype(bool))
    dom = o["banksy_domain"].values.astype(str)[qc]
    is_l1 = o["banksy_is_l1"].values.astype(bool)[qc]
    depth_old = o["predicted_norm_depth_prev"].values[qc]
    coords = adata.obsm["spatial"][qc][:, :2]
    lay = assign_discrete_layers(depth_old)
    lay[dom == "Vascular"] = "Vascular"
    sm = smooth_layers_spatial(coords=coords, layers=lay, domains=dom,
                               is_l1_banksy=is_l1, depths=depth_old, verbose=False)
    full = np.full(adata.n_obs, "Unassigned", dtype=object)
    full[np.where(qc)[0]] = sm
    return full


def main():
    # ---- quantitative roughness across ALL samples ----
    print(f"{'sample':10s} {'rough_old':>10s} {'rough_new':>10s} {'Δ%':>7s}")
    rough_rows = []
    for f in sorted(os.listdir(H5AD_DIR)):
        if not f.endswith("_annotated.h5ad"):
            continue
        sid = f.replace("_annotated.h5ad", "")
        a = ad.read_h5ad(os.path.join(H5AD_DIR, f), backed="r")
        if "predicted_norm_depth_prev" not in a.obs:
            continue
        coords = np.asarray(a.obsm["spatial"])[:, :2]
        ro = roughness(coords, a.obs["predicted_norm_depth_prev"].values.astype(float))
        rn = roughness(coords, a.obs["predicted_norm_depth"].values.astype(float))
        rough_rows.append((sid, ro, rn))
        print(f"{sid:10s} {ro:10.4f} {rn:10.4f} {100*(rn-ro)/ro:+6.1f}%")
    if rough_rows:
        ros = np.array([r[1] for r in rough_rows]); rns = np.array([r[2] for r in rough_rows])
        print(f"{'MEAN':10s} {ros.mean():10.4f} {rns.mean():10.4f} {100*(rns.mean()-ros.mean())/ros.mean():+6.1f}%")

    # ---- spatial visual comparison for representative samples ----
    for sid in PLOT_SAMPLES:
        p = os.path.join(H5AD_DIR, f"{sid}_annotated.h5ad")
        if not os.path.exists(p):
            continue
        a = ad.read_h5ad(p)
        coords = a.obsm["spatial"][:, :2]
        new_layer = a.obs["layer"].values.astype(str)
        old_layer = recon_old_layers(a)
        dnew = a.obs["predicted_norm_depth"].values.astype(float)
        dold = a.obs["predicted_norm_depth_prev"].values.astype(float)

        fig, ax = plt.subplots(2, 2, figsize=(13, 13), facecolor="black")
        s = 1.2
        # depth row
        for j, (dd, ttl) in enumerate([(dold, "Depth — OLD (K=50, 24-donor)"),
                                       (dnew, "Depth — NEW (K=100, 6 low-CPS)")]):
            m = np.isfinite(dd)
            sc = ax[0, j].scatter(coords[m, 0], coords[m, 1], c=dd[m], s=s, cmap="turbo",
                                  vmin=0, vmax=1, rasterized=True, linewidths=0)
            ax[0, j].set_title(ttl, color="white", fontsize=15)
            cb = fig.colorbar(sc, ax=ax[0, j], fraction=0.035); cb.set_label("norm depth", color="white")
            cb.ax.yaxis.set_tick_params(color="white"); plt.setp(cb.ax.get_yticklabels(), color="white")
        # layer row
        for j, (ll, ttl) in enumerate([(old_layer, "Layers — OLD"), (new_layer, "Layers — NEW")]):
            for lay in LAYER_ORDER:
                mm = ll == lay
                col = LAYER_COLORS.get(lay, (0.5, 0.5, 0.5))
                ax[1, j].scatter(coords[mm, 0], coords[mm, 1], s=s, color=col,
                                 rasterized=True, linewidths=0, label=f"{lay} ({mm.sum():,})")
            ax[1, j].set_title(ttl, color="white", fontsize=15)
            leg = ax[1, j].legend(markerscale=6, fontsize=8, loc="upper right",
                                  facecolor="black", labelcolor="white", framealpha=0.4)
        for a_ in ax.ravel():
            a_.set_facecolor("black"); a_.set_aspect("equal"); a_.set_xticks([]); a_.set_yticks([])
        fig.suptitle(f"{sid}: layer annotation comparison (pia→WM)", color="white", fontsize=17)
        plt.tight_layout()
        op = os.path.join(OUT, f"{sid}_layer_compare.png")
        fig.savefig(op, dpi=140, facecolor="black", bbox_inches="tight")
        plt.close(fig)
        print(f"  saved {op}")


if __name__ == "__main__":
    main()
