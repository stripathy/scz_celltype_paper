#!/usr/bin/env python3
"""Full-section spatial maps coloured by the NEW deployed layer labels
(bottom-third K=100 + k=30 smoothing). A few representative sections in a grid."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import anndata as ad

_HERE = os.path.dirname(os.path.abspath(__file__))
_SPATIAL = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_SPATIAL, "code", "analysis"))
sys.path.insert(0, os.path.join(_SPATIAL, "code"))
from config import H5AD_DIR
from modules.depth_model import LAYER_COLORS

OUT = os.path.join(_SPATIAL, "output/depth_validation/layer_compare")
SAMPLES = ["Br8667", "Br6032", "Br5588", "Br5400", "Br2719", "Br5622"]
LAYER_ORDER = ["L1", "L2/3", "L4", "L5", "L6", "WM", "Vascular"]


def dx_of(obs):
    for c in ["diagnosis", "Diagnosis", "dx", "Dx", "group", "Group", "condition"]:
        if c in obs.columns:
            v = str(obs[c].iloc[0])
            return v
    return ""


def main():
    ncol = 3
    nrow = int(np.ceil(len(SAMPLES) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(6 * ncol, 6 * nrow), facecolor="black")
    axes = np.atleast_1d(axes).ravel()
    for ax, sid in zip(axes, SAMPLES):
        p = os.path.join(H5AD_DIR, f"{sid}_annotated.h5ad")
        if not os.path.exists(p):
            ax.axis("off"); continue
        a = ad.read_h5ad(p)
        xy = a.obsm["spatial"][:, :2]
        lay = a.obs["layer"].values.astype(str)
        dx = dx_of(a.obs)
        for L in LAYER_ORDER:
            m = lay == L
            if m.sum():
                ax.scatter(xy[m, 0], xy[m, 1], s=1.0, color=LAYER_COLORS.get(L, (.5, .5, .5)),
                           rasterized=True, linewidths=0)
        ttl = f"{sid}" + (f"  ({dx})" if dx else "")
        ax.set_title(ttl, color="white", fontsize=16)
        ax.set_facecolor("black"); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.invert_yaxis()
    for ax in axes[len(SAMPLES):]:
        ax.axis("off")
    handles = [Patch(facecolor=LAYER_COLORS.get(L, (.5, .5, .5)), label=L) for L in LAYER_ORDER]
    fig.legend(handles=handles, loc="lower center", ncol=len(LAYER_ORDER), fontsize=14,
               facecolor="black", labelcolor="white", framealpha=0.3, markerscale=2)
    fig.suptitle("New deployed cortical layer labels (bottom-third K=100 + k=30 smoothing)",
                 color="white", fontsize=20)
    plt.tight_layout(rect=[0, 0.04, 1, 0.97])
    op = os.path.join(OUT, "new_layers_full_sections.png")
    fig.savefig(op, dpi=150, facecolor="black", bbox_inches="tight")
    plt.close(fig)
    print(f"saved {op}")


if __name__ == "__main__":
    main()
