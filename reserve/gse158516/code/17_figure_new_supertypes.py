#!/usr/bin/env python3
"""
Which supertypes change status with the 8th dataset, and where they sit in cortex.

The laminar panel is the point: the paper's claim is that Sst depletion is concentrated
in upper-layer supertypes, and the depth-vs-effect gradient gets *stronger* with
GSE158516 added -- which is not what adding noise to a spurious gradient would do.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

plt.rcParams.update({"font.size": 15, "axes.labelsize": 17, "axes.titlesize": 18,
                     "xtick.labelsize": 13, "ytick.labelsize": 13, "legend.fontsize": 13})

REPO = "/Users/shreejoy/Github/scz_celltype_paper"
OUTD = f"{REPO}/reserve/gse158516/output"

r = pd.read_csv(f"{OUTD}/meta8_composition.csv")
d = pd.read_csv(f"{REPO}/spatial/output/depth_platform/supertype_depth_platform_summary.csv")
d["depth"] = d[["median_MERFISH", "median_Xenium"]].mean(1)
# left join: the depth table covers neuronal supertypes only, and some status changes
# are non-neuronal (Astro_3, Endo_1), which must still appear in panel (b)
m = r.merge(d[["supertype", "subclass", "depth"]], left_on="cell_type",
            right_on="supertype", how="left")

GAIN10 = ["Sst_20", "L5/6 NP_4"]
GAIN20 = ["Sst_11", "L5/6 NP_2", "Pvalb_6"]
LOST = ["Sst_3", "Pvalb_14", "Astro_3", "Endo_1", "L2/3 IT_7"]

fig, axes = plt.subplots(1, 2, figsize=(20, 8.5))

# ---- (a) Sst depth vs compositional effect, 7 vs 8 datasets -----------------
ax = axes[0]
sst = m[m["subclass"] == "Sst"].sort_values("depth")
for col, c, lab in [("beta7", "0.45", "7 datasets"), ("beta_fe", "#C44E52", "8 datasets")]:
    rho, p = stats.spearmanr(sst["depth"], sst[col])
    ax.scatter(sst["depth"], sst[col], s=150, color=c, edgecolor="k", zorder=3,
               label=f"{lab}  ($\\rho$ = {rho:+.2f}, $p$ = {p:.4f})")
    z = np.polyfit(sst["depth"], sst[col], 1)
    xs = np.linspace(sst["depth"].min(), sst["depth"].max(), 10)
    ax.plot(xs, np.polyval(z, xs), color=c, lw=2.5, alpha=0.85)
for _, x in sst.iterrows():
    if x["fdr_fe"] < 0.20 or x["cell_type"] in GAIN20:
        ax.annotate(x["cell_type"], (x["depth"], x["beta_fe"]), fontsize=12,
                    xytext=(7, -4), textcoords="offset points")
ax.axhline(0, color="0.7", lw=1.2)
ax.axvline(0.35, color="#4C72B0", lw=1.5, ls="--")
ax.text(0.355, ax.get_ylim()[1] * 0.62, "upper layer\n(depth < 0.35)", color="#4C72B0",
        ha="left", va="top", fontsize=13)
ax.set_xlabel("cortical depth (0 = pia, 1 = white matter)")
ax.set_ylabel("compositional effect (CLR $\\beta$)")
ax.set_title("Sst depletion is laminar — and the gradient sharpens\nwith the 8th dataset (n = 16 Sst supertypes)")
ax.legend(loc="upper left")

# ---- (b) status changes ------------------------------------------------------
ax = axes[1]
ch = m[m["cell_type"].isin(GAIN10 + GAIN20 + LOST)].copy()
ch["grp"] = np.where(ch["cell_type"].isin(GAIN10), "gained FDR<0.10",
                     np.where(ch["cell_type"].isin(GAIN20), "gained FDR<0.20", "lost"))
ch = ch.sort_values(["grp", "beta_fe"])
y = np.arange(len(ch))[::-1]
ax.errorbar(ch["beta7"], y + 0.19, xerr=1.96 * ch["se7"], fmt="o", color="0.45",
            ms=9, lw=2, capsize=3, label="7 datasets", zorder=3)
ax.errorbar(ch["beta_fe"], y - 0.19, xerr=1.96 * ch["se_fe"], fmt="o", color="#C44E52",
            ms=9, lw=2, capsize=3, label="8 datasets", zorder=3)
for i, (_, x) in enumerate(ch.iterrows()):
    ax.scatter(x["beta8_new"], y[i], marker="^", s=110, color="#4C72B0",
               edgecolor="k", zorder=4, label="GSE158516 alone" if i == 0 else None)
ax.axvline(0, color="k", lw=1.2)
ax.set_yticks(y)
ax.set_yticklabels([f"{c}  (d={dd:.2f})" if np.isfinite(dd) else f"{c}  (non-neuronal)"
                    for c, dd in zip(ch["cell_type"], ch["depth"])])
ax.set_ylim(-0.8, len(ch) - 0.2)
for g, c in [("gained FDR<0.10", "#55A868"), ("gained FDR<0.20", "#8172B2"), ("lost", "#DD8452")]:
    idx = np.where(ch["grp"].to_numpy() == g)[0]
    if len(idx):
        ax.axhspan(y[idx].min() - 0.5, y[idx].max() + 0.5, color=c, alpha=0.11, zorder=0)
        ax.text(0.985, y[idx].mean(), g, transform=ax.get_yaxis_transform(), ha="right",
                va="center", color=c, fontsize=13, weight="bold")
ax.set_xlabel("compositional effect (CLR $\\beta$)")
ax.set_title("Supertypes whose significance status changes\n(d = cortical depth)")
ax.legend(loc="lower left")

fig.suptitle("New and lost compositional signals after adding GSE158516", fontsize=23)
fig.tight_layout()
fig.savefig(f"{OUTD}/fig_new_supertypes.png", dpi=110)
print("wrote fig_new_supertypes.png")

ch[["cell_type", "subclass", "depth", "beta7", "fdr7", "beta8_new", "p8_new",
    "beta_fe", "fdr_fe", "fdr_dl", "concordant", "grp"]].to_csv(
    f"{OUTD}/status_changes.csv", index=False)
print("wrote status_changes.csv")
