#!/usr/bin/env python3
"""
What adding GSE158516 does to the supertype compositional meta-analysis.

Three views: how far the pooled estimates move, how the significance ranking changes,
and a forest of the paper's eight headline supertypes.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.size": 15, "axes.labelsize": 17, "axes.titlesize": 18,
                     "xtick.labelsize": 13, "ytick.labelsize": 13, "legend.fontsize": 13})

BASE = "/Users/shreejoy/Github/scz_celltype_paper/reserve/gse158516"
OUTD = f"{BASE}/output"
HEADLINE = ["Sst_2", "Sst_3", "Sst_20", "Sst_22", "Sst_25", "L6b_1", "L6b_2", "L6b_4"]
C7, C8, CNEW = "0.35", "#C44E52", "#4C72B0"

r = pd.read_csv(f"{OUTD}/meta8_composition.csv")
fig, axes = plt.subplots(1, 3, figsize=(25, 8.5))

# ---- (a) how far do the pooled estimates move? ------------------------------
ax = axes[0]
is_h = r["cell_type"].isin(HEADLINE)
ax.scatter(r.loc[~is_h, "beta7"], r.loc[~is_h, "beta_fe"], s=45, color="0.75",
           edgecolor="0.5", zorder=2, label="other supertypes")
ax.scatter(r.loc[is_h, "beta7"], r.loc[is_h, "beta_fe"], s=170, color=CNEW,
           edgecolor="k", zorder=3, label="headline supertypes")
for _, x in r[is_h].iterrows():
    ax.annotate(x["cell_type"], (x["beta7"], x["beta_fe"]), fontsize=12,
                xytext=(7, -3), textcoords="offset points")
lim = 1.1 * max(r[["beta7", "beta_fe"]].abs().max())
ax.plot([-lim, lim], [-lim, lim], "k--", lw=1.5, alpha=0.6)
ax.axhline(0, color="0.8", lw=1); ax.axvline(0, color="0.8", lw=1)
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
ax.set_xlabel("7-dataset pooled $\\beta$"); ax.set_ylabel("8-dataset pooled $\\beta$")
ax.set_title("Pooled effect sizes barely move")
ax.text(0.04, 0.95, f"median |shift| = {r['shift'].abs().median():.3f}\n"
                    f"GSE158516 weight = {100*r['weight8'].median():.1f}%",
        transform=ax.transAxes, va="top", fontsize=15)
ax.legend(loc="lower right")

# ---- (b) significance before vs after ---------------------------------------
ax = axes[1]
x, y = -np.log10(r["fdr7"]), -np.log10(r["fdr_fe"])
ax.scatter(x[~is_h], y[~is_h], s=45, color="0.75", edgecolor="0.5", zorder=2)
ax.scatter(x[is_h], y[is_h], s=170, color=CNEW, edgecolor="k", zorder=3)
for _, xx in r[is_h].iterrows():
    ax.annotate(xx["cell_type"], (-np.log10(xx["fdr7"]), -np.log10(xx["fdr_fe"])),
                fontsize=12, xytext=(7, -3), textcoords="offset points")
lim = 1.08 * max(x.max(), y.max())
ax.plot([0, lim], [0, lim], "k--", lw=1.5, alpha=0.6)
for t, c in [(0.05, "#C44E52"), (0.10, "#DD8452")]:
    ax.axhline(-np.log10(t), color=c, lw=1.3, ls=":")
    ax.axvline(-np.log10(t), color=c, lw=1.3, ls=":")
    ax.text(lim * 0.985, -np.log10(t) + 0.05, f"FDR {t:g}", color=c, ha="right", fontsize=12)
ax.set_xlim(0, lim); ax.set_ylim(0, lim)
ax.set_xlabel("$-\\log_{10}$ FDR, 7 datasets")
ax.set_ylabel("$-\\log_{10}$ FDR, 8 datasets")
ax.set_title("Above the diagonal = strengthened")

# ---- (c) forest of the headline supertypes ----------------------------------
ax = axes[2]
h = r.set_index("cell_type").loc[HEADLINE].reset_index()
y = np.arange(len(h))[::-1]
for off, bcol, scol, c, lab in [(0.26, "beta7", "se7", C7, "7 datasets"),
                                (0.00, "beta8_new", "se8_new", C8, "GSE158516 alone"),
                                (-0.26, "beta_fe", "se_fe", CNEW, "8 datasets")]:
    ax.errorbar(h[bcol], y + off, xerr=1.96 * h[scol], fmt="o", color=c, ms=10,
                lw=2.2, capsize=4, label=lab, zorder=3)
ax.axvline(0, color="k", lw=1.2)
ax.set_yticks(y); ax.set_yticklabels(h["cell_type"])
ax.set_ylim(-0.7, len(h) - 0.3)
ax.set_xlabel("compositional effect (CLR $\\beta$, SCZ vs control)")
ax.set_title("Headline supertypes: 95% CIs")
ax.legend(loc="upper left")
# reserve a clear strip on the right for the 8-dataset FDRs
lo, hi = ax.get_xlim()
ax.set_xlim(lo, hi + 0.42 * (hi - lo))
for i, (_, x_) in enumerate(h.iterrows()):
    ax.text(0.995, y[i] - 0.26, f"FDR {x_['fdr_fe']:.3g}", transform=ax.get_yaxis_transform(),
            ha="right", va="center", fontsize=12, color=CNEW)

fig.suptitle("Adding GSE158516 as an 8th dataset to the supertype compositional meta-analysis",
             fontsize=23)
fig.tight_layout()
fig.savefig(f"{OUTD}/fig_meta8.png", dpi=110)
print("wrote fig_meta8.png")
