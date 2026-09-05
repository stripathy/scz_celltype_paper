#!/usr/bin/env python3
"""
Result figures: GSE158516 versus the paper's snRNA-seq meta-analysis.

Panel layout mirrors the paper's own cross-platform comparisons (Fig. 2j for DE,
Fig. 3f for composition) so the two are read the same way.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

plt.rcParams.update({"font.size": 15, "axes.labelsize": 17, "axes.titlesize": 18,
                     "xtick.labelsize": 13, "ytick.labelsize": 13, "legend.fontsize": 13})

BASE = "/Users/shreejoy/Github/scz_celltype_paper/external_validation/gse158516"
OUTD = f"{BASE}/output"
DX = {"Control": "#3B75AF", "Schizophrenia": "#D1615D"}

fig, axes = plt.subplots(2, 3, figsize=(24, 15))

# ---- (a) transcriptome-wide DE concordance ---------------------------------
j = pd.read_csv(f"{OUTD}/compare_de_joined.csv.gz")
sig = j[j["meta_FDR"] < 0.10]
ax = axes[0, 0]
ax.scatter(j["meta_logFC"], j["logFC"], s=2, alpha=0.12, color="0.7", linewidths=0,
           rasterized=True, label=f"all pairs (n={len(j):,})")
ax.scatter(sig["meta_logFC"], sig["logFC"], s=9, alpha=0.5, color="#C44E52",
           linewidths=0, rasterized=True, label=f"meta FDR<0.10 (n={len(sig):,})")
lim = np.percentile(np.abs(np.r_[sig["meta_logFC"], sig["logFC"]]), 99.5)
xs = np.linspace(-lim, lim, 10)
ax.plot(xs, xs, "k--", lw=1.5, alpha=0.6)
b = np.polyfit(sig["meta_logFC"], sig["logFC"], 1)
ax.plot(xs, np.polyval(b, xs), color="#C44E52", lw=2.5)
r, p = stats.pearsonr(sig["meta_logFC"], sig["logFC"])
same = np.mean(np.sign(sig["meta_logFC"]) == np.sign(sig["logFC"]))
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
ax.axhline(0, color="0.8", lw=1); ax.axvline(0, color="0.8", lw=1)
ax.set_xlabel("meta-analytic log$_2$FC (7 datasets)")
ax.set_ylabel("GSE158516 log$_2$FC")
ax.set_title("Cell-type-specific DE")
ax.text(0.04, 0.95, f"$r$ = {r:.2f}, slope = {b[0]:.2f}\n{same:.0%} sign-concordant",
        transform=ax.transAxes, va="top", fontsize=15)
ax.legend(loc="lower right", markerscale=3)

# ---- (b) the named genes ----------------------------------------------------
nt = pd.read_csv(f"{OUTD}/compare_named_genes.csv")
ax = axes[0, 1]
for ct, c in [("Sst", "#55A868"), ("Pvalb", "#C44E52"), ("Chandelier", "#8172B2")]:
    s = nt[nt["cell_type"] == ct]
    ax.scatter(s["meta_logFC"], s["gse_logFC"], s=160, color=c, edgecolor="k",
               zorder=3, label=ct)
    for _, r_ in s.iterrows():
        ax.annotate(r_["gene"], (r_["meta_logFC"], r_["gse_logFC"]),
                    fontsize=11, xytext=(6, 4), textcoords="offset points")
lim = 1.15 * max(nt[["meta_logFC", "gse_logFC"]].abs().max())
ax.plot([-lim, lim], [-lim, lim], "k--", lw=1.5, alpha=0.6)
ax.axhline(0, color="0.8", lw=1); ax.axvline(0, color="0.8", lw=1)
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
ax.set_xlabel("meta-analytic log$_2$FC"); ax.set_ylabel("GSE158516 log$_2$FC")
ax.set_title(f"Genes named in the paper ({int(nt.same_sign.sum())}/{len(nt)} concordant)")
ax.legend(loc="lower right")

# ---- (c) SST / PVALB --------------------------------------------------------
mk = pd.read_csv(f"{OUTD}/compare_markers.csv")
ax = axes[0, 2]
y = np.arange(len(mk))[::-1]
ax.barh(y - 0.17, mk["meta_logFC"], height=0.32, color="0.35", label="meta-analysis (7 datasets)")
ax.barh(y + 0.17, mk["gse_logFC"], height=0.32, color="#D1615D", label="GSE158516")
for i, (_, r_) in enumerate(mk.iterrows()):
    ax.text(r_["meta_logFC"] - 0.02, y[i] - 0.17, f"FDR={r_['meta_FDR']:.3g}",
            ha="right", va="center", fontsize=11)
    ax.text(r_["gse_logFC"] - 0.02, y[i] + 0.17, f"p={r_['gse_P']:.2g}",
            ha="right", va="center", fontsize=11)
ax.set_yticks(y)
ax.set_yticklabels([f"{r_['gene']}\nin {r_['cell_type']}" for _, r_ in mk.iterrows()])
ax.axvline(0, color="k", lw=1)
ax.set_xlabel("log$_2$FC (SCZ vs control)")
ax.set_ylim(-0.6, len(mk) - 0.4)
ax.set_title("Canonical interneuron markers")
ax.legend(loc="lower left")

# ---- (d) compositional concordance -----------------------------------------
for ax, level, ttl in [(axes[1, 0], "subclass", "Composition — subclass"),
                       (axes[1, 1], "supertype", "Composition — supertype")]:
    x = pd.read_csv(f"{OUTD}/compare_composition_{level}.csv")
    for comp, c, mk_ in [("neuronal", "#8172B2", "o"), ("non-neuronal", "#CCB974", "s")]:
        s = x[x["compartment"] == comp]
        ax.scatter(s["meta_beta"], s["logFC"], s=110 if level == "subclass" else 55,
                   color=c, marker=mk_, edgecolor="k", alpha=0.85, label=comp, zorder=3)
    n = x[x["compartment"] == "neuronal"]
    r, p = stats.pearsonr(n["meta_beta"], n["logFC"])
    lim = 1.15 * max(x[["meta_beta", "logFC"]].abs().max())
    ax.plot([-lim, lim], [-lim, lim], "k--", lw=1.5, alpha=0.5)
    ax.axhline(0, color="0.8", lw=1); ax.axvline(0, color="0.8", lw=1)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("meta-analytic crumblr $\\beta$"); ax.set_ylabel("GSE158516 crumblr logFC")
    ax.set_title(ttl)
    ax.text(0.04, 0.95, f"neuronal $r$ = {r:.2f}\n$P$ = {p:.1e}  (n = {len(n)})",
            transform=ax.transAxes, va="top", fontsize=15)
    if level == "subclass":
        for _, r_ in x.iterrows():
            if abs(r_["meta_beta"]) > 0.1 or abs(r_["logFC"]) > 0.2:
                ax.annotate(r_["cell_type"], (r_["meta_beta"], r_["logFC"]),
                            fontsize=11, xytext=(6, 4), textcoords="offset points")
    ax.legend(loc="lower right")

# ---- (e) the 8 headline supertypes -----------------------------------------
h = pd.read_csv(f"{OUTD}/compare_headline_supertypes.csv")
ax = axes[1, 2]
y = np.arange(len(h))[::-1]
ax.barh(y - 0.17, h["meta_beta"], height=0.32, color="0.35", label="meta-analysis")
ax.barh(y + 0.17, h["gse_logFC"], height=0.32, color="#D1615D", label="GSE158516")
ax.set_yticks(y); ax.set_yticklabels(h["supertype"])
ax.axvline(0, color="k", lw=1)
ax.set_xlabel("compositional effect (SCZ vs control)")
ax.set_title(f"Headline supertypes ({int(h.same_sign.sum())}/{len(h)} concordant)")
ax.legend(loc="lower right")
ax.set_ylim(-0.6, len(h) - 0.4)
xr = ax.get_xlim()
for i, (_, r_) in enumerate(h.iterrows()):
    off = 0.012 * (xr[1] - xr[0])
    neg = r_["gse_logFC"] < 0
    ax.text(r_["gse_logFC"] + (off if neg else -off), y[i] + 0.17,
            f"p={r_['gse_P']:.2f}", ha="left" if neg else "right",
            va="center", fontsize=11)

fig.suptitle("GSE158516 (Reiner et al., n=26) vs the 7-dataset snRNA-seq meta-analysis",
             fontsize=24)
fig.tight_layout()
fig.savefig(f"{OUTD}/fig_results.png", dpi=110)
print("wrote fig_results.png")

# ---- per-donor proportions of the vulnerable Sst set ------------------------
counts = pd.read_csv(f"{OUTD}/counts_supertype.csv", index_col=0)
meta = pd.read_csv(f"{BASE}/data/sample_metadata.csv").set_index("Deidentified ID")
comp = (pd.read_csv(f"{OUTD}/celltype_compartment.csv")
        .drop_duplicates("cell_type").set_index("cell_type")["compartment"])
neuronal = [c for c in counts.columns if comp.get(c) == "neuronal"]
meta = meta.loc[counts.index]

fig, axes = plt.subplots(1, 4, figsize=(21, 6))
sets = [("Sst_25", ["Sst_25"]), ("Sst_2", ["Sst_2"]),
        ("vulnerable Sst\n(Sst_2/3/20/22/25)", ["Sst_2", "Sst_3", "Sst_20", "Sst_22", "Sst_25"]),
        ("L6b (all)", [c for c in counts.columns if c.startswith("L6b_")])]
for ax, (name, types) in zip(axes, sets):
    types = [t for t in types if t in counts.columns]
    prop = 100 * counts[types].sum(1) / counts[neuronal].sum(1)
    for i, dx in enumerate(["Control", "Schizophrenia"]):
        v = prop[meta["Diagnosis"] == dx]
        ax.scatter(np.random.normal(i, 0.07, len(v)), v, s=130, color=DX[dx],
                   edgecolor="k", alpha=0.9, zorder=3)
        ax.hlines(v.mean(), i - 0.25, i + 0.25, color="k", lw=3, zorder=4)
    a = prop[meta["Diagnosis"] == "Control"]; b = prop[meta["Diagnosis"] == "Schizophrenia"]
    t, p = stats.ttest_ind(a, b, equal_var=False)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Control", "SCZ"])
    ax.set_ylabel("% of neurons")
    ax.set_title(f"{name}\nWelch $P$ = {p:.2f}")
fig.suptitle("GSE158516 per-donor proportions (n = 14 control / 12 SCZ)", fontsize=22)
fig.tight_layout()
fig.savefig(f"{OUTD}/fig_proportions.png", dpi=110)
print("wrote fig_proportions.png")
