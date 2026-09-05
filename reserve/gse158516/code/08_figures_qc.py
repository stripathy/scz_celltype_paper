#!/usr/bin/env python3
"""
QC and annotation-validation figures for GSE158516.

Nothing downstream is trustworthy unless these look right, so they get drawn before
any disease comparison: does QC attrition track diagnosis (it must not), did Harmony
actually mix the samples, do the transferred labels agree with unsupervised structure,
and do the canonical markers land on the cells that carry the label.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from matplotlib.colors import LogNorm

plt.rcParams.update({"font.size": 15, "axes.labelsize": 17, "axes.titlesize": 19,
                     "xtick.labelsize": 13, "ytick.labelsize": 13, "legend.fontsize": 13})

BASE = "/Users/shreejoy/Github/scz_celltype_paper/reserve/gse158516"
DATA = "/Users/shreejoy/Github/shared_data/GSE158516"
OUTD = f"{BASE}/output"

DX_COLOR = {"Control": "#3B75AF", "Schizophrenia": "#D1615D", "not_in_supplement": "#999999"}

obs = pd.read_parquet(f"{DATA}/annotated/query_obs.parquet")
umap = np.load(f"{DATA}/annotated/query_umap.npy")
qc = pd.read_csv(f"{OUTD}/qc_per_sample.csv")
meta = pd.read_csv(f"{BASE}/data/sample_metadata.csv").rename(
    columns={"Deidentified ID": "sample"})
qc = qc.merge(meta[["sample", "Diagnosis", "Age", "PMI"]], on="sample", how="left")
qc["Diagnosis"] = qc["Diagnosis"].fillna("not_in_supplement")
obs["diagnosis"] = obs["diagnosis"].astype(str)

# ============================ Figure 1: QC ==================================
fig, axes = plt.subplots(2, 3, figsize=(22, 12))

ax = axes[0, 0]
w = 0.4
for i, (_, r) in enumerate(qc.sort_values("sample").iterrows()):
    ax.bar(i, r["n_called"], color="0.8", width=w * 2)
    ax.bar(i, r["n_pass"], color=DX_COLOR[r["Diagnosis"]], width=w * 2)
ax.set_xticks(range(len(qc)))
ax.set_xticklabels(qc.sort_values("sample")["sample"], rotation=90, fontsize=10)
ax.set_ylabel("nuclei")
ax.set_title("Nuclei called (grey) and passing QC (colour)")
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in DX_COLOR.values()]
ax.legend(handles, DX_COLOR.keys(), loc="upper right")

for ax, col, lab in [(axes[0, 1], "frac_removed", "fraction of called nuclei removed"),
                     (axes[0, 2], "frac_fail_mito", "fraction with mito >= 5%")]:
    for i, dx in enumerate(["Control", "Schizophrenia"]):
        v = qc[qc["Diagnosis"] == dx][col]
        ax.scatter(np.random.normal(i, 0.06, len(v)), v, s=90, alpha=0.85,
                   color=DX_COLOR[dx], edgecolor="k", zorder=3)
        ax.hlines(v.mean(), i - 0.2, i + 0.2, color="k", lw=3, zorder=4)
    v = qc[qc["Diagnosis"] == "not_in_supplement"][col]
    ax.scatter(np.random.normal(2, 0.06, len(v)), v, s=90, alpha=0.85,
               color=DX_COLOR["not_in_supplement"], edgecolor="k", zorder=3)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["Control", "SCZ", "dropped\nby authors"])
    ax.set_ylabel(lab)
    ax.set_title(lab.capitalize())

ax = axes[1, 0]
# xscale/yscale must be given to hexbin: it bins before any axis transform, so setting
# a log scale afterwards leaves every point in one giant bin.
hb = ax.hexbin(obs["num_counts"], obs["num_genes"], gridsize=90, bins="log",
               xscale="log", yscale="log", cmap="magma", mincnt=1)
plt.colorbar(hb, ax=ax, label="nuclei (log)")
ax.set_xlabel("UMIs per nucleus"); ax.set_ylabel("genes per nucleus")
ax.set_title("Post-QC depth")

ax = axes[1, 1]
ax.hist(obs["mito_fraction"] * 100, bins=60, color="0.3")
ax.set_xlabel("mitochondrial %"); ax.set_ylabel("nuclei")
ax.set_title("Post-QC mitochondrial fraction")

ax = axes[1, 2]
for lv, c in [("class", "#4C72B0"), ("subclass", "#55A868"), ("supertype", "#C44E52")]:
    ax.hist(obs[f"{lv}_confidence"], bins=40, histtype="step", lw=3, label=lv, color=c)
ax.set_xlabel("label-transfer confidence"); ax.set_ylabel("nuclei")
ax.set_title("Confidence by taxonomic level")
ax.legend()

fig.suptitle("GSE158516 (Reiner et al.) — quality control across 32 samples", fontsize=24)
fig.tight_layout()
fig.savefig(f"{OUTD}/fig_qc.png", dpi=110)
print("wrote fig_qc.png")

# ====================== Figure 2: annotation ================================
sub_order = obs["subclass"].value_counts().index.tolist()
palette = plt.cm.tab20(np.linspace(0, 1, 20)).tolist() + plt.cm.tab20b(np.linspace(0, 1, 20)).tolist()
sub_color = {s: palette[i % len(palette)] for i, s in enumerate(sub_order)}

fig, axes = plt.subplots(2, 2, figsize=(22, 20))
rng = np.random.default_rng(0)
o = rng.permutation(len(obs))                     # avoid draw-order artefacts

ax = axes[0, 0]
ax.scatter(umap[o, 0], umap[o, 1], s=0.6, alpha=0.5, linewidths=0,
           c=[sub_color[s] for s in obs["subclass"].to_numpy()[o]], rasterized=True)
for s in sub_order:
    m = obs["subclass"] == s
    if m.sum() > 300:
        ax.text(np.median(umap[m.to_numpy(), 0]), np.median(umap[m.to_numpy(), 1]), s,
                fontsize=13, weight="bold", ha="center",
                bbox=dict(fc="white", alpha=0.65, ec="none", pad=1.2))
ax.set_title("Transferred SEA-AD subclass")

ax = axes[0, 1]
samples = sorted(obs["sample_id"].unique())
scol = {s: plt.cm.gist_ncar(i / len(samples)) for i, s in enumerate(samples)}
ax.scatter(umap[o, 0], umap[o, 1], s=0.6, alpha=0.5, linewidths=0,
           c=[scol[s] for s in obs["sample_id"].to_numpy()[o]], rasterized=True)
ax.set_title("Sample (Harmony batch mixing check)")

ax = axes[1, 0]
ax.scatter(umap[o, 0], umap[o, 1], s=0.6, alpha=0.5, linewidths=0,
           c=[DX_COLOR.get(d, "#999999") for d in obs["diagnosis"].to_numpy()[o]],
           rasterized=True)
ax.set_title("Diagnosis")
ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=c, label=k, ms=12)
                   for k, c in DX_COLOR.items()], loc="best")

ax = axes[1, 1]
if "cluster" in obs.columns:
    ct = pd.crosstab(obs["cluster"], obs["subclass"], normalize="index")
    ct = ct.loc[:, [c for c in sub_order if c in ct.columns]]
    im = ax.imshow(ct.values, aspect="auto", cmap="magma", vmin=0, vmax=1)
    ax.set_xticks(range(ct.shape[1]))
    ax.set_xticklabels(ct.columns, rotation=90, fontsize=11)
    ax.set_yticks(range(ct.shape[0]))
    ax.set_yticklabels(ct.index, fontsize=8)
    ax.set_ylabel("de novo Leiden cluster")
    ax.set_title("Leiden cluster composition by transferred subclass")
    plt.colorbar(im, ax=ax, label="fraction of cluster")
    purity = ct.max(1)
    print(f"cluster purity: median {purity.median():.3f}, "
          f"{(purity > 0.8).mean():.1%} of clusters >80% one subclass")

for ax in axes.ravel()[:3]:
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
fig.suptitle("GSE158516 — SEA-AD label transfer and unsupervised structure", fontsize=24)
fig.tight_layout()
fig.savefig(f"{OUTD}/fig_annotation.png", dpi=110)
print("wrote fig_annotation.png")

# ==================== Figure 3: marker validation ===========================
MARKERS = {
    "Sst": ["SST"], "Pvalb": ["PVALB"], "Vip": ["VIP"], "Lamp5": ["LAMP5"],
    "Sncg": ["SNCG"], "Chandelier": ["UNC5B"], "L2/3 IT": ["CUX2", "RORB"],
    "L4 IT": ["RORB"], "L5 IT": ["THEMIS"], "L6 CT": ["SYT6"], "L6b": ["CTGF", "NR4A2"],
    "L5/6 NP": ["TSHZ2"], "Astrocyte": ["AQP4", "GFAP"], "Oligodendrocyte": ["PLP1", "MBP"],
    "OPC": ["PDGFRA"], "Microglia-PVM": ["P2RY12", "CSF1R"], "Endothelial": ["CLDN5"],
    "VLMC": ["COL1A2"],
}
genes = sorted({g for v in MARKERS.values() for g in v})

import anndata as ad
import os
files = sorted(f for f in os.listdir(f"{DATA}/per_sample_qc") if f.endswith(".h5ad"))
_v = ad.read_h5ad(f"{DATA}/per_sample_qc/{files[0]}", backed="r").var_names
present = [g for g in genes if g in _v]        # same gene space in every sample
mean_expr, frac_expr = {}, {}
for f in files:
    a = ad.read_h5ad(f"{DATA}/per_sample_qc/{f}")
    X = a[:, present].X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    tot = np.asarray(a.X.sum(1)).ravel()
    tot[tot == 0] = 1
    cp10k = X / tot[:, None] * 1e4
    lab = obs.set_index("_index").loc[a.obs_names, "subclass"].to_numpy()
    for s in np.unique(lab):
        m = lab == s
        mean_expr.setdefault(s, []).append((np.log1p(cp10k[m]).sum(0), m.sum()))
        frac_expr.setdefault(s, []).append(((X[m] > 0).sum(0), m.sum()))
    del a

rows_m, rows_f = {}, {}
for s in mean_expr:
    n = sum(c for _, c in mean_expr[s])
    rows_m[s] = sum(v for v, _ in mean_expr[s]) / n
    rows_f[s] = sum(v for v, _ in frac_expr[s]) / n
M = pd.DataFrame(rows_m, index=present).T
F = pd.DataFrame(rows_f, index=present).T
order = [s for s in MARKERS if s in M.index] + [s for s in M.index if s not in MARKERS]
M, F = M.loc[order], F.loc[order]
Mz = M.div(M.max(0).replace(0, 1), axis=1)

fig, ax = plt.subplots(figsize=(max(12, 0.75 * len(present)), 0.55 * len(order) + 3))
xs, ys = np.meshgrid(np.arange(M.shape[1]), np.arange(M.shape[0]))
sc = ax.scatter(xs.ravel(), ys.ravel(), s=(F.values.ravel() * 320) + 4,
                c=Mz.values.ravel(), cmap="Reds", vmin=0, vmax=1, edgecolor="0.6", lw=0.4)
ax.set_xticks(range(M.shape[1])); ax.set_xticklabels(M.columns, rotation=90)
ax.set_yticks(range(M.shape[0])); ax.set_yticklabels(M.index)
ax.set_ylim(-0.8, M.shape[0] - 0.2); ax.invert_yaxis()
plt.colorbar(sc, ax=ax, label="mean log1p CP10K\n(scaled per gene)")
ax.set_title("Canonical markers by transferred subclass\n(dot size = fraction of nuclei expressing)")
fig.tight_layout()
fig.savefig(f"{OUTD}/fig_markers.png", dpi=110)
print("wrote fig_markers.png")
M.to_csv(f"{OUTD}/marker_mean_expression.csv")
