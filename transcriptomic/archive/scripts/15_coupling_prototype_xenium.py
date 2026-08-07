"""
15_coupling_prototype_xenium.py — PROTOTYPE of the DE x composition coupling plan
(notes/plan_de_composition_coupling.md), run on the Xenium data (24 donors, 12/12)
because the per-donor matrices are on hand. This validates the scoring machinery;
it is NOT a powered result (n = 12 SCZ -> wide CIs on the within-SCZ correlation).

Per donor we build two leave-one-donor-out (LOO) disease-alignment scores:
  S_DE   = mean over cell types of cosine( donor's pseudobulk deviation from
           controls , SCZ-vs-control DE direction )   [cell-type-specific, aggregated]
  S_comp = cosine( donor's CLR composition deviation from controls ,
           SCZ-vs-control compositional direction )
Both directions and the control centroid are recomputed leaving the scored donor
out (no circularity). Then: within-SCZ corr(S_DE, S_comp) is the question; within
controls is the negative control; pooled is reported for context (trivially +).

Inputs (per-donor, same canonical corr cells as crumblr + edgeR DE):
  ~/Github/SCZ_Xenium/output/de/pseudobulk_subclass.csv          celltype,donor,gene,count
  ~/Github/SCZ_Xenium/output/de/pseudobulk_subclass_samples.csv  celltype,donor,total,diagnosis,...
  ../spatial/output/crumblr/crumblr_input_subclass_corr.csv      celltype,donor,count,total,diagnosis
Outputs:
  results/tables/coupling_prototype_xenium_scores.csv
  results/figures/coupling_prototype_xenium.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

XEN = os.path.expanduser("~/Github/SCZ_Xenium/output/de")
CRUMBLR = "../spatial/output/crumblr/crumblr_input_subclass_corr.csv"
OUT_T, OUT_F = "results/tables", "results/figures"


def cosine(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(np.dot(a, b) / (na * nb)) if na > 0 and nb > 0 else np.nan


# ---- load ----
pb = pd.read_csv(f"{XEN}/pseudobulk_subclass.csv")                 # celltype, donor, gene, count
samp = pd.read_csv(f"{XEN}/pseudobulk_subclass_samples.csv")       # celltype, donor, total, diagnosis ...
cr = pd.read_csv(CRUMBLR)                                          # celltype, donor, count, total, diagnosis

dxmap = samp.drop_duplicates("donor").set_index("donor")["diagnosis"].to_dict()
cr_donors = set(cr["donor"].unique())
donors = sorted(d for d in dxmap if d in cr_donors)               # donors present in BOTH modalities
dropped = sorted(set(dxmap) - cr_donors)
SCZ = [d for d in donors if dxmap[d] == "SCZ"]
CTRL = [d for d in donors if dxmap[d] == "Control"]
tot = samp.set_index(["celltype", "donor"])["total"].to_dict()    # pseudobulk lib size per (type,donor)
print(f"{len(donors)} donors with both expr+composition ({len(SCZ)} SCZ / {len(CTRL)} Control)"
      + (f"; dropped (no composition): {dropped}" if dropped else ""))

# ---- expression: per-celltype logCPM matrix (donors x genes) ----
pb = pb[pb.apply(lambda r: (r.celltype, r.donor) in tot, axis=1)].copy()
pb["logcpm"] = np.log2(pb["count"] / pb.apply(lambda r: tot[(r.celltype, r.donor)], axis=1) * 1e6 + 1)
logcpm = {}                                                       # celltype -> (donors x genes) df
for c, g in pb.groupby("celltype"):
    M = g.pivot_table(index="donor", columns="gene", values="logcpm", fill_value=0.0)
    keep = (M > 1).mean(0) >= 0.5                                 # genes expressed in >=50% donors
    M = M.loc[:, keep]
    if M.shape[1] >= 10 and M.shape[0] >= 18:                     # enough genes + donors for a direction
        logcpm[c] = M
print(f"cell types used for S_DE: {len(logcpm)}  (>=10 expressed genes, >=18 donors)")

# ---- composition: CLR (donors x celltypes) ----
cw = cr.pivot_table(index="donor", columns="celltype", values="count", fill_value=0.0)
prop = (cw + 0.5) / (cw + 0.5).sum(1).values[:, None]             # +0.5 pseudocount, renormalise
clr = np.log(prop) - np.log(prop).mean(1).values[:, None]        # CLR

# ---- LOO scores ----
def loo_dir_center(M, d):
    """control centroid and SCZ-vs-control direction over M (donors x feat), excluding donor d."""
    co = [x for x in CTRL if x != d and x in M.index]
    so = [x for x in SCZ if x != d and x in M.index]
    mu = M.loc[co].mean(0)
    return mu, M.loc[so].mean(0) - mu

rows = []
for d in donors:
    # S_DE: mean cosine over cell types
    cs = []
    for c, M in logcpm.items():
        if d not in M.index:
            continue
        mu, beta = loo_dir_center(M, d)
        cs.append(cosine(M.loc[d] - mu, beta))
    s_de = float(np.nanmean(cs)) if cs else np.nan
    # S_comp: single cosine over cell types
    mu, gamma = loo_dir_center(clr, d)
    s_comp = cosine(clr.loc[d] - mu, gamma)
    rows.append(dict(donor=d, diagnosis=dxmap[d], n_celltypes=len(cs), S_DE=s_de, S_comp=s_comp))

sc = pd.DataFrame(rows)
os.makedirs(OUT_T, exist_ok=True)
sc.to_csv(f"{OUT_T}/coupling_prototype_xenium_scores.csv", index=False)

# ---- stats ----
def report(label, sub):
    if len(sub) < 4:
        return
    r, p = stats.pearsonr(sub.S_DE, sub.S_comp)
    rs, ps = stats.spearmanr(sub.S_DE, sub.S_comp)
    print(f"  {label:16s} n={len(sub):2d}  Pearson r={r:+.2f} (p={p:.2f})  Spearman={rs:+.2f} (p={ps:.2f})")

print("\n=== sanity: score means by group (expect SCZ > Control on both) ===")
print(sc.groupby("diagnosis")[["S_DE", "S_comp"]].mean().round(3).to_string())
print("\n=== coupling: corr(S_DE, S_comp) ===")
report("within SCZ", sc[sc.diagnosis == "SCZ"])
report("within Control", sc[sc.diagnosis == "Control"])
report("pooled (context)", sc)

# ---- scatter ----
fig, ax = plt.subplots(figsize=(4.6, 4.4))
col = {"Control": "#0072B2", "SCZ": "#D55E00"}
for dx_, g in sc.groupby("diagnosis"):
    ax.scatter(g.S_DE, g.S_comp, c=col[dx_], s=46, edgecolor="grey25" if False else "0.25",
               linewidth=0.4, alpha=0.9, label=dx_)
scz = sc[sc.diagnosis == "SCZ"]
if len(scz) >= 4:                                                 # within-SCZ fit (the question)
    b1, b0 = np.polyfit(scz.S_DE, scz.S_comp, 1)
    xs = np.linspace(scz.S_DE.min(), scz.S_DE.max(), 50)
    ax.plot(xs, b0 + b1 * xs, color=col["SCZ"], lw=1.3, ls="--")
    r, p = stats.pearsonr(scz.S_DE, scz.S_comp)
    ax.text(0.04, 0.96, f"within SCZ: r = {r:+.2f}, p = {p:.2f}  (n={len(scz)})",
            transform=ax.transAxes, va="top", ha="left", fontsize=9, color=col["SCZ"])
ax.axhline(0, color="0.8", lw=0.6); ax.axvline(0, color="0.8", lw=0.6)
ax.set_xlabel("transcriptional SCZ-alignment  (S_DE)", fontsize=10)
ax.set_ylabel("compositional SCZ-alignment  (S_comp)", fontsize=10)
ax.set_title(f"Xenium prototype: DE vs composition coupling\n"
             f"(per donor, LOO cosine scores; n={len(SCZ)} SCZ — methods check, underpowered)",
             fontsize=9)
ax.legend(frameon=False, fontsize=9, loc="lower right")
ax.tick_params(labelsize=8)
fig.tight_layout()
fig.savefig(f"{OUT_F}/coupling_prototype_xenium.png", dpi=200)
print(f"\nSaved {OUT_T}/coupling_prototype_xenium_scores.csv + {OUT_F}/coupling_prototype_xenium.png")
