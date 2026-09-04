#!/usr/bin/env python3
"""
APPROXIMATE "Xenium-realistic" degradation of the panel resolvability benchmark.

Exploratory (not for supplement): the clean snRNA panel F1 is optimistic because snRNA
lacks Xenium's spatial spillover/segmentation contamination. Calibration (see
_xenium_panel_depth.csv) showed Xenium actually has MORE panel-gene counts/cell than
snRNA (median 871 vs 482) — so depth is NOT the degrader; the Xenium-specific noise is
neighbor-transcript SPILLOVER. We model:

    observed_g ~ Poisson( D * [ (1-c)*p_true_g + c*p_ambient_g ] )

  - D            = per-cell total panel counts, drawn from the matched REAL Xenium
                   per-subclass depth distribution (calibrated) -> realistic counting noise
  - p_true       = the cell's own snRNA panel composition
  - p_ambient    = global mean panel composition (the spillover/ambient pool, dominated by
                   abundant neighbors) -> approximates segmentation bleed
  - c            = spillover fraction, SWEPT 0..0.3 (the key uncertain knob; ~0.1-0.3 is a
                   plausible Xenium range). c=0 isolates depth+Poisson (should ~= clean panel).

Then re-run the same correlation classifier + leave-one-donor-out CV; report subclass(24)
and Sst-supertype(16) F1 vs c, against the clean-panel and transcriptome references.

Caveats (flag, don't hide): ambient is GLOBAL not neighbor-structured (real spillover is
local, often same-layer); per-gene probe efficiency not modeled; for types where Xenium is
deeper than snRNA, Poisson "upsamples" a profile estimated from fewer counts. Approximation.

Outputs (output/depth_validation/lieber_layers/):
  report_degraded_sweep.csv, report_degraded_per_celltype.csv, report_degraded.png/.pdf
"""
import os, sys
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from sklearn.model_selection import GroupKFold
from sklearn.metrics import classification_report
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Helvetica","Arial","DejaVu Sans"],
    "font.size":11,"axes.titlesize":12.5,"axes.labelsize":11,"xtick.labelsize":10,"ytick.labelsize":10,
    "legend.fontsize":9,"axes.linewidth":0.6,"pdf.fonttype":42,"ps.fonttype":42})

BASE=os.path.expanduser("~/Github/SCZ_Xenium"); NICOLE=os.path.expanduser("~/Github/shared_data/nicole_sea_ad_snrnaseq_reference.h5ad")
_SP=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); OUT=os.path.join(_SP,"output/depth_validation/lieber_layers")
sys.path.insert(0,os.path.join(_SP,"code","modules")); from constants import SUBCLASS_TO_CLASS
CS=[0.0,0.1,0.2,0.3]; SEED=0
def dense(x): return x.toarray() if sp.issparse(x) else np.asarray(x)
def corr_dense(Xtr,ytr,Xte,cand=None):
    types=[t for t in (cand if cand is not None else sorted(set(ytr))) if (ytr==t).sum()>0]
    C=np.vstack([Xtr[ytr==t].mean(0) for t in types]); Cc=C-C.mean(1,keepdims=True); Cc/=(np.linalg.norm(Cc,axis=1,keepdims=True)+1e-9)
    Xc=Xte-Xte.mean(1,keepdims=True); Xc/=(np.linalg.norm(Xc,axis=1,keepdims=True)+1e-9)
    return np.array(types)[(Xc@Cc.T).argmax(1)]

# ── load snRNA raw panel counts ──
xh=os.path.join(BASE,"output/h5ad"); panel=list(ad.read_h5ad(os.path.join(xh,[f for f in sorted(os.listdir(xh)) if f.endswith("_annotated.h5ad")][0]),backed="r").var_names)
print("loading reference...",flush=True); a=ad.read_h5ad(NICOLE)
SUB=a.obs["Subclass"].astype(str).values; SUP=a.obs["Supertype"].astype(str).values; donor=a.obs["donor_id"].astype(str).values
pg=[g for g in panel if g in set(a.var_names)]; raw=dense(a[:,pg].X).astype(np.float32); del a
N=len(SUB); nd=len(set(donor)); subclasses=sorted(set(SUB)); sst=sorted([s for s in set(SUP) if s.startswith("Sst_")],key=lambda s:int(s.split("_")[1]))
tot=raw.sum(1); p_true=raw/np.where(tot[:,None]==0,1,tot[:,None]); ambient=p_true.mean(0)   # global ambient composition
print(f"N={N:,}, {nd} donors, panel={len(pg)}g, {len(subclasses)} subclasses, {len(sst)} Sst",flush=True)

# ── Xenium per-subclass depth targets (calibrated) ──
xd=pd.read_csv(os.path.join(OUT,"_xenium_panel_depth.csv")); xdep={s:g.total.values for s,g in xd.groupby("subclass")}
xall=xd.total.values; rng=np.random.default_rng(SEED)
Dtarget=np.array([rng.choice(xdep.get(SUB[i],xall)) for i in range(N)],dtype=np.float32)   # per-cell Xenium-like depth

def cv_f1(X):
    ps=np.empty(N,object); pt=np.empty(N,object)
    for tr,te in GroupKFold(nd).split(np.empty((N,1)),SUB,donor):
        ps[te]=corr_dense(X[tr],SUB[tr],X[te]); pt[te]=corr_dense(X[tr],SUP[tr],X[te])
    rs=classification_report(SUB,ps,labels=subclasses,output_dict=True,zero_division=0)
    rt=classification_report(SUP,pt,labels=sst,output_dict=True,zero_division=0)
    f1s=pd.Series({s:rs[s]["f1-score"] for s in subclasses}); f1t=pd.Series({s:rt[s]["f1-score"] for s in sst})
    return f1s,f1t
def lognorm(C):
    t=C.sum(1); X=C/np.where(t[:,None]==0,1,t[:,None])*1e4; return np.log1p(X).astype(np.float32)

# ── clean panel baseline (undegraded) ──
f1s0,f1t0=cv_f1(lognorm(raw)); print(f"clean panel: subclass med F1={f1s0.median():.3f}, Sst med={f1t0.median():.3f}, Sst_25={f1t0['Sst_25']:.3f}",flush=True)

# ── Xenium-degraded sweep over spillover fraction c ──
rows=[{"scenario":"clean panel","c":-1,"sub_med":f1s0.median(),"sst_med":f1t0.median(),"sst25":f1t0["Sst_25"]}]
per_ct={"clean":(f1s0,f1t0)}
for c in CS:
    mixed=(1-c)*p_true + c*ambient[None,:]
    obs=rng.poisson(Dtarget[:,None]*mixed).astype(np.float32)
    f1s,f1t=cv_f1(lognorm(obs))
    rows.append({"scenario":f"Xenium-degraded c={c}","c":c,"sub_med":f1s.median(),"sst_med":f1t.median(),"sst25":f1t["Sst_25"]})
    per_ct[f"c{c}"]=(f1s,f1t)
    print(f"c={c}: subclass med F1={f1s.median():.3f}, Sst med={f1t.median():.3f}, Sst_25={f1t['Sst_25']:.3f}",flush=True)
sweep=pd.DataFrame(rows); sweep.round(4).to_csv(os.path.join(OUT,"report_degraded_sweep.csv"),index=False)

# transcriptome ceiling refs (from the reporting tables resolvability_report.py writes)
_RES=os.path.join(_SP,"output/celltyping_supplement/data")
try:
    sc_all=pd.read_csv(os.path.join(_RES,"resolvability_subclass.csv")).set_index("cell_type")["F1_all"].median()
    _ss=pd.read_csv(os.path.join(_RES,"resolvability_sst_supertype.csv")).set_index("cell_type")
    ss_all=_ss["F1_all"].median(); ss25_all=_ss.loc["Sst_25","F1_all"]
except Exception: sc_all=ss_all=ss25_all=np.nan

# per-celltype table at a central spillover (c=0.2): clean vs degraded
cc="c0.2"; sub_tab=pd.DataFrame({"F1_clean_panel":f1s0,"F1_xenium_deg":per_ct[cc][0]}); sub_tab["drop"]=sub_tab.F1_clean_panel-sub_tab.F1_xenium_deg
sub_tab["class"]=[SUBCLASS_TO_CLASS.get(s,"Non-neuronal") for s in sub_tab.index]
sst_tab=pd.DataFrame({"F1_clean_panel":f1t0,"F1_xenium_deg":per_ct[cc][1]}); sst_tab["drop"]=sst_tab.F1_clean_panel-sst_tab.F1_xenium_deg
sub_tab.sort_values("F1_xenium_deg").round(4).to_csv(os.path.join(OUT,"report_degraded_per_celltype.csv"))
print("\n=== sweep ===\n"+sweep.round(3).to_string(index=False))
print(f"\ntranscriptome-ceiling medians: subclass={sc_all:.3f}, Sst={ss_all:.3f}, Sst_25={ss25_all:.3f}")
print("\n=== subclass F1 at c=0.2 (clean panel -> Xenium-degraded), most-degraded first ===")
print(sub_tab.sort_values("drop",ascending=False).head(10)[["class","F1_clean_panel","F1_xenium_deg","drop"]].round(3).to_string())
print("\n=== Sst supertype F1 at c=0.2 ===")
print(sst_tab.sort_values("F1_clean_panel",ascending=False)[["F1_clean_panel","F1_xenium_deg","drop"]].round(3).to_string())

# ── figure ──
fig,ax=plt.subplots(1,2,figsize=(13,5.4))
xs=["clean\npanel"]+[f"c={c}" for c in CS]
sub=[f1s0.median()]+[sweep[sweep.c==c].sub_med.values[0] for c in CS]
sstm=[f1t0.median()]+[sweep[sweep.c==c].sst_med.values[0] for c in CS]
s25=[f1t0["Sst_25"]]+[sweep[sweep.c==c].sst25.values[0] for c in CS]
xi=np.arange(len(xs))
ax[0].plot(xi,sub,"-o",color="#1565C0",label="subclass (median of 24)")
ax[0].axhline(sc_all,ls="--",color="#1565C0",alpha=.5,lw=1); ax[0].text(len(xs)-1,sc_all+.005,"transcriptome",color="#1565C0",fontsize=7.5,ha="right")
ax[0].set_xticks(xi); ax[0].set_xticklabels(xs); ax[0].set_ylabel("classification F1 (donor-CV)"); ax[0].set_ylim(0,1.02)
ax[0].set_title("Subclass identity is robust to\nXenium-like spillover"); ax[0].legend(loc="lower left",frameon=False); ax[0].spines[["top","right"]].set_visible(False)
ax[0].axvspan(0.5,len(xs)-0.5,color="grey",alpha=0.05)
ax[1].plot(xi,sstm,"-o",color="#2E7D32",label="Sst supertype (median of 16)")
ax[1].plot(xi,s25,"-o",color="#b8860b",label="Sst_25")
ax[1].axhline(ss_all,ls="--",color="#2E7D32",alpha=.5,lw=1); ax[1].axhline(ss25_all,ls="--",color="#b8860b",alpha=.5,lw=1)
ax[1].text(len(xs)-1,ss_all+.005,"transcriptome (median)",color="#2E7D32",fontsize=7.5,ha="right")
ax[1].set_xticks(xi); ax[1].set_xticklabels(xs); ax[1].set_ylabel("classification F1 (donor-CV)"); ax[1].set_ylim(0,1.02)
ax[1].set_title("Sst supertype F1 erodes with spillover\n(Sst_25 stays the most robust)"); ax[1].legend(loc="lower left",frameon=False); ax[1].spines[["top","right"]].set_visible(False)
ax[1].axvspan(0.5,len(xs)-0.5,color="grey",alpha=0.05)
fig.suptitle("Panel resolvability under approximate Xenium-like degradation (Poisson at calibrated Xenium depth + ambient spillover c)",fontsize=11.5,y=1.02)
fig.tight_layout()
for ext in ("png","pdf"): fig.savefig(os.path.join(OUT,f"report_degraded.{ext}"),dpi=200,bbox_inches="tight")
print("\nsaved report_degraded.png/.pdf + report_degraded_{sweep,per_celltype}.csv")
