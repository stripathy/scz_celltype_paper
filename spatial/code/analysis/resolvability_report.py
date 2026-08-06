#!/usr/bin/env python3
"""
Publication-ready reporting tables: per-cell-type classification F1 on the Xenium
PANEL vs the full TRANSCRIPTOME, at two resolutions:
  - subclass (24-way classification)
  - Sst supertype (16 Sst types; full-multiclass over all 137 supertypes, reported for Sst)

Same machinery throughout: nearest-centroid Pearson-correlation classifier,
leave-one-donor-out CV on the SEA-AD neurotypical snRNA reference (5 donors).
Feature sets (each cell normalized in its own space, mirroring how it is measured):
  - panel : 300 Xenium-panel genes, normalized over the panel (as deployed on Xenium)
  - all   : all ~36k genes, full-library normalized (literal full transcriptome)
  - hvg   : 2000 HVG, full-library normalized (de-noised transcriptome ceiling; ~= all)

Caveat baked into the report text: snRNA is cleaner than Xenium, so absolute F1 is an
upper bound; this measures the INFORMATION the panel carries (separability ceiling),
not deployed Xenium accuracy. The panel-vs-transcriptome GAP is the citable quantity.

Outputs:
  output/celltyping_supplement/data/resolvability_subclass.csv
  output/celltyping_supplement/data/resolvability_sst_supertype.csv
      -- the two tables the R figure scripts read (S2 panels c/e). These used to be
      written into lieber_layers/ under report_*.csv names and copied here by hand,
      which meant re-running this script updated nothing the figure actually read.
  output/depth_validation/lieber_layers/report_resolvability.png/.pdf
      -- matplotlib preview, not used by the submission figure.
"""
import os, sys
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, scipy.sparse as sp
from sklearn.model_selection import GroupKFold
from sklearn.metrics import classification_report
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

matplotlib.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Helvetica","Arial","DejaVu Sans"],
    "font.size":11,"axes.titlesize":12.5,"axes.labelsize":11,"xtick.labelsize":9.5,"ytick.labelsize":9.5,
    "legend.fontsize":9,"axes.linewidth":0.6,"pdf.fonttype":42,"ps.fonttype":42})

BASE=os.path.expanduser("~/Github/SCZ_Xenium")
NICOLE=os.path.expanduser("~/Github/shared_data/nicole_sea_ad_snrnaseq_reference.h5ad")
_SP=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT=os.path.join(_SP,"output/depth_validation/lieber_layers")
# The figure-input tables live with the other S2/S3 inputs, under the names the R reads.
RES_OUT=os.path.join(_SP,"output/celltyping_supplement/data"); os.makedirs(RES_OUT,exist_ok=True)
sys.path.insert(0,os.path.join(_SP,"code","modules")); from constants import SUBCLASS_TO_CLASS
CLASS_COLORS={"Glutamatergic":"#E65100","GABAergic":"#2E7D32","Non-neuronal":"#1565C0"}
N_HVG=2000

def corr_dense(Xtr,ytr,Xte,cand=None):
    types=[t for t in (cand if cand is not None else sorted(set(ytr))) if (ytr==t).sum()>0]
    C=np.vstack([Xtr[ytr==t].mean(0) for t in types])
    Cc=C-C.mean(1,keepdims=True); Cc/=(np.linalg.norm(Cc,axis=1,keepdims=True)+1e-9)
    Xc=Xte-Xte.mean(1,keepdims=True); Xc/=(np.linalg.norm(Xc,axis=1,keepdims=True)+1e-9)
    return np.array(types)[(Xc@Cc.T).argmax(1)]
def corr_sparse(Xtr,ytr,Xte,batch=4000):
    types=sorted(set(ytr)); C=np.vstack([np.asarray(Xtr[ytr==t].mean(0)).ravel() for t in types])
    Cc=C-C.mean(1,keepdims=True); Cc/=(np.linalg.norm(Cc,axis=1,keepdims=True)+1e-9)
    out=np.empty(Xte.shape[0],dtype=object)
    for i in range(0,Xte.shape[0],batch):
        xb=Xte[i:i+batch].toarray(); xb=xb-xb.mean(1,keepdims=True); xb/=(np.linalg.norm(xb,axis=1,keepdims=True)+1e-9)
        out[i:i+batch]=np.array(types)[(xb@Cc.T).argmax(1)]
    return out
def dense(x): return x.toarray() if sp.issparse(x) else np.asarray(x)

# ── load reference, three feature spaces ──
xh=os.path.join(BASE,"output/h5ad")
panel=list(ad.read_h5ad(os.path.join(xh,[f for f in sorted(os.listdir(xh)) if f.endswith("_annotated.h5ad")][0]),backed="r").var_names)
print("loading reference...",flush=True); a=ad.read_h5ad(NICOLE)
SUB=a.obs["Subclass"].astype(str).values; SUP=a.obs["Supertype"].astype(str).values
donor=a.obs["donor_id"].astype(str).values; N=len(SUB); nd=len(set(donor))
pg=[g for g in panel if g in set(a.var_names)]
ap=a[:,pg].copy(); sc.pp.normalize_total(ap,target_sum=1e4); sc.pp.log1p(ap); Xp=dense(ap.X); del ap
sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a); sc.pp.highly_variable_genes(a,n_top_genes=N_HVG)
hv=list(a.var_names[a.var["highly_variable"].values]); Xh=dense(a[:,hv].X); Xa=a.X.tocsr(); del a
subclasses=sorted(set(SUB)); supertypes=sorted(set(SUP)); sst=sorted([s for s in supertypes if s.startswith("Sst_")],key=lambda s:int(s.split("_")[1]))
print(f"N={N:,}, {nd} donors | panel={len(pg)}g hvg={len(hv)}g all={Xa.shape[1]}g | {len(subclasses)} subclasses, {len(sst)} Sst supertypes",flush=True)

# ── CV: predict subclass(24-way) and supertype(full) for each feature set ──
P={k:np.empty(N,object) for k in ["sub_panel","sup_panel","sub_hvg","sup_hvg","sub_all","sup_all"]}
for fi,(tr,te) in enumerate(GroupKFold(nd).split(np.empty((N,1)),SUB,donor)):
    print(f"  fold {fi+1}/{nd}...",flush=True)
    P["sub_panel"][te]=corr_dense(Xp[tr],SUB[tr],Xp[te]); P["sup_panel"][te]=corr_dense(Xp[tr],SUP[tr],Xp[te])
    P["sub_hvg"][te]=corr_dense(Xh[tr],SUB[tr],Xh[te]);   P["sup_hvg"][te]=corr_dense(Xh[tr],SUP[tr],Xh[te])
    P["sub_all"][te]=corr_sparse(Xa[tr],SUB[tr],Xa[te]);  P["sup_all"][te]=corr_sparse(Xa[tr],SUP[tr],Xa[te])

def prf(gt,pred,labels,tag):
    r=classification_report(gt,pred,labels=labels,output_dict=True,zero_division=0)
    return pd.DataFrame([{"cell_type":s,f"prec_{tag}":r[s]["precision"],f"recall_{tag}":r[s]["recall"],
                          f"F1_{tag}":r[s]["f1-score"]} for s in labels]).set_index("cell_type")

# subclass table
sc_df=prf(SUB,P["sub_panel"],subclasses,"panel").join(prf(SUB,P["sub_hvg"],subclasses,"hvg")).join(prf(SUB,P["sub_all"],subclasses,"all"))
sc_df["n"]=[int((SUB==s).sum()) for s in subclasses]; sc_df["class"]=[SUBCLASS_TO_CLASS.get(s,"Non-neuronal") for s in subclasses]
sc_df["gap_all"]=sc_df.F1_all-sc_df.F1_panel
sc_df=sc_df.sort_values(["class","F1_panel"],ascending=[True,False])
sc_df.round(4).to_csv(os.path.join(RES_OUT,"resolvability_subclass.csv"))

# Sst supertype table
ss_df=prf(SUP,P["sup_panel"],sst,"panel").join(prf(SUP,P["sup_hvg"],sst,"hvg")).join(prf(SUP,P["sup_all"],sst,"all"))
ss_df["n"]=[int((SUP==s).sum()) for s in sst]; ss_df["gap_all"]=ss_df.F1_all-ss_df.F1_panel
ss_df=ss_df.sort_values("F1_panel",ascending=False)
ss_df.round(4).to_csv(os.path.join(RES_OUT,"resolvability_sst_supertype.csv"))

pd.set_option("display.width",200)
print("\n================ SUBCLASS (24-way) — F1 panel vs transcriptome ================")
print(sc_df[["n","class","F1_panel","F1_hvg","F1_all","gap_all"]].round(3).to_string())
print(f"\nsubclass median F1: panel={sc_df.F1_panel.median():.3f}  all-genes={sc_df.F1_all.median():.3f}  (mean gap {sc_df.gap_all.mean():.3f})")
print("\n================ Sst SUPERTYPE (16) — F1 panel vs transcriptome ================")
print(ss_df[["n","prec_panel","recall_panel","F1_panel","F1_hvg","F1_all","gap_all"]].round(3).to_string())
print(f"\nSst median F1: panel={ss_df.F1_panel.median():.3f}  all-genes={ss_df.F1_all.median():.3f}  (mean gap {ss_df.gap_all.mean():.3f})")

# ── figure: dumbbell panel -> all-genes, subclass (left) + Sst supertype (right) ──
fig,ax=plt.subplots(1,2,figsize=(13.5,6.4),gridspec_kw={"width_ratios":[1.05,1.0]})
def dumbbell(axx,d,order,ttl,color_by_class=False,hi=None):
    yi=np.arange(len(order))
    for i,s in enumerate(order):
        axx.plot([d.loc[s,"F1_panel"],d.loc[s,"F1_all"]],[i,i],"-",color="#cfcfcf",lw=2,zorder=1)
    pc=[CLASS_COLORS[d.loc[s,"class"]] for s in order] if color_by_class else "#4393c3"
    axx.scatter(d.loc[order,"F1_panel"],yi,s=34,c=pc,zorder=3,edgecolor="white",linewidth=.4,label="panel (300 g)")
    axx.scatter(d.loc[order,"F1_all"],yi,s=34,c="#222222",zorder=3,edgecolor="white",linewidth=.4,label="all genes (~36k)")
    axx.set_yticks(yi); axx.set_yticklabels(order)
    if hi and hi in order:
        axx.get_yticklabels()[order.index(hi)].set_fontweight("bold"); axx.get_yticklabels()[order.index(hi)].set_color("#b8860b")
    axx.set_xlabel("classification F1 (leave-one-donor-out CV)"); axx.set_xlim(0,1.0); axx.set_title(ttl)
    axx.spines[["top","right"]].set_visible(False)
order_sc=sc_df.sort_values("F1_panel").index.tolist()
dumbbell(ax[0],sc_df,order_sc,"Subclass identity (24-way)",color_by_class=True)
ax[0].legend(loc="lower right",frameon=False,handletextpad=.2)
ax[0].legend(handles=[Line2D([0],[0],marker="o",ls="",mfc=CLASS_COLORS[k],mec="white",label=k) for k in CLASS_COLORS]
             +[Line2D([0],[0],marker="o",ls="",mfc="#222",mec="white",label="all genes")],loc="lower right",frameon=False,fontsize=8)
order_ss=ss_df.sort_values("F1_panel").index.tolist()
dumbbell(ax[1],ss_df,order_ss,"Sst supertype resolution (16)",hi="Sst_25"); ax[1].legend(loc="lower right",frameon=False,handletextpad=.2)
fig.suptitle("Cell-type resolvability on the Xenium panel vs the full transcriptome (SEA-AD snRNA, donor-CV)",fontsize=12.5,y=1.01)
fig.tight_layout()
for ext in ("png","pdf"): fig.savefig(os.path.join(OUT,f"report_resolvability.{ext}"),dpi=200,bbox_inches="tight")
print("\nsaved report_resolvability.png/.pdf + report_{subclass,sst_supertype}_resolvability.csv")
