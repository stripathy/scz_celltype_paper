#!/usr/bin/env python3
"""
Systematic per-supertype resolvability on the Xenium panel.

Generalizes sst25_resolvability.py to the whole SEA-AD taxonomy: how cleanly can
each supertype (and subclass) be recovered from the 300 Xenium panel genes alone?

Substrate: SEA-AD neurotypical-MTG snRNAseq (whole-transcriptome, ground-truth
Class/Subclass/Supertype), restricted to the Xenium panel genes and normalized over
that panel (mirrors what Xenium actually measures). Classifier: nearest-centroid by
Pearson correlation (faithful, simplified replica of the pipeline's correlation
classifier). CV: leave-one-donor-out (5 donors) — no donor leakage.

Three scorings, per cell, cross-validated:
  - full-multiclass : classify among ALL supertypes (absolute resolvability; matches
                      the Sst_25 analysis).
  - within-subclass : classify among the supertypes of the cell's TRUE subclass
                      (supertype resolvability given the subclass — faithful to the
                      deployed two-stage classifier's stage 2). NB singleton subclasses
                      score trivially 1.0 here (flagged via n_supertypes).
  - subclass-level  : classify among the 24 subclass centroids (context).

Caveat: snRNAseq is cleaner than Xenium (more counts/cell, no spatial spillover), so
absolute F1 is an UPPER BOUND; the RANKING across types is the robust read-out.

Outputs (output/depth_validation/lieber_layers/):
  supertype_resolvability.csv   (137 rows: recall/prec/F1 full + within, n, subclass, class)
  subclass_resolvability.csv    (24 rows)
  supertype_resolvability.png/.pdf  (landscape + SST family + full-vs-within)
"""
import os, sys
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, scipy.sparse as sp
from sklearn.model_selection import GroupKFold
from sklearn.metrics import classification_report
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

matplotlib.rcParams.update({
    "font.family":"sans-serif","font.sans-serif":["Helvetica","Arial","DejaVu Sans"],
    "font.size":11,"axes.titlesize":12.5,"axes.labelsize":11,"xtick.labelsize":9,
    "ytick.labelsize":9,"legend.fontsize":9,"axes.linewidth":0.6,"pdf.fonttype":42,"ps.fonttype":42})

BASE=os.path.expanduser("~/Github/SCZ_Xenium")
SNRNA=os.path.expanduser("~/Github/shared_data/nicole_sea_ad_snrnaseq_reference.h5ad")
_HERE=os.path.dirname(os.path.abspath(__file__)); _SP=os.path.dirname(os.path.dirname(_HERE))
OUT=os.path.join(_SP,"output/depth_validation/lieber_layers")
sys.path.insert(0,os.path.join(_SP,"code","modules"))
from constants import SUBCLASS_TO_CLASS
CLASS_COLORS={"Glutamatergic":"#E65100","GABAergic":"#2E7D32","Non-neuronal":"#1565C0"}
CLASS_RANK={"GABAergic":0,"Glutamatergic":1,"Non-neuronal":2}

def corr_classify(Xtr,ytr,Xte,cand=None):
    types=[t for t in (cand if cand is not None else sorted(set(ytr))) if (ytr==t).sum()>0]
    C=np.vstack([Xtr[ytr==t].mean(0) for t in types])
    Cc=C-C.mean(1,keepdims=True); Cc/=(np.linalg.norm(Cc,axis=1,keepdims=True)+1e-9)
    Xc=Xte-Xte.mean(1,keepdims=True); Xc/=(np.linalg.norm(Xc,axis=1,keepdims=True)+1e-9)
    return np.array(types)[(Xc@Cc.T).argmax(1)]

# ── load snRNAseq on the Xenium panel ──
xh=os.path.join(BASE,"output/h5ad")
panel=list(ad.read_h5ad(os.path.join(xh,[f for f in sorted(os.listdir(xh)) if f.endswith("_annotated.h5ad")][0]),
                        backed="r").var_names)
ref=ad.read_h5ad(SNRNA,backed="r"); genes=[g for g in panel if g in set(ref.var_names)]
sub=ref[:,genes].to_memory(); sc.pp.normalize_total(sub,target_sum=1e4); sc.pp.log1p(sub)
X=sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
ST=sub.obs["Supertype"].astype(str).values; SUB=sub.obs["Subclass"].astype(str).values
donor=sub.obs["donor_id"].astype(str).values; nd=len(set(donor))
MIN_N=25   # min reference cells to estimate F1 (rarer types can't be CV-assessed, not "unresolvable")
supertypes=sorted(set(ST)); subclasses=sorted(set(SUB))
CLS={s:SUBCLASS_TO_CLASS.get(s,"Non-neuronal") for s in subclasses}
ST2SUB={t:SUB[ST==t][0] for t in supertypes}
print(f"snRNAseq: {len(genes)}/{len(panel)} panel genes, {nd} donors, "
      f"{len(ST):,} cells, {len(supertypes)} supertypes, {len(subclasses)} subclasses",flush=True)

# ── cross-validated predictions ──
def cv(mode):
    pred=np.empty(len(ST),dtype=object); pred[:]=""
    for tr,te in GroupKFold(nd).split(X,ST,donor):
        if mode=="full":
            pred[te]=corr_classify(X[tr],ST[tr],X[te])
        elif mode=="within":
            for s in np.unique(SUB[te]):
                cand=sorted(set(ST[tr][SUB[tr]==s]))
                if not cand: continue
                idx=te[SUB[te]==s]; mtr=SUB[tr]==s
                pred[idx]=corr_classify(X[tr][mtr],ST[tr][mtr],X[idx],cand=cand)
        elif mode=="subclass":
            pred[te]=corr_classify(X[tr],SUB[tr],X[te])
    return pred
print("CV: full-multiclass...",flush=True); pf=cv("full")
print("CV: within-subclass...",flush=True); pw=cv("within")
print("CV: subclass-level...",flush=True); psub=cv("subclass")

def rows(gt,pred,labels,tag):
    rep=classification_report(gt,pred,labels=labels,output_dict=True,zero_division=0)
    return pd.DataFrame([{"name":s,f"recall_{tag}":rep[s]["recall"],f"prec_{tag}":rep[s]["precision"],
                          f"F1_{tag}":rep[s]["f1-score"]} for s in labels]).set_index("name")

df=rows(ST,pf,supertypes,"full").join(rows(ST,pw,supertypes,"within"))
df["n"]=[int((ST==s).sum()) for s in supertypes]
df["subclass"]=[ST2SUB[s] for s in supertypes]; df["class"]=[CLS[ST2SUB[s]] for s in supertypes]
df["reliable"]=df.n>=MIN_N   # F1 only estimable with enough reference cells across donors
nsup=df.groupby("subclass").size(); df["n_supertypes_in_subclass"]=df["subclass"].map(nsup)
dfsub=rows(SUB,psub,subclasses,"subclass"); dfsub["n"]=[int((SUB==s).sum()) for s in subclasses]
dfsub["class"]=[CLS[s] for s in subclasses]; dfsub["n_supertypes"]=dfsub.index.map(nsup)
df=df.sort_values(["class","subclass","F1_full"],ascending=[True,True,False])
df.round(4).to_csv(os.path.join(OUT,"supertype_resolvability.csv"))
dfsub.sort_values(["class","F1_subclass"],ascending=[True,False]).round(4).to_csv(os.path.join(OUT,"subclass_resolvability.csv"))

# ── summary (reliable subset only: F1 is meaningless for n<MIN_N) ──
dr=df[df.reliable]
print(f"\nReliable supertypes (n>={MIN_N}): {len(dr)}/{len(df)}  "
      f"({(~df.reliable).sum()} excluded for too few reference cells, mostly rare -SEAAD states)")
print("\n=== median F1 by class (reliable only) ===")
print(dr.groupby("class").agg(n_supertypes=("F1_full","size"),
      F1_full_median=("F1_full","median"),F1_within_median=("F1_within","median")).round(3).to_string())
print("\n=== 8 best-resolved supertypes (full-multiclass F1, reliable) ===")
print(dr.sort_values("F1_full",ascending=False).head(8)[["subclass","class","n","F1_full","F1_within"]].to_string())
print("\n=== 8 worst-resolved supertypes (full-multiclass F1, reliable) ===")
print(dr.sort_values("F1_full").head(8)[["subclass","class","n","F1_full","F1_within"]].to_string())
sst=df[df.subclass=="Sst"].sort_values("F1_full",ascending=False)
print(f"\n=== SST family (16), by full-multiclass F1 ===")
print(sst[["n","recall_full","prec_full","F1_full","F1_within"]].round(3).to_string())
from scipy.stats import spearmanr
print(f"\nF1_full vs log10(n) [reliable only]: Spearman {spearmanr(np.log10(dr.n),dr.F1_full)[0]:+.2f}")

# ════ figure (landscape + scatter use the reliable subset; SST all reliable) ════
fig=plt.figure(figsize=(16,5.4)); gs=fig.add_gridspec(1,3,width_ratios=[1.55,1.0,1.0],wspace=0.32)
# (a) landscape: full F1 by subclass, points=supertypes, colored by class
axA=fig.add_subplot(gs[0,0])
medF=dr.groupby("subclass")["F1_full"].median()
sub_order=sorted([s for s in subclasses if (dr.subclass==s).any()],key=lambda s:(CLASS_RANK[CLS[s]],-medF[s]))
rng=np.random.default_rng(0)
for i,s in enumerate(sub_order):
    d=dr[dr.subclass==s]; xj=i+rng.uniform(-0.22,0.22,len(d))
    axA.scatter(xj,d.F1_full,s=14,c=CLASS_COLORS[CLS[s]],edgecolor="white",linewidth=.3,alpha=.9,zorder=3)
    axA.plot([i-0.3,i+0.3],[medF[s]]*2,color="black",lw=1.4,zorder=4)
axA.set_xticks(range(len(sub_order))); axA.set_xticklabels(sub_order,rotation=90)
axA.set_ylabel("full-multiclass F1  (Xenium panel, donor-CV)"); axA.set_ylim(0,1.0)
axA.set_title("Supertype resolvability across the taxonomy\n(point = supertype, bar = subclass median)")
axA.text(0.99,0.02,f"n≥{MIN_N} cells; {(~df.reliable).sum()} rare types excluded",transform=axA.transAxes,
         ha="right",va="bottom",fontsize=8,fontstyle="italic",color="grey")
axA.spines[["top","right"]].set_visible(False)
axA.legend(handles=[Line2D([0],[0],marker="o",ls="",mfc=c,mec="white",label=k) for k,c in CLASS_COLORS.items()],
           loc="lower right",frameon=False)
# (b) SST family bars (full F1), Sst_25 highlighted
axB=fig.add_subplot(gs[0,1]); so=sst.sort_values("F1_full").index.tolist(); yi=np.arange(len(so))
cols=["#b8860b" if s=="Sst_25" else CLASS_COLORS["GABAergic"] for s in so]
axB.barh(yi,sst.loc[so,"F1_full"],color=cols,edgecolor="black",linewidth=.4)
axB.set_yticks(yi); axB.set_yticklabels(so,fontsize=8.5)
axB.get_yticklabels()[so.index("Sst_25")].set_fontweight("bold"); axB.get_yticklabels()[so.index("Sst_25")].set_color("#b8860b")
for i,s in enumerate(so): axB.text(sst.loc[s,"F1_full"]+.01,i,f"{sst.loc[s,'F1_full']:.2f}",va="center",fontsize=7.5)
axB.set_xlabel("full-multiclass F1"); axB.set_xlim(0,.95); axB.set_title("SST family\n(Sst_25 = best)")
axB.spines[["top","right"]].set_visible(False)
# (c) full vs within-subclass F1, all supertypes
axC=fig.add_subplot(gs[0,2])
axC.plot([0,1],[0,1],"--",color="grey",lw=.8,zorder=1)
for k,c in CLASS_COLORS.items():
    d=dr[dr["class"]==k]; axC.scatter(d.F1_full,d.F1_within,s=14,c=c,edgecolor="white",linewidth=.3,alpha=.85,label=k,zorder=3)
axC.scatter([sst.loc["Sst_25","F1_full"]],[sst.loc["Sst_25","F1_within"]],s=55,c="#b8860b",edgecolor="black",linewidth=.8,zorder=5)
axC.annotate("Sst_25",(sst.loc["Sst_25","F1_full"],sst.loc["Sst_25","F1_within"]),xytext=(5,-9),
             textcoords="offset points",fontweight="bold",color="#b8860b",fontsize=9)
axC.set_xlabel("full-multiclass F1"); axC.set_ylabel("within-subclass F1")
axC.set_title("Subclass leakage vs within-subclass\nambiguity (below line = subclass leak)")
axC.set_xlim(0,1.0); axC.set_ylim(0,1.02); axC.spines[["top","right"]].set_visible(False)
fig.tight_layout()
for ext in ("png","pdf"): fig.savefig(os.path.join(OUT,f"supertype_resolvability.{ext}"),dpi=200,bbox_inches="tight")
print("\nsaved supertype_resolvability.png/.pdf + supertype_resolvability.csv + subclass_resolvability.csv")
