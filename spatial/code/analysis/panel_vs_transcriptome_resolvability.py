#!/usr/bin/env python3
"""
Panel-vs-transcriptome control: is a supertype's poor resolvability the PANEL's
fault (discriminating genes exist but aren't on the 300-gene Xenium panel) or the
BIOLOGY's (the supertype boundary is intrinsically fuzzy — continuous even with the
whole transcriptome)?

Holds the classifier (nearest-centroid Pearson correlation) and CV (leave-one-donor-
out, 5 donors) FIXED; only the feature set changes:
  - panel  : 300 Xenium panel genes, normalized over the panel (mirrors Xenium as deployed)
  - hvg    : 2000 highly-variable genes, full-library normalized (the transcriptome ceiling)
  - allgene: all ~36k genes (batched), full-library normalized (literal full transcriptome,
             included to rebut "you didn't use the WHOLE transcriptome"; for a correlation
             classifier the extra noise genes can only dilute, so hvg is the fair ceiling)

Diagnostic per supertype: gap = F1(ceiling) - F1(panel).
  well-resolved   : panel F1 >= 0.6
  panel-limited   : panel F1 < 0.6 AND gap >= 0.15   (rescuable by a better panel)
  biology-limited : panel F1 < 0.6 AND gap < 0.15    (intrinsically fuzzy; no panel fixes it)

Reliability filter n>=25 reference cells (rarer types can't be CV-estimated).
Caveat: snRNA is cleaner than Xenium so absolute F1 is an upper bound; and the ceiling
is "best a CORRELATION classifier does" — a type fuzzy under correlation is very likely
genuinely hard, but a stronger classifier could rescue a few.

Outputs (output/depth_validation/lieber_layers/):
  panel_vs_transcriptome_resolvability.csv  +  .png/.pdf
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
    "font.size":11,"axes.titlesize":12.5,"axes.labelsize":11,"xtick.labelsize":10,
    "ytick.labelsize":10,"legend.fontsize":9,"axes.linewidth":0.6,"pdf.fonttype":42,"ps.fonttype":42})

BASE=os.path.expanduser("~/Github/SCZ_Xenium")
SNRNA=os.path.expanduser("~/Github/shared_data/nicole_sea_ad_snrnaseq_reference.h5ad")
_HERE=os.path.dirname(os.path.abspath(__file__)); _SP=os.path.dirname(os.path.dirname(_HERE))
OUT=os.path.join(_SP,"output/depth_validation/lieber_layers")
sys.path.insert(0,os.path.join(_SP,"code","modules")); from constants import SUBCLASS_TO_CLASS
MIN_N=25; N_HVG=2000; GAP_HI=0.15; PANEL_OK=0.6
CAT_COL={"well-resolved":"#2E7D32","panel-limited":"#E65100","biology-limited":"#6A1B9A"}

def corr_dense(Xtr,ytr,Xte,cand=None):
    types=[t for t in (cand if cand is not None else sorted(set(ytr))) if (ytr==t).sum()>0]
    C=np.vstack([Xtr[ytr==t].mean(0) for t in types])
    Cc=C-C.mean(1,keepdims=True); Cc/=(np.linalg.norm(Cc,axis=1,keepdims=True)+1e-9)
    Xc=Xte-Xte.mean(1,keepdims=True); Xc/=(np.linalg.norm(Xc,axis=1,keepdims=True)+1e-9)
    return np.array(types)[(Xc@Cc.T).argmax(1)]

def corr_sparse(Xtr,ytr,Xte,batch=4000):
    types=sorted(set(ytr))
    C=np.vstack([np.asarray(Xtr[ytr==t].mean(0)).ravel() for t in types])
    Cc=C-C.mean(1,keepdims=True); Cc/=(np.linalg.norm(Cc,axis=1,keepdims=True)+1e-9)
    out=np.empty(Xte.shape[0],dtype=object)
    for i in range(0,Xte.shape[0],batch):
        xb=Xte[i:i+batch].toarray(); xb=xb-xb.mean(1,keepdims=True); xb/=(np.linalg.norm(xb,axis=1,keepdims=True)+1e-9)
        out[i:i+batch]=np.array(types)[(xb@Cc.T).argmax(1)]
    return out

# ── load + two normalizations (panel-over-panel; transcriptome-over-all) ──
xh=os.path.join(BASE,"output/h5ad")
panel=list(ad.read_h5ad(os.path.join(xh,[f for f in sorted(os.listdir(xh)) if f.endswith("_annotated.h5ad")][0]),
                        backed="r").var_names)
print("loading full snRNAseq reference to memory...",flush=True)
a=ad.read_h5ad(SNRNA); ST=a.obs["Supertype"].astype(str).values; SUB=a.obs["Subclass"].astype(str).values
donor=a.obs["donor_id"].astype(str).values; N=len(ST); nd=len(set(donor))
supertypes=sorted(set(ST)); subclasses=sorted(set(SUB))
CLS={s:SUBCLASS_TO_CLASS.get(s,"Non-neuronal") for s in subclasses}; ST2SUB={t:SUB[ST==t][0] for t in supertypes}
pg=[g for g in panel if g in set(a.var_names)]

ap=a[:,pg].copy(); sc.pp.normalize_total(ap,target_sum=1e4); sc.pp.log1p(ap)   # panel-normalized
Xp=ap.X.toarray() if sp.issparse(ap.X) else np.asarray(ap.X); del ap
sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)                         # transcriptome-normalized (in place)
sc.pp.highly_variable_genes(a,n_top_genes=N_HVG)
hv=a.var["highly_variable"].values; Xh=a[:,hv].X.toarray() if sp.issparse(a.X) else np.asarray(a[:,hv].X)
Xa=a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
print(f"N={N:,} cells, {nd} donors, {len(supertypes)} supertypes | panel={Xp.shape[1]}g, hvg={Xh.shape[1]}g, all={Xa.shape[1]}g",flush=True)

# ── CV ──
def cv(kind,X,within=False):
    pred=np.empty(N,dtype=object); pred[:]=""
    for tr,te in GroupKFold(nd).split(np.empty((N,1)),ST,donor):
        if within:
            for s in np.unique(SUB[te]):
                cand=sorted(set(ST[tr][SUB[tr]==s]))
                if not cand: continue
                idx=te[SUB[te]==s]; m=SUB[tr]==s
                pred[idx]=corr_dense(X[tr][m],ST[tr][m],X[idx],cand)
        elif kind=="sparse":
            pred[te]=corr_sparse(X[tr],ST[tr],X[te])
        else:
            pred[te]=corr_dense(X[tr],ST[tr],X[te])
    return pred

def f1(pred):
    r=classification_report(ST,pred,labels=supertypes,output_dict=True,zero_division=0)
    return pd.Series({s:r[s]["f1-score"] for s in supertypes})

print("CV panel (full, within)...",flush=True); pf=f1(cv("dense",Xp)); pw=f1(cv("dense",Xp,within=True))
print("CV hvg (full, within)...",flush=True);   hf=f1(cv("dense",Xh)); hw=f1(cv("dense",Xh,within=True))
print("CV all-genes (full, batched)...",flush=True); af=f1(cv("sparse",Xa))

df=pd.DataFrame({"panel_full":pf,"panel_within":pw,"hvg_full":hf,"hvg_within":hw,"allgene_full":af})
df["n"]=[int((ST==s).sum()) for s in supertypes]; df["subclass"]=[ST2SUB[s] for s in supertypes]
df["class"]=[CLS[ST2SUB[s]] for s in supertypes]; df["reliable"]=df.n>=MIN_N
df["gap"]=df.hvg_full-df.panel_full
def cat(r):
    if not r.reliable: return "untestable"
    if r.panel_full>=PANEL_OK: return "well-resolved"
    return "panel-limited" if r.gap>=GAP_HI else "biology-limited"
df["category"]=df.apply(cat,axis=1)
df.sort_values(["class","subclass","panel_full"],ascending=[True,True,False]).round(4).to_csv(
    os.path.join(OUT,"panel_vs_transcriptome_resolvability.csv"))

dr=df[df.reliable]
print(f"\nReliable supertypes n>={MIN_N}: {len(dr)}/{len(df)}")
print("\n=== median F1: panel vs ceiling (reliable) ===")
print(dr.groupby("class")[["panel_full","hvg_full","allgene_full"]].median().round(3).to_string())
print(f"\nALL-GENES vs HVG ceiling (median full F1): all={dr.allgene_full.median():.3f}, hvg={dr.hvg_full.median():.3f} "
      f"(close → ceiling robust to gene count)")
print("\n=== category counts (reliable) ===")
print(dr.category.value_counts().to_string())
print("\n=== panel-LIMITED: poor on panel, rescued by transcriptome (biggest gaps) ===")
print(dr[dr.category=="panel-limited"].sort_values("gap",ascending=False).head(10)[
    ["subclass","class","n","panel_full","hvg_full","gap"]].to_string())
print("\n=== biology-LIMITED: fuzzy even with the full transcriptome ===")
print(dr[dr.category=="biology-limited"].sort_values("panel_full").head(10)[
    ["subclass","class","n","panel_full","hvg_full","gap"]].to_string())
sst=dr[dr.subclass=="Sst"].sort_values("panel_full",ascending=False)
print("\n=== SST family: panel vs ceiling ===")
print(sst[["n","panel_full","hvg_full","gap","category"]].round(3).to_string())

# ════ figure ════
fig,ax=plt.subplots(1,2,figsize=(14.5,6.2),gridspec_kw={"width_ratios":[1.25,1.0]})
# (A) panel vs ceiling scatter
ax[0].plot([0,1],[0,1],"--",color="grey",lw=.8,zorder=1)
ax[0].axhspan(0,1,xmax=0,alpha=0)  # noop keeps axis
for c,col in CAT_COL.items():
    d=dr[dr.category==c]; ax[0].scatter(d.panel_full,d.hvg_full,s=22,c=col,edgecolor="white",linewidth=.3,
                                         alpha=.9,label=c,zorder=3)
# label the biggest-gap (panel-limited) types + all Sst
lab=set(dr[dr.category=="panel-limited"].sort_values("gap",ascending=False).head(8).index)|set(sst.index[:3])|{"Sst_25"}
for s in lab:
    if s in dr.index: ax[0].annotate(s,(dr.loc[s,"panel_full"],dr.loc[s,"hvg_full"]),xytext=(3,3),
                                      textcoords="offset points",fontsize=7.5,color="#333")
ax[0].scatter([dr.loc["Sst_25","panel_full"]],[dr.loc["Sst_25","hvg_full"]],s=60,facecolor="none",
              edgecolor="#b8860b",linewidth=1.6,zorder=5)
ax[0].set_xlabel("F1 on Xenium panel (300 genes, as deployed)"); ax[0].set_ylabel("F1 ceiling (transcriptome, 2000 HVG)")
ax[0].set_xlim(0,1.0); ax[0].set_ylim(0,1.0); ax[0].set_aspect("equal")
ax[0].text(.30,.96,"above line = panel-limited\n(genes exist, off-panel)",fontsize=8.5,va="top",color="#E65100",style="italic")
ax[0].text(.62,.10,"low + on line =\nbiology-limited",fontsize=8.5,va="bottom",color="#6A1B9A",style="italic")
ax[0].set_title("Is poor resolution the panel's fault or the biology's?"); ax[0].legend(loc="lower right",frameon=False)
ax[0].spines[["top","right"]].set_visible(False)
# (B) SST family dumbbell: panel -> ceiling
so=sst.sort_values("panel_full").index.tolist(); yi=np.arange(len(so))
for i,s in enumerate(so):
    ax[1].plot([sst.loc[s,"panel_full"],sst.loc[s,"hvg_full"]],[i,i],"-",color="#cccccc",lw=2,zorder=1)
ax[1].scatter(sst.loc[so,"panel_full"],yi,s=34,c="#4393c3",label="panel (300g)",zorder=3,edgecolor="white",linewidth=.4)
ax[1].scatter(sst.loc[so,"hvg_full"],yi,s=34,c="#222222",label="ceiling (HVG)",zorder=3,edgecolor="white",linewidth=.4)
ax[1].set_yticks(yi); ax[1].set_yticklabels(so,fontsize=9)
ax[1].get_yticklabels()[so.index("Sst_25")].set_fontweight("bold"); ax[1].get_yticklabels()[so.index("Sst_25")].set_color("#b8860b")
ax[1].set_xlabel("full-multiclass F1"); ax[1].set_xlim(0,1.0)
ax[1].set_title("SST: ceilings are high, so hard Sst types are\npanel-limited (rescuable); Sst_25 nearest its ceiling")
ax[1].legend(loc="lower right",frameon=False); ax[1].spines[["top","right"]].set_visible(False)
fig.tight_layout()
for ext in ("png","pdf"): fig.savefig(os.path.join(OUT,f"panel_vs_transcriptome_resolvability.{ext}"),dpi=200,bbox_inches="tight")
print("\nsaved panel_vs_transcriptome_resolvability.png/.pdf + .csv")
