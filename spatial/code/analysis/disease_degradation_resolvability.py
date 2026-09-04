#!/usr/bin/env python3
"""
Disease stress-test: does AD pathology degrade panel-restricted cell-type
classification MORE than full-transcriptome classification? I.e. when you drop to
the Xenium panel, are calls fragile to disease-driven expression shifts?

(Not SCZ subjects — AD — but a clean probe for how disease status interacts with
marker loss, since the SEA-AD snRNAseq spans the full neuropathology spectrum.)

Design (faithful to deployment: neurotypical reference -> diseased query):
  - Centroids built from the 5-donor NEUROTYPICAL reference (nicole h5ad), in two
    feature spaces: panel (300 genes, panel-normalized) and ceiling (2000 HVG,
    full-library normalized). Subclass + supertype centroids each.
  - Query = disease cells (ADNC Not AD->High, CPS 0.15-0.93, 84 donors) from the full
    SEA-AD MTG snRNAseq (a representative contiguous 250k-row cross-section, capped per
    donor). NO donor overlap with the reference.
  - Classify query with reference centroids; score vs SEA-AD ground-truth labels.
  - Per donor: subclass & supertype accuracy under panel vs ceiling; regress on CPS.
    KEY: if the panel accuracy falls with CPS faster than the ceiling (the gap widens),
    dropping markers makes classification vulnerable to disease.

Outputs (output/depth_validation/lieber_layers/):
  disease_degradation_per_donor.csv  +  disease_degradation_resolvability.png/.pdf
"""
import os, sys, time
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, scipy.sparse as sp, h5py
from scipy.stats import linregress
from sklearn.metrics import classification_report
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

matplotlib.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Helvetica","Arial","DejaVu Sans"],
    "font.size":11,"axes.titlesize":12.5,"axes.labelsize":11,"xtick.labelsize":10,"ytick.labelsize":10,
    "legend.fontsize":9,"axes.linewidth":0.6,"pdf.fonttype":42,"ps.fonttype":42})

BASE=os.path.expanduser("~/Github/SCZ_Xenium")
NICOLE=os.path.expanduser("~/Github/shared_data/nicole_sea_ad_snrnaseq_reference.h5ad")
FULL=os.path.join(BASE,"data/reference/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad")
_SP=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT=os.path.join(_SP,"output/depth_validation/lieber_layers")
N_HVG=2000; SPAN=250000; CHUNK=50000; CAP=1000   # read SPAN rows in CHUNK-sized pieces (memory-bounded), cap cells/donor
t0=time.time(); tic=lambda m:print(f"[{time.time()-t0:5.0f}s] {m}",flush=True)

def centroids(X,labels):
    types=sorted(set(labels)); C=np.vstack([X[labels==t].mean(0) for t in types])
    Cc=C-C.mean(1,keepdims=True); Cc/=(np.linalg.norm(Cc,axis=1,keepdims=True)+1e-9)
    return np.array(types),Cc
def classify(Xq,types,Cc):
    Xc=Xq-Xq.mean(1,keepdims=True); Xc/=(np.linalg.norm(Xc,axis=1,keepdims=True)+1e-9)
    return types[(Xc@Cc.T).argmax(1)]
def dense(x): return x.toarray() if sp.issparse(x) else np.asarray(x)

# ════ 1. reference centroids (neurotypical) ════
xh=os.path.join(BASE,"output/h5ad")
panel=list(ad.read_h5ad(os.path.join(xh,[f for f in sorted(os.listdir(xh)) if f.endswith("_annotated.h5ad")][0]),
                        backed="r").var_names)
tic("loading neurotypical reference..."); a=ad.read_h5ad(NICOLE)
pg=[g for g in panel if g in set(a.var_names)]
sub=a.obs["Subclass"].astype(str).values; sup=a.obs["Supertype"].astype(str).values
ap=a[:,pg].copy(); sc.pp.normalize_total(ap,target_sum=1e4); sc.pp.log1p(ap); Xpr=dense(ap.X); del ap
sub_tp,Cp_sub=centroids(Xpr,sub); sup_tp,Cp_sup=centroids(Xpr,sup); del Xpr
sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a); sc.pp.highly_variable_genes(a,n_top_genes=N_HVG)
hv=list(a.var_names[a.var["highly_variable"].values]); Xhr=dense(a[:,hv].X)
sub_th,Ch_sub=centroids(Xhr,sub); sup_th,Ch_sup=centroids(Xhr,sup); del Xhr,a
tic(f"reference centroids built: {len(sup_tp)} supertypes, {len(sub_tp)} subclasses, panel={len(pg)}g hvg={len(hv)}g")

# ════ 2. query via raw h5py (no full-obs load): contiguous X slabs + per-chunk obs slices ════
def obs_col(f,name,cs,ce):                       # old-format categoricals live in obs['__categories']
    if name in f["obs"]["__categories"]:
        codes=f["obs"][name][cs:ce]; cats=np.array([c.decode() if isinstance(c,bytes) else str(c) for c in f["obs"]["__categories"][name][:]])
        return np.where(codes>=0,cats[np.clip(codes,0,None)],"nan")
    return f["obs"][name][cs:ce]                  # numeric (e.g. CPS float)
def read_csr_rows(g,cs,ce,ng):
    ip=g["indptr"][cs:ce+1]; p0,p1=int(ip[0]),int(ip[-1])
    return sp.csr_matrix((g["data"][p0:p1],g["indices"][p0:p1],ip-p0),shape=(ce-cs,ng))
def lognorm(Mc,lib,target=1e4):
    lib=np.asarray(lib,dtype=np.float32).ravel().copy(); lib[lib==0]=1
    D=Mc.toarray().astype(np.float32); D/=lib[:,None]; D*=target; np.log1p(D,out=D); return D

f=h5py.File(FULL,"r"); Ug=f["layers"]["UMIs"]; ng=int(f["X"].attrs["shape"][1])   # UMIs = raw counts (X is pre-lognormalized)
vi=f["var"].attrs.get("_index","_index"); gfull=np.array([g.decode() if isinstance(g,bytes) else str(g) for g in f["var"][vi][:]])
gpos={g:i for i,g in enumerate(gfull)}
miss=[g for g in pg+hv if g not in gpos]; assert not miss, f"genes missing from full dataset: {miss[:5]}"
panel_idx=np.array([gpos[g] for g in pg]); hv_idx=np.array([gpos[g] for g in hv])
kept={}; Xqp_p=[]; Xqh_p=[]; meta_p=[]
for cs in range(0,SPAN,CHUNK):
    ce=min(cs+CHUNK,SPAN); tic(f"chunk rows {cs:,}-{ce:,}...")
    dser=obs_col(f,"Donor ID",cs,ce); adncc=obs_col(f,"Overall AD neuropathological Change",cs,ce)
    didx=np.where(adncc!="Reference")[0]; rng=np.random.default_rng(cs); rng.shuffle(didx)
    take=[i for i in didx if kept.get(dser[i],0)<CAP]
    for i in take: kept[dser[i]]=kept.get(dser[i],0)+1
    if not take: continue
    loc=np.sort(np.array(take)); M=read_csr_rows(Ug,cs,ce,ng)[loc]
    Xqp_p.append(lognorm(M[:,panel_idx],np.asarray(M[:,panel_idx].sum(1)).ravel()))
    Xqh_p.append(lognorm(M[:,hv_idx],np.asarray(M.sum(1)).ravel())); del M
    supc=obs_col(f,"Supertype",cs,ce); subc=obs_col(f,"Subclass",cs,ce); cpsc=obs_col(f,"Continuous Pseudo-progression Score",cs,ce)
    meta_p.append(pd.DataFrame({"donor":dser[loc],"CPS":pd.to_numeric(pd.Series(cpsc[loc]),errors="coerce").values,
        "ADNC":adncc[loc],"gt_sub":subc[loc],"gt_sup":supc[loc]}))
f.close()
Xqp=np.vstack(Xqp_p); Xqh=np.vstack(Xqh_p); meta=pd.concat(meta_p,ignore_index=True)
gt_sub=meta.gt_sub.values; gt_sup=meta.gt_sup.values; donor=meta.donor.values; CPS=meta.CPS.values; ADNC=meta.ADNC.values
tic(f"query: {len(meta):,} disease cells, {len(set(donor))} donors; classifying...")

P={"sub_panel":classify(Xqp,sub_tp,Cp_sub),"sub_ceil":classify(Xqh,sub_th,Ch_sub),
   "sup_panel":classify(Xqp,sup_tp,Cp_sup),"sup_ceil":classify(Xqh,sup_th,Ch_sup)}
df=pd.DataFrame({"donor":donor,"CPS":CPS,"ADNC":ADNC,"gt_sub":gt_sub,"gt_sup":gt_sup,
    "c_sub_panel":P["sub_panel"]==gt_sub,"c_sub_ceil":P["sub_ceil"]==gt_sub,
    "c_sup_panel":P["sup_panel"]==gt_sup,"c_sup_ceil":P["sup_ceil"]==gt_sup})
tic("classified.")

P={"sub_panel":classify(Xqp,sub_tp,Cp_sub),"sub_ceil":classify(Xqh,sub_th,Ch_sub),
   "sup_panel":classify(Xqp,sup_tp,Cp_sup),"sup_ceil":classify(Xqh,sup_th,Ch_sup)}
df=pd.DataFrame({"donor":donor,"CPS":CPS,"ADNC":ADNC,"gt_sub":gt_sub,"gt_sup":gt_sup,
    "c_sub_panel":P["sub_panel"]==gt_sub,"c_sub_ceil":P["sub_ceil"]==gt_sub,
    "c_sup_panel":P["sup_panel"]==gt_sup,"c_sup_ceil":P["sup_ceil"]==gt_sup})
tic("classified.")

# ════ 3. per-donor aggregation + CPS regression ════
g=df.groupby("donor").agg(n=("CPS","size"),CPS=("CPS","first"),ADNC=("ADNC","first"),
    sub_panel=("c_sub_panel","mean"),sub_ceil=("c_sub_ceil","mean"),
    sup_panel=("c_sup_panel","mean"),sup_ceil=("c_sup_ceil","mean"))
g["sub_gap"]=g.sub_ceil-g.sub_panel; g["sup_gap"]=g.sup_ceil-g.sup_panel
g.sort_values("CPS").round(4).to_csv(os.path.join(OUT,"disease_degradation_per_donor.csv"))

def reg(y): r=linregress(g.CPS,g[y]); return r.slope,r.pvalue,r.rvalue
print("\n=== per-donor accuracy ~ CPS (slope = change in accuracy per unit CPS; CPS spans ~0.8) ===")
for y in ["sub_panel","sub_ceil","sup_panel","sup_ceil","sub_gap","sup_gap"]:
    s,p,r=reg(y); print(f"  {y:10s}: slope={s:+.3f}  p={p:.1e}  r={r:+.2f}  (mean {g[y].mean():.3f})")
print("\n=== accuracy by ADNC (pooled cells) ===")
ord_=["Not AD","Low","Intermediate","High"]
ad_=df.groupby("ADNC").agg(sub_panel=("c_sub_panel","mean"),sub_ceil=("c_sub_ceil","mean"),
    sup_panel=("c_sup_panel","mean"),sup_ceil=("c_sup_ceil","mean"),n=("CPS","size")).reindex(ord_)
print(ad_.round(3).to_string())
ad_.round(4).to_csv(os.path.join(OUT,"disease_degradation_by_adnc.csv"))

# per-class supertype degradation (low vs high CPS quartile)
lo=df[df.CPS<=g.CPS.quantile(.25)]; hi=df[df.CPS>=g.CPS.quantile(.75)]
print(f"\n=== supertype accuracy: low-CPS vs high-CPS (panel vs ceiling) ===")
print(f"  PANEL : low={lo.c_sup_panel.mean():.3f}  high={hi.c_sup_panel.mean():.3f}  drop={lo.c_sup_panel.mean()-hi.c_sup_panel.mean():+.3f}")
print(f"  CEIL  : low={lo.c_sup_ceil.mean():.3f}  high={hi.c_sup_ceil.mean():.3f}  drop={lo.c_sup_ceil.mean()-hi.c_sup_ceil.mean():+.3f}")

# ════ 4. figure ════
fig,ax=plt.subplots(1,3,figsize=(15.5,5.2))
def scat_reg(axx,yp,yc,ylab,ttl):
    for y,c,lab in [(yp,"#4393c3","panel (300g)"),(yc,"#222222","ceiling (HVG)")]:
        axx.scatter(g.CPS,g[y],s=22,c=c,alpha=.8,edgecolor="white",linewidth=.3,zorder=3,label=lab)
        r=linregress(g.CPS,g[y]); xs=np.array([g.CPS.min(),g.CPS.max()])
        axx.plot(xs,r.intercept+r.slope*xs,"-",color=c,lw=1.6,zorder=4)
        axx.text(.04,.10 if c=="#4393c3" else .03,f"{lab.split()[0]} slope={r.slope:+.2f} (p={r.pvalue:.0e})",
                 transform=axx.transAxes,fontsize=8,color=c)
    axx.set_xlabel("donor CPS (AD pseudo-progression)"); axx.set_ylabel(ylab); axx.set_title(ttl)
    axx.legend(loc="upper right",frameon=False); axx.spines[["top","right"]].set_visible(False)
scat_reg(ax[0],"sup_panel","sup_ceil","per-donor supertype accuracy","Supertype accuracy falls with AD pathology")
# gap vs CPS
r=linregress(g.CPS,g.sup_gap)
ax[1].scatter(g.CPS,g.sup_gap,s=24,c="#E65100",alpha=.85,edgecolor="white",linewidth=.3,zorder=3)
xs=np.array([g.CPS.min(),g.CPS.max()]); ax[1].plot(xs,r.intercept+r.slope*xs,"-",color="#E65100",lw=1.8,zorder=4)
ax[1].text(.04,.93,f"gap slope={r.slope:+.2f}\np={r.pvalue:.1e}, r={r.rvalue:+.2f}",transform=ax[1].transAxes,fontsize=9,va="top",color="#E65100")
ax[1].set_xlabel("donor CPS"); ax[1].set_ylabel("ceiling − panel accuracy (gap)")
ax[1].set_title("Does dropping markers hurt MORE\nin disease? (gap vs pathology)"); ax[1].spines[["top","right"]].set_visible(False)
# ADNC bars
x=np.arange(len(ord_)); w=0.38
ax[2].bar(x-w/2,ad_["sup_panel"],w,color="#4393c3",label="panel",edgecolor="black",linewidth=.4)
ax[2].bar(x+w/2,ad_["sup_ceil"],w,color="#222222",label="ceiling",edgecolor="black",linewidth=.4)
ax[2].set_xticks(x); ax[2].set_xticklabels(ord_,rotation=20,ha="right"); ax[2].set_ylim(0,1)
ax[2].set_ylabel("supertype accuracy"); ax[2].set_title("Accuracy by neuropathology stage"); ax[2].legend(frameon=False)
ax[2].spines[["top","right"]].set_visible(False)
fig.tight_layout()
for ext in ("png","pdf"): fig.savefig(os.path.join(OUT,f"disease_degradation_resolvability.{ext}"),dpi=200,bbox_inches="tight")
tic("saved disease_degradation_resolvability.png/.pdf + per-donor/adnc CSVs")
