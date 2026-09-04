#!/usr/bin/env python3
"""
Is Sst_25 the most resolvable Sst supertype — ON THE XENIUM PANEL specifically?

Motivation: Sst_25 has a distinctive Xenium-panel signature (notably GAD2-low
while GAD1/SST/SLC32A1 stay normal — see check_sst25_gad2_snrnaseq.py). A
correlation classifier keys on POSITIVE markers, so the real question is whether
Sst_25's overall profile is more separable than the other 15 Sst supertypes when
restricted to the genes Xenium actually measures.

Two ground-truth tests (true SEA-AD Supertype labels in both):
  (A) MERFISH benchmark — the DEPLOYED correlation classifier on 1.89M real
      spatial cells (merfish_reclassified.h5ad). CAVEAT: the 180-gene MERFISH
      panel LACKS GAD1 and SST (the pan-Sst anchors); only 23 genes overlap
      Xenium. So it is the WRONG panel for this question, and indeed Sst_25
      resolves poorly there (leaks out of the Sst subclass).
  (B) snRNAseq restricted to the 300 Xenium panel genes — correlation-centroid
      classifier, leave-one-donor-out CV. The Xenium-faithful test. Includes a
      within-Sst 16-way variant and a GAD2-drop ablation to isolate GAD2's role.

Outputs: sst25_resolvability.csv + sst25_resolvability.png/.pdf in
output/depth_validation/lieber_layers/.
"""
import os
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, scipy.sparse as sp
from sklearn.model_selection import GroupKFold
from sklearn.metrics import classification_report
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

matplotlib.rcParams.update({
    "font.family":"sans-serif","font.sans-serif":["Helvetica","Arial","DejaVu Sans"],
    "font.size":11,"axes.titlesize":12,"axes.labelsize":11,"xtick.labelsize":10,
    "ytick.labelsize":10,"legend.fontsize":9,"axes.linewidth":0.6,
    "pdf.fonttype":42,"ps.fonttype":42})

BASE=os.path.expanduser("~/Github/SCZ_Xenium")
MERF=os.path.join(BASE,"output/merfish_benchmark/merfish_reclassified.h5ad")
SNRNA=os.path.expanduser("~/Github/shared_data/nicole_sea_ad_snrnaseq_reference.h5ad")
OUT=os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                 "output/depth_validation/lieber_layers")
HI="Sst_25"; GOLD="#b8860b"

def corr_classify(Xtr,ytr,Xte,cand=None):
    """Nearest-centroid by Pearson correlation across genes (mirrors the pipeline classifier)."""
    types=[t for t in (cand if cand is not None else sorted(set(ytr))) if (ytr==t).sum()>0]
    C=np.vstack([Xtr[ytr==t].mean(0) for t in types])
    Cc=C-C.mean(1,keepdims=True); Cc/=(np.linalg.norm(Cc,axis=1,keepdims=True)+1e-9)
    Xc=Xte-Xte.mean(1,keepdims=True); Xc/=(np.linalg.norm(Xc,axis=1,keepdims=True)+1e-9)
    return np.array(types)[(Xc@Cc.T).argmax(1)]

def report_rows(gt,pred,types,tag):
    rep=classification_report(gt,pred,labels=types,output_dict=True,zero_division=0)
    return pd.DataFrame([{"supertype":s,f"recall_{tag}":rep[s]["recall"],
        f"prec_{tag}":rep[s]["precision"],f"F1_{tag}":rep[s]["f1-score"]} for s in types]).set_index("supertype")

# ── Xenium panel ──
xh=os.path.join(BASE,"output/h5ad")
panel=list(ad.read_h5ad(os.path.join(xh,[f for f in sorted(os.listdir(xh)) if f.endswith("_annotated.h5ad")][0]),
                        backed="r").var_names)

# ════ (A) MERFISH benchmark — deployed classifier, real spatial, wrong panel ════
mo=ad.read_h5ad(MERF,backed="r").obs
mgt=mo["Supertype"].astype(str).values; mpred=mo["corr_supertype"].astype(str).values
sst=sorted([s for s in set(mgt) if s.startswith("Sst_")],key=lambda s:int(s.split("_")[1]))
dM=report_rows(mgt,mpred,sst,"merf")
for s in sst: dM.loc[s,"stay_merf"]=np.isin(mpred[mgt==s],sst).mean()
print(f"(A) MERFISH: {len(sst)} Sst supertypes, {np.isin(mgt,sst).sum():,} GT Sst cells. "
      f"Sst_25 F1={dM.loc[HI,'F1_merf']:.3f} (rank {int(dM['F1_merf'].rank(ascending=False)[HI])}/16), "
      f"stay_in_Sst={dM.loc[HI,'stay_merf']:.2f}")

# ════ (B) snRNAseq on the Xenium panel — the faithful test ════
ref=ad.read_h5ad(SNRNA,backed="r")
genes=[g for g in panel if g in set(ref.var_names)]
sub=ref[:,genes].to_memory(); sc.pp.normalize_total(sub,target_sum=1e4); sc.pp.log1p(sub)
X=sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
ST=sub.obs["Supertype"].astype(str).values; donor=sub.obs["donor_id"].astype(str).values
gidx={g:i for i,g in enumerate(genes)}; keep_noG2=[i for g,i in gidx.items() if g!="GAD2"]
print(f"(B) snRNAseq: {len(genes)}/{len(panel)} panel genes present, {len(set(donor))} donors, "
      f"{np.isin(ST,sst).sum():,} Sst cells")

def cv_pred(Xmat,within_sst=False):
    pred=np.empty(len(ST),dtype=object); pred[:]=""
    for tr,te in GroupKFold(len(set(donor))).split(Xmat,ST,donor):
        if within_sst:
            mtr=np.isin(ST[tr],sst); mte=np.isin(ST[te],sst)
            pred[te[mte]]=corr_classify(Xmat[tr][mtr],ST[tr][mtr],Xmat[te[mte]],cand=sst)
        else:
            pred[te]=corr_classify(Xmat[tr],ST[tr],Xmat[te])
    return pred

pf=cv_pred(X); pn=cv_pred(X[:,keep_noG2]); pw=cv_pred(X,within_sst=True)
dX=report_rows(ST,pf,sst,"full").join(report_rows(ST,pn,sst,"noG2"))
mask_sst=np.isin(ST,sst)
dX=dX.join(report_rows(ST[mask_sst],pw[mask_sst],sst,"wsst")[["recall_wsst","F1_wsst"]])
dX["n_gt"]=[int((ST==s).sum()) for s in sst]
for s in sst: dX.loc[s,"stay_xen"]=np.isin(pf[ST==s],sst).mean()
dX["dF1_GAD2"]=dX["F1_full"]-dX["F1_noG2"]; dX["dPrec_GAD2"]=dX["prec_full"]-dX["prec_noG2"]

df=dX.join(dM); df.round(4).to_csv(os.path.join(OUT,"sst25_resolvability.csv"))
print(f"\nSst_25 on Xenium panel: F1={df.loc[HI,'F1_full']:.3f} (rank {int(df['F1_full'].rank(ascending=False)[HI])}/16), "
      f"precision={df.loc[HI,'prec_full']:.3f} (rank {int(df['prec_full'].rank(ascending=False)[HI])}/16), "
      f"recall={df.loc[HI,'recall_full']:.3f} (rank {int(df['recall_full'].rank(ascending=False)[HI])}/16)")
print(f"GAD2 ablation on Sst_25: ΔF1={df.loc[HI,'dF1_GAD2']:+.3f}, ΔPrec={df.loc[HI,'dPrec_GAD2']:+.3f} "
      f"(others median ΔF1={df.drop(HI)['dF1_GAD2'].median():+.3f})")

# Sst_25 top positive markers on the panel (z across the 16 Sst), + GAD2 for contrast
mean=pd.DataFrame({s:X[ST==s].mean(0) for s in sst},index=genes).T
z=(mean-mean.mean(0))/(mean.std(0)+1e-9)
top=z.loc[HI].sort_values(ascending=False).head(6)
print(f"\nSst_25 top positive panel markers (enrichment z across 16 Sst): "
      + ", ".join(f"{g}({z.loc[HI,g]:+.1f})" for g in top.index)
      + f"  | GAD2 z={z.loc[HI,'GAD2']:+.1f}")

# ════ figure ════
order=df.sort_values("F1_full").index.tolist(); yi=np.arange(len(order))
fig,ax=plt.subplots(1,3,figsize=(15.5,5.2))
# A: F1 ranking (Xenium panel)
cols=[GOLD if s==HI else "#9ecae1" for s in order]
ax[0].barh(yi,df.loc[order,"F1_full"],color=cols,edgecolor="black",linewidth=0.4)
ax[0].set_yticks(yi); ax[0].set_yticklabels(order,fontsize=9)
ax[0].get_yticklabels()[order.index(HI)].set_fontweight("bold"); ax[0].get_yticklabels()[order.index(HI)].set_color(GOLD)
for i,s in enumerate(order): ax[0].text(df.loc[s,"F1_full"]+.01,i,f"{df.loc[s,'F1_full']:.2f}",va="center",fontsize=8)
ax[0].set_xlabel("F1 (Xenium panel, leave-one-donor-out CV)"); ax[0].set_xlim(0,.95)
ax[0].set_title("Sst_25 is the best-resolved Sst\nsupertype on the Xenium panel"); ax[0].spines[["top","right"]].set_visible(False)
# B: precision vs recall (Xenium)
ax[1].scatter(df["recall_full"],df["prec_full"],s=20+df["n_gt"]/15,c="#9ecae1",edgecolor="black",linewidth=.5,zorder=3)
ax[1].scatter([df.loc[HI,"recall_full"]],[df.loc[HI,"prec_full"]],s=20+df.loc[HI,"n_gt"]/15,c=GOLD,edgecolor="black",linewidth=.8,zorder=4)
ax[1].annotate(HI,(df.loc[HI,"recall_full"],df.loc[HI,"prec_full"]),xytext=(6,-2),textcoords="offset points",
               fontweight="bold",color=GOLD,fontsize=10)
ax[1].set_xlabel("recall (sensitivity)"); ax[1].set_ylabel("precision (call purity)")
ax[1].set_title("Resolvability is precision-driven:\nSst_25 calls are the purest"); ax[1].spines[["top","right"]].set_visible(False)
ax[1].set_xlim(.35,.95); ax[1].set_ylim(.3,1.0)
# C: panel-dependence — stay-in-Sst MERFISH vs Xenium (paired)
for s in sst:
    c=GOLD if s==HI else "#bbbbbb"; lw=2.2 if s==HI else 0.8; z_=5 if s==HI else 2
    ax[2].plot([0,1],[df.loc[s,"stay_merf"],df.loc[s,"stay_xen"]],"-o",color=c,lw=lw,ms=4,zorder=z_)
ax[2].annotate(HI,(0,df.loc[HI,"stay_merf"]),xytext=(-6,0),textcoords="offset points",ha="right",fontweight="bold",color=GOLD,fontsize=10)
ax[2].set_xticks([0,1]); ax[2].set_xticklabels(["MERFISH panel\n(no GAD1/SST)","Xenium panel\n(GAD1+SST)"])
ax[2].set_ylabel("fraction of true cells kept in Sst subclass"); ax[2].set_xlim(-.35,1.25); ax[2].set_ylim(0,1.02)
ax[2].set_title("Why the panel matters: GAD1/SST\nanchor Sst_25 as Sst"); ax[2].spines[["top","right"]].set_visible(False)
fig.tight_layout()
for ext in ("png","pdf"): fig.savefig(os.path.join(OUT,f"sst25_resolvability.{ext}"),dpi=200,bbox_inches="tight")
print("\nsaved sst25_resolvability.png/.pdf + .csv")
