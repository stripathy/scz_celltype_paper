#!/usr/bin/env python3
"""Visualize our layers vs the Lieber/Kwon spatial-domain layers."""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import anndata as ad

_HERE=os.path.dirname(os.path.abspath(__file__)); _SP=os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0,os.path.join(_SP,"code","analysis"))
from config import H5AD_DIR
LAB=os.path.join(_SP,"output/depth_validation/lieber_layers")
SPD2LAYER={"spd07":"L1/M","spd06":"L2/3","spd02":"L3/4","spd05":"L5","spd03":"L6","spd01":"WMtz","spd04":"WM"}
THEIR_ORDER=["L1/M","L2/3","L3/4","L5","L6","WMtz","WM"]          # their pia->WM
OUR_ORDER=["L1","L2/3","L4","L5","L6","WM"]
LCOL={"L1/M":"#D81B60","L2/3":"#1E88E5","L3/4":"#00ACC1","L5":"#43A047","L6":"#FB8C00","WMtz":"#8E24AA","WM":"#757575"}

j=pd.read_csv(os.path.join(LAB,"joined_cells.csv"))
j["their_layer"]=j["their_spd"].map(SPD2LAYER)

fig=plt.figure(figsize=(17,11),facecolor="white")
# A: our depth per their layer
axA=fig.add_subplot(2,2,1)
data=[j.loc[j.their_layer==L,"our_depth"].dropna().values for L in THEIR_ORDER]
bp=axA.violinplot(data,showmedians=True,vert=False)
for b,L in zip(bp['bodies'],THEIR_ORDER): b.set_facecolor(LCOL[L]); b.set_alpha(.8)
axA.set_yticks(range(1,len(THEIR_ORDER)+1)); axA.set_yticklabels(THEIR_ORDER)
axA.invert_yaxis(); axA.set_xlabel("our predicted depth (pia 0 -> 1 WM)",fontsize=12)
axA.set_title("A. Our depth within each Lieber layer\n(L1/M should be ~0 but isn't)",fontsize=13)
axA.axvline(0.12,ls=":",c="grey"); axA.grid(axis="x",alpha=.3)

# B: confusion their_layer x our_layer (row-normalized %)
axB=fig.add_subplot(2,2,2)
cm=j[j.our_layer.isin(OUR_ORDER) & j.their_layer.notna()]
ct=pd.crosstab(cm.their_layer,cm.our_layer).reindex(index=THEIR_ORDER,columns=OUR_ORDER).fillna(0)
pct=ct.div(ct.sum(1),axis=0)*100
im=axB.imshow(pct.values,cmap="viridis",aspect="auto",vmin=0,vmax=100)
axB.set_xticks(range(len(OUR_ORDER))); axB.set_xticklabels(OUR_ORDER)
axB.set_yticks(range(len(THEIR_ORDER))); axB.set_yticklabels(THEIR_ORDER)
axB.set_xlabel("OUR layer",fontsize=12); axB.set_ylabel("Lieber layer",fontsize=12)
for i in range(len(THEIR_ORDER)):
    for k in range(len(OUR_ORDER)):
        v=pct.values[i,k]
        if v>=8: axB.text(k,i,f"{v:.0f}",ha="center",va="center",color="white" if v<60 else "black",fontsize=9)
axB.set_title("B. Lieber layer vs our layer (row %)",fontsize=13); fig.colorbar(im,ax=axB,fraction=.04,label="% of Lieber layer")

# C/D: spatial for one sample, our layer vs their layer
sid="Br8667"; a=ad.read_h5ad(os.path.join(H5AD_DIR,f"{sid}_annotated.h5ad"))
xy=a.obsm["spatial"][:,:2]; n=a.n_obs
ids=pd.Series([f"{sid}_{i+1}" for i in range(n)])
their=pd.read_csv(os.path.join(LAB,"label_transfer_N24_k50_smoothed_labels.csv")); their.columns=["cellid","spd"]
sp=ids.map(their.set_index("cellid")["spd"]).map(SPD2LAYER).values
ourL=a.obs["layer"].values.astype(str)
LAYER_COLORS={"L1":(.9,.3,.3),"L2/3":(.3,.8,.3),"L4":(.3,.3,.9),"L5":(.9,.6,.1),
              "L6":(.7,.3,.8),"WM":(.5,.5,.5),"Vascular":(.85,.85,.2)}
axC=fig.add_subplot(2,2,3); axD=fig.add_subplot(2,2,4)
for L in OUR_ORDER+["Vascular"]:
    m=ourL==L
    if m.sum(): axC.scatter(xy[m,0],xy[m,1],s=1,color=LAYER_COLORS.get(L,(.5,.5,.5)),rasterized=True,lw=0,label=L)
axC.set_title(f"C. {sid} — OUR layers",fontsize=13); axC.legend(markerscale=5,fontsize=7,loc="upper right")
for L in THEIR_ORDER:
    m=sp==L
    if m.sum(): axD.scatter(xy[m,0],xy[m,1],s=1,color=LCOL[L],rasterized=True,lw=0,label=L)
axD.set_title(f"D. {sid} — LIEBER layers (note L1/M (pink) is scattered, not at pia)",fontsize=12)
axD.legend(markerscale=5,fontsize=7,loc="upper right")
for ax in (axC,axD): ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.invert_yaxis()
plt.tight_layout()
op=os.path.join(LAB,"lieber_vs_ours_comparison.png"); fig.savefig(op,dpi=140,bbox_inches="tight"); print("saved",op)
