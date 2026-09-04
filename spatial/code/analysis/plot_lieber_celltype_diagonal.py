#!/usr/bin/env python3
"""Subclass x Kwon-cell-type concordance as a blocked-diagonal matrix.
SEA-AD subclasses on x (diagonal-ordered), Kwon et al. cell types on y
(biological groups, excitatory at top). Cell values are the % of each subclass's
cells assigned to each Kwon type (columns sum to ~100).

Styled to the manuscript convention: 7.1 in canvas, 7 pt base text, no in-plot
title, bold reserved for panel labels, neutral separators."""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.colors import LinearSegmentedColormap

_HERE=os.path.dirname(os.path.abspath(__file__)); _SP=os.path.dirname(os.path.dirname(_HERE))
LAB=os.path.join(_SP,"output/depth_validation/lieber_layers")
ct=pd.read_csv(os.path.join(LAB,"subclass_vs_celltype.csv"),index_col=0)
ct=ct.drop(index=[i for i in ["Unassigned"] if i in ct.index])           # drop our untyped
ct=ct.loc[:, ct.sum(0)>0]                                                # drop empty cols (NA)

# biological Lieber order: Ex upper->deep, interneurons, glia/vascular, ambiguous
COL_ORDER=["L2/3 Ex","L4/5 Ex","L5 Ex","L6 Ex","MGE","CGE","Ast","Oligo",
           "Ambig/In/Endo","Mic","Endo","Ambig/Oligo","NA"]
lieber=[c for c in COL_ORDER if c in ct.columns]
ct=ct.reindex(columns=lieber)
pct=ct.div(ct.sum(1),axis=0)*100                                          # % of our subclass

# order our subclasses by dominant Lieber type (then share) -> diagonal
argmax=pct.values.argmax(1); maxval=pct.values.max(1)
sub_order=np.lexsort((-maxval, argmax))
pct=pct.iloc[sub_order]
subs=list(pct.index)                                                      # x: our subclasses (diagonal order)
M=pct.values.T                                                            # rows=Lieber (y), cols=our subclass (x)

GROUPS=[("Excitatory",["L2/3 Ex","L4/5 Ex","L5 Ex","L6 Ex"]),
        ("Inhibitory",["MGE","CGE"]),
        ("Glia / vascular",["Ast","Oligo","Ambig/In/Endo","Mic","Endo"]),
        ("Ambiguous",["Ambig/Oligo","NA"])]
grp_of={c:gi for gi,(_,g) in enumerate(GROUPS) for c in g}

# ColorBrewer Reds-9 with the lowest stop (#fff5f0, a pale pink) swapped for pure
# white, so an off-diagonal zero is indistinguishable from the page background.
CMAP=LinearSegmentedColormap.from_list("Reds_from_white",
        ["#ffffff","#fee0d2","#fcbba1","#fc9272","#fb6a4a",
         "#ef3b2c","#cb181d","#a50f15","#67000d"])

BASE=7.0; AXIS_TITLE=BASE-0.5; AXIS_TEXT=BASE-1.0
matplotlib.rcParams.update({"font.size":BASE,"pdf.fonttype":42,
                            "axes.linewidth":0.4,"xtick.major.width":0.4,
                            "ytick.major.width":0.4})
fig,ax=plt.subplots(figsize=(7.1,4.0))
im=ax.imshow(M,cmap=CMAP,vmin=0,vmax=100,aspect="auto")
ax.set_xticks(range(len(subs))); ax.set_xticklabels(subs,rotation=45,ha="right",fontsize=AXIS_TEXT)
ax.set_yticks(range(len(lieber))); ax.set_yticklabels(lieber,fontsize=AXIS_TEXT)
ax.set_xlabel("SEA-AD subclass",fontsize=AXIS_TITLE,labelpad=4)
ax.set_ylabel("Kwon et al. cell type",fontsize=AXIS_TITLE)
# place explicitly: the rotated group headers occupy x = -0.175 axes-fraction,
# so the axis title has to sit outside them rather than rely on labelpad
ax.yaxis.set_label_coords(-0.235, 0.5)
for i in range(len(lieber)):
    for j in range(len(subs)):
        v=M[i,j]
        if v>=10: ax.text(j,i,f"{v:.0f}",ha="center",va="center",
                          color="black" if v<60 else "white",fontsize=BASE-2.5)
# horizontal separators between Lieber groups (y)
acc=0; row_bounds=[]
for _,g in GROUPS:
    acc+=sum(c in lieber for c in g)
    if 0<acc<len(lieber): row_bounds.append(acc)
for b in row_bounds: ax.axhline(b-0.5,color="grey20" if False else "#444",lw=0.5)
# vertical separators where our-subclass dominant Lieber group changes (x)
col_grp=[grp_of.get(lieber[a],99) for a in pct.values.argmax(1)]
for j in range(1,len(subs)):
    if col_grp[j]!=col_grp[j-1]: ax.axvline(j-0.5,color="#444",lw=0.5)
# Lieber group headers on the left (rotated), blended transform (x in axes frac, y in data)
trans=mtransforms.blended_transform_factory(ax.transAxes,ax.transData)
acc=0
for name,g in GROUPS:
    n=sum(c in lieber for c in g)
    if n: ax.text(-0.175,acc+n/2-0.5,name,transform=trans,rotation=90,ha="center",va="center",
                  fontsize=AXIS_TEXT,color="#333"); acc+=n
cb=fig.colorbar(im,fraction=0.022,pad=0.012)
cb.set_label("% of subclass",fontsize=AXIS_TITLE)
cb.ax.tick_params(labelsize=AXIS_TEXT,width=0.4)
cb.outline.set_linewidth(0.4)
FIGDIR=os.path.join(_SP,"supplemental_figures"); os.makedirs(FIGDIR,exist_ok=True)   # consolidated supplemental-figure home
matplotlib.rcParams["pdf.fonttype"]=42                                               # editable PDF text
op=os.path.join(FIGDIR,"celltype_comparison_diagonal")
fig.savefig(op+".png",dpi=400,bbox_inches="tight"); fig.savefig(op+".pdf",bbox_inches="tight")
print("saved",op+".png/.pdf  shape",M.shape)
