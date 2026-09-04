#!/usr/bin/env python3
"""
Validation side-check: does Sst_25 genuinely express GAD2 at very low levels
relative to the other Sst supertypes, in clean snRNAseq?

Motivation: in the Xenium SST supertype marker dot plot (panel b), Sst_25 showed
conspicuously low GAD2. Xenium is a 300-gene targeted assay, so low signal could
be a detection/classification artifact. The SEA-AD neurotypical-MTG snRNAseq
reference (whole-transcriptome, high-quality) is the ground truth.

For the 16 Sst supertypes (same set as the Xenium panel) we compute per-supertype
mean normalized expression and % expressing for GAD1, GAD2, SLC32A1. If Sst_25 is
low for GAD2 *specifically* (and normal for GAD1 / SLC32A1 in the same cells), the
effect is real biology and not a global quality/depth confound.
"""
import os
import numpy as np, pandas as pd
import anndata as ad, scanpy as sc, scipy.sparse as sp
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

matplotlib.rcParams.update({
    "font.family":"sans-serif","font.sans-serif":["Helvetica","Arial","DejaVu Sans"],
    "font.size":12,"axes.titlesize":15,"axes.labelsize":13,"xtick.labelsize":11,
    "ytick.labelsize":12,"legend.fontsize":11,"axes.linewidth":0.6,
    "pdf.fonttype":42,"ps.fonttype":42})

REF=os.path.expanduser("~/Github/shared_data/nicole_sea_ad_snrnaseq_reference.h5ad")
_HERE=os.path.dirname(os.path.abspath(__file__)); _SP=os.path.dirname(os.path.dirname(_HERE))
OUT=os.path.join(_SP,"output/depth_validation/lieber_layers")
GENES=["GAD1","GAD2","SLC32A1"]; GCOL={"GAD1":"#1b9e77","GAD2":"#d95f02","SLC32A1":"#7570b3"}
HILITE="Sst_25"

# ---- load Sst cells, normalize ----
a=ad.read_h5ad(REF, backed="r")
st=a.obs["Supertype"].astype(str)
sub=a[st.str.match(r"^Sst_\d+$").values].to_memory()
sc.pp.normalize_total(sub,target_sum=1e4); sc.pp.log1p(sub)
order=sorted(sub.obs["Supertype"].astype(str).unique(), key=lambda s:int(s.split("_")[1]))
def col(cells,g):
    x=cells[:,g].X; return (x.toarray().ravel() if sp.issparse(x) else np.asarray(x).ravel())

rows=[]
for s in order:
    cells=sub[sub.obs["Supertype"].astype(str)==s]; rec={"supertype":s,"n":cells.n_obs}
    for g in GENES:
        x=col(cells,g); rec[f"{g}_mean"]=float(x.mean()); rec[f"{g}_pct"]=float((x>0).mean()*100)
    rows.append(rec)
df=pd.DataFrame(rows).set_index("supertype")
df.round(3).to_csv(os.path.join(OUT,"sst25_gad2_snrnaseq_check.csv"))

# order x-axis by GAD2 mean (desc) so the Sst_25 outlier sits at the far right
xo=df.sort_values("GAD2_mean",ascending=False).index.tolist()
print("SEA-AD snRNAseq — GAD1/GAD2/SLC32A1 by Sst supertype:\n",df.round(3).to_string())
oth=df.drop(HILITE)
for g in GENES:
    print(f"  {g}: {HILITE} mean={df.loc[HILITE,g+'_mean']:.3f} ({df.loc[HILITE,g+'_pct']:.0f}%)  | "
          f"other median mean={oth[g+'_mean'].median():.3f} ({oth[g+'_pct'].median():.0f}%)")

# ---- figure: grouped bars, mean (top) + %expressing (bottom); Sst_25 highlighted ----
fig,axes=plt.subplots(2,1,figsize=(11,7.2),sharex=True)
x=np.arange(len(xo)); w=0.27
for ax,metric,ylab in [(axes[0],"mean","mean norm. expr.\n(log1p CP10K)"),
                       (axes[1],"pct","% cells expressing")]:
    for j,g in enumerate(GENES):
        vals=df.loc[xo,f"{g}_{metric}"].values
        ax.bar(x+(j-1)*w,vals,w,color=GCOL[g],label=g if ax is axes[0] else None,
               edgecolor="black",linewidth=0.4)
    ax.set_ylabel(ylab); ax.spines[["top","right"]].set_visible(False)
    ax.axvspan(list(xo).index(HILITE)-0.5,list(xo).index(HILITE)+0.5,color="gold",alpha=0.22,zorder=0)
axes[0].legend(handles=[Patch(facecolor=GCOL[g],edgecolor="black",label=g) for g in GENES],
               title=None,frameon=False,ncol=3,loc="upper right")
axes[0].set_title("Sst_25 specifically downregulates GAD2 (SEA-AD snRNAseq, n=12,393 Sst cells)")
axes[1].set_xticks(x); axes[1].set_xticklabels(xo,rotation=90)
fig.tight_layout(); fig.canvas.draw()   # populate tick-label artists before styling
hi=list(xo).index(HILITE)
lab=axes[1].get_xticklabels()[hi]; lab.set_fontweight("bold"); lab.set_color("#b8860b")
fig.savefig(os.path.join(OUT,"sst25_gad2_snrnaseq_check.png"),dpi=200,bbox_inches="tight")
fig.savefig(os.path.join(OUT,"sst25_gad2_snrnaseq_check.pdf"),bbox_inches="tight")
print("saved",os.path.join(OUT,"sst25_gad2_snrnaseq_check.png/.pdf"))
