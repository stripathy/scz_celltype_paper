#!/usr/bin/env python3
"""Why do some supertypes (esp. non-neuronal) have poor Xenium-vs-MERFISH depth concordance?
Test: concordance error vs (a) within-supertype depth spread (non-laminar?), (b) cell count, (c) class."""
import os, sys
import numpy as np, pandas as pd
import anndata as ad
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

_HERE=os.path.dirname(os.path.abspath(__file__)); _SP=os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0,os.path.join(_SP,"code","analysis"))
from config import MERFISH_PATH, CORTICAL_LAYERS
LAB=os.path.join(_SP,"output/depth_validation/lieber_layers")
DD=os.path.join(_SP,"output/celltyping_supplement/data_bottomthird")
HELDOUT=["H21.33.011","H21.33.015","H21.33.028"]   # the 3 held-out MERFISH donors used in the comparison

df=pd.read_csv(os.path.join(DD,"depth_supertype.csv"))   # celltype, merfish_depth, xenium_depth, n_merfish_cells, n_xenium_cells, klass
df["abs_resid"]=(df.xenium_depth-df.merfish_depth).abs()
df["min_n"]=df[["n_merfish_cells","n_xenium_cells"]].min(1)

# within-supertype depth SPREAD in the MERFISH reference (held-out donors, cortical) = "how laminar"
a=ad.read_h5ad(MERFISH_PATH,backed="r"); o=a.obs
mk=(o["Donor ID"].astype(str).isin(HELDOUT)) & (o["Normalized depth from pia"].notna()) \
   & (o["Layer annotation"].astype(str).isin(CORTICAL_LAYERS))
m=pd.DataFrame({"st":o.loc[mk,"Supertype"].astype(str).values,
                "d":o.loc[mk,"Normalized depth from pia"].astype(float).values})
spread=m.groupby("st")["d"].agg(depth_sd="std",
        depth_iqr=lambda x:x.quantile(.75)-x.quantile(.25))
df=df.merge(spread,left_on="celltype",right_index=True,how="left")

print(f"=== 12 worst-concordance supertypes (|Xenium-MERFISH median depth|) ===")
cols=["celltype","klass","abs_resid","merfish_depth","xenium_depth","n_merfish_cells","n_xenium_cells","depth_sd"]
print(df.sort_values("abs_resid",ascending=False).head(12)[cols].to_string(index=False))

print("\n=== drivers ===")
for v in ["depth_sd","depth_iqr"]:
    r=pearsonr(df[v].fillna(df[v].mean()),df.abs_resid)[0]
    print(f"  abs_resid vs {v} (non-laminar spread): r={r:+.2f}")
r=pearsonr(np.log10(df.min_n),df.abs_resid)[0]
print(f"  abs_resid vs log10(min cell count): r={r:+.2f}")
print("\n  mean abs_resid by class:")
print(df.groupby("klass").agg(mean_abs_resid=("abs_resid","mean"),
      median_depth_sd=("depth_sd","median"),median_min_n=("min_n","median"),n_supertypes=("celltype","size")).round(3).to_string())

# figure
fig,ax=plt.subplots(1,2,figsize=(14,6))
CC={"GABAergic":"#2E7D32","Glutamatergic":"#E65100","Non-neuronal":"#1565C0"}
for kl,g in df.groupby("klass"):
    ax[0].scatter(g.depth_sd,g.abs_resid,c=CC.get(kl,"grey"),label=kl,s=30,alpha=.8)
    ax[1].scatter(g.min_n,g.abs_resid,c=CC.get(kl,"grey"),label=kl,s=30,alpha=.8)
for _,r in df.sort_values("abs_resid",ascending=False).head(8).iterrows():
    ax[0].annotate(r.celltype,(r.depth_sd,r.abs_resid),fontsize=7)
    ax[1].annotate(r.celltype,(r.min_n,r.abs_resid),fontsize=7)
ax[0].set_xlabel("within-supertype depth SD (MERFISH) — higher = less laminar"); ax[0].set_ylabel("|Xenium - MERFISH| median depth")
ax[1].set_xlabel("min cell count (Xenium/MERFISH)"); ax[1].set_xscale("log"); ax[1].set_ylabel("|Xenium - MERFISH| median depth")
ax[0].legend(); ax[0].set_title("Concordance error vs laminarity"); ax[1].set_title("Concordance error vs cell count")
plt.tight_layout(); op=os.path.join(LAB,"supertype_depth_concordance_drivers.png")
fig.savefig(op,dpi=140,bbox_inches="tight"); print(f"\nsaved {op}")
