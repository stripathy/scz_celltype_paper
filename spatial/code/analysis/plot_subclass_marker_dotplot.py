#!/usr/bin/env python3
"""
Block-diagonal marker dot plot: Xenium-panel genes that best distinguish each
SEA-AD subclass, in the Xenium data itself (corr_subclass, cortical QC-pass cells).

Markers are picked greedily by enrichment (z-scored mean expression across
subclasses), mostly unique per subclass; subclasses and their markers are ordered
biologically (excitatory by layer -> inhibitory -> non-neuronal) so the high-
expression dots fall on the diagonal. Dot size = % expressing, color = mean
expression scaled 0-1 per gene.
"""
import os, sys
import numpy as np, pandas as pd
import anndata as ad, scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- house style (Fig-09 spirit): legible sized sans fonts, editable PDF text ----
matplotlib.rcParams.update({
    "font.family":"sans-serif","font.sans-serif":["Helvetica","Arial","DejaVu Sans"],
    "font.size":8,"axes.titlesize":10,"axes.labelsize":9,"xtick.labelsize":7,
    "ytick.labelsize":8,"legend.fontsize":7,"legend.title_fontsize":8,
    "axes.linewidth":0.5,"pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none"})

def save_styled(fig, stem):
    """italicize gene (x-tick) labels and write png+pdf."""
    try:
        gax=max(fig.axes,key=lambda a:len(a.get_xticklabels()))
        for t in gax.get_xticklabels(): t.set_fontstyle("italic")
    except Exception: pass
    fig.savefig(stem+".png",dpi=300,bbox_inches="tight")
    fig.savefig(stem+".pdf",bbox_inches="tight")
    print("saved",stem+".png/.pdf")

_HERE=os.path.dirname(os.path.abspath(__file__)); _SP=os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0,os.path.join(_SP,"code","analysis"))
from config import H5AD_DIR, EXCLUDE_SAMPLES
OUT=os.path.join(_SP,"output/depth_validation/lieber_layers")   # reuse fig dir

EXC=["L2/3 IT","L4 IT","L5 IT","L5 ET","L5/6 NP","L6 IT","L6 IT Car3","L6 CT","L6b"]      # excitatory
INH=["Lamp5","Lamp5 Lhx6","Sncg","Vip","Pax6","Sst","Sst Chodl","Pvalb","Chandelier"]    # inhibitory
NN =["Astrocyte","Oligodendrocyte","OPC","Microglia-PVM","Endothelial","VLMC"]            # non-neuronal
ORDER=INH+EXC+NN   # inhibitory first, then excitatory, then non-neuronal
PAN_INH=["GAD1","GAD2","SLC32A1"]    # pan-inhibitory (GABAergic) anchor; pan-excitatory is data-driven (no SLC17A7/SATB2 in panel)
K=2            # markers per subclass
MIN_FRAC=0.33  # prioritize markers expressed in > 33% of the expected subclass
MIN_CELLS=50

# ---- load cortical QC-pass Xenium cells (expression + corr_subclass) ----
parts=[]
for f in sorted(os.listdir(H5AD_DIR)):
    if not f.endswith("_annotated.h5ad"): continue
    sid=f.replace("_annotated.h5ad","")
    if sid in EXCLUDE_SAMPLES: continue
    a=ad.read_h5ad(os.path.join(H5AD_DIR,f))
    m=(a.obs.get("corr_qc_pass",pd.Series(True,index=a.obs.index)).astype(bool)
       & (a.obs["spatial_domain"].astype(str)=="Cortical"))
    parts.append(a[m.values,:].copy())
adata=ad.concat(parts); adata.obs["corr_subclass"]=adata.obs["corr_subclass"].astype(str)
sc.pp.normalize_total(adata,target_sum=1e4); sc.pp.log1p(adata)
print(f"cells={adata.n_obs:,}  genes={adata.n_vars}",flush=True)

subs=[s for s in ORDER if (adata.obs.corr_subclass==s).sum()>=MIN_CELLS]
# mean log-norm expression per subclass x gene
X=adata.X
mean=pd.DataFrame(index=subs,columns=adata.var_names,dtype=float)
frac=pd.DataFrame(index=subs,columns=adata.var_names,dtype=float)
for s in subs:
    sub=adata[adata.obs.corr_subclass==s]
    xs=sub.X
    mean.loc[s]=np.asarray(xs.mean(0)).ravel()
    frac.loc[s]=np.asarray((xs>0).mean(0)).ravel()
# z-score each gene across subclasses -> enrichment
z=(mean-mean.mean(0))/(mean.std(0)+1e-9)
# exclude non-marker genes (mitochondrial / pseudo-mito) — high everywhere, not cell-type specific
bad=[g for g in adata.var_names if g.startswith(("MT-","MTRNR","MT_"))]
z.loc[:, bad]= -1e9
print(f"excluded {len(bad)} mito-type genes from marker candidates: {bad}",flush=True)

# greedy selection: each subclass takes its top-K most-enriched UNUSED genes,
# PRIORITIZING genes detected in > MIN_FRAC of that subclass; fall back to the
# best-enriched gene only if fewer than K clear the threshold.
# pan-inhibitory anchor only; the panel has no clean pan-excitatory gene (no SLC17A7/SLC17A6/SATB2)
PAN_INH_p=[g for g in PAN_INH if g in adata.var_names]
print(f"pan-Inh anchor: {PAN_INH_p}",flush=True)
used=set(PAN_INH_p); markers={}   # anchor shown separately, not re-picked per subclass
for s in subs:
    ranked=[g for g in z.loc[s].sort_values(ascending=False).index if g not in used]
    picks=[g for g in ranked if frac.loc[s,g]>MIN_FRAC][:K]        # prioritize >33% expressing
    if len(picks)<K:
        picks+=[g for g in ranked if g not in picks][:K-len(picks)]  # fill with best-enriched
    if not picks: picks=[z.loc[s].sort_values(ascending=False).index[0]]
    for g in picks: used.add(g)
    markers[s]=picks
pd.DataFrame([(s,g,round(z.loc[s,g],2),round(frac.loc[s,g],3)) for s in subs for g in markers[s]],
            columns=["subclass","gene","enrichment_z","pct_expressing"]).to_csv(os.path.join(OUT,"subclass_marker_dotplot.csv"),index=False)
print("markers/subclass:",{s:markers[s] for s in subs},flush=True)

# ---- dot plot (scanpy), diagonal ordering ----
adata.obs["corr_subclass"]=pd.Categorical(adata.obs["corr_subclass"],categories=subs,ordered=True)
# insert pan-inhibitory block at the start of the inhibitory section (no pan-excitatory block)
markers_full={}
for s in subs:
    if s in INH and "pan-Inh" not in markers_full: markers_full["pan-Inh"]=PAN_INH_p
    markers_full[s]=markers[s]

# ---- export tidy long data for the ggplot/cowplot figure (R renders) ----
SECTION={"pan-Inh":"pan-Inh", **{s:"Inhibitory" for s in INH}, **{s:"Excitatory" for s in EXC}, **{s:"Non-neuronal" for s in NN}}
_genes=[g for grp in markers_full for g in markers_full[grp]]
_sm=mean[_genes]; _sm=(_sm-_sm.min(0))/(_sm.max(0)-_sm.min(0)+1e-9)   # per-gene 0-1 scaled mean
_rows=[]; _gi=0
for grp in markers_full:
    for g in markers_full[grp]:
        for ci,ct in enumerate(subs):
            _rows.append(dict(cell_type=ct,gene=g,section=SECTION.get(grp,grp),gene_owner=grp,
                              pct=float(frac.loc[ct,g]),scaled_mean=float(_sm.loc[ct,g]),ct_idx=ci,gene_idx=_gi))
        _gi+=1
pd.DataFrame(_rows).to_csv(os.path.join(OUT,"subclass_dotplot_long.csv"),index=False)
print(f"exported subclass_dotplot_long.csv ({len(_rows)} rows, {_gi} genes, {len(subs)} subclasses)",flush=True)

ngenes=sum(len(v) for v in markers_full.values())
sc.pl.dotplot(adata, markers_full, groupby="corr_subclass", categories_order=subs,
              standard_scale="var", cmap="Reds", show=False,
              figsize=(0.27*ngenes+3, 0.34*len(subs)+1.5),
              dot_max=0.9, colorbar_title="scaled mean\nexpression", size_title="% expressing")
save_styled(plt.gcf(), os.path.join(OUT,"subclass_marker_dotplot"))
