#!/usr/bin/env python3
"""
Traditional grouped marker dot plot for SST (+ Sst Chodl) supertypes, with
supertypes ordered by median predicted cortical depth (pia -> WM).

For each supertype the top-enriched Xenium-panel genes that are detected in
>=10% of that supertype's cells are shown (1-3 each). No cross-subclass or
enrichment-floor filters are applied: genes that genuinely co-vary with a
supertype are kept even if they are more highly expressed in another cell type,
since such signal can still aid classification. (Only mitochondrial QC genes are
excluded.) A pan-SST / GABAergic anchor block is shown first.

Dot size = % expressing; color = mean expression scaled 0-1 per gene.
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
OUT=os.path.join(_SP,"output/depth_validation/lieber_layers")
MIN_CELLS=150; KMAX=2; MIN_FRAC=0.10   # up to 2 markers/type, each detected in >=10% of that supertype
PAN_GENES=["SST","CALB1","GAD1","GAD2","SLC32A1","PVALB"]   # pan-SST / GABAergic anchors

# ---- load Sst + Sst Chodl cortical cells ----
parts=[]
for f in sorted(os.listdir(H5AD_DIR)):
    if not f.endswith("_annotated.h5ad"): continue
    sid=f.replace("_annotated.h5ad","")
    if sid in EXCLUDE_SAMPLES: continue
    a=ad.read_h5ad(os.path.join(H5AD_DIR,f))
    m=(a.obs.get("corr_qc_pass",pd.Series(True,index=a.obs.index)).astype(bool)
       & (a.obs["spatial_domain"].astype(str)=="Cortical")
       & (a.obs["corr_subclass"].astype(str)=="Sst")
       & a.obs["predicted_norm_depth"].notna())
    if m.any(): parts.append(a[m.values,:].copy())
adata=ad.concat(parts); adata.obs["st"]=adata.obs["corr_supertype"].astype(str)
sc.pp.normalize_total(adata,target_sum=1e4); sc.pp.log1p(adata)
print(f"Sst cells={adata.n_obs:,}  supertypes={adata.obs.st.nunique()}",flush=True)

# ---- per-supertype median depth, mean + frac expression; keep well-populated, sort pia->WM ----
keep=[s for s in adata.obs.st.unique() if (adata.obs.st==s).sum()>=MIN_CELLS]
mean=pd.DataFrame(index=keep,columns=adata.var_names,dtype=float); frac=mean.copy(); meta=[]
for s in keep:
    sub=adata[adata.obs.st==s]
    mean.loc[s]=np.asarray(sub.X.mean(0)).ravel()
    frac.loc[s]=np.asarray((sub.X>0).mean(0)).ravel()
    meta.append({"st":s,"n":sub.n_obs,"med_depth":float(np.median(sub.obs.predicted_norm_depth))})
meta=pd.DataFrame(meta).set_index("st").sort_values("med_depth")
# Sst supertypes depth-ordered (pia->WM); Sst Chodl (distinct branch) appended at the far right/bottom
order=([s for s in meta.index if not s.startswith("Sst Chodl")]
       + [s for s in meta.index if s.startswith("Sst Chodl")])
print(f"kept {len(order)} supertypes (n>={MIN_CELLS}); depth {meta.med_depth.min():.2f}-{meta.med_depth.max():.2f}",flush=True)

# ---- markers: enrichment z, exclude only mito, greedy top-K with >=10% expressing (depth order) ----
z=(mean-mean.mean(0))/(mean.std(0)+1e-9)
bad=[g for g in adata.var_names if g.startswith(("MT-","MTRNR","MT_"))]; z.loc[:,bad]=-1e9
PAN=[g for g in PAN_GENES if g in adata.var_names]
print(f"anchor (pan-SST / GABAergic) markers in panel: {PAN}",flush=True)
used=set(PAN); markers={}
for s in order:
    ranked=[g for g in z.loc[s].sort_values(ascending=False).index if g not in used]
    passing=[g for g in ranked if frac.loc[s,g]>=MIN_FRAC][:KMAX]   # >=10% expressing, up to 3
    if not passing: passing=[ranked[0]]                            # fallback: best available
    for g in passing: used.add(g)
    markers[s]=passing
pd.DataFrame([(s,g,round(z.loc[s,g],2),round(frac.loc[s,g],3),round(meta.med_depth[s],3),int(meta.n[s]))
              for s in order for g in markers[s]],
            columns=["supertype","gene","z","pct_expressing","supertype_med_depth","n"]
            ).to_csv(os.path.join(OUT,"sst_supertype_markers.csv"),index=False)

# ---- row labels carry the depth so the pia->WM ordering is explicit ----
lab={s:f"{s}  (depth {meta.med_depth[s]:.2f})" for s in order}
adata.obs["st_lab"]=pd.Categorical(adata.obs["st"].map(lab),categories=[lab[s] for s in order],ordered=True)
markers_full={"pan-SST / GABA": PAN, **markers}

# ---- export tidy long data for the ggplot/cowplot figure (R renders) ----
_genes=[g for grp in markers_full for g in markers_full[grp]]
_sm=mean.loc[order][_genes]; _sm=(_sm-_sm.min(0))/(_sm.max(0)-_sm.min(0)+1e-9)
_rows=[]; _gi=0
for grp in markers_full:
    for g in markers_full[grp]:
        for si,s in enumerate(order):
            _rows.append(dict(cell_type=s,gene=g,gene_owner=grp,
                              section=("anchor" if grp=="pan-SST / GABA" else ("Sst Chodl" if s.startswith("Sst Chodl") else "Sst")),
                              pct=float(frac.loc[s,g]),scaled_mean=float(_sm.loc[s,g]),
                              ct_idx=si,gene_idx=_gi,depth=round(float(meta.med_depth[s]),3)))
        _gi+=1
pd.DataFrame(_rows).to_csv(os.path.join(OUT,"sst_dotplot_long.csv"),index=False)
print(f"exported sst_dotplot_long.csv ({len(_rows)} rows, {_gi} genes, {len(order)} supertypes)",flush=True)

ngenes=sum(len(v) for v in markers_full.values())
sc.pl.dotplot(adata, markers_full, groupby="st_lab", categories_order=[lab[s] for s in order],
              standard_scale="var", cmap="Reds", show=False,
              figsize=(0.30*ngenes+3, 0.34*len(order)+1.6),
              dot_max=0.9, colorbar_title="scaled mean\nexpression", size_title="% expressing")
save_styled(plt.gcf(), os.path.join(OUT,"sst_supertype_depth_dotplot"))
