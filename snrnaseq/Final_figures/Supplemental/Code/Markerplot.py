#mamba activate brisc
from brisc import SingleCell,concat_obs
import polars as pl,numpy as np,matplotlib.pyplot as plt,os
from pathlib import Path
from scipy import sparse
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.gridspec import GridSpec,GridSpecFromSubplotSpec

os.chdir("P1_SCZ_paper/Final_figures/Supplemental"); OUT=Path("."); N=3; HI="Sst_25"; HI_COL="#B8860B"

# -------------------- LOAD --------------------
HBCC1=SingleCell("/scratch/nendresz/PsychAD/Data/HBCC_1_symbols.rds").with_columns_obs(pl.col("Age").cast(pl.Float64,strict=False)).filter_obs((pl.col("Age")>20)&(pl.col("Age")<70)).qc(allow_float=True,max_mito_fraction=None)
HBCC2=SingleCell("/scratch/nendresz/PsychAD/Data/HBCC_2_symbols.rds").with_columns_obs(pl.col("Age").cast(pl.Float64,strict=False)).filter_obs((pl.col("Age")>20)&(pl.col("Age")<70)).qc(allow_float=True,max_mito_fraction=None); HBCC=concat_obs([HBCC1,HBCC2],flexible=True)

MSSM1=SingleCell("/scratch/nendresz/PsychAD/Data/MSSM_1_symbols.rds").with_columns_obs(pl.col("Age").cast(pl.Float64,strict=False)).filter_obs(pl.col("Age")<70).qc(allow_float=True,max_mito_fraction=None)
MSSM2=SingleCell("/scratch/nendresz/PsychAD/Data/MSSM_2_symbols.rds").with_columns_obs(pl.col("Age").cast(pl.Float64,strict=False)).filter_obs(pl.col("Age")<70).qc(allow_float=True,max_mito_fraction=None)
MSSM3=SingleCell("/scratch/nendresz/PsychAD/Data/MSSM_3_symbols.rds").with_columns_obs(pl.col("Age").cast(pl.Float64,strict=False)).filter_obs(pl.col("Age")<70).qc(allow_float=True,max_mito_fraction=None)
MSSM4=SingleCell("/scratch/nendresz/PsychAD/Data/MSSM_4_symbols.rds").with_columns_obs(pl.col("Age").cast(pl.Float64,strict=False)).filter_obs(pl.col("Age")<70).qc(allow_float=True,max_mito_fraction=None); MSSM=concat_obs([MSSM1,MSSM2,MSSM3,MSSM4],flexible=True)

OFC=SingleCell("/project/rrg-shreejoy/nendresz/Supertypes/OFC_updated.rds").filter_obs(pl.col("Age")<70).qc(allow_float=True)
Bat=SingleCell("/project/rrg-shreejoy/nendresz/Supertypes/Bat_updated.rds").filter_obs(pl.col("Age")<70).qc(allow_float=True)
Ruz=SingleCell("/project/rrg-shreejoy/nendresz/Supertypes/Ruz_updated.rds").filter_obs(pl.col("Age")<70).qc(allow_float=True)
Ruz_McLean=Ruz.filter_obs(pl.col("Cohort")=="McLean"); Ruz_MtSinai=Ruz.filter_obs(pl.col("Cohort")=="MtSinai")
Multi=SingleCell("/scratch/nendresz/P1_Brain_scope/Files/Multiome.rds").with_columns_obs(pl.col("Age_death").cast(pl.Float64,strict=False),pl.col("Disorder").alias("Diagnosis")).filter_obs((pl.col("Age_death")<70)&pl.col("Disorder").is_in(["control","Schizophrenia"])).qc(allow_float=True,max_mito_fraction=None)

MSSM=MSSM.with_columns_obs(pl.lit("MSSM 2").alias("Cohort")); HBCC=HBCC.with_columns_obs(pl.lit("HBCC").alias("Cohort")); OFC=OFC.with_columns_obs(pl.lit("Fröhlich").alias("Cohort")); Bat=Bat.with_columns_obs(pl.lit("Batiuk").alias("Cohort")); Ruz_McLean=Ruz_McLean.with_columns_obs(pl.lit("McLean").alias("Cohort")); Ruz_MtSinai=Ruz_MtSinai.with_columns_obs(pl.lit("MSSM 1").alias("Cohort")); Multi=Multi.with_columns_obs(pl.lit("Multiome").alias("Cohort"))
scs=[MSSM,HBCC,OFC,Bat,Ruz_McLean,Ruz_MtSinai,Multi]; labels=["MSSM 2","HBCC","Fröhlich","Batiuk","McLean","MSSM 1","Multiome"]

# -------------------- PREP COHORTS + FIND COHORT MARKERS --------------------
scs=[sc.filter_obs(pl.col("Donor").is_in(sc.obs.group_by("Donor").len().filter(pl.col("len")>=500)["Donor"])) for sc in scs]
scs=[sc.with_columns_obs(pl.when(pl.col("Diagnosis").cast(pl.String).str.to_lowercase()=="control").then(pl.lit("Control")).when(pl.col("Diagnosis").cast(pl.String).str.to_lowercase()=="schizophrenia").then(pl.lit("Schizophrenia")).otherwise(pl.col("Diagnosis").cast(pl.String)).alias("Diagnosis")) for sc in scs]
scs=[sc.with_columns_obs(pl.col("predicted.id").cast(pl.String).str.replace(r"_[0-9].*$","").alias("Subclass")) for sc in scs]
scs=[sc.with_columns_obs(pl.col("Subclass").replace({"Astro":"Astrocyte","Endo":"Endothelial","Lamp5_Lhx6":"Lamp5 Lhx6","Micro-PVM":"Microglia-PVM","Oligo":"Oligodendrocyte"}).alias("Subclass")) for sc in scs]
scs=[sc.normalize() for sc in scs]; combined=concat_obs(scs,dataset_column="batch",dataset_labels=labels,flexible=True)
markers=combined.find_markers("Subclass"); markers.write_csv(OUT/"Files/cohort_subclass_markers.csv")

# -------------------- REFERENCE --------------------
ref=SingleCell("/scratch/nendresz/1_SCZ_project/1_Reference/Data/MTG_ref.rds").qc(allow_float=True,max_mito_fraction=None)
ref=ref.with_columns_obs(pl.col("Supertype").cast(pl.String).str.replace(r"_[0-9].*$","").alias("Subclass"))
ref=ref.with_columns_obs(pl.col("Subclass").replace({"Astro":"Astrocyte","Endo":"Endothelial","Lamp5_Lhx6":"Lamp5 Lhx6","Micro-PVM":"Microglia-PVM","Oligo":"Oligodendrocyte"}).alias("Subclass")).normalize()

# -------------------- SST + FIND COHORT SST MARKERS --------------------
scs_sst=[sc.filter_obs(pl.col("Subclass")=="Sst").with_columns_obs(pl.col("predicted.id").cast(pl.String).alias("Sst_supertype")) for sc in scs]
combined_sst=concat_obs(scs_sst,dataset_column="batch",dataset_labels=labels,flexible=True)
sst_markers=combined_sst.find_markers("Sst_supertype"); sst_markers.write_csv(OUT/"Files/cohort_sst_supertype_markers.csv")
ref_sst=ref.filter_obs(pl.col("Subclass")=="Sst").with_columns_obs(pl.col("Supertype").cast(pl.String).alias("Sst_supertype"))

# -------------------- SETTINGS --------------------
GABA=["GAD1","GAD2","SLC32A1"]; GABA_SST=["SST","CALB1","GAD1","GAD2","SLC32A1"]
SUB_ORDER=["Lamp5","Sncg","Vip","Pax6","Sst","Sst Chodl","Pvalb","Chandelier","L2/3 IT","L4 IT","L5 IT","L5 ET","L5/6 NP","L6 IT","L6 IT Car3","L6 CT","L6b","Astrocyte","Oligodendrocyte","OPC","Microglia-PVM","Endothelial","VLMC"]
SST_ORDER=["Sst_23","Sst_25","Sst_11","Sst_22","Sst_20","Sst_2","Sst_3","Sst_19","Sst_13","Sst_10","Sst_9","Sst_5","Sst_4","Sst_12","Sst_1","Sst_7"]
INH=SUB_ORDER[:8]; EXC=SUB_ORDER[8:17]; NON=SUB_ORDER[17:]

# -------------------- TOP 3 COHORT MARKERS --------------------
def top_genes(m,types,n=N):
    x=(m.filter(pl.col("cell_type").is_in(types)&pl.col("gene").is_not_null()&pl.col("fold_change").is_not_null()&(pl.col("fold_change")>0)).sort(["cell_type","fold_change"],descending=[False,True]).group_by("cell_type",maintain_order=True).head(n))
    out=[]
    for ct in types:
        for g in x.filter(pl.col("cell_type")==ct)["gene"].to_list():
            if g not in out: out.append(g)
    return out

def subclass_sections():
    out={"GABA":GABA}; used=set(GABA)
    for sec,types in [("Inhibitory",INH),("Excitatory",EXC),("Non-neuronal",NON)]:
        g=[x for x in top_genes(markers,types) if x not in used]; out[sec]=g; used.update(g)
    return out

def sst_sections(): return {"GABA/SST":GABA_SST,"Sst":[x for x in top_genes(sst_markers,SST_ORDER) if x not in GABA_SST]}

SUB_SECTIONS=subclass_sections(); SST_SECTIONS=sst_sections()
print("\nCOHORT SUBCLASS MARKERS:",SUB_SECTIONS); print("\nCOHORT SST MARKERS:",SST_SECTIONS)

# -------------------- EXPRESSION --------------------
def gene_names(sc):
    for a in ["var_names","gene_names"]:
        if hasattr(sc,a):
            x=getattr(sc,a); x=x() if callable(x) else x
            if x is not None and len(x)==sc.X.shape[1]: return list(map(str,x))
    for c in ["gene","Gene","symbol","gene_symbol","feature_name","_index","index"]:
        if c in sc.var.columns:
            x=sc.var[c].cast(pl.String).to_list()
            if len(x)==sc.X.shape[1]: return x
    if sc.var.width==1:return sc.var[sc.var.columns[0]].cast(pl.String).to_list()
    raise ValueError(f"Cannot locate gene names; var columns={sc.var.columns}")

def summarise(sc,group,sections,order):
    names=gene_names(sc); lookup={g:i for i,g in enumerate(names)}
    sections={s:[g for g in gs if g in lookup] for s,gs in sections.items()}; genes=list(dict.fromkeys(g for gs in sections.values() for g in gs))
    groups=np.asarray(sc.obs[group].cast(pl.String).to_list()); order=[ct for ct in order if np.any(groups==ct)]; X=sc.X[:,[lookup[g] for g in genes]]; rows=[]
    for ct in order:
        Xi=X[groups==ct]
        if sparse.issparse(Xi):
            mean=np.asarray(Xi.mean(axis=0)).ravel(); nnz=Xi.getnnz(axis=0) if hasattr(Xi,"getnnz") else Xi.count_nonzero(axis=0); pct=np.asarray(nnz).ravel()/Xi.shape[0]
        else: Xi=np.asarray(Xi); mean=Xi.mean(axis=0); pct=(Xi>0).mean(axis=0)
        rows.extend((ct,g,float(mu),float(p)) for g,mu,p in zip(genes,mean,pct))
    d=pl.DataFrame(rows,schema=["cell_type","gene","mean","pct"],orient="row").with_columns(pl.when(pl.col("mean").max().over("gene")>pl.col("mean").min().over("gene")).then((pl.col("mean")-pl.col("mean").min().over("gene"))/(pl.col("mean").max().over("gene")-pl.col("mean").min().over("gene"))).otherwise(0).alias("scaled_mean"))
    return d,sections,order

def add_strip(ax,label):
    ax.add_patch(Rectangle((0,1.01),1,.065,transform=ax.transAxes,clip_on=False,facecolor="#E6E6E6",edgecolor="none")); ax.text(.5,1.042,label,transform=ax.transAxes,ha="center",va="center",fontsize=10)

def panel(fig,spec,sc,group,sections,order,title="",is_sst=False,show_title=True,title_offset=.030):
    d,sections,order=summarise(sc,group,sections,order); gs=GridSpecFromSubplotSpec(1,len(sections),subplot_spec=spec,width_ratios=[max(len(g),1) for g in sections.values()],wspace=.025); axes=[]
    for i,(section,genes) in enumerate(sections.items()):
        ax=fig.add_subplot(gs[0,i],sharey=axes[0] if axes else None); axes.append(ax); z=d.filter(pl.col("gene").is_in(genes)).to_pandas(); xm={g:j for j,g in enumerate(genes)}; ym={ct:j for j,ct in enumerate(order)}
        ax.scatter(z["gene"].map(xm),z["cell_type"].map(ym),s=2+52*z["pct"],c=z["scaled_mean"],cmap="Reds",vmin=0,vmax=1,edgecolors="none")
        ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes,rotation=90,ha="center",fontsize=10); ax.set_yticks(range(len(order)))
        if i==0:
            ax.set_yticklabels(order,fontsize=11)
            if is_sst:
                for lab in ax.get_yticklabels():
                    if lab.get_text()==HI: lab.set_color(HI_COL); lab.set_fontweight("bold")
        else: ax.tick_params(axis="y",left=False,labelleft=False)
        ax.set_xlim(-.55,len(genes)-.45); ax.set_ylim(len(order)-.5,-.5); ax.tick_params(axis="x",length=2,width=.4,pad=1); ax.tick_params(axis="y",length=2,width=.4,pad=2)
        for sp in ax.spines.values(): sp.set_linewidth(.45)
        add_strip(ax,section)
    p0,p1=axes[0].get_position(),axes[-1].get_position()
    if show_title: fig.text((p0.x0+p1.x1)/2,p0.y1+title_offset,title,ha="center",va="bottom",fontsize=12)
    return axes,d

# -------------------- FIGURE --------------------
fig=plt.figure(figsize=(22,14)); outer=GridSpec(2,2,figure=fig,left=.065,right=.985,top=.93,bottom=.16,wspace=.15,hspace=.31,height_ratios=[1.10,.90])
a,d1=panel(fig,outer[0,0],combined,"Subclass",SUB_SECTIONS,SUB_ORDER,"Seven-cohort marker expression")
b,d2=panel(fig,outer[0,1],ref,"Subclass",SUB_SECTIONS,SUB_ORDER,"SEA-AD reference marker expression")
c,d3=panel(fig,outer[1,0],combined_sst,"Sst_supertype",SST_SECTIONS,SST_ORDER,"",True,False,.025)
d,d4=panel(fig,outer[1,1],ref_sst,"Sst_supertype",SST_SECTIONS,SST_ORDER,"",True,False,.025)

for axes,label in zip([a,b,c,d],"abcd"): axes[0].text(-.10,1.11,label,transform=axes[0].transAxes,fontsize=15,fontweight="bold",ha="right",va="bottom")

# -------------------- LEGENDS + SAVE --------------------
size_handles=[Line2D([],[],marker="o",linestyle="",markerfacecolor="black",markeredgecolor="none",markersize=np.sqrt(2+52*p)*1.25) for p in [.1,.25,.5,.9]]
fig.legend(size_handles,["10%","25%","50%","90%"],title="% expressing",loc="lower left",bbox_to_anchor=(.09,.050),ncol=4,frameon=False,handletextpad=.7,columnspacing=2,handlelength=1.6,fontsize=10,title_fontsize=11)
cax=fig.add_axes([.53,.050,.26,.014]); cb=fig.colorbar(ScalarMappable(norm=Normalize(0,1),cmap="Reds"),cax=cax,orientation="horizontal"); cb.set_ticks([0,.5,1]); cb.ax.tick_params(labelsize=9,length=2,pad=2); cb.outline.set_linewidth(.5); fig.text(.66,.073,"Scaled mean expression",ha="center",va="bottom",fontsize=11)

for ext in ["png"]:
    kw={"dpi":600} if ext=="png" else {}; fig.savefig(OUT/f"Figures/FigureS1_four_compiled_cohort_markers.{ext}",bbox_inches="tight",facecolor="white",**kw)
d1.write_csv(OUT/"Files/cohort_subclass_dotplot.csv"); d2.write_csv(OUT/"Files/reference_cohort_subclass_markers_dotplot.csv"); d3.write_csv(OUT/"Files/cohort_sst_dotplot.csv"); d4.write_csv(OUT/"Files/reference_cohort_sst_markers_dotplot.csv")
