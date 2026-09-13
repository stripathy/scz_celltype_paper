# Figure 1
# mamba activate brisc

from brisc import SingleCell,concat_obs
from matplotlib.lines import Line2D
import polars as pl
import pandas as pd
import matplotlib.pyplot as plt
import os

os.chdir("scz_celltype_paper/snrnaseq/Final_figures")

# -------------------- LOAD DATA --------------------

HBCC1=SingleCell("/scratch/nendresz/PsychAD/Data/HBCC_1_symbols.rds").with_columns_obs(pl.col("Age").cast(pl.Float64,strict=False)).filter_obs((pl.col("Age")>20)&(pl.col("Age")<70)).qc(allow_float=True,max_mito_fraction=None)
HBCC2=SingleCell("/scratch/nendresz/PsychAD/Data/HBCC_2_symbols.rds").with_columns_obs(pl.col("Age").cast(pl.Float64,strict=False)).filter_obs((pl.col("Age")>20)&(pl.col("Age")<70)).qc(allow_float=True,max_mito_fraction=None)
HBCC=concat_obs([HBCC1,HBCC2],flexible=True)

MSSM1=SingleCell("/scratch/nendresz/PsychAD/Data/MSSM_1_symbols.rds").with_columns_obs(pl.col("Age").cast(pl.Float64,strict=False)).filter_obs(pl.col("Age")<70).qc(allow_float=True,max_mito_fraction=None)
MSSM2=SingleCell("/scratch/nendresz/PsychAD/Data/MSSM_2_symbols.rds").with_columns_obs(pl.col("Age").cast(pl.Float64,strict=False)).filter_obs(pl.col("Age")<70).qc(allow_float=True,max_mito_fraction=None)
MSSM3=SingleCell("/scratch/nendresz/PsychAD/Data/MSSM_3_symbols.rds").with_columns_obs(pl.col("Age").cast(pl.Float64,strict=False)).filter_obs(pl.col("Age")<70).qc(allow_float=True,max_mito_fraction=None)
MSSM4=SingleCell("/scratch/nendresz/PsychAD/Data/MSSM_4_symbols.rds").with_columns_obs(pl.col("Age").cast(pl.Float64,strict=False)).filter_obs(pl.col("Age")<70).qc(allow_float=True,max_mito_fraction=None)
MSSM=concat_obs([MSSM1,MSSM2,MSSM3,MSSM4],flexible=True)

OFC=SingleCell("/project/rrg-shreejoy/nendresz/Supertypes/OFC_updated.rds").filter_obs(pl.col("Age")<70).qc(allow_float=True)
Bat=SingleCell("/project/rrg-shreejoy/nendresz/Supertypes/Bat_updated.rds").filter_obs(pl.col("Age")<70).qc(allow_float=True)

Ruz=SingleCell("/project/rrg-shreejoy/nendresz/Supertypes/Ruz_updated.rds").filter_obs(pl.col("Age")<70).qc(allow_float=True)
Ruz_McLean=Ruz.filter_obs(pl.col("Cohort")=="McLean")
Ruz_MtSinai=Ruz.filter_obs(pl.col("Cohort")=="MtSinai")

Multi=SingleCell("/scratch/nendresz/P1_Brain_scope/Files/Multiome.rds").with_columns_obs(
    pl.col("Age_death").cast(pl.Float64,strict=False),
    pl.col("Disorder").alias("Diagnosis")
).filter_obs(
    (pl.col("Age_death")<70)&pl.col("Disorder").is_in(["control","Schizophrenia"])
).qc(allow_float=True,max_mito_fraction=None)

# -------------------- COHORT LABELS --------------------

MSSM=MSSM.with_columns_obs(pl.lit("MSSM 2").alias("Cohort"))
HBCC=HBCC.with_columns_obs(pl.lit("HBCC").alias("Cohort"))
OFC=OFC.with_columns_obs(pl.lit("Fröhlich").alias("Cohort"))
Bat=Bat.with_columns_obs(pl.lit("Batiuk").alias("Cohort"))
Ruz_McLean=Ruz_McLean.with_columns_obs(pl.lit("McLean").alias("Cohort"))
Ruz_MtSinai=Ruz_MtSinai.with_columns_obs(pl.lit("MSSM 1").alias("Cohort"))
Multi=Multi.with_columns_obs(pl.lit("Multiome").alias("Cohort"))

scs=[MSSM,HBCC,OFC,Bat,Ruz_McLean,Ruz_MtSinai,Multi]
labels=["MSSM 2","HBCC","Fröhlich","Batiuk","McLean","MSSM 1","Multiome"]

# Keep donors with at least 500 cells
scs=[sc.filter_obs(pl.col("Donor").is_in(sc.obs.group_by("Donor").len().filter(pl.col("len")>=500)["Donor"])) for sc in scs]

# Standardize diagnosis
scs=[sc.with_columns_obs(
    pl.when(pl.col("Diagnosis").cast(pl.String).str.to_lowercase()=="control").then(pl.lit("Control"))
      .when(pl.col("Diagnosis").cast(pl.String).str.to_lowercase()=="schizophrenia").then(pl.lit("Schizophrenia"))
      .otherwise(pl.col("Diagnosis").cast(pl.String)).alias("Diagnosis")
) for sc in scs]

# Make subclass
scs=[sc.with_columns_obs(pl.col("predicted.id").cast(pl.String).str.replace(r"_[0-9].*$","").alias("Subclass")) for sc in scs]
scs=[sc.with_columns_obs(pl.col("Subclass").replace({
    "Astro":"Astrocyte",
    "Endo":"Endothelial",
    "Lamp5_Lhx6":"Lamp5 Lhx6",
    "Micro-PVM":"Microglia-PVM",
    "Oligo":"Oligodendrocyte"
}).alias("Subclass")) for sc in scs]

# -------------------- FULL DATA UMAP --------------------

scs_all=list(scs[0].hvg(*scs[1:],batch_column="Donor"))
scs_all=[sc.normalize() for sc in scs_all]
scs_all=list(scs_all[0].pca(*scs_all[1:]))
scs_all=list(scs_all[0].harmonize(*scs_all[1:],theta=6,max_iterations=20, overwrite=True))

combined=concat_obs(scs_all,dataset_column="batch",dataset_labels=labels,flexible=True)
combined=combined.neighbors(PC_key="harmony").shared_neighbors()
combined=combined.umap(PC_key="harmony",hogwild=True)

# -------------------- SST-ONLY OBJECT AND UMAP --------------------

sst_scs=[sc.filter_obs(pl.col("Subclass")=="Sst") for sc in scs]
sst_scs=list(sst_scs[0].hvg(*sst_scs[1:],batch_column=None,span=0.8,overwrite=True))
sst_scs=[sc.normalize() for sc in sst_scs]
sst_scs=list(sst_scs[0].pca(*sst_scs[1:]))

combined_sst=concat_obs(sst_scs,dataset_column="batch",dataset_labels=labels,flexible=True)
combined_sst=combined_sst.harmonize(batch_column="batch",theta=6,max_iterations=20,overwrite=True)
combined_sst=combined_sst.neighbors(PC_key="harmony",overwrite=True).shared_neighbors(overwrite=True).umap(PC_key="harmony",hogwild=True,overwrite=True)

# -------------------- COLOURS --------------------

colours=pd.read_csv("/scratch/nendresz/scz_celltype_paper/snrnaseq/Compositional_analysis/Files/cluster_order_and_colors.csv")

subclass_cols=colours.drop_duplicates("subclass_label").set_index("subclass_label")["subclass_color"].to_dict()
sst_cols=colours[colours["cluster_label"].str.startswith("Sst_")].set_index("cluster_label")["cluster_color"].to_dict()

cohort_cols={
    "Batiuk":"#1f93c6",
    "Fröhlich":"#ff4500",
    "HBCC":"#8ac926",
    "McLean":"#d100ff",
    "MSSM 1":"#11823b",
    "MSSM 2":"#2ec4b6",
    "Multiome":"#a87400"
}

diagnosis_cols={"Control":"#0a7ad0ff","Schizophrenia":"#d73027"}
layer_cols={"L1":"#D73027","L2/3":"#56B4E9","L4":"#009E73","L5":"#E69F00","L6":"#D55E00","WM":"#CC79A7","Vascular":"#7F7F7F"}

# Save separate SST UMAP
fig,ax=plt.subplots(figsize=(6,6))
combined_sst.plot_umap("predicted.id",ax=ax,colormap=sst_cols,first_color=None,stride=None,legend=False,label=True)
ax.axis("off")
plt.savefig("Figures/combined_sst_umap_predicted_id.png",dpi=300,bbox_inches="tight",pad_inches=.01)
plt.close()

# -------------------- SPATIAL DATA --------------------

meta=pd.read_csv("Data/xenium_metadata.csv")
meta=meta[meta.qc_pass.astype(str)=="True"]

def plot_spatial(ax,df,col,cmap,title="",mirror=False):
    df=df.copy()
    if mirror: df["x"]=df["x"].max()-df["x"]
    ax.scatter(df.x,df.y,c=df[col].map(cmap).fillna("lightgrey"),s=.35,lw=0,rasterized=True)
    ax.set_aspect("equal",adjustable="box"); ax.set_anchor("N"); ax.invert_yaxis(); ax.axis("off")
    if title: ax.set_title(title,fontweight="bold",fontsize=14,pad=1)


# -------------------- FINAL COMPOSITE --------------------
cohort_order=["MSSM 2","HBCC","Fröhlich","MSSM 1","McLean","Batiuk","Multiome"]
panel_size=20; umap_label_size=12; legend_text_size=12; legend_title_size=14; legend_marker_size=9

fig=plt.figure(figsize=(20,8))
outer=fig.add_gridspec(2,1,height_ratios=[2,1],hspace=-.015)

# UMAP row
gs_umap=outer[0].subgridspec(1,4,width_ratios=[1,.72,1.15,1.15],wspace=.005)
axs_umap=[fig.add_subplot(gs_umap[0,i]) for i in range(4)]

combined.plot_umap("Subclass",ax=axs_umap[0],colormap=subclass_cols,first_color=None,stride=None,legend=False,label=True,label_kwargs=dict(fontsize=umap_label_size,fontweight="normal"))
combined_sst.plot_umap("predicted.id",ax=axs_umap[1],colormap=sst_cols,first_color=None,stride=None,legend=False,label=True,label_kwargs=dict(fontsize=umap_label_size,fontweight="normal"))
combined.plot_umap("Cohort",ax=axs_umap[2],colormap=cohort_cols,first_color=None,stride=None,legend=False,label=False)
combined.plot_umap("Diagnosis",ax=axs_umap[3],colormap=diagnosis_cols,first_color=None,stride=None,legend=False,label=False)

axs_umap[2].legend(handles=[Line2D([0],[0],marker="o",ls="",mfc=cohort_cols[x],mec="none",ms=legend_marker_size,label=x) for x in cohort_order],title="Cohort",loc="upper left",bbox_to_anchor=(.76,1),frameon=True,fontsize=legend_text_size,title_fontsize=legend_title_size,labelspacing=.25,handletextpad=.35)
axs_umap[3].legend(handles=[Line2D([0],[0],marker="o",ls="",mfc=c,mec="none",ms=legend_marker_size,label=x) for x,c in diagnosis_cols.items()],title="Diagnosis",loc="upper left",bbox_to_anchor=(.72,1),frameon=True,fontsize=legend_text_size,title_fontsize=legend_title_size,labelspacing=.25,handletextpad=.35)

for ax in axs_umap: ax.axis("off"); ax.set_anchor("N")

# Spatial row — e has no legend; f retains layer legend
gs_sp=outer[1].subgridspec(1,3,width_ratios=[2.5,2,.65],wspace=.02)
gs_sub=gs_sp[0,0].subgridspec(1,2,wspace=-.05); axs_sub=[fig.add_subplot(gs_sub[0,i]) for i in range(2)]
gs_layer=gs_sp[0,1].subgridspec(1,2,wspace=-.05); axs_layer=[fig.add_subplot(gs_layer[0,i]) for i in range(2)]
ax_layerleg=fig.add_subplot(gs_sp[0,2]); ax_layerleg.axis("off")

for ax,sid,mirror,title in zip(axs_sub,["Br5931","Br1139"],[False,True],["Control","Schizophrenia"]): plot_spatial(ax,meta[meta.sample_id==sid],"subclass",subclass_cols,title,mirror)
for ax,sid,mirror,title in zip(axs_layer,["Br5931","Br1139"],[False,True],["Control","Schizophrenia"]): plot_spatial(ax,meta[meta.sample_id==sid],"layer",layer_cols,title,mirror)

axs_sub[0].set_anchor("E"); axs_sub[1].set_anchor("W"); axs_layer[0].set_anchor("E"); axs_layer[1].set_anchor("W")

lay_handles=[Line2D([0],[0],marker="o",ls="",mfc=c,mec="none",ms=legend_marker_size,label=k) for k,c in layer_cols.items() if k in meta.layer.values]
ax_layerleg.legend(handles=lay_handles,title="Layer",loc="center left",bbox_to_anchor=(0,.5),frameon=False,ncol=1,fontsize=legend_text_size,title_fontsize=legend_title_size,labelspacing=.25,handletextpad=.3)

# Panel labels
axs_umap[0].text(.03,.98,"b",transform=axs_umap[0].transAxes,ha="left",va="top",fontsize=panel_size,fontweight="bold")
axs_umap[2].text(.03,.98,"c",transform=axs_umap[2].transAxes,ha="left",va="top",fontsize=panel_size,fontweight="bold")
axs_umap[3].text(.03,.98,"d",transform=axs_umap[3].transAxes,ha="left",va="top",fontsize=panel_size,fontweight="bold")
axs_sub[0].text(-.03,1.02,"e",transform=axs_sub[0].transAxes,ha="left",va="top",fontsize=panel_size,fontweight="bold")
axs_layer[0].text(-.03,1.02,"f",transform=axs_layer[0].transAxes,ha="left",va="top",fontsize=panel_size,fontweight="bold")

# Save
fig.subplots_adjust(left=.008,right=.992,top=.998,bottom=.008)
plt.savefig("Figures/Figure1b_f.png",dpi=300,bbox_inches="tight",pad_inches=.01)
