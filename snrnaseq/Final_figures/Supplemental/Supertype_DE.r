library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(grid)

BASE <- 13
AXIS_TITLE <- 12
AXIS_TEXT <- 11
LEGEND_TEXT <- 11
TEXT_SIZE <- 3.5
PANEL_LABEL <- 16

# PANEL A 
all_meta_counts <- read.csv("/scratch/nendresz/P1_Compositional_analysis/Files/7_cohorts_metadata_names.csv",row.names=1,check.names=FALSE)
de_all_supertype <- read.csv("/scratch/nendresz/P1_SCZ_DE_fresh/Files/DE_all_cohorts_meta_supertype.csv")

meta_cols <- c("Cohort","Donor","Age","Sex","Diagnosis","PMI")
cell_cols <- setdiff(colnames(all_meta_counts),meta_cols)
EXC <- c("L2_3 IT","L4 IT","L5 IT","L5 ET","L5_6 NP","L6 CT","L6 IT","L6 IT Car3","L6b")
INH <- c("Lamp5","Lamp5-Lhx6","Pax6","Pvalb","Sncg","Sst","Sst Chodl","Vip","Chandelier")
GLI <- c("Astro","Oligo","OPC","Micro-PVM","Endo","VLMC")
CLASS_COL <- c(Excitatory="#117733",Inhibitory="#882255",Glia="#DDCC77")
class_of <- function(ct) case_when(ct%in%EXC~"Excitatory",ct%in%INH~"Inhibitory",ct%in%GLI~"Glia",TRUE~NA_character_)

all_meta_counts$total_cells <- rowSums(all_meta_counts[,cell_cols],na.rm=TRUE)

abundance <- all_meta_counts |>
  select(all_of(cell_cols),total_cells) |>
  mutate(across(all_of(cell_cols),~.x/total_cells)) |>
  summarise(across(all_of(cell_cols),~mean(.x,na.rm=TRUE))) |>
  pivot_longer(everything(),names_to="cell_type",values_to="mean_prop") |>
  filter(cell_type!="total_cells") |>
  mutate(cell_type=gsub("_([0-9]+)$","-\\1",cell_type))

de_counts <- de_all_supertype |>
  filter(!is.na(FDR),FDR<.1) |>
  group_by(cell_type) |>
  summarise(n_DE=n_distinct(gene),.groups="drop")

plot_df <- abundance |>
  left_join(de_counts,by="cell_type") |>
  mutate(n_DE=coalesce(n_DE,0L),subclass=sub("-[0-9]+$","",cell_type),class=class_of(subclass))

plot_df_nonzero <- plot_df |> filter(n_DE>0,!is.na(class))
ct <- cor.test(plot_df_nonzero$mean_prop,plot_df_nonzero$n_DE,method="spearman",exact=FALSE)
rho <- unname(ct$estimate)

lab <- bind_rows(
  plot_df_nonzero |> slice_max(n_DE,n=3),
  plot_df_nonzero |> slice_max(mean_prop,n=3),
  plot_df_nonzero |> slice_min(mean_prop,n=3)
) |> distinct(cell_type,.keep_all=TRUE)

p_a <- ggplot(plot_df_nonzero,aes(mean_prop,n_DE)) +
  geom_smooth(method="lm",se=FALSE,colour=scales::alpha("grey25",.3),linewidth=.5,formula=y~x) +
  geom_point(aes(colour=class),size=.9,alpha=.9) +
  geom_text_repel(data=lab,aes(label=cell_type),size=TEXT_SIZE,min.segment.length=0,segment.size=.25,
                  segment.colour="grey55",box.padding=.5,point.padding=.4,force=5,force_pull=.5,
                  max.overlaps=Inf,max.time=2,seed=3,colour="grey15") +
  annotate("text",x=Inf,y=-Inf,label=sprintf("rho == %.2f",rho),parse=TRUE,
           hjust=1.05,vjust=-.5,size=TEXT_SIZE,colour="grey25") +
  scale_colour_manual(values=CLASS_COL,name=NULL) +
  scale_x_log10(breaks=c(.001,.01,.1),labels=c("0.1%","1%","10%")) +
  labs(x="Cell proportion (%)",y="DE genes (#)") +
  theme_cowplot(font_size=BASE) +
  theme(text=element_text(size=BASE),axis.title=element_text(size=AXIS_TITLE),
        axis.text=element_text(size=AXIS_TEXT),legend.text=element_text(size=LEGEND_TEXT),
        legend.position="right",legend.key.size=unit(.3,"cm"),
        axis.line=element_line(linewidth=.25),axis.ticks=element_line(linewidth=.25),
        plot.background=element_blank(),plot.margin=margin(5.5,5.5,5.5,5.5))

# PANEL B 
sup <- read.csv("/scratch/nendresz/P1_SCZ_DE_fresh/Files/DE_all_cohorts_meta_supertype.csv")
meta <- read.csv("/scratch/nendresz/P1_SCZ_DE_fresh/Files/DE_genes_all_cells_scz.csv")

sst_order <- c("Sst_23","Sst_25","Sst_11","Sst_22","Sst_20","Sst_2","Sst_3","Sst_19",
               "Sst_13","Sst_10","Sst_9","Sst_5","Sst_4","Sst_12","Sst_1","Sst_7")

sst_super <- sup %>%
  filter(gene=="SST",cohort=="Meta-analysis",grepl("^Sst[-_]",cell_type)) %>%
  transmute(label=gsub("-","_",cell_type),estimate=logFC,pval=PValue,padj=FDR,
            se=abs(logFC/qnorm(1-PValue/2)),source="Supertype") %>%
  mutate(ci.lb=estimate-1.96*se,ci.ub=estimate+1.96*se,
         sig=case_when(padj<.01~"***",padj<.05~"**",padj<.1~"*",padj<.2~"+",pval<.05~"•",TRUE~""))

sst_sub <- meta %>%
  filter(cell_type=="Sst",genes=="SST") %>%
  transmute(label="Sst (subclass)",estimate,se,pval,padj,ci.lb,ci.ub,source="Subclass",
            sig=case_when(padj<.01~"***",padj<.05~"**",padj<.1~"*",padj<.2~"+",pval<.05~"•",TRUE~""))

sst_forest <- bind_rows(sst_super,sst_sub) %>%
  mutate(label=factor(label,levels=c("Sst (subclass)",rev(sst_order))))

p_b <- ggplot(sst_forest,aes(estimate,label)) +
  geom_vline(xintercept=0,lty=2,color="grey60") +
  geom_errorbarh(aes(xmin=ci.lb,xmax=ci.ub),height=.18,linewidth=.6) +
  geom_point(aes(shape=source),size=3) +
  geom_text(aes(x=ci.lb-.03,label=sig),hjust=1,size=TEXT_SIZE,fontface="bold") +
  scale_shape_manual(values=c(Supertype=16,Subclass=18)) +
  labs(x=expression("SCZ log"[2]*"FC (95% CI)"),y=NULL,shape=NULL) +
  theme_classic(base_size=BASE) +
  theme(text=element_text(size=BASE),axis.title=element_text(size=AXIS_TITLE),
        axis.text=element_text(size=AXIS_TEXT),legend.position="none",
        plot.margin=margin(5.5,20,5.5,5.5))

# COMBINE 
p_combined <- plot_grid(p_a,p_b,labels=c("a","b"),label_size=PANEL_LABEL,label_fontface="bold",
                        nrow=1,rel_widths=c(1.3,1),align="h",axis="tb")

p_combined


ggsave("/scratch/nendresz/FINAL_FIGS/Paper/DE_abundance_and_SST_forest.png",
       p_combined,width=10,height=5,dpi=600,bg="white")

ggsave("/scratch/nendresz/FINAL_FIGS/Paper/DE_abundance_and_SST_forest.svg",
       p_combined,width=10,height=5,dpi=600,bg="white")

       ggsave("FINAL_FIGS/SST_forest_supertypes_subclass.png",p_b,width=6,height=5,dpi=600)

       ggsave("/scratch/nendresz/FINAL_FIGS/Paper/supertype_DE_vs_abundance.png",p_b,width=4.5,height=3.5,dpi=600,bg="white")
