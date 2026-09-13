# Figure 3A) Cell composition barchart
library(ggplot2)
library(dplyr)
library(ggtext) 
setwd("scz_celltype_paper/snrnaseq")

BASE <- 23

AXIS_TITLE <- BASE - 0.5
AXIS_TEXT  <- BASE - 1
GEOM_TEXT  <- BASE * 0.30

final_results <- read.csv("Compositional_analysis/Files/Meta_Neurons.csv")

colours <- read.csv("Compositional_analysis/Files/cluster_order_and_colors.csv")

final_results$FDR <- final_results$padj

plot_df <- final_results %>%
  left_join(colours, by = c("CellType" = "cluster_label")) %>%
  mutate(signif_label = case_when(FDR < 0.01 ~ "***", FDR < 0.05 ~ "**", FDR < 0.1 ~ "*", FDR >= 0.1 & FDR < 0.2 ~ "+", TRUE ~ ""),
         CellType = factor(CellType, levels = colours$cluster_label))

celltype_labels <- setNames(case_when(
  plot_df$FDR < 0.1 ~ paste0("<span style='color:red;font-size:16pt'><b>", plot_df$CellType, "</b></span>"),
  plot_df$FDR < 0.2 ~ paste0("<span style='color:red'><i>", plot_df$CellType, "</i></span>"),
  TRUE ~ paste0("<span style='color:black'>", plot_df$CellType, "</span>")
), as.character(plot_df$CellType))

p3a <- ggplot(plot_df, aes(x = CellType, y = estimate, fill = cluster_color)) +
  geom_col(color = "black", width = 0.8) +
  geom_errorbar(aes(ymin = estimate - se, ymax = estimate + se), width = 0.5) +
  geom_text(aes(label = signif_label, y = estimate + sign(estimate) * (se + 0.02)), fontface = "bold",
            vjust = ifelse(plot_df$estimate >= 0, 0, 1), size = GEOM_TEXT, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_identity() +
  scale_x_discrete(labels = celltype_labels) +
  annotate("text", x = -Inf, y = 0.6, label = "Increased abundance in SCZ", hjust = -0.1, color = "black", size = GEOM_TEXT) +
  annotate("text", x = -Inf, y = -0.4, label = "Decreased abundance in SCZ", hjust = -0.1, color = "black", size = GEOM_TEXT) +
  theme_classic(base_size = BASE) +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = AXIS_TEXT), axis.title.y = element_text(size = AXIS_TITLE),
         axis.ticks.y = element_blank(), panel.grid = element_blank(),
        axis.title.x = element_blank(), plot.caption = element_blank()) +
  labs(y = "SCZ abundance change (β ± SE)")

ggsave("Final_figures/Figures/Barchart_SCZ_meta_FINAL.png", p3a, width = 21, height = 8,  dpi = 600)


# Figure 3B
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggsignif)

merged_all <- read.csv("Final_figures/Data/neuron_props_7_cohorts.csv")
merged_all$cohort <- recode(merged_all$cohort, "Frölich" = "Fröhlich", "MtSinai" = "MSSM 1", "MSSM" = "MSSM 2")

res <- read.csv("Compositional_analysis/Files/crumblr_results_neurons.csv")
res$Cohort <- recode(res$Cohort, "Bat" = "Batiuk", "OFC" = "Fröhlich", "Multi" = "Multiome",
                     "Ruz_McLean" = "McLean", "Ruz_MtSinai" = "MSSM 1", "MSSM" = "MSSM 2")
res_25 <- res %>% filter(CellType == "Sst_25")

df_xen <- read.csv("Final_figures/Data/xen_Sst_proportions.csv")
res_xen <- read.csv("Final_figures/Data/xenium_crumblr_results_supertype_neuronal.csv") %>% filter(celltype == "Sst_25")
res_xen$CellType <- res_xen$celltype

plot_df <- merged_all %>%
  group_by(cohort, Donor, Diagnosis) %>%
  summarise(Sst_25 = mean(Sst_25, na.rm = TRUE), .groups = "drop") %>%
  mutate(Cohort = cohort, Diagnosis = recode(Diagnosis, "Control" = "CON", "Schizophrenia" = "SCZ")) %>%
  select(Cohort, Donor, Diagnosis, Sst_25)

xen_df <- df_xen %>%
  filter(subtype == "Sst_25") %>%
  mutate(Cohort = "Xenium", Donor = paste0("Xenium_", row_number()),
         Diagnosis = recode(diagnosis, "Control" = "CON", "Schizophrenia" = "SCZ"),
         Sst_25 = proportion) %>%
  select(Cohort, Donor, Diagnosis, Sst_25)

plot_df <- bind_rows(plot_df, xen_df) %>% filter(Sst_25 > 0)

annot_df <- bind_rows(
  res_25 %>% select(Cohort, P.Value),
  res_xen %>% transmute(Cohort = "Xenium", P.Value)
) %>%
  mutate(label = paste0("p = ", signif(P.Value, 2)), xmin = 1, xmax = 2, y_position = 5)

cohort_order <- c("Batiuk", "MSSM 1", "Fröhlich", "McLean", "HBCC", "Multiome", "MSSM 2", "Xenium")
plot_df$Cohort <- factor(plot_df$Cohort, levels = cohort_order)
annot_df$Cohort <- factor(annot_df$Cohort, levels = cohort_order)

p3b <- ggplot(plot_df, aes(x = Diagnosis, y = Sst_25 * 100)) +
  geom_boxplot(aes(fill = Diagnosis), alpha = 0.35, outlier.shape = NA, width = 0.55, linewidth = 0.4) +
  geom_jitter(aes(color = Diagnosis), width = 0.12, size = 1, alpha = 0.8) +
  geom_signif(data = annot_df, aes(xmin = xmin, xmax = xmax, annotations = label, y_position = y_position),
              manual = TRUE, inherit.aes = FALSE, textsize = GEOM_TEXT, tip_length = 0.01) +
  scale_y_continuous(limits = c(0, 5.5), expand = expansion(mult = c(0, 0))) +
  ylab("Sst_25 proportion\n(% of all neurons)") +
  facet_wrap(~Cohort, ncol = 8) +
  scale_fill_manual(values = c("CON" = "#0a7ad0ff", "SCZ" = "#6e0505ff"), guide = "none") +
  scale_color_manual(values = c("CON" = "#0a7ad0ff", "SCZ" = "#6e0505ff"), guide = "none") +
  theme_classic(base_size = BASE) +
  theme(panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank(), axis.text.x = element_text(size = AXIS_TEXT), axis.text.y = element_text(size = AXIS_TEXT),
        axis.title.x = element_blank(), axis.title.y = element_text(size = AXIS_TITLE),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        strip.text = element_text(size = AXIS_TEXT, face = "bold"),
        strip.background = element_rect(fill = "white", color = "black"))

ggsave("Final_figures/Figures/Boxplots3b_with_Xenium.png", p3b, width = 14, height = 8, dpi = 600)

# Figure 3C
library(dplyr)
library(ggplot2)
plot_data <- read.csv("Final_figures/Data/plotdata.csv")
res_xen <- read.csv("Final_figures/Data/xenium_crumblr_results_supertype_neuronal.csv") %>% filter(celltype == "Sst_25") 
res_xen$CellType <- res_xen$celltype
xen_n <- read.csv("Final_figures/Data/xen_Sst_proportions.csv") %>% filter(subtype=="Sst_25") %>% nrow()

df_sst25 <- plot_data %>%
  filter(CellType=="Sst_25") %>%
  mutate(Cohort=recode(Cohort,"MtSinai"="MSSM 1","MSSM"="MSSM 2")) %>%
  bind_rows(res_xen %>% transmute(CellType="Sst_25",Cohort="Xenium",estimate=logFC,se=abs(logFC/t),n=xen_n,padj=adj.P.Val)) %>%
  mutate(ci_low=estimate-1.96*se,ci_high=estimate+1.96*se,
         signif_label=case_when(padj<0.01~"***",padj<0.05~"**",padj<0.1~"*", pval <0.05~".",TRUE~""))

meta_n <- df_sst25 %>% filter(!Cohort %in% c("Meta-analysis","Xenium")) %>% summarise(n=sum(n,na.rm=TRUE)) %>% pull(n)
df_sst25$n[df_sst25$Cohort=="Meta-analysis"] <- meta_n

plot_order <- c("Xenium","Meta-analysis",
                df_sst25 %>% filter(!Cohort %in% c("Meta-analysis","Xenium")) %>% arrange(estimate) %>% pull(Cohort) %>% as.character())
df_sst25$Cohort <- factor(df_sst25$Cohort, levels=plot_order)

p3c <- ggplot(df_sst25,aes(Cohort,estimate)) +
  geom_hline(yintercept=0,linetype="dashed") +
  geom_vline(xintercept=1.5,linetype="dotted",color="grey85") +
geom_errorbar(data=filter(df_sst25,!Cohort %in% c("Meta-analysis","Xenium")),aes(ymin=ci_low,ymax=ci_high),width=0,linewidth=.8,color="grey35") +
geom_errorbar(data=filter(df_sst25,Cohort=="Meta-analysis"),aes(ymin=ci_low,ymax=ci_high),width=0,linewidth=.8,color="black") +
geom_errorbar(data=filter(df_sst25,Cohort=="Xenium"),aes(ymin=ci_low,ymax=ci_high),width=0,linewidth=.8,color="#1b9e77") +
  geom_point(data=filter(df_sst25,!Cohort %in% c("Meta-analysis","Xenium")),aes(size=n),shape=21,fill="grey35",color="grey35",stroke=.5) +
  geom_point(data=filter(df_sst25,Cohort=="Meta-analysis"),shape=23,size=8,fill="black",color="black") +
  geom_point(data=filter(df_sst25,Cohort=="Xenium"),aes(size=n),shape=24,fill="#1b9e77",color="#1b9e77",show.legend=FALSE) +
geom_point(data=filter(df_sst25,signif_label=="."),
           aes(y=ci_low-.08),shape=16,size=1.8) +
geom_text(data=filter(df_sst25,signif_label!="" & signif_label!="."),
          aes(y=ci_low-.12,label=signif_label),
          fontface="bold",size=GEOM_TEXT)+
  scale_size_continuous(range=c(2.5,9),breaks=c(50,200,400),limits=c(0,meta_n),name = "n") +
  scale_x_discrete(limits=plot_order) +
  coord_flip() +
  labs(x=NULL, y = expression(atop("Sst_25 abundance change", "(" * beta * " ± 95% CI)"))) +
  theme_classic(base_size=BASE) +
theme(axis.text=element_text(size=AXIS_TEXT),axis.title.x=element_text(size=AXIS_TITLE),
      axis.title.y=element_blank(),axis.ticks=element_blank(),panel.grid=element_blank(),
      legend.position = "none")

ggsave("Final_figures/Figures/Forestplot3c.png",p3c,width=10,height=5,dpi=600)
######

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)
seu <- readRDS("Final_figures/Data/Xenium_SCZ_R.rds")

Depleted <- c("Sst_2","Sst_25","Sst_22","Sst_20","Sst_3")
Unaffected <- c("Sst_11","Sst_23","Sst_9","Sst_13","Sst_19","Sst_5","Sst_10","Sst_4","Sst_1","Sst_12","Sst_7")
label_sst <- function(x) ifelse(x %in% Depleted,"Depleted",ifelse(x %in% Unaffected,"Unaffected","Other"))
density_cols <- c(Depleted="#f90202ff",Unaffected="#000042")
scale_factor <- 1.15

rotate_coords <- function(x,y,angle,center=NULL){
  a <- angle*pi/180
  if(is.null(center)) center <- c(mean(x,na.rm=TRUE),mean(y,na.rm=TRUE))
  xc <- x-center[1]; yc <- y-center[2]
  data.frame(x=xc*cos(a)-yc*sin(a)+center[1],y=xc*sin(a)+yc*cos(a)+center[2])
}

scale_coords <- function(x,y,sf,center){
  data.frame(x=center[1]+(x-center[1])*sf,y=center[2]+(y-center[2])*sf)
}

scale_limits <- function(lim,sf){
  m <- mean(lim)
  m+(lim-m)*sf
}

prepare_sample <- function(seu,sample,xlim,ylim,angle,sf){
  obj <- subset(seu,subset=sample_id==sample)
  xy <- FetchData(obj,vars=c("x","y")); xy$y <- -xy$y
  xy <- rotate_coords(xy$x,xy$y,angle)
  obj$x <- xy$x; obj$y <- xy$y

  keep <- with(xy,x>=xlim[1] & x<=xlim[2] & y>=ylim[1] & y<=ylim[2])
  obj <- obj[,keep]

  xy <- scale_coords(obj$x,obj$y,sf,c(mean(xlim),mean(ylim)))
  obj$x <- xy$x; obj$y <- xy$y

  df <- obj@meta.data
  df$`SST type` <- factor(label_sst(df$supertype),levels=c("Depleted","Unaffected","Other"))
  df
}

xlim_con <- c(200,2300); ylim_con <- c(-5200,-2400)
xlim_scz <- c(8700,10500); ylim_scz <- c(-5900,-3500)

con_df <- prepare_sample(seu,"Br5931",xlim_con,ylim_con,-4,scale_factor)
scz_df <- prepare_sample(seu,"Br1139",xlim_scz,ylim_scz,12,scale_factor)

plot_xlim_con <- scale_limits(xlim_con,scale_factor)
plot_ylim_con <- scale_limits(ylim_con,scale_factor)
plot_xlim_scz <- scale_limits(xlim_scz,scale_factor)
plot_ylim_scz <- scale_limits(ylim_scz,scale_factor)

plot_theme <- theme_classic() +
  theme(axis.title=element_blank(),axis.text=element_blank(),
        axis.ticks=element_blank(),axis.line=element_blank())

make_spatial_density <- function(df,xlim,ylim,title,density_width=.15,density_gap=.025){
  xr <- diff(xlim)
  dens_base <- xlim[2]+density_gap*xr
  dens_width <- density_width*xr
  sst_df <- filter(df,`SST type` %in% c("Depleted","Unaffected"))

  dens <- sst_df %>%
    group_by(`SST type`) %>%
    group_modify(~{
      if(nrow(.x)<2) return(tibble(y=numeric(),density=numeric()))
      d <- density(.x$y,from=ylim[1],to=ylim[2],cut=0,n=512)
      tibble(y=d$x,density=d$y)
    }) %>%
    ungroup() %>%
    mutate(x=dens_base+density/max(density,na.rm=TRUE)*dens_width)

  dens_poly <- dens %>%
    group_by(`SST type`) %>%
    group_modify(~bind_rows(
      tibble(y=.x$y[1],x=dens_base),
      select(.x,y,x),
      tibble(y=.x$y[nrow(.x)],x=dens_base)
    )) %>%
    ungroup()

  layer_pos <- df %>%
    filter(!is.na(layer),layer!="",!grepl("Vascular",layer,ignore.case=TRUE)) %>%
    group_by(layer) %>%
    summarise(y=median(y,na.rm=TRUE),.groups="drop")

  ggplot() +
    geom_point(data=filter(df,!`SST type` %in% c("Depleted","Unaffected")),
               aes(x,y),color="grey85",size=2,alpha=.6) +
    geom_point(data=filter(df,`SST type`=="Unaffected"),
               aes(x,y,color=`SST type`),size=3) +
    geom_point(data=filter(df,`SST type`=="Depleted"),
               aes(x,y,color=`SST type`),size=3) +
    geom_point(data=sst_df,aes(x,y),shape=1,color="black",size=1,stroke=.1) +
    geom_polygon(data=dens_poly,
                 aes(x=x,y=y,fill=`SST type`,group=`SST type`),
                 alpha=.5,color=NA) +
    geom_vline(xintercept=dens_base,color="grey70",linewidth=.25) +
    geom_text(data=layer_pos,
              aes(x=xlim[1]+.015*xr,y=y,label=layer),
              inherit.aes=FALSE,hjust=0,size=GEOM_TEXT,color="black") +
scale_color_manual(values=density_cols,breaks=c("Depleted","Unaffected"),
                   labels=c("Depleted","Not depleted"),name=NULL) +
scale_fill_manual(values=density_cols,guide="none") +
guides(color=guide_legend(nrow=1,byrow=TRUE))+
theme(
  legend.position="bottom",
  legend.direction="horizontal",
  legend.spacing.x=unit(10,"pt"),
  legend.key.width=unit(18,"pt"),
  legend.text=element_text(size=AXIS_TITLE)
)+
coord_cartesian(xlim=c(xlim[1],dens_base+dens_width),ylim=ylim,expand=FALSE,clip="on") +
ggtitle(title) + plot_theme +
theme(
  legend.position = "bottom",
  legend.title = element_blank(),
  legend.text = element_text(size = AXIS_TITLE),
  plot.title = element_text(hjust = .5, size = AXIS_TITLE, face = "bold"),
  plot.margin = margin(0,2,2,2)
)}
p_con <- make_spatial_density(con_df,plot_xlim_con,plot_ylim_con,"Control")

p_scz <- make_spatial_density(scz_df,plot_xlim_scz,plot_ylim_scz,"Schizophrenia") +
  theme(legend.position = "none")

legend_d <- ggplot(
  data.frame(
    x = c(.38, .55),
    label = c("Depleted", "Not depleted"),
    col = density_cols
  ),
  aes(x, 0)
) +
  geom_point(aes(color=col),size=3) +
  geom_text(aes(label=label),hjust=0,nudge_x=.02,size=GEOM_TEXT) +
  scale_color_identity() +
  coord_cartesian(xlim=c(0,1),ylim=c(-.1,.1),clip="off") +
  theme_void()

  p3d_plots <- plot_grid(
  p_con+theme(legend.position="none"),
  p_scz+theme(legend.position="none"),
  ncol=2,rel_widths=c(1,1)
)

p3d <- plot_grid(p3d_plots,legend_d,ncol=1,rel_heights=c(1,.06))

ggsave("Final_figures/Figures/Spatial3d.png",p3d,width=16,height=8,dpi=600,bg="white")

#CONCORDANCE
library(dplyr)
library(ggplot2)
library(ggrepel)

res <- read.csv("Final_figures/Data/xenium_crumblr_results_supertype_neuronal.csv")
res$CellType <- res$celltype
final_results <- read.csv("Compositional_analysis/Files/Meta_Neurons.csv")
types <- read.csv("Compositional_analysis/Files/cluster_order_and_colors.csv")

df_concord <- res %>%
  select(CellType, logFC) %>%
  inner_join(final_results %>% select(CellType, estimate, padj), by = "CellType") %>%
  inner_join(types %>% select(cluster_label, class_label), by = c("CellType" = "cluster_label")) %>%
  mutate(sig = padj < 0.2)

spearman_test <- cor.test(df_concord$logFC, df_concord$estimate, method = "spearman", exact = FALSE)
rho_label <- paste0("\u03C1 = ", round(unname(spearman_test$estimate), 2), ", p < 0.001")

p3f <- ggplot(df_concord,aes(estimate,logFC)) +
geom_point(data=filter(df_concord,!sig),aes(fill=class_label),shape=21,size=2,color="black",stroke=.3) +
geom_point(data=filter(df_concord,sig),aes(fill=class_label),shape=21,size=4,color="black",stroke=.8) +
  geom_smooth(method="lm",se=FALSE,linetype="dashed",linewidth=.4,color="grey30") +
  geom_vline(xintercept=0,linetype="dotted",linewidth=.3) +
  geom_hline(yintercept=0,linetype="dotted",linewidth=.3) +
  geom_text_repel(data=filter(df_concord,sig),aes(label=CellType),color="black",size=GEOM_TEXT,box.padding=.6,point.padding=.3,force=4,force_pull=.2,min.segment.length=0,segment.color="grey50",segment.size=.25,max.overlaps=Inf,seed=3,show.legend=FALSE)+
  annotate("text",x=min(df_concord$estimate,na.rm=T)+.01,y=max(df_concord$logFC,na.rm=T)-.02,label=rho_label,hjust=0,vjust=1,size=GEOM_TEXT) +
  scale_fill_manual(values=c("Neuronal: GABAergic"="#8E4585","Neuronal: Glutamatergic"="#2E8B57")) +
 labs (x = expression("Abundance change (meta-analysis, " * beta * ")"),
  y = expression("Abundance change (Xenium, " * beta * ")")) +
  theme_classic(base_size=BASE) +
  theme(panel.grid=element_blank(),axis.text=element_text(size=AXIS_TEXT,color="black"),
        axis.title=element_text(size=AXIS_TITLE,color="black"),axis.ticks=element_blank(),
        legend.position="none")

ggsave("Final_figures/Figures/concordanceplot.png", p3f, width=9.5, height=8, dpi=600)


#
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)

seu <- readRDS("Final_figures/Data/Xenium_SCZ_R.rds")
final_results <- read.csv("Compositional_analysis/Files/Meta_Neurons.csv")
colours <- read.csv("Compositional_analysis/Files/cluster_order_and_colors.csv")

depth_effect <- seu@meta.data %>%
  filter(qc_pass %in% c(TRUE,"True"),subclass=="Sst",!is.na(predicted_norm_depth)) %>%
  group_by(CellType=supertype) %>%
  summarise(mean_depth=mean(predicted_norm_depth),.groups="drop") %>%
  inner_join(final_results %>% select(CellType,estimate,padj),by="CellType") %>%
  left_join(colours %>% select(cluster_label,cluster_color),by=c("CellType"="cluster_label"))

spearman_test <- cor.test(depth_effect$estimate,depth_effect$mean_depth,method="spearman",exact=FALSE)
rho_label <- paste0("\u03C1 = ",round(unname(spearman_test$estimate),2),", p = ",format.pval(spearman_test$p.value,digits=2,eps=0.001))

p3h <- ggplot(depth_effect,aes(estimate,mean_depth)) +
  geom_smooth(method="lm",se=FALSE,linetype="dashed",linewidth=.4,color="grey30") +
  geom_vline(xintercept=0,linetype="dotted",linewidth=.3,color="grey40") +
  geom_hline(yintercept=mean(depth_effect$mean_depth),linetype="dotted",linewidth=.3,color="grey40") +
  geom_point(data=filter(depth_effect,padj>=.2),aes(fill=cluster_color),shape=21,color="black",size=5,stroke=.4) +
  geom_point(data=filter(depth_effect,padj<.2),aes(fill=cluster_color),shape=21,color="black",size=5,stroke=1) +
  geom_text_repel(aes(label=CellType),color="black",size=GEOM_TEXT,box.padding=.7,point.padding=.4,force=4,force_pull=.5,min.segment.length=0,segment.color="grey50",segment.size=.25,max.overlaps=Inf,seed=1,show.legend=FALSE)+
  annotate("text",x=min(depth_effect$estimate)+.01,y=max(depth_effect$mean_depth)-.02,label=rho_label,hjust=0,vjust=1,size=GEOM_TEXT) +
  scale_fill_identity() + scale_y_reverse() +
  labs(x=expression("Abundance change (meta-analysis, "*beta*")"),y="Mean cortical depth") +
  theme_classic(base_size=BASE) +
  theme(axis.text=element_text(size=AXIS_TEXT,color="black"),axis.title=element_text(size=AXIS_TITLE,color="black"),
        axis.ticks=element_blank(),panel.grid=element_blank(),legend.position="none")

ggsave("Final_figures/Figures/p3h_Sst_depth_effect.png",p3h,width=8,height=6,dpi=600)

library(cowplot)

top_row <- plot_grid(p3a, labels="a", label_size=26, label_fontface="bold")

middle_row <- plot_grid(p3c, p3b, ncol=2, labels=c("b","c"), rel_widths=c(0.75,2),
                        label_size=26, label_fontface="bold")

bottom_row <- plot_grid(p3d,plot_grid(p3h,p3f,ncol=1,labels=c("e","f"),
align="v",axis="lr",label_size=26,label_fontface="bold"),ncol=2,rel_widths=c(1.5,.8),
labels=c("d",""),label_size=26,label_fontface="bold")

full <- plot_grid(top_row,middle_row,bottom_row,ncol=1,rel_heights=c(1,.75,1.8))

ggsave("Final_figures/Figures/Figure3_composite.png", full, width=22, height=24, dpi=900, bg="white")
