#mamba activate r_plot
library(dplyr)
library(ggplot2)
library(scales)
library(readr)
library(patchwork)
setwd("FINAL_FIGS")


library(dplyr); library(tidyr); library(ggplot2); library(scales)
library(readr); library(patchwork); library(ggpattern)

BASE <- 16; AXIS_TITLE <- 18; AXIS_TEXT <- 16; BAR_TEXT <- 5
HEADER <- 6; LEGEND_TITLE <- 13; LEGEND_TEXT <- 16; HEADER_X <- 92

all_meta_counts <- read.csv("/scratch/nendresz/P1_Compositional_analysis/Files/7_cohorts_metadata_names.csv", row.names=1,check.names=FALSE)

all_meta_counts$Cohort <- recode(all_meta_counts$Cohort,"OFC"="Fröhlich","MtSinai"="MSSM 1","MSSM"="MSSM 2")

meta_xen <- read_csv("/scratch/nendresz/Xenium/metadata.csv",show_col_types=FALSE) %>% filter(qc_pass)
if(!"Diagnosis" %in% names(meta_xen) && "diagnosis" %in% names(meta_xen)) meta_xen <- rename(meta_xen,Diagnosis=diagnosis)

clean_dx <- function(x) case_when(
  grepl("control|ctrl|^con$",x,ignore.case=TRUE)~"Control",
  grepl("schizophrenia|scz|case",x,ignore.case=TRUE)~"Schizophrenia",
  TRUE~as.character(x)
)

all_meta_counts <- mutate(all_meta_counts,Diagnosis=clean_dx(Diagnosis))
meta_xen <- mutate(meta_xen,Diagnosis=clean_dx(Diagnosis))

cohort_cols <- c("Batiuk"="#1f93c6","Fröhlich"="#ff4500","HBCC"="#8ac926","McLean"="#d100ff",
                 "MSSM 1"="#11823b","MSSM 2"="#2ec4b6","Multiome"="#a87400","Xenium"="#1b9e77")

cohort_donors <- bind_rows(
  all_meta_counts %>% transmute(Cohort,Donor=as.character(Donor),Diagnosis) %>% distinct(),
  meta_xen %>% transmute(Cohort="Xenium",Donor=as.character(sample_id),Diagnosis) %>% distinct()
) %>%
  filter(Cohort %in% names(cohort_cols),Diagnosis %in% c("Control","Schizophrenia")) %>%
  count(Cohort,Diagnosis,name="n") %>%
  complete(Cohort,Diagnosis=c("Control","Schizophrenia"),fill=list(n=0)) %>%
  pivot_wider(names_from=Diagnosis,values_from=n) %>%
  rename(n_control=Control,n_case=Schizophrenia) %>%
  mutate(n_donors=n_control+n_case)

sn_order <- cohort_donors %>% filter(Cohort!="Xenium") %>% arrange(n_donors) %>% pull(Cohort)

cohort_donors <- cohort_donors %>% mutate(
  y=case_when(
    Cohort=="Xenium"~1.42,
    Cohort==sn_order[1]~1.78,Cohort==sn_order[2]~1.97,Cohort==sn_order[3]~2.16,
    Cohort==sn_order[4]~2.35,Cohort==sn_order[5]~2.54,Cohort==sn_order[6]~2.73,
    Cohort==sn_order[7]~2.92
  ),
  xmin_control=0,xmax_control=n_control,
  xmin_case=n_control,xmax_case=n_control+n_case,
  donor_label=paste0(n_control,"/",n_case)
)
legend_df <- tibble(Diagnosis=factor(c("Control","Schizophrenia"),levels=c("Control","Schizophrenia")),
                    xmin=NA_real_,xmax=NA_real_,ymin=NA_real_,ymax=NA_real_)

p_cohort_donors <- ggplot(cohort_donors) +
  geom_rect(aes(xmin=xmin_control,xmax=xmax_control,ymin=y-.09,ymax=y+.09,fill=Cohort),colour=NA) +
  geom_rect_pattern(
    aes(xmin=xmin_case,xmax=xmax_case,ymin=y-.09,ymax=y+.09,fill=Cohort),
    pattern="stripe",pattern_fill="white",pattern_colour="white",
    pattern_angle=45,pattern_density=.35,pattern_spacing=.035,
    colour="black",linewidth=.3,show.legend=FALSE
  ) +
  geom_rect_pattern(
    data=legend_df,
    aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,pattern=Diagnosis),
    fill="grey55",colour="black",linewidth=.3,
    pattern_fill="white",pattern_colour="white",
    pattern_angle=45,pattern_density=.35,pattern_spacing=.035,
    inherit.aes=FALSE,show.legend=TRUE
  ) +
  geom_text(aes(n_donors,y,label=donor_label),hjust=-.15,size=BAR_TEXT) +
  annotate("segment",x=0,xend=55,y=3.10,yend=3.10,colour="grey60",linewidth=.8) +
  annotate("segment",x=129,xend=210,y=3.10,yend=3.10,colour="grey60",linewidth=.8) +
  annotate("text",x=HEADER_X,y=3.10,label="snRNA-seq",fontface="bold",size=HEADER) +
  annotate("segment",x=0,xend=66,y=1.60,yend=1.60,colour="grey60",linewidth=.8) +
  annotate("segment",x=118,xend=210,y=1.60,yend=1.60,colour="grey60",linewidth=.8) +
  annotate("text",x=HEADER_X,y=1.60,label="Spatial",fontface="bold",size=HEADER) +
  scale_fill_manual(values=cohort_cols,guide="none") +
  scale_pattern_manual(
    name=NULL,
    values=c("Control"="none","Schizophrenia"="stripe")
  ) +
  guides(pattern=guide_legend(
    direction="vertical",byrow=TRUE,
    override.aes=list(
      fill="grey55",colour="black",
      pattern_fill="white",pattern_colour="white"
    )
  )) +
  scale_y_continuous(
    breaks=cohort_donors$y,labels=cohort_donors$Cohort,
    limits=c(1.30,3.18),expand=c(0,0)
  ) +
  scale_x_continuous(
    breaks=c(0,50,100,150,200),limits=c(0,240),
    expand=expansion(mult=c(.01,0))
  ) +
labs(x="Number of donors", y=NULL, caption=NULL) +
  theme_classic(base_size=BASE) +
  theme(
    axis.title=element_text(size=AXIS_TITLE),
    axis.text=element_text(size=AXIS_TEXT),
    axis.ticks.y=element_blank(),
    legend.position=c(.64,.29),
    legend.justification=c(.5,.5),
    legend.direction="vertical",
    legend.text=element_text(size=LEGEND_TEXT),
    legend.key.height=unit(.35,"cm"),
    legend.key.width=unit(.55,"cm"),
    legend.spacing.y=unit(.02,"cm"),
    legend.background=element_rect(fill="white",colour=NA),
    plot.caption=element_text(size=10,hjust=0),
    plot.margin=margin(0,6,0,6)
  )

ggsave("Figures/cohort_donor_counts.png", p_cohort_donors,width=6,height=4,dpi=300,bg="white")




# Vertical total cell-count bars 
sn_order <- c("Multiome","Batiuk","McLean","MSSM 1","Fröhlich","HBCC","MSSM 2")

cell_cols <- setdiff(
  colnames(all_meta_counts),
  c("Cohort","Donor","Age","Sex","Diagnosis","PMI"))

plot_cells <- bind_rows(
  all_meta_counts %>%
    filter(Cohort %in% sn_order) %>%
    mutate(total_cells=rowSums(across(all_of(cell_cols)),na.rm=TRUE)) %>%
    group_by(Cohort) %>%
    summarise(Cells=sum(total_cells),.groups="drop") %>%
    mutate(Technology="snRNA-seq"),
  tibble(Cohort="Xenium",Cells=nrow(meta_xen),Technology="Xenium")
) %>%
  mutate(
    Cohort=factor(Cohort,levels=c(sn_order,"Xenium")),
    Technology=factor(Technology,levels=c("snRNA-seq","Xenium")))

p_cells <- ggplot(plot_cells,aes(x=Technology,y=Cells,fill=Cohort)) +
  geom_col(width=.65) +
  scale_fill_manual(
    values=cohort_cols[c(sn_order,"Xenium")],
    breaks=c(sn_order,"Xenium")
  ) +
  scale_y_continuous(
    labels=comma,
    expand=expansion(mult=c(0,.05))
  ) +
  labs(x=NULL,y="Number of cells") +
  theme_classic(base_size=BASE) +
  theme(
    axis.title=element_text(size=AXIS_TITLE),
    axis.text=element_text(size=AXIS_TEXT),
    axis.text.x=element_text(size=AXIS_TEXT,angle=45,hjust=1,vjust=1),
    axis.ticks.x=element_blank(),
    legend.title=element_text(size=LEGEND_TITLE,face="bold"),
    legend.text=element_text(size=LEGEND_TEXT),
    legend.position="none"
  )

ggsave("Figures/snrnaseq_xenium_cell_counts.png",  p_cells,width=4,height=5,dpi=300,bg="white")

# Combined panel
tog <- p_cells | free(p_cohort_donors, side="b")
tog <- tog + plot_layout(widths=c(.5,1.6))

ggsave("Figures/1a.png",tog,width=11,height=5,dpi=300,bg="white")
ggsave("Figures/1a.svg",tog,width=11,height=5,dpi=300,bg="white")
