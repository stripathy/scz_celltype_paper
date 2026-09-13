library(dplyr)

setwd("scz_celltype_paper/snrnaseq")

cohort <- read.csv("Compositional_analysis/Files/crumblr_results_neurons.csv")
meta <- read.csv("Compositional_analysis/Files/Meta_Neurons.csv")

meta_plot <- meta %>%
  mutate(Cohort="Meta-analysis", neglogP=-log10(pval),
         sig_label=case_when(padj<.01~"***",padj<.05~"**",padj<.1~"*",padj<.2~"+",pval<.05~".",TRUE~"")) %>%
  select(CellType,estimate,padj,neglogP,sig_label,Cohort,pval,se)

cohort_plot <- cohort %>%
  mutate(Cohort=recode(Cohort,"Bat"="Batiuk","OFC"="Fröhlich","Multi"="Multiome",
                       "Ruz_McLean"="McLean","Ruz_MtSinai"="MSSM 1","MSSM"="MSSM 2"),
         estimate=logFC,se=abs(logFC/t),padj=as.numeric(adj.P.Val),pval=as.numeric(P.Value),
         neglogP=-log10(pval),
         sig_label=case_when(padj<.01~"***",padj<.05~"**",padj<.1~"*",padj<.2~"+",pval<.05~".",TRUE~"")) %>%
  select(CellType,estimate,padj,neglogP,sig_label,Cohort,pval,se)

plot_data <- bind_rows(meta_plot,cohort_plot) %>%
  mutate(n=case_when(Cohort=="Batiuk"~15,Cohort=="HBCC"~130,Cohort=="MSSM 2"~185,
                     Cohort=="McLean"~32,Cohort=="MSSM 1"~35,Cohort=="Fröhlich"~61,
                     Cohort=="Multiome"~11,Cohort=="Meta-analysis"~469,TRUE~NA_real_))

write.csv(plot_data,"Final_figures/Data/plotdata.csv",row.names=FALSE)