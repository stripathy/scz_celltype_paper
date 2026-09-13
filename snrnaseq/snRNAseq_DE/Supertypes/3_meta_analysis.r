# conda activate r_env_meta_analysis

setwd("scz_celltype_paper/snrnaseq/snRNAseq_DE/Supertypes")
library(dplyr)
library(metafor)

cohorts <- c("Batiuk","Frohlich","MSSM1","McLean","MSSM2","HBCC","Multiome")

# Get all supertypes from DE result files
files <- list.files("Files", pattern="^DE_results_.*\\.rds$")
cell_types <- basename(files)

for(cohort in cohorts)
  cell_types <- sub(paste0("^DE_results_",cohort,"_"),"",cell_types)

cell_types <- unique(sub("\\.rds$","",cell_types))

# Load cohort DE results
de_list <- list()

for(ct in cell_types){
  de_list[[ct]] <- list()

  for(cohort in cohorts){
    f <- paste0("Files/DE_results_",cohort,"_",ct,".rds")

    if(!file.exists(f)){
      message("Skipping missing: ",f)
      next
    }

    de_list[[ct]][[cohort]] <- readRDS(f)
  }
}

# Meta-analysis
meta_all <- list()

for(ct in cell_types){

  if(length(de_list[[ct]])==0) next
  message("Meta-analysis for ",ct)

  all_data <- bind_rows(lapply(names(de_list[[ct]]),function(cohort)
    de_list[[ct]][[cohort]] %>%
      mutate(cohort=cohort,SE=abs(logFC/t))))

  meta_results <- list()

  for(g in unique(all_data$genes)){

    df <- all_data %>%
      filter(genes==g,!is.na(SE),SE!=0,!is.na(logFC))

    if(nrow(df)>4){

      res <- tryCatch({

        m <- rma(yi=logFC,sei=SE,data=df,method="REML")

        data.frame(
          cell_type=ct,
          genes=g,
          estimate=as.numeric(m$b),
          se=m$se,
          pval=m$pval,
          ci.lb=m$ci.lb,
          ci.ub=m$ci.ub,
          k=m$k,
          tau2=m$tau2,
          I2=m$I2
        )

      },error=function(e) NULL)

      if(!is.null(res))
        meta_results[[g]] <- res
    }
  }

  final_results <- bind_rows(meta_results)

  if(nrow(final_results)>0)
    final_results$padj <- p.adjust(final_results$pval,method="fdr")

  meta_all[[ct]] <- final_results

  write.csv(
    final_results,
    paste0("Files/meta_results_",ct,".csv"),
    row.names=FALSE
  )
}

meta_supertype <- bind_rows(meta_all)

write.csv(meta_supertype,"Files/DE_meta_results_supertype.csv",row.names=FALSE)