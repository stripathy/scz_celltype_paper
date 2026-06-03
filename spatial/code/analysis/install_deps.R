if (!require("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
BiocManager::install(c("limma", "edgeR", "variancePartition", "zellkonverter"), update=FALSE, ask=FALSE)
