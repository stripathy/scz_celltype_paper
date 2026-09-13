library(Seurat)
library(AnnotationDbi)
library(org.Hs.eg.db)

update_genes <- function(seu) {
  genes <- rownames(seu)
  genes_clean <- sub("\\..*$", "", genes)  # remove Ensembl version suffixes if present

  gene_symbols <- mapIds(
    org.Hs.eg.db,
    keys = genes_clean,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

  gene_symbols <- unname(gene_symbols)
  gene_symbols[is.na(gene_symbols) | gene_symbols == ""] <- genes[is.na(gene_symbols) | gene_symbols == ""]
  gene_symbols <- make.unique(gene_symbols)

  rownames(seu) <- gene_symbols
  seu}


HBCC1 <- readRDS("/scratch/nendresz/PsychAD/Data/HBCC_1.rds")
HBCC2 <- readRDS("/scratch/nendresz/PsychAD/Data/HBCC_2.rds")
HBCC1$Age <- as.numeric(as.character(HBCC1$Age))
HBCC2$Age <- as.numeric(as.character(HBCC2$Age))


HBCC1 <- update_genes(HBCC1)
HBCC2 <- update_genes(HBCC2)

saveRDS(HBCC1, "/scratch/nendresz/PsychAD/Data/HBCC_1_symbols.rds")
saveRDS(HBCC2, "/scratch/nendresz/PsychAD/Data/HBCC_2_symbols.rds")


MSSM1 <- readRDS("/scratch/nendresz/PsychAD/Data/MSSM_1.rds")
MSSM2 <- readRDS("/scratch/nendresz/PsychAD/Data/MSSM_2.rds")
MSSM3 <- readRDS("/scratch/nendresz/PsychAD/Data/MSSM_3.rds")
MSSM4 <- readRDS("/scratch/nendresz/PsychAD/Data/MSSM_4.rds")

MSSM1$Age <- as.numeric(as.character(MSSM1$Age))
MSSM2$Age <- as.numeric(as.character(MSSM2$Age))
MSSM3$Age <- as.numeric(as.character(MSSM3$Age))
MSSM4$Age <- as.numeric(as.character(MSSM4$Age))

MSSM1 <- update_genes(MSSM1)
MSSM2 <- update_genes(MSSM2)
MSSM3 <- update_genes(MSSM3)
MSSM4 <- update_genes(MSSM4)


saveRDS(MSSM1, "/scratch/nendresz/PsychAD/Data/MSSM_1_symbols.rds")
saveRDS(MSSM2, "/scratch/nendresz/PsychAD/Data/MSSM_2_symbols.rds")
saveRDS(MSSM3, "/scratch/nendresz/PsychAD/Data/MSSM_3_symbols.rds")
saveRDS(MSSM4, "/scratch/nendresz/PsychAD/Data/MSSM_4_symbols.rds")