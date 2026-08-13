#!/usr/bin/env Rscript
#
# crumblr compositional analysis of the SEA-AD DLPFC (A9) cohort along the
# Alzheimer's continuous pseudo-progression score (CPS). Produces the AD axis
# of Figure 4j.
#
# Run build_seaad_cps_input.py first.
#
# Model
# -----
#   ~ scale(CPS) + sex + scale(age) + scale(PMI)
#
# The reported coefficient is scale(CPS), so logFC is the compositional slope
# per standard deviation of CPS; negative = the supertype declines as AD
# advances. Gabitto et al. (2024) model SEA-AD A9 composition with
# sex + age + race + CPS and no PMI -- SEA-AD is a rapid-autopsy cohort with PMI
# restricted to 3.2-11.4 h, so there is little variance to model. We add PMI for
# consistency with the schizophrenia compositional analyses in this paper; it
# leaves the panel-4j correlation unchanged (rho = 0.868 either way). Race is
# omitted: among these 80 donors it is 74 White / 3 Asian / 3 multiple, which is
# too imbalanced to estimate.
#
# Supertypes present in fewer than 50% of donors are dropped, matching the
# filter used in the schizophrenia crumblr analyses.
#
# Output
# ------
#   results/intermediates/seaad_dlpfc/crumblr_results_supertype_neurons.csv

suppressPackageStartupMessages({
  library(crumblr); library(dreamlet); library(variancePartition); library(limma)
})

repo <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=",
          commandArgs(trailingOnly = FALSE), value = TRUE)[1])), "..", ".."))
io <- file.path(repo, "crossdisorder", "results")

in_file <- file.path(io, "crumblr_input_mtg_supertype_neurons.csv")
if (!file.exists(in_file))
  stop("Missing ", in_file, " -- run build_seaad_cps_input.py first")

df <- read.csv(in_file, stringsAsFactors = FALSE)
cat(sprintf("%d rows, %d donors, %d supertypes\n",
            nrow(df), length(unique(df$donor)), length(unique(df$celltype))))

# ---- donor x celltype count matrix --------------------------------------
wide <- reshape(df[, c("donor", "celltype", "count")], idvar = "donor",
                timevar = "celltype", direction = "wide")
rownames(wide) <- wide$donor; wide$donor <- NULL
colnames(wide) <- sub("^count\\.", "", colnames(wide))
wide[is.na(wide)] <- 0

keep <- colMeans(wide > 0) >= 0.5
cat(sprintf("keeping %d / %d supertypes (>=50%% donor presence)\n",
            sum(keep), length(keep)))
wide <- wide[, keep, drop = FALSE]

# ---- donor covariates ----------------------------------------------------
meta <- unique(df[, c("donor", "CPS", "sex", "age", "pmi")])
rownames(meta) <- meta$donor
meta <- meta[rownames(wide), ]
meta$sex     <- factor(meta$sex)
meta$CPS_num <- as.numeric(meta$CPS)
meta$age_num <- as.numeric(meta$age)
meta$pmi_num <- as.numeric(meta$pmi)
cat(sprintf("%d donors; CPS [%.3f, %.3f]; PMI [%.1f, %.1f] h\n", nrow(meta),
            min(meta$CPS_num), max(meta$CPS_num),
            min(meta$pmi_num), max(meta$pmi_num)))

# ---- fit -----------------------------------------------------------------
cobj <- crumblr(as.matrix(wide))
fit  <- eBayes(dream(cobj, ~ scale(CPS_num) + sex + scale(age_num) + scale(pmi_num),
                     meta))

res <- topTable(fit, coef = "scale(CPS_num)", number = Inf, sort.by = "none")
res$celltype <- rownames(res)
res$SE  <- res$logFC / res$t            # dream returns t, not SE directly
res$FDR <- p.adjust(res$P.Value, method = "BH")
res <- res[order(res$P.Value), ]

out_file <- file.path(io, "crumblr_results_mtg_supertype_neurons.csv")
write.csv(res, out_file, row.names = FALSE)

cat(sprintf("FDR<0.05: %d | FDR<0.10: %d | nominal p<0.05: %d\n",
            sum(res$FDR < 0.05), sum(res$FDR < 0.10), sum(res$P.Value < 0.05)))
sst <- res[grepl("^Sst_[0-9]+$", res$celltype), c("celltype", "logFC", "SE", "P.Value", "FDR")]
cat("\nSst supertypes (negative logFC = declines as AD advances):\n")
print(sst[order(sst$logFC), ], row.names = FALSE, digits = 3)
cat(sprintf("\n-> %s\n", out_file))
