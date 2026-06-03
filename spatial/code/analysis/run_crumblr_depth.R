#!/usr/bin/env Rscript
#
# CLR compositional analysis with depth x diagnosis interaction
#
# Tests whether cell type composition changes with cortical depth
# differ between SCZ and Control.
#
# Approach: manually compute CLR (centered log-ratio) transform of
# count data binned by depth, then fit linear mixed models via lme4.
#
# Model per cell type:
#   CLR(proportion) ~ depth * diagnosis + sex + scale(age) + scale(pmi) + (1|donor)
#
# This avoids crumblr+dream's rank-deficiency issues by fitting one
# model per cell type rather than using the joint framework.
#
# Input:  output/crumblr/crumblr_depth_input_subclass.csv
# Output: output/crumblr/crumblr_depth_results_subclass.csv
#         output/crumblr/crumblr_depth_interaction_subclass.csv

library(lme4)
library(lmerTest)  # for p-values via Satterthwaite

cat("Libraries loaded\n")

# ── Paths ──────────────────────────────────────────────────────────
base_dir <- path.expand("~/Github/SCZ_Xenium")
in_dir <- file.path(base_dir, "output", "crumblr")
out_dir <- in_dir

# ── Load data ──────────────────────────────────────────────────────
fpath <- file.path(in_dir, "crumblr_depth_input_subclass.csv")
cat(sprintf("Reading: %s\n", fpath))
df <- read.csv(fpath, stringsAsFactors = FALSE)

n_donors <- length(unique(df$donor))
n_types <- length(unique(df$celltype))
n_bins <- length(unique(df$depth_bin))
cat(sprintf("  %d rows, %d donors, %d cell types, %d depth bins\n",
            nrow(df), n_donors, n_types, n_bins))

# ── Compute CLR transform per observation (donor x depth_bin) ──────
# CLR_i = log(count_i / geometric_mean(counts))
# Use pseudocount of 0.5 for zeros
cat("\nComputing CLR transform...\n")

df$obs_id <- paste(df$donor, df$depth_bin, sep = ":")

# Pivot to wide for CLR computation
wide <- reshape(df[, c("obs_id", "celltype", "count")],
                idvar = "obs_id", timevar = "celltype",
                direction = "wide")
rownames(wide) <- wide$obs_id
wide$obs_id <- NULL
colnames(wide) <- sub("^count\\.", "", colnames(wide))
wide[is.na(wide)] <- 0

# Add pseudocount and compute CLR
mat <- as.matrix(wide) + 0.5
log_mat <- log(mat)
geom_mean <- rowMeans(log_mat)  # log of geometric mean
clr_mat <- log_mat - geom_mean  # CLR transform

cat(sprintf("  CLR matrix: %d obs x %d cell types\n",
            nrow(clr_mat), ncol(clr_mat)))

# ── Build metadata per observation ─────────────────────────────────
meta <- df[!duplicated(df$obs_id), c("obs_id", "donor", "depth_bin",
                                      "depth_midpoint", "diagnosis",
                                      "sex", "age", "pmi")]
rownames(meta) <- meta$obs_id
meta <- meta[rownames(clr_mat), ]

meta$diagnosis <- factor(meta$diagnosis, levels = c("Control", "SCZ"))
meta$sex <- factor(meta$sex)
meta$depth <- as.numeric(meta$depth_midpoint)

cat(sprintf("  %d Control donors, %d SCZ donors\n",
            sum(!duplicated(meta$donor[meta$diagnosis == "Control"])),
            sum(!duplicated(meta$donor[meta$diagnosis == "SCZ"]))))

# ── Filter: cell types present in >=50% of donors ─────────────────
# Check from original long-format data (simpler, avoids wide matrix indexing)
donor_type_counts <- aggregate(count ~ donor + celltype, data = df,
                               FUN = sum, na.rm = TRUE)
donor_type_counts <- donor_type_counts[donor_type_counts$count > 0, ]
type_donor_n <- table(donor_type_counts$celltype)
donor_presence <- type_donor_n / n_donors
keep_types <- names(donor_presence[donor_presence >= 0.5])
# Also ensure they're in the CLR matrix
celltypes <- intersect(keep_types, colnames(clr_mat))
cat(sprintf("  Testing %d / %d types (>=50%% donor presence)\n",
            length(celltypes), ncol(clr_mat)))

# ── Fit mixed models per cell type ─────────────────────────────────
cat("\nFitting mixed models per cell type...\n")
cat(sprintf("  Formula: CLR ~ depth * diagnosis + sex + scale(age) + scale(pmi) + (1|donor)\n\n"))

all_results <- data.frame(celltype=character(), coefficient=character(),
                          logFC=numeric(), SE=numeric(), t=numeric(),
                          P.Value=numeric(), FDR=numeric(),
                          stringsAsFactors=FALSE)

for (ct in celltypes) {
  # Build per-celltype data frame
  mdf <- data.frame(
    clr = clr_mat[meta$obs_id, ct],
    depth = meta$depth,
    diagnosis = meta$diagnosis,
    sex = meta$sex,
    age = as.numeric(meta$age),
    pmi = as.numeric(meta$pmi),
    donor = factor(meta$donor)
  )

  # Fit mixed model
  fit <- tryCatch({
    lmer(clr ~ depth * diagnosis + sex + scale(age) + scale(pmi) + (1|donor),
         data = mdf, REML = TRUE)
  }, error = function(e) {
    cat(sprintf("  %s: model failed (%s)\n", ct, e$message))
    return(NULL)
  })

  if (is.null(fit)) next

  # Extract fixed effects with Satterthwaite p-values
  coefs <- summary(fit)$coefficients

  # Get all coefficients
  for (coef_name in rownames(coefs)) {
    if (coef_name == "(Intercept)") next

    # Map coefficient names to clean labels
    clean_name <- coef_name
    if (grepl("diagnosisSCZ$", coef_name) && !grepl(":", coef_name)) {
      clean_name <- "diagnosisSCZ"
    } else if (coef_name == "depth") {
      clean_name <- "depth"
    } else if (grepl("depth:diagnosisSCZ", coef_name)) {
      clean_name <- "depth:diagnosisSCZ"
    } else {
      next  # skip sex, age, pmi coefficients
    }

    est <- coefs[coef_name, "Estimate"]
    se <- coefs[coef_name, "Std. Error"]
    tval <- coefs[coef_name, "t value"]
    # Satterthwaite p-value (from lmerTest)
    pval <- coefs[coef_name, "Pr(>|t|)"]

    all_results <- rbind(all_results, data.frame(
      celltype = ct,
      coefficient = clean_name,
      logFC = est,
      SE = se,
      t = tval,
      P.Value = pval,
      stringsAsFactors = FALSE
    ))
  }
}

# ── FDR correction (per coefficient, across cell types) ────────────
for (coef_name in unique(all_results$coefficient)) {
  mask <- all_results$coefficient == coef_name
  all_results$FDR[mask] <- p.adjust(all_results$P.Value[mask], method = "BH")
}

# Sort
all_results <- all_results[order(all_results$coefficient, all_results$P.Value), ]

# ── Save all results ───────────────────────────────────────────────
out_all <- file.path(out_dir, "crumblr_depth_results_subclass.csv")
write.csv(all_results, out_all, row.names = FALSE)
cat(sprintf("\nSaved: %s (%d rows)\n", basename(out_all), nrow(all_results)))

# ── Save interaction results separately ────────────────────────────
interact <- all_results[all_results$coefficient == "depth:diagnosisSCZ", ]
out_int <- file.path(out_dir, "crumblr_depth_interaction_subclass.csv")
write.csv(interact, out_int, row.names = FALSE)
cat(sprintf("Saved: %s (%d rows)\n", basename(out_int), nrow(interact)))

# ── Print summaries ────────────────────────────────────────────────
for (coef_name in c("diagnosisSCZ", "depth", "depth:diagnosisSCZ")) {
  sub <- all_results[all_results$coefficient == coef_name, ]
  if (nrow(sub) == 0) next

  n_fdr05 <- sum(sub$FDR < 0.05, na.rm = TRUE)
  n_fdr10 <- sum(sub$FDR < 0.10, na.rm = TRUE)
  n_nom05 <- sum(sub$P.Value < 0.05, na.rm = TRUE)

  cat(sprintf("\n══ %s ══\n", coef_name))
  cat(sprintf("  FDR < 0.05: %d | FDR < 0.10: %d | nom p < 0.05: %d / %d\n",
              n_fdr05, n_fdr10, n_nom05, nrow(sub)))

  # Print all results sorted by p-value
  sub <- sub[order(sub$P.Value), ]
  cat(sprintf("\n  %-25s  %8s  %8s  %10s  %8s\n",
              "Celltype", "logFC", "t", "P.Value", "FDR"))
  cat(paste0("  ", paste(rep("-", 65), collapse = ""), "\n"))
  for (i in 1:nrow(sub)) {
    r <- sub[i, ]
    star <- ifelse(r$FDR < 0.05, "***",
            ifelse(r$FDR < 0.10, "**",
            ifelse(r$P.Value < 0.05, "*", "")))
    cat(sprintf("  %-25s  %+8.4f  %8.2f  %10.6f  %8.4f  %s\n",
                r$celltype, r$logFC, r$t, r$P.Value, r$FDR, star))
  }
}

cat("\nDone!\n")
