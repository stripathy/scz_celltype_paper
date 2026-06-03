#!/usr/bin/env Rscript
#
# crumblr compositional analysis: SEA-AD pseudoprogression
#
# Runs crumblr + dream on cell type proportions from three SEA-AD datasets:
#   - snRNAseq (89 donors)
#   - MERFISH original labels (27 donors)
#   - MERFISH reclassified by our pipeline (27 donors)
#
# Tests: association of cell type proportions with Continuous Pseudo-progression
# Score (CPS), controlling for age and sex (same covariates as SEA-AD paper).
#
# Formula: ~ scale(CPS) + sex + scale(age)
#
# Input:  output/crumblr_seaad/crumblr_input_{source}_{level}.csv
# Output: output/crumblr_seaad/crumblr_results_{source}_{level}.csv
#         output/crumblr_seaad/crumblr_results_all.csv

library(crumblr)
library(dreamlet)
library(variancePartition)
library(limma)

cat("Libraries loaded\n")

# ── Paths ──────────────────────────────────────────────────────────
base_dir <- path.expand("~/Github/SCZ_Xenium")
io_dir <- file.path(base_dir, "output", "crumblr_seaad")

# ── Discover input files ───────────────────────────────────────────
input_files <- Sys.glob(file.path(io_dir, "crumblr_input_*.csv"))
cat(sprintf("Found %d input files:\n", length(input_files)))
for (f in input_files) cat(sprintf("  %s\n", basename(f)))

all_results <- list()
idx <- 1

for (fpath in input_files) {
  # Parse source and level from filename: crumblr_input_{source}_{level}.csv
  # Level can be: subclass, supertype, subclass_neurons, supertype_neurons
  fname <- tools::file_path_sans_ext(basename(fpath))
  parts <- sub("crumblr_input_", "", fname)
  # Match known level suffixes (longest first to avoid partial matches)
  level_patterns <- c("_subclass_neurons", "_supertype_neurons",
                       "_subclass", "_supertype")
  level <- ""
  source <- parts
  for (lp in level_patterns) {
    if (grepl(paste0(lp, "$"), parts)) {
      level <- sub("^_", "", lp)
      source <- sub(paste0(lp, "$"), "", parts)
      break
    }
  }
  label <- sprintf("%s / %s", source, level)
  cat(sprintf("\n══ Processing: %s ══\n", label))

  # ── Load data ──────────────────────────────────────────────────
  df <- read.csv(fpath, stringsAsFactors = FALSE)
  cat(sprintf("  %d rows, %d donors, %d cell types\n",
              nrow(df), length(unique(df$donor)), length(unique(df$celltype))))

  # ── Pivot to wide count matrix (donors × cell types) ──────────
  count_wide <- reshape(df[, c("donor", "celltype", "count")],
                        idvar = "donor", timevar = "celltype",
                        direction = "wide")
  rownames(count_wide) <- count_wide$donor
  count_wide$donor <- NULL
  colnames(count_wide) <- sub("^count\\.", "", colnames(count_wide))
  count_wide[is.na(count_wide)] <- 0

  # ── Filter: cell types present in ≥50% of samples ─────────────
  presence <- colMeans(count_wide > 0)
  keep <- presence >= 0.5
  cat(sprintf("  Keeping %d / %d types (≥50%% presence)\n",
              sum(keep), length(keep)))
  count_wide <- count_wide[, keep, drop = FALSE]

  if (ncol(count_wide) < 2) {
    cat("  Skipping: fewer than 2 cell types\n")
    next
  }

  # ── Build metadata ─────────────────────────────────────────────
  meta <- unique(df[, c("donor", "CPS", "sex", "age")])
  rownames(meta) <- meta$donor
  meta <- meta[rownames(count_wide), ]
  meta$sex <- factor(meta$sex)
  meta$CPS_num <- as.numeric(meta$CPS)
  meta$age_num <- as.numeric(meta$age)

  cat(sprintf("  CPS range: [%.3f, %.3f], n=%d donors\n",
              min(meta$CPS_num, na.rm = TRUE),
              max(meta$CPS_num, na.rm = TRUE),
              nrow(meta)))

  # ── Run crumblr ────────────────────────────────────────────────
  count_mat <- as.matrix(count_wide)

  cobj <- tryCatch(
    crumblr(count_mat),
    error = function(e) {
      cat(sprintf("  crumblr error: %s\n", e$message))
      return(NULL)
    }
  )
  if (is.null(cobj)) next

  # ── Fit dream model ────────────────────────────────────────────
  # Same covariates as SEA-AD paper: CPS + sex + age
  form <- ~ scale(CPS_num) + sex + scale(age_num)

  fit <- tryCatch({
    f <- dream(cobj, form, meta)
    eBayes(f)
  }, error = function(e) {
    cat(sprintf("  dream error: %s\n", e$message))
    return(NULL)
  })
  if (is.null(fit)) next

  # ── Extract results for CPS effect ────────────────────────────
  res <- topTable(fit, coef = "scale(CPS_num)", number = Inf, sort.by = "none")
  res$celltype <- rownames(res)
  res$level <- level
  res$source <- source

  # Standard error: SE = logFC / t
  res$SE <- res$logFC / res$t

  # Per-level FDR
  res$FDR <- p.adjust(res$P.Value, method = "BH")

  # Save individual results
  res_sorted <- res[order(res$P.Value), ]
  out_file <- file.path(io_dir, sprintf("crumblr_results_%s_%s.csv", source, level))
  write.csv(res_sorted, out_file, row.names = FALSE)
  cat(sprintf("  Saved: %s (%d types)\n", basename(out_file), nrow(res)))

  all_results[[idx]] <- res
  idx <- idx + 1

  # ── Print summary ──────────────────────────────────────────────
  n_fdr05 <- sum(res$FDR < 0.05)
  n_fdr10 <- sum(res$FDR < 0.10)
  n_nom05 <- sum(res$P.Value < 0.05)
  cat(sprintf("  FDR < 0.05: %d | FDR < 0.10: %d | nom p < 0.05: %d\n",
              n_fdr05, n_fdr10, n_nom05))

  # Show significant hits
  sig <- res[res$FDR < 0.10, ]
  sig <- sig[order(sig$P.Value), ]
  if (nrow(sig) > 0) {
    cat("  FDR < 0.10 hits:\n")
    for (i in 1:nrow(sig)) {
      r <- sig[i, ]
      d <- ifelse(r$logFC > 0, "↑CPS", "↓CPS")
      cat(sprintf("    %-30s %s logFC=%+.4f p=%.5f FDR=%.4f\n",
                  r$celltype, d, r$logFC, r$P.Value, r$FDR))
    }
  } else {
    cat("  No FDR < 0.10 hits\n")
  }

  # Also show top 5 nominal hits
  top5 <- res_sorted[1:min(5, nrow(res_sorted)), ]
  cat("  Top 5 nominal hits:\n")
  for (i in 1:nrow(top5)) {
    r <- top5[i, ]
    d <- ifelse(r$logFC > 0, "↑CPS", "↓CPS")
    cat(sprintf("    %-30s %s logFC=%+.4f p=%.5f FDR=%.4f\n",
                r$celltype, d, r$logFC, r$P.Value, r$FDR))
  }
}

# ── Combine all results ────────────────────────────────────────────
if (length(all_results) > 0) {
  combined <- do.call(rbind, all_results)
  combined <- combined[order(combined$source, combined$level, combined$P.Value), ]
  out_file <- file.path(io_dir, "crumblr_results_all.csv")
  write.csv(combined, out_file, row.names = FALSE)
  cat(sprintf("\nSaved combined: %s (%d total rows)\n",
              basename(out_file), nrow(combined)))
} else {
  cat("\nNo results to combine!\n")
}

cat("\nDone!\n")
