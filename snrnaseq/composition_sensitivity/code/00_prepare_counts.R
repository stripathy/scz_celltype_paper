# Build the analysis table from Nicole's raw per-donor cell-count export.
#
# Raw input (not tracked; ask Nicole for the current version):
#   7_cohorts_metadata_names.csv - one row per donor, columns Cohort, Donor, Age, Sex,
#   Diagnosis, PMI followed by one column per SEA-AD supertype holding nucleus counts.
# Harmonization applied here: Sex and Diagnosis case-folded (Schizophrenia -> SCZ),
# dataset names mapped to the ones used in the paper (MSSM -> MSSM 2, MtSinai -> MSSM 1,
# OFC -> Fröhlich), and the 109 neuronal supertypes of Fig. 3a selected. PMI is in hours
# in every dataset; the Multiome dataset has none, and one MSSM 2 donor is coded 0 h,
# which 01 treats as missing.
#
# Usage (from anywhere):
#   Rscript code/00_prepare_counts.R [path/to/7_cohorts_metadata_names.csv]   # verify only
#   Rscript code/00_prepare_counts.R [path] --write                          # rewrite data/
# Without --write the script only checks that it reproduces the tracked
# data/neuron_counts_469donors.csv, so a new export can be diffed against the old one.
suppressPackageStartupMessages({library(readr); library(dplyr)})
args <- commandArgs(trailingOnly = FALSE)
fp   <- sub("^--file=", "", args[grep("^--file=", args)])
HERE <- if (length(fp)) normalizePath(file.path(dirname(fp), "..")) else normalizePath(getwd())
opts <- commandArgs(trailingOnly = TRUE)
write_out <- "--write" %in% opts
raw_path  <- setdiff(opts, "--write")
raw_path  <- if (length(raw_path)) raw_path[1] else "~/Downloads/7_cohorts_metadata_names.csv"
stopifnot(file.exists(path.expand(raw_path)))

sup <- read_csv(file.path(HERE, "data/neuronal_supertypes_109.csv"), show_col_types = FALSE)$supertype
raw <- read_csv(path.expand(raw_path), show_col_types = FALSE)
raw <- raw[, !grepl("^\\.\\.\\.|^Unnamed", names(raw))]
missing <- setdiff(sup, names(raw))
if (length(missing)) stop("supertypes missing from the raw export: ", paste(missing, collapse = ", "))

out <- raw |>
  transmute(dataset = recode(Cohort, MSSM = "MSSM 2", MtSinai = "MSSM 1", OFC = "Fröhlich"),
            donor = Donor,
            dx  = recode(tools::toTitleCase(tolower(Diagnosis)), Schizophrenia = "SCZ"),
            age = Age,
            sex = tools::toTitleCase(tolower(Sex)),
            PMI = PMI) |>
  bind_cols(raw[, sup])
stopifnot(all(out$dx %in% c("Control", "SCZ")), all(out$sex %in% c("Male", "Female")))
cat(sprintf("%d donors, %d datasets, %d neuronal supertypes\n", nrow(out), n_distinct(out$dataset), length(sup)))
print(table(out$dataset, out$dx))
cat("PMI missing by dataset: ");
print(out |> group_by(dataset) |> summarise(n_missing = sum(is.na(PMI) | PMI == 0), .groups = "drop") |> filter(n_missing > 0))

dest <- file.path(HERE, "data/neuron_counts_469donors.csv")
if (write_out) {
  write_csv(out, dest); cat(sprintf("wrote %s\n", dest))
} else if (file.exists(dest)) {
  old <- read_csv(dest, show_col_types = FALSE)
  same <- identical(dim(old), dim(out)) && isTRUE(all.equal(as.data.frame(old), as.data.frame(out)))
  cat(if (same) "matches the tracked data/neuron_counts_469donors.csv\n"
      else "DIFFERS from the tracked data/neuron_counts_469donors.csv (re-run with --write to update)\n")
} else cat("no tracked copy to compare against; re-run with --write\n")
