#!/usr/bin/env Rscript
# ============================================================================
# _xenium_exclusions.R — cell types omitted from Xenium figures as too rare
# ============================================================================
# Some SEA-AD types are recovered in numbers too small to support any per-donor
# estimate in this Xenium dataset. Showing them anyway implies a measurement the
# data cannot support, so they are dropped from every Xenium figure.
#
# Lamp5 Lhx6: 16 cells in the 742,103-cell cortical analysis set, spread over
# 11 of 24 donors (1-3 each) — 37x rarer than the next-rarest subclass, Sst
# Chodl (586). Only 0.7% of its quality-passing calls fall in cortex at all
# (2,549 of 2,598 are white matter), the lowest of any subclass, and those cells
# match their centroid poorly (median r 0.502 vs 0.794 for cortical cells). In
# this dataset the label absorbs poorly matched white-matter cells rather than
# marking the cortical interneuron type it names.
#
# This is a cell-recovery limit, not a panel limit — the 300-gene panel resolves
# Lamp5 Lhx6 at F1 0.98 in the snRNA-seq benchmark, which is exactly why it must
# be dropped explicitly: the resolvability panel would otherwise show it near the
# top of the figure as though it were a well-measured type.
#
# The disease analyses exclude it through their own thresholds (crumblr's
# >=50%-of-donors presence filter; DE's >=10-cells-per-donor pseudobulk floor),
# so this file only governs what is *plotted*. Documented in Supplementary
# Methods SM1, "Cell types too rare to analyse in the Xenium data".
# ============================================================================

XENIUM_TOO_RARE_SUBCLASS  <- c("Lamp5 Lhx6")
XENIUM_TOO_RARE_SUPERTYPE <- c("Lamp5_Lhx6_1")

# Drop the excluded types from a data frame, announcing what went so the removal
# is never silent.
drop_too_rare <- function(df, col = "cell_type", level = c("subclass", "supertype")) {
  level <- match.arg(level)
  drop <- if (level == "subclass") XENIUM_TOO_RARE_SUBCLASS else XENIUM_TOO_RARE_SUPERTYPE
  hit <- df[[col]] %in% drop
  if (any(hit))
    message(sprintf("  excluded as too rare in Xenium (SM1): %s",
                    paste(unique(df[[col]][hit]), collapse = ", ")))
  df[!hit, , drop = FALSE]
}
