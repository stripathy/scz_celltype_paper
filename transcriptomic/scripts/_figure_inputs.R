#!/usr/bin/env Rscript
# ============================================================================
# _figure_inputs.R — Figure 2's binding to the repo-wide snapshot guard
# ============================================================================
# The mechanism lives in ../shared/figure_inputs.R and is shared with Figure 4;
# this file only fixes Figure 2's snapshot directory and refresh command so the
# scripts can call fig_input("x.csv") without repeating either.
#
# Refresh with:  Rscript scripts/00_refresh_figure_inputs.R
# ============================================================================

source("../shared/figure_inputs.R")

FIG_INPUTS   <- "data/figure_inputs"
FIG_MANIFEST <- file.path(FIG_INPUTS, "MANIFEST.tsv")
FIG_REFRESH  <- "Rscript scripts/00_refresh_figure_inputs.R"

fig_input <- function(name, strict = TRUE)
  fi_path(name, FIG_INPUTS, strict = strict, refresh_cmd = FIG_REFRESH)

# verify the whole set at once (used by the composite before it draws anything)
fig_check <- function(strict = TRUE)
  fi_check(FIG_INPUTS, strict = strict, refresh_cmd = FIG_REFRESH)
