#!/usr/bin/env Rscript
# ============================================================================
# shared/figure_inputs.R — staleness guard for committed figure-input snapshots
# ============================================================================
# Every figure in this paper is drawn from a *snapshot* of an analysis output,
# not from the analysis output itself: transcriptomic/data/figure_inputs/ for
# Fig. 2, genetics/results/figures/r_panels/ for Fig. 4. Snapshots are what make
# a figure reproducible from a clone, but they rot silently — rerun the upstream
# analysis and the figure keeps rendering the old numbers with no error. That is
# not hypothetical: it is how Fig. 2 came to show "76% concordant" while the
# manuscript said 72%.
#
# The fix is a MANIFEST.tsv beside each snapshot set, recording for every file
# the canonical source it came from and that source's checksum at capture time.
# Before drawing, the figure script re-checksums the sources and refuses to run
# if any has moved on.
#
# Required MANIFEST.tsv columns: file, source, source_md5, refreshed_at.
#
# Usage from a figure script:
#     source("../shared/figure_inputs.R")          # path relative to component
#     fi_check("data/figure_inputs")               # verify the whole set, or
#     p <- fi_path("de_results_subclass.csv", "data/figure_inputs")
#
# When a canonical source is absent (fresh clone, collaborator machine) its
# check is skipped and the snapshot is trusted — self-contained reproduction
# must keep working without the upstream repos.
# ============================================================================

fi_manifest <- function(dir) {
  p <- file.path(dir, "MANIFEST.tsv")
  if (!file.exists(p)) return(NULL)
  utils::read.delim(p, stringsAsFactors = FALSE)
}

# Must match figure_inputs.py::fingerprint. Two kinds exist: an md5, and a cheap
# "size-mtime:<bytes>:<epoch>" tag used for very large or bulk inputs, where
# hashing gigabytes on every render buys nothing (any real rewrite moves both).
#
# Which kind to compute is decided by the RECORDED value's format, never by the
# file's current size. Deciding by size means that the moment a writer and a
# reader disagree about the threshold — or a file grows past it — every check
# fails with a spurious "source has changed", which is worse than no check at
# all because it trains you to ignore it.
FI_HASH_LIMIT <- 256 * 1024 * 1024

fi_cheap_fp <- function(path)
  sprintf("size-mtime:%.0f:%.0f", file.size(path), trunc(as.numeric(file.mtime(path))))

fi_fingerprint <- function(path, like = NULL) {
  if (!file.exists(path)) return("")
  cheap <- if (!is.null(like) && nzchar(like)) startsWith(like, "size-mtime:")
           else file.size(path) > FI_HASH_LIMIT
  if (cheap) fi_cheap_fp(path) else unname(tools::md5sum(path))
}

# Returns a character vector of human-readable problems (empty if all good).
fi_problems <- function(dir) {
  man <- fi_manifest(dir)
  if (is.null(man)) return(character(0))
  out <- character(0)
  for (i in seq_len(nrow(man))) {
    row  <- man[i, ]
    snap <- file.path(dir, row$file)
    if (!file.exists(snap)) {
      out <- c(out, sprintf("  %s — snapshot file is missing", row$file)); next
    }
    src <- path.expand(row$source)
    if (!nzchar(row$source_md5) || !file.exists(src)) next   # unreachable source: skip
    if (!identical(fi_fingerprint(src, like = row$source_md5), row$source_md5))
      out <- c(out, sprintf("  %s — source has changed since %s\n      source: %s",
                            row$file, row$refreshed_at, row$source))
  }
  out
}

# Verify an entire snapshot directory. strict = TRUE stops; FALSE warns.
fi_check <- function(dir, strict = TRUE, refresh_cmd = NULL) {
  probs <- fi_problems(dir)
  if (!length(probs)) return(invisible(TRUE))
  if (is.null(refresh_cmd)) refresh_cmd <- "the component's refresh script"
  msg <- paste0(
    sprintf("STALE FIGURE INPUTS in %s\n", dir),
    paste(probs, collapse = "\n"), "\n",
    "  The figure would render numbers that no longer match the analysis.\n",
    sprintf("  Fix: %s", refresh_cmd))
  if (isTRUE(strict)) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  invisible(FALSE)
}

# Resolve one snapshot file, checking just that entry.
fi_path <- function(name, dir, strict = TRUE, refresh_cmd = NULL) {
  p <- file.path(dir, name)
  if (!file.exists(p))
    stop(sprintf("Missing figure input '%s' in %s.\n  Fix: %s", name, dir,
                 if (is.null(refresh_cmd)) "run the component's refresh script" else refresh_cmd),
         call. = FALSE)
  man <- fi_manifest(dir)
  if (!is.null(man) && name %in% man$file) {
    row <- man[man$file == name, ][1, ]
    src <- path.expand(row$source)
    if (nzchar(row$source_md5) && file.exists(src) &&
        !identical(fi_fingerprint(src, like = row$source_md5), row$source_md5)) {
      msg <- sprintf(paste0(
        "STALE FIGURE INPUT: '%s'\n",
        "  snapshot was taken from : %s\n",
        "  that source has CHANGED since %s\n",
        "  the figure would render numbers that no longer match the analysis.\n",
        "  Fix: %s"),
        name, row$source, row$refreshed_at,
        if (is.null(refresh_cmd)) "run the component's refresh script" else refresh_cmd)
      if (isTRUE(strict)) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
    }
  }
  p
}
