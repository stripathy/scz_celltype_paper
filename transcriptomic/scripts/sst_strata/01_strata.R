#!/usr/bin/env Rscript
# Step 1 | Define the three Sst depletion strata.
#
# depleted     crumblr estimate < 0 and FDR < 0.20   (Sst_2, Sst_3, Sst_20, Sst_22, Sst_25)
# intermediate estimate < 0, not significant
# non_depleted estimate >= 0
#
# The FDR < 0.20 cut is a deliberate co-author decision: it admits Sst_20, whose
# depletion is trend-level in the meta-analysis (FDR = 0.19) but which replicates
# the laminar and genetic pattern of the other four.
#
# Writes the single source of truth for stratum membership; every downstream step
# reads it rather than re-listing supertypes.
source("transcriptomic/scripts/sst_strata/_common.R")

crum <- read_csv(P$crumblr, show_col_types = FALSE) |>
  filter(str_starts(CellType, "Sst"), !str_detect(CellType, "Chodl")) |>
  mutate(stratum = case_when(estimate < 0 & padj < 0.20 ~ "depleted",
                             estimate < 0               ~ "intermediate",
                             TRUE                       ~ "non_depleted"))

depth <- read_csv(P$depth, show_col_types = FALSE) |>
  filter(subclass == "Sst") |>
  select(supertype, depth_xenium = median_Xenium, depth_merfish = median_MERFISH)

strata <- crum |> left_join(depth, by = c(CellType = "supertype")) |> arrange(estimate)
write_csv(strata, file.path(P$out, "strata_definition.csv"))

say("strata (n = %d supertypes):", nrow(strata))
print(strata |> select(CellType, estimate, se, padj, stratum, depth_xenium), n = 20)
print(strata |> group_by(stratum) |>
        summarise(n = n(), mean_depth = round(mean(depth_xenium, na.rm = TRUE), 2),
                  members = paste(CellType, collapse = ", ")))
ct <- cor.test(strata$estimate, strata$depth_xenium, method = "spearman")
say("Spearman(crumblr estimate, Xenium median depth) = %.2f, P = %.4f, n = %d",
    ct$estimate, ct$p.value, nrow(strata))
