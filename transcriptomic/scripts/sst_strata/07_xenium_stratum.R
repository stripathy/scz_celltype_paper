#!/usr/bin/env Rscript
# Step 7 | Xenium cross-platform check at stratum level.
#
# PROVENANCE NOTE. The results text quotes Xenium stratum-level z values for SST
# and VGF and the panel-coverage counts for the graded modules. Those numbers
# previously lived in a CSV with no producer script in the repo (written ad hoc in
# an earlier session). This script reconstructs that file exactly -- the recipe
# was recovered by testing candidate aggregations against the orphaned CSV, and
# reproduces it to 0.0000 on both logFC and z:
#   per stratum, unweighted MEAN of the per-supertype Xenium logFC, with the
#   per-supertype z combined by Stouffer's method (sum(z)/sqrt(n)).
# Unweighted rather than inverse-variance because the 300-gene panel resolves the
# supertypes unequally; weighting by precision would let the best-resolved
# supertype dominate a stratum it only partly represents.
#
# Xenium is a whole-cell platform on an independent cohort (LIBD, 24 DLPFC
# sections), so agreement here is genuine cross-platform replication -- but only
# for the ~144 panel genes, which is why the graded modules cannot be tested.
source("transcriptomic/scripts/sst_strata/_common.R")

st_of <- stratum_of()
sig <- read_csv(file.path(P$pb, "stratum_gene_signatures.csv"), show_col_types = FALSE)

xen <- read_csv(P$spatial_de, show_col_types = FALSE) |>
  filter(celltype %in% names(st_of)) |>
  mutate(stratum = st_of[celltype],
         z_i = sign(logFC) * abs(qnorm(PValue / 2)))
xen_st <- xen |> group_by(stratum, gene) |>
  summarise(xen_logFC = mean(logFC), xen_z = sum(z_i) / sqrt(n()),
            n_supertypes = n(), .groups = "drop")

conc <- sig |> select(stratum, gene, sn_est = estimate, sn_se = se, sn_z = z) |>
  inner_join(xen_st, by = c("stratum", "gene"))
write_csv(conc, file.path(P$out, "xenium_stratum_concordance.csv"))
say("panel genes with stratum-level estimates on both platforms: %d (%d rows)",
    n_distinct(conc$gene), nrow(conc))

say("=== SST / VGF across strata, both platforms ===")
print(conc |> filter(gene %in% c("SST", "VGF")) |>
        transmute(gene, stratum, sn_z = round(sn_z, 1), xen_z = round(xen_z, 1)) |>
        arrange(gene, match(stratum, STRATA_LEVELS)))

say("=== direction agreement and correlation, per stratum ===")
print(conc |> group_by(stratum) |>
        summarise(n = n(), r = round(cor(sn_z, xen_z), 2),
                  pct_same_sign = round(100 * mean(sign(sn_z) == sign(xen_z))),
                  .groups = "drop") |> arrange(match(stratum, STRATA_LEVELS)))

# Module coverage: why the graded modules are untestable on this panel
gsea <- read_csv(file.path(P$pb, "gsea_all_signatures.csv"), show_col_types = FALSE)
mods <- modules(gsea)
panel <- unique(conc$gene)
say("=== module coverage on the Xenium panel ===")
for (k in names(mods))
  say("  %-12s %d of %d leading-edge genes on the panel", k,
      sum(mods[[k]] %in% panel), length(mods[[k]]))
