# Purpose: shared bulk egene enrichment.

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(scales)
  library(sandwich)
  library(lmtest)
})

# M10
# Enrichment of ieGenes/cseGenes among bin-shared bulk eGenes

bins <- c("bin1", "bin2", "bin3", "bin4")

bin_labels <- c(
  bin1 = "Prenatal",
  bin2 = "Early postnatal",
  bin3 = "Adolescent",
  bin4 = "Adult"
)

outdir <- "M10_results"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

standardize_cell_names <- function(x) {
  cell_name_map <- c(
    adip  = "adipocyte",
    cap   = "capillary_endothelial_cell",
    fap   = "fibro_adipogenic_progenitor_cell",
    glial = "peripheral_glial",
    lec   = "lymphatic_endothelial_cell",
    mac   = "macrophage",
    mf1   = "type_i_myonuclei",
    mf2ab = "type_ii_a_b_myonuclei",
    mf2x  = "type_ii_x_myonuclei",
    msc   = "muscle_stem_cell",
    peri  = "pericyte",
    tc    = "t_cell",
    teno  = "tenocyte"
  )

  x <- as.character(x)
  mapped <- unname(cell_name_map[x])
  mapped[is.na(mapped)] <- x[is.na(mapped)]
  mapped
}

to_double_array <- function(x) {
  out <- array(
    as.numeric(x),
    dim = dim(x),
    dimnames = dimnames(x)
  )
  storage.mode(out) <- "double"
  out
}

# 1. Build bulk-eGene sharing classes
bulk_hits <- as.data.frame(fread(cfg$bulk_hits), stringsAsFactors = FALSE)

required_bulk_cols <- c("pheno_id", "group")
missing_bulk_cols <- setdiff(required_bulk_cols, colnames(bulk_hits))
if (length(missing_bulk_cols) > 0) {
  stop("bulk_hits is missing required columns: ", paste(missing_bulk_cols, collapse = ", "))
}

bulk_gene_bin <- bulk_hits %>%
  transmute(
    gene_id = as.character(pheno_id),
    bin = as.character(group)
  ) %>%
  filter(!is.na(gene_id), bin %in% bins) %>%
  distinct(gene_id, bin)

bulk_sharing <- bulk_gene_bin %>%
  count(gene_id, name = "n_sig_bins") %>%
  mutate(
    bulk_sharing_class = if_else(
      n_sig_bins >= 2L,
      "bin-shared",
      "bin-specific"
    )
  )

bulk_gene_bin <- bulk_gene_bin %>%
  left_join(bulk_sharing, by = "gene_id")

# 2. Restrict ieQTL/cseQTL discoveries to eligible bin-cell pairs
required_frac_cols <- c("bin", "cell", "eligible")
missing_frac_cols <- setdiff(required_frac_cols, colnames(frac_summary))
if (length(missing_frac_cols) > 0) {
  stop("frac_summary is missing required columns: ", paste(missing_frac_cols, collapse = ", "))
}

eligible_pairs <- frac_summary %>%
  transmute(
    bin = as.character(bin),
    cell = standardize_cell_names(cell),
    eligible = as.logical(eligible)
  ) %>%
  filter(bin %in% bins, eligible) %>%
  distinct(bin, cell)

if (nrow(eligible_pairs) == 0L) {
  stop("No eligible bin-cell combinations remain after filtering.")
}

required_hit_cols <- c("gene_id", "bin", "cell")
missing_ie_cols <- setdiff(required_hit_cols, colnames(ie_hits))
missing_cse_cols <- setdiff(required_hit_cols, colnames(cse_hits))

if (length(missing_ie_cols) > 0) {
  stop("ie_hits is missing required columns: ", paste(missing_ie_cols, collapse = ", "))
}
if (length(missing_cse_cols) > 0) {
  stop("cse_hits is missing required columns: ", paste(missing_cse_cols, collapse = ", "))
}

ie_gene_bin <- ie_hits %>%
  transmute(
    gene_id = as.character(gene_id),
    bin = as.character(bin),
    cell = standardize_cell_names(cell)
  ) %>%
  filter(!is.na(gene_id), bin %in% bins) %>%
  semi_join(eligible_pairs, by = c("bin", "cell")) %>%
  distinct(gene_id, bin) %>%
  mutate(ieGene = 1L)

cse_gene_bin <- cse_hits %>%
  transmute(
    gene_id = as.character(gene_id),
    bin = as.character(bin),
    cell = standardize_cell_names(cell)
  ) %>%
  filter(!is.na(gene_id), bin %in% bins) %>%
  semi_join(eligible_pairs, by = c("bin", "cell")) %>%
  distinct(gene_id, bin) %>%
  mutate(cseGene = 1L)

# 3. Build analysis table
m10_df <- bulk_gene_bin %>%
  left_join(ie_gene_bin, by = c("gene_id", "bin")) %>%
  left_join(cse_gene_bin, by = c("gene_id", "bin")) %>%
  mutate(
    ieGene = replace_na(ieGene, 0L),
    cseGene = replace_na(cseGene, 0L),
    bin = factor(bin, levels = bins),
    bulk_sharing_class = factor(
      bulk_sharing_class,
      levels = c("bin-specific", "bin-shared")
    )
  )

stopifnot(!anyDuplicated(m10_df[c("gene_id", "bin")]))

fwrite(m10_df, file.path(outdir, "M10_analysis_table.tsv"), sep = "\t")

# 4. Per-bin Fisher exact tests
run_fisher_by_bin <- function(df, outcome) {
  split(df, df$bin, drop = TRUE) %>%
    imap_dfr(function(dat, bin_name) {
      tab <- table(
        factor(dat$bulk_sharing_class, levels = c("bin-specific", "bin-shared")),
        factor(dat[[outcome]], levels = c(0, 1))
      )

      n_specific <- sum(tab["bin-specific", ])
      n_shared <- sum(tab["bin-shared", ])
      prop_specific <- if (n_specific > 0) tab["bin-specific", "1"] / n_specific else NA_real_
      prop_shared <- if (n_shared > 0) tab["bin-shared", "1"] / n_shared else NA_real_

      if (any(rowSums(tab) == 0L) || any(colSums(tab) == 0L)) {
        return(tibble(
          bin = bin_name,
          outcome = outcome,
          n_specific = n_specific,
          n_shared = n_shared,
          outcome_specific = unname(tab["bin-specific", "1"]),
          outcome_shared = unname(tab["bin-shared", "1"]),
          prop_specific = prop_specific,
          prop_shared = prop_shared,
          OR = NA_real_,
          conf_low = NA_real_,
          conf_high = NA_real_,
          p_value = NA_real_
        ))
      }

      ft <- fisher.test(tab)

      tibble(
        bin = bin_name,
        outcome = outcome,
        n_specific = n_specific,
        n_shared = n_shared,
        outcome_specific = unname(tab["bin-specific", "1"]),
        outcome_shared = unname(tab["bin-shared", "1"]),
        prop_specific = prop_specific,
        prop_shared = prop_shared,
        OR = unname(ft$estimate),
        conf_low = ft$conf.int[1],
        conf_high = ft$conf.int[2],
        p_value = ft$p.value
      )
    })
}

fisher_res <- bind_rows(
  run_fisher_by_bin(m10_df, "ieGene"),
  run_fisher_by_bin(m10_df, "cseGene")
) %>%
  group_by(outcome) %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  ungroup()

fwrite(fisher_res, file.path(outdir, "M10_fisher_by_bin.tsv"), sep = "\t")

# 5. CMH sensitivity analysis
run_cmh <- function(df, outcome) {
  form <- reformulate(c("bulk_sharing_class", outcome, "bin"))
  arr <- xtabs(form, data = df)
  arr <- to_double_array(arr)

  fit <- mantelhaen.test(arr, correct = FALSE)

  tibble(
    outcome = outcome,
    method = "CMH",
    common_OR = unname(fit$estimate),
    conf_low = fit$conf.int[1],
    conf_high = fit$conf.int[2],
    p_value = fit$p.value
  )
}

cmh_res <- bind_rows(
  run_cmh(m10_df, "ieGene"),
  run_cmh(m10_df, "cseGene")
)

fwrite(cmh_res, file.path(outdir, "M10_CMH_results.tsv"), sep = "\t")

# 6. Preferred pooled inference: bin-adjusted logistic regression
#    with gene-clustered robust standard errors
run_cluster_robust_logistic <- function(df, outcome) {
  dat <- df %>%
    mutate(
      outcome_value = as.integer(.data[[outcome]]),
      shared = as.integer(bulk_sharing_class == "bin-shared")
    )

  fit <- glm(
    outcome_value ~ shared + bin,
    family = binomial(),
    data = dat
  )

  robust_vcov <- sandwich::vcovCL(
    fit,
    cluster = dat$gene_id,
    type = "HC0"
  )

  robust_test <- lmtest::coeftest(fit, vcov. = robust_vcov)

  beta <- unname(robust_test["shared", "Estimate"])
  se <- unname(robust_test["shared", "Std. Error"])
  p <- unname(robust_test["shared", "Pr(>|z|)"])

  tibble(
    outcome = outcome,
    method = "Cluster-robust logistic regression",
    common_OR = exp(beta),
    conf_low = exp(beta - 1.96 * se),
    conf_high = exp(beta + 1.96 * se),
    p_value = p
  )
}

robust_res <- bind_rows(
  run_cluster_robust_logistic(m10_df, "ieGene"),
  run_cluster_robust_logistic(m10_df, "cseGene")
)

pooled_res <- bind_rows(robust_res, cmh_res) %>%
  mutate(
    outcome_label = recode(
      outcome,
      ieGene = "ieGenes",
      cseGene = "cseGenes"
    )
  )

fwrite(pooled_res, file.path(outdir, "M10_pooled_results.tsv"), sep = "\t")

# 7. Publication-ready forest plot
plot_bin <- fisher_res %>%
  transmute(
    outcome,
    outcome_label = recode(outcome, ieGene = "ieGenes", cseGene = "cseGenes"),
    row_label = unname(bin_labels[as.character(bin)]),
    OR,
    conf_low,
    conf_high,
    p_value,
    fdr,
    estimate_type = "Per-bin Fisher"
  )

plot_pooled <- robust_res %>%
  transmute(
    outcome,
    outcome_label = recode(outcome, ieGene = "ieGenes", cseGene = "cseGenes"),
    row_label = "Overall",
    OR = common_OR,
    conf_low,
    conf_high,
    p_value,
    fdr = NA_real_,
    estimate_type = "Bin-adjusted pooled"
  )

plot_df <- bind_rows(plot_bin, plot_pooled) %>%
  mutate(
    row_label = factor(
      row_label,
      levels = rev(c("Prenatal", "Early postnatal", "Adolescent", "Adult", "Overall"))
    ),
    is_overall = estimate_type == "Bin-adjusted pooled",
    sig_label = case_when(
      is.na(p_value) ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

label_col <- c("Adolescent"="#35B779", "Prenatal"="#440154", "Early postnatal"="#31688E", "Adult"="#FDE725", "Overall"="black")

p_forest <- ggplot(plot_df, aes(x = OR, y = row_label, color = row_label)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.45) +
  geom_errorbarh(
    aes(xmin = conf_low, xmax = conf_high, linewidth = is_overall),
    height = 0.16,
    na.rm = TRUE
  ) +
  geom_point(aes(size = is_overall), na.rm = TRUE) +
  geom_text(aes(x = conf_high, label = sig_label), hjust = -0.45, size = 3.4, na.rm = TRUE) +
  facet_wrap(~outcome_label, nrow = 1) +
  scale_x_log10() +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 16), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 1.6, `TRUE` = 1.6), guide = "none") +
  scale_linewidth_manual(values = c(`FALSE` = 0.55, `TRUE` = 0.55), guide = "none") +
  scale_color_manual(values = label_col) +
  labs(
    x = "Odds ratio for bin-shared versus bin-specific bulk eGenes",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(colour = "black"),
    panel.spacing.x = grid::unit(1.1, "lines"),
    plot.margin = margin(5.5, 18, 5.5, 5.5)
  )

print(p_forest)

ggsave(
  file.path(outdir, "M10_forest_publication.pdf"),
  p_forest,
  width = 7.2,
  height = 4.6,
  useDingbats = FALSE
)

ggsave(
  file.path(outdir, "M10_forest_publication.png"),
  p_forest,
  width = 7.2,
  height = 4.6,
  dpi = 600
)

# 8. Companion proportion plot
prop_df <- fisher_res %>%
  select(bin, outcome, prop_specific, prop_shared) %>%
  pivot_longer(
    cols = c(prop_specific, prop_shared),
    names_to = "bulk_class",
    values_to = "proportion"
  ) %>%
  mutate(
    bulk_class = recode(
      bulk_class,
      prop_specific = "Bin-specific",
      prop_shared = "Bin-shared"
    ),
    outcome = recode(outcome, ieGene = "ieGenes", cseGene = "cseGenes"),
    bin = factor(bin, levels = bins)
  )

message(
  "M10 analysis completed. Results written to: ",
  normalizePath(outdir, mustWork = FALSE)
)
