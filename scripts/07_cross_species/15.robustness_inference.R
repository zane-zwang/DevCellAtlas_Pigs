# Purpose: robustness inference.

# Input: pairwise human–pig and human–mouse correlation tables and summary objects
# created by 16.similarity.R. Output: balanced and bootstrap uncertainty tables.
# Delta is median human–pig correlation minus median human–mouse correlation.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

set.seed(12345)

# Parameters
N_BALANCE <- 1000
N_BOOT    <- 1000

out_dir <- "./human_centered_robustness"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Check required objects
required_objects <- c(
  "celltype_cor_pairwise",
  "celltype_cor_sample",
  "identity_stats",
  "human_spec_cor_pairwise",
  "human_spec_cor_sample",
  "human_spec_stats"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0) {
  stop(
    "Missing required object(s): ",
    paste(missing_objects, collapse = ", "),
    "\nRun 16.similarity.R before this analysis."
  )
}

# Helpers

detect_comparator_col <- function(df) {

  if (!"human_sample" %in% colnames(df)) {
    stop("Pairwise table must contain `human_sample`.")
  }

  sample_cols <- grep(
    "sample",
    colnames(df),
    value = TRUE,
    ignore.case = TRUE
  )

  candidates <- setdiff(sample_cols, "human_sample")

  preferred <- c(
    "other_sample",
    "comparator_sample",
    "target_sample",
    "nonhuman_sample",
    "species_sample",
    "sample2",
    "pig_mouse_sample",
    "comparison_sample"
  )

  hit <- preferred[preferred %in% candidates]

  if (length(hit) >= 1) {
    return(hit[1])
  }

  if (length(candidates) == 1) {
    return(candidates)
  }

  stop(
    "Could not uniquely detect comparator-sample column. Candidates: ",
    paste(candidates, collapse = ", ")
  )
}


standardize_pairwise <- function(df) {

  comparator_col <- detect_comparator_col(df)

  df %>%
    mutate(
      comparator_sample = .data[[comparator_col]]
    ) %>%
    select(
      celltype,
      comparison,
      human_sample,
      comparator_sample,
      cor,
      everything()
    )
}


make_cor_matrices <- function(df_ct) {

  hp <- df_ct %>%
    filter(comparison == "human_pig") %>%
    select(human_sample, comparator_sample, cor) %>%
    distinct()

  hm <- df_ct %>%
    filter(comparison == "human_mouse") %>%
    select(human_sample, comparator_sample, cor) %>%
    distinct()

  human_use <- intersect(
    unique(hp$human_sample),
    unique(hm$human_sample)
  )

  pig_use <- unique(hp$comparator_sample)
  mouse_use <- unique(hm$comparator_sample)

  hp_mat <- matrix(
    NA_real_,
    nrow = length(human_use),
    ncol = length(pig_use),
    dimnames = list(human_use, pig_use)
  )

  hm_mat <- matrix(
    NA_real_,
    nrow = length(human_use),
    ncol = length(mouse_use),
    dimnames = list(human_use, mouse_use)
  )

  hp_i <- match(hp$human_sample, human_use)
  hp_j <- match(hp$comparator_sample, pig_use)
  ok_hp <- !is.na(hp_i) & !is.na(hp_j)
  hp_mat[cbind(hp_i[ok_hp], hp_j[ok_hp])] <- hp$cor[ok_hp]

  hm_i <- match(hm$human_sample, human_use)
  hm_j <- match(hm$comparator_sample, mouse_use)
  ok_hm <- !is.na(hm_i) & !is.na(hm_j)
  hm_mat[cbind(hm_i[ok_hm], hm_j[ok_hm])] <- hm$cor[ok_hm]

  list(
    hp = hp_mat,
    hm = hm_mat,
    n_human = length(human_use),
    n_pig = length(pig_use),
    n_mouse = length(mouse_use)
  )
}


row_median_na <- function(x) {
  apply(x, 1, median, na.rm = TRUE)
}


observed_from_mats <- function(mats) {

  hp <- row_median_na(mats$hp)
  hm <- row_median_na(mats$hm)

  delta <- hp - hm
  complete <- is.finite(delta)

  tibble(
    n_human = sum(complete),
    n_pig = mats$n_pig,
    n_mouse = mats$n_mouse,
    median_HP = median(hp[complete], na.rm = TRUE),
    median_HM = median(hm[complete], na.rm = TRUE),
    median_delta = median(delta[complete], na.rm = TRUE),
    prop_pig_closer = mean(delta[complete] > 0, na.rm = TRUE)
  )
}


leave_one_out_from_mats <- function(mats) {

  hp0 <- mats$hp
  hm0 <- mats$hm

  out <- list()
  k <- 1

  # remove one pig comparator
  if (ncol(hp0) > 1) {
    for (j in seq_len(ncol(hp0))) {

      hp <- row_median_na(hp0[, -j, drop = FALSE])
      hm <- row_median_na(hm0)

      delta <- hp - hm
      ok <- is.finite(delta)

      out[[k]] <- tibble(
        omitted_comparison = "human_pig",
        omitted_sample = colnames(hp0)[j],
        n_complete = sum(ok),
        median_delta = median(delta[ok], na.rm = TRUE),
        prop_pig_closer = mean(delta[ok] > 0, na.rm = TRUE)
      )
      k <- k + 1
    }
  }

  # remove one mouse comparator
  if (ncol(hm0) > 1) {
    for (j in seq_len(ncol(hm0))) {

      hp <- row_median_na(hp0)
      hm <- row_median_na(hm0[, -j, drop = FALSE])

      delta <- hp - hm
      ok <- is.finite(delta)

      out[[k]] <- tibble(
        omitted_comparison = "human_mouse",
        omitted_sample = colnames(hm0)[j],
        n_complete = sum(ok),
        median_delta = median(delta[ok], na.rm = TRUE),
        prop_pig_closer = mean(delta[ok] > 0, na.rm = TRUE)
      )
      k <- k + 1
    }
  }

  bind_rows(out)
}


balance_from_mats <- function(mats, n_iter = 5000) {

  hp0 <- mats$hp
  hm0 <- mats$hm

  n_target <- min(ncol(hp0), ncol(hm0))

  if (n_target < 2) {
    return(
      tibble(
        iteration = seq_len(n_iter),
        n_target = n_target,
        delta = NA_real_,
        prop_pig_closer = NA_real_
      )
    )
  }

  map_dfr(
    seq_len(n_iter),
    function(i) {

      hp_idx <- sample(seq_len(ncol(hp0)), n_target, replace = FALSE)
      hm_idx <- sample(seq_len(ncol(hm0)), n_target, replace = FALSE)

      hp <- row_median_na(hp0[, hp_idx, drop = FALSE])
      hm <- row_median_na(hm0[, hm_idx, drop = FALSE])

      delta_h <- hp - hm
      ok <- is.finite(delta_h)

      tibble(
        iteration = i,
        n_target = n_target,
        delta = median(delta_h[ok], na.rm = TRUE),
        prop_pig_closer = mean(delta_h[ok] > 0, na.rm = TRUE)
      )
    }
  )
}


bootstrap_from_mats <- function(mats, n_iter = 5000) {

  hp0 <- mats$hp
  hm0 <- mats$hm

  nh <- nrow(hp0)
  np <- ncol(hp0)
  nm <- ncol(hm0)

  map_dfr(
    seq_len(n_iter),
    function(i) {

      # Same resampled human rows are used for HP and HM.
      h_idx <- sample(seq_len(nh), nh, replace = TRUE)
      p_idx <- sample(seq_len(np), np, replace = TRUE)
      m_idx <- sample(seq_len(nm), nm, replace = TRUE)

      hp_boot <- hp0[h_idx, p_idx, drop = FALSE]
      hm_boot <- hm0[h_idx, m_idx, drop = FALSE]

      hp <- row_median_na(hp_boot)
      hm <- row_median_na(hm_boot)

      delta_h <- hp - hm
      ok <- is.finite(delta_h)

      tibble(
        iteration = i,
        delta = median(delta_h[ok], na.rm = TRUE),
        prop_pig_closer = mean(delta_h[ok] > 0, na.rm = TRUE)
      )
    }
  )
}


summarise_distribution <- function(df) {

  x <- df$delta

  tibble(
    n_iteration = sum(is.finite(x)),
    median_delta = median(x, na.rm = TRUE),
    ci_low = unname(quantile(x, 0.025, na.rm = TRUE)),
    ci_high = unname(quantile(x, 0.975, na.rm = TRUE)),
    prop_delta_gt0 = mean(x > 0, na.rm = TRUE),
    median_prop_pig_closer =
      median(df$prop_pig_closer, na.rm = TRUE)
  )
}


run_robustness <- function(pairwise_df,
                           analysis_name,
                           n_balance = N_BALANCE,
                           n_boot = N_BOOT) {

  message("===== ", analysis_name, " =====")

  pair_df <- standardize_pairwise(pairwise_df)

  celltypes <- unique(pair_df$celltype)

  observed_all <- list()
  loo_all <- list()
  balance_all <- list()
  boot_all <- list()

  for (ct in celltypes) {

    message("  ", ct)

    df_ct <- pair_df %>%
      filter(celltype == ct)

    mats <- make_cor_matrices(df_ct)

    observed_all[[ct]] <- observed_from_mats(mats) %>%
      mutate(celltype = ct, .before = 1)

    loo_all[[ct]] <- leave_one_out_from_mats(mats) %>%
      mutate(celltype = ct, .before = 1)

    balance_all[[ct]] <- balance_from_mats(
      mats,
      n_iter = n_balance
    ) %>%
      mutate(celltype = ct, .before = 1)

    boot_all[[ct]] <- bootstrap_from_mats(
      mats,
      n_iter = n_boot
    ) %>%
      mutate(celltype = ct, .before = 1)
  }

  observed <- bind_rows(observed_all)
  loo <- bind_rows(loo_all)
  balance <- bind_rows(balance_all)
  boot <- bind_rows(boot_all)

  loo_summary <- loo %>%
    group_by(celltype) %>%
    summarise(
      min_median_delta = min(median_delta, na.rm = TRUE),
      max_median_delta = max(median_delta, na.rm = TRUE),
      min_prop_pig_closer = min(prop_pig_closer, na.rm = TRUE),
      max_prop_pig_closer = max(prop_pig_closer, na.rm = TRUE),
      any_sign_flip = any(median_delta < 0, na.rm = TRUE),
      .groups = "drop"
    )

  balance_summary <- balance %>%
    group_by(celltype) %>%
    group_modify(~summarise_distribution(.x)) %>%
    ungroup() %>%
    left_join(
      balance %>%
        group_by(celltype) %>%
        summarise(
          n_target = first(n_target),
          .groups = "drop"
        ),
      by = "celltype"
    )

  boot_summary <- boot %>%
    group_by(celltype) %>%
    group_modify(~summarise_distribution(.x)) %>%
    ungroup()

  final_summary <- observed %>%
    left_join(
      loo_summary,
      by = "celltype"
    ) %>%
    left_join(
      balance_summary %>%
        transmute(
          celltype,
          balance_n_target = n_target,
          balance_median_delta = median_delta,
          balance_ci_low = ci_low,
          balance_ci_high = ci_high,
          balance_prop_delta_gt0 = prop_delta_gt0,
          balance_median_prop_pig_closer =
            median_prop_pig_closer
        ),
      by = "celltype"
    ) %>%
    left_join(
      boot_summary %>%
        transmute(
          celltype,
          bootstrap_median_delta = median_delta,
          bootstrap_ci_low = ci_low,
          bootstrap_ci_high = ci_high,
          bootstrap_prop_delta_gt0 = prop_delta_gt0,
          bootstrap_median_prop_pig_closer =
            median_prop_pig_closer
        ),
      by = "celltype"
    ) %>%
    mutate(
      observed_delta_positive = median_delta > 0,
      loo_no_sign_flip = !any_sign_flip,
      balance_ci_above_zero = balance_ci_low > 0,
      bootstrap_ci_above_zero = bootstrap_ci_low > 0
    )

  prefix <- gsub("[^A-Za-z0-9]+", "_", analysis_name)

  write.csv(
    observed,
    file.path(out_dir, paste0(prefix, "_01_observed.csv")),
    row.names = FALSE
  )

  write.csv(
    loo,
    file.path(out_dir, paste0(prefix, "_02_leave_one_out_iterations.csv")),
    row.names = FALSE
  )

  write.csv(
    loo_summary,
    file.path(out_dir, paste0(prefix, "_03_leave_one_out_summary.csv")),
    row.names = FALSE
  )

  write.csv(
    balance,
    file.path(out_dir, paste0(prefix, "_04_balance_iterations.csv")),
    row.names = FALSE
  )

  write.csv(
    balance_summary,
    file.path(out_dir, paste0(prefix, "_05_balance_summary.csv")),
    row.names = FALSE
  )

  write.csv(
    boot,
    file.path(out_dir, paste0(prefix, "_06_bootstrap_iterations.csv")),
    row.names = FALSE
  )

  write.csv(
    boot_summary,
    file.path(out_dir, paste0(prefix, "_07_bootstrap_summary.csv")),
    row.names = FALSE
  )

  write.csv(
    final_summary,
    file.path(out_dir, paste0(prefix, "_08_final_summary.csv")),
    row.names = FALSE
  )

  list(
    observed = observed,
    leave_one_out = loo,
    leave_one_out_summary = loo_summary,
    balance = balance,
    balance_summary = balance_summary,
    bootstrap = boot,
    bootstrap_summary = boot_summary,
    final_summary = final_summary
  )
}

# 1. Cell-type identity conservation

identity_robustness <- run_robustness(
  pairwise_df = celltype_cor_pairwise,
  analysis_name = "identity",
  n_balance = N_BALANCE,
  n_boot = N_BOOT
)

identity_final <- identity_stats %>%
  rename(
    source_summary_n_human_sample = n_human_sample,
    source_summary_median_cor_HP = median_cor_HP,
    source_summary_median_cor_HM = median_cor_HM,
    source_summary_delta_HP_minus_HM = delta_HP_minus_HM,
    source_summary_prop_pig_closer = prop_pig_closer,
    source_summary_wilcox_p = p_value,
    source_summary_wilcox_FDR = p_adj
  ) %>%
  left_join(
    identity_robustness$final_summary,
    by = "celltype"
  )

write.csv(
  identity_final,
  file.path(out_dir, "identity_combined_summary.csv"),
  row.names = FALSE
)

# 2. Human-specific program retention

human_spec_robustness <- run_robustness(
  pairwise_df = human_spec_cor_pairwise,
  analysis_name = "human_specific_program",
  n_balance = N_BALANCE,
  n_boot = N_BOOT
)

human_spec_final <- human_spec_stats %>%
  rename(
    source_summary_n_human_sample = n_human_sample,
    source_summary_median_cor_HP = median_cor_HP,
    source_summary_median_cor_HM = median_cor_HM,
    source_summary_delta_HP_minus_HM = delta_HP_minus_HM,
    source_summary_prop_pig_closer = prop_pig_closer,
    source_summary_wilcox_p = p_value,
    source_summary_wilcox_FDR = p_adj
  ) %>%
  left_join(
    human_spec_robustness$final_summary,
    by = "celltype"
  )

write.csv(
  human_spec_final,
  file.path(out_dir, "human_specific_program_combined_summary.csv"),
  row.names = FALSE
)

# 3. Minimal manuscript-facing tables
# These deliberately emphasize effect sizes and bootstrap uncertainty.
# Paired Wilcoxon results are descriptive; balancing and bootstrap results
# provide the robustness estimates.

identity_manuscript_table <- identity_final %>%
  transmute(
    celltype,
    n_human = n_human,
    n_pig = n_pig,
    n_mouse = n_mouse,
    median_cor_human_pig = median_HP,
    median_cor_human_mouse = median_HM,
    observed_delta = median_delta,
    observed_prop_pig_closer = prop_pig_closer,
    leave_one_out_sign_flip = any_sign_flip,
    balanced_delta = balance_median_delta,
    balanced_ci_low = balance_ci_low,
    balanced_ci_high = balance_ci_high,
    balanced_prop_delta_gt0 = balance_prop_delta_gt0,
    bootstrap_delta = bootstrap_median_delta,
    bootstrap_ci_low = bootstrap_ci_low,
    bootstrap_ci_high = bootstrap_ci_high,
    bootstrap_prop_delta_gt0 = bootstrap_prop_delta_gt0
  )

write.csv(
  identity_manuscript_table,
  file.path(out_dir, "identity_manuscript_table.csv"),
  row.names = FALSE
)

human_spec_manuscript_table <- human_spec_final %>%
  transmute(
    celltype,
    n_human = n_human,
    n_pig = n_pig,
    n_mouse = n_mouse,
    median_cor_human_pig = median_HP,
    median_cor_human_mouse = median_HM,
    observed_delta = median_delta,
    observed_prop_pig_closer = prop_pig_closer,
    leave_one_out_sign_flip = any_sign_flip,
    balanced_delta = balance_median_delta,
    balanced_ci_low = balance_ci_low,
    balanced_ci_high = balance_ci_high,
    balanced_prop_delta_gt0 = balance_prop_delta_gt0,
    bootstrap_delta = bootstrap_median_delta,
    bootstrap_ci_low = bootstrap_ci_low,
    bootstrap_ci_high = bootstrap_ci_high,
    bootstrap_prop_delta_gt0 = bootstrap_prop_delta_gt0
  )

write.csv(
  human_spec_manuscript_table,
  file.path(out_dir, "human_specific_program_manuscript_table.csv"),
  row.names = FALSE
)

message("Human-centered robustness analysis complete.")
message("Output directory: ", normalizePath(out_dir))

message("\nPrimary files:")
message("  identity_combined_summary.csv")
message("  identity_manuscript_table.csv")
message("  human_specific_program_combined_summary.csv")
message("  human_specific_program_manuscript_table.csv")

