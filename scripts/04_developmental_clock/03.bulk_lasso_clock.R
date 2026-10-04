# Purpose: bulk lasso clock.

suppressPackageStartupMessages({
  library(glmnet)
  library(dplyr)
  library(ggplot2)
})

# Bulk skeletal-muscle developmental clock
# Input:
#   pseudobulk_data/bulk_tpm.csv
#   pseudobulk_data/ages
#       one developmental age per sample, in the same order as
#       columns of bulk_tpm.csv; prenatal ages are negative,
#       postnatal ages are positive, birth = 0
# Model:
#   repeated age-stratified train/test split
#   training-only feature selection
#   LASSO regression (alpha = 1)
#   5-fold CV for lambda
#   stable clock genes = selected in >=5/10 repeats

set.seed(1234)

# 1. Parameters
n_repeats <- 10
train_prop <- 0.80

mean_tpm_cutoff <- 5
variance_quantile <- 0.25

spearman_cutoff <- 0.30
quadratic_p_cutoff <- 0.05

inner_folds <- 5
stable_min_repeats <- ceiling(n_repeats / 2)

out_dir <- "bulk_clock_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Read input data
tpm_gene_by_sample <- read.csv(
  "pseudobulk_data/bulk_tpm.csv",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

age_df <- read.table(
  "pseudobulk_data/ages",
  header = FALSE
)

ages <- as.numeric(age_df[[1]])

stopifnot(
  ncol(tpm_gene_by_sample) == length(ages),
  !anyNA(ages)
)

# samples x genes
x_all <- t(as.matrix(tpm_gene_by_sample))
storage.mode(x_all) <- "double"

if (is.null(rownames(x_all))) {
  rownames(x_all) <- paste0("sample_", seq_len(nrow(x_all)))
}

cat("Samples:", nrow(x_all), "\n")
cat("Genes:", ncol(x_all), "\n")
cat("Developmental age range:", paste(range(ages), collapse = " to "), "days\n")

# 3. Helper functions

# Age-stratified outer split.
# For ages with >=2 replicates, at least one sample is retained
# for testing and at least one for training.
make_age_stratified_split <- function(age, train_prop = 0.80) {
  train_idx <- integer(0)

  for (a in unique(age)) {
    ids <- which(age == a)
    n <- length(ids)

    if (n == 1) {
      # Cannot split a singleton age across training and test.
      train_idx <- c(train_idx, ids)
    } else {
      n_train <- round(train_prop * n)
      n_train <- max(1, min(n - 1, n_train))
      train_idx <- c(train_idx, sample(ids, n_train))
    }
  }

  train_idx <- sort(unique(train_idx))
  test_idx <- setdiff(seq_along(age), train_idx)

  if (length(test_idx) == 0) {
    stop("No samples were assigned to the outer test set.")
  }

  list(train = train_idx, test = test_idx)
}


# Training-only feature selection.
select_candidate_genes <- function(
  x_train,
  y_train,
  mean_tpm_cutoff = 5,
  variance_quantile = 0.25,
  spearman_cutoff = 0.30,
  quadratic_p_cutoff = 0.05
) {
  # 1) Mean-expression filter
  gene_mean <- colMeans(x_train, na.rm = TRUE)
  keep <- is.finite(gene_mean) & gene_mean > mean_tpm_cutoff
  x <- x_train[, keep, drop = FALSE]

  if (ncol(x) < 2) {
    stop("Too few genes remain after mean-expression filtering.")
  }

  # 2) Variance filter
  gene_var <- apply(x, 2, var, na.rm = TRUE)
  var_cut <- quantile(
    gene_var,
    probs = variance_quantile,
    na.rm = TRUE,
    names = FALSE
  )

  keep <- is.finite(gene_var) & gene_var > var_cut
  x <- x[, keep, drop = FALSE]

  if (ncol(x) < 2) {
    stop("Too few genes remain after variance filtering.")
  }

  # 3) Age-associated filter
  spearman_rho <- apply(
    x,
    2,
    function(v) {
      suppressWarnings(
        cor(v, y_train, method = "spearman", use = "complete.obs")
      )
    }
  )

  quadratic_p <- apply(
    x,
    2,
    function(v) {
      fit <- tryCatch(
        lm(v ~ y_train + I(y_train^2)),
        error = function(e) NULL
      )

      if (is.null(fit)) {
        return(NA_real_)
      }

      tab <- summary(fit)$coefficients

      if (!"I(y_train^2)" %in% rownames(tab)) {
        return(NA_real_)
      }

      tab["I(y_train^2)", "Pr(>|t|)"]
    }
  )

  keep_age <- (
    abs(spearman_rho) > spearman_cutoff |
      quadratic_p < quadratic_p_cutoff
  )
  keep_age[is.na(keep_age)] <- FALSE

  genes <- colnames(x)[keep_age]

  if (length(genes) < 2) {
    stop("Too few genes remain after age-associated filtering.")
  }

  genes
}


# Training-set scaling; the same parameters are applied to test data.
scale_train_test <- function(x_train, x_test) {
  center <- colMeans(x_train, na.rm = TRUE)
  scale_sd <- apply(x_train, 2, sd, na.rm = TRUE)

  valid <- is.finite(scale_sd) & scale_sd > 0

  x_train <- x_train[, valid, drop = FALSE]
  x_test <- x_test[, valid, drop = FALSE]
  center <- center[valid]
  scale_sd <- scale_sd[valid]

  x_train_scaled <- sweep(x_train, 2, center, "-")
  x_train_scaled <- sweep(x_train_scaled, 2, scale_sd, "/")

  x_test_scaled <- sweep(x_test, 2, center, "-")
  x_test_scaled <- sweep(x_test_scaled, 2, scale_sd, "/")

  list(
    train = x_train_scaled,
    test = x_test_scaled,
    center = center,
    scale = scale_sd
  )
}


# Balanced random folds for inner CV.
make_inner_folds <- function(n, k = 5) {
  k <- min(k, n)
  sample(rep(seq_len(k), length.out = n))
}


# 4. Repeated outer train/test evaluation
prediction_list <- vector("list", n_repeats)
performance_list <- vector("list", n_repeats)
clock_gene_list_bulk <- vector("list", n_repeats)
model_info_list <- vector("list", n_repeats)

for (rep in seq_len(n_repeats)) {

  cat("\n===== Repeat", rep, "of", n_repeats, "=====\n")
  set.seed(1234 + rep)

  # Outer split
  split <- make_age_stratified_split(
    age = ages,
    train_prop = train_prop
  )

  train_idx <- split$train
  test_idx <- split$test

  x_train_raw <- x_all[train_idx, , drop = FALSE]
  x_test_raw <- x_all[test_idx, , drop = FALSE]

  y_train <- ages[train_idx]
  y_test <- ages[test_idx]

  # Feature selection using training samples only
  candidate_genes <- select_candidate_genes(
    x_train = x_train_raw,
    y_train = y_train,
    mean_tpm_cutoff = mean_tpm_cutoff,
    variance_quantile = variance_quantile,
    spearman_cutoff = spearman_cutoff,
    quadratic_p_cutoff = quadratic_p_cutoff
  )

  cat("Training samples:", length(train_idx), "\n")
  cat("Test samples:", length(test_idx), "\n")
  cat("Selected candidate genes:", length(candidate_genes), "\n")

  x_train <- x_train_raw[, candidate_genes, drop = FALSE]
  x_test <- x_test_raw[, candidate_genes, drop = FALSE]

  # Scale using training-set parameters only
  scaled <- scale_train_test(x_train, x_test)

  x_train_scaled <- scaled$train
  x_test_scaled <- scaled$test

  stopifnot(
    nrow(x_train_scaled) == length(y_train),
    nrow(x_test_scaled) == length(y_test),
    all(is.finite(x_train_scaled)),
    all(is.finite(x_test_scaled))
  )

  # LASSO: alpha fixed at 1; lambda selected by 5-fold CV
  foldid <- make_inner_folds(
    n = nrow(x_train_scaled),
    k = inner_folds
  )

  cv_fit <- cv.glmnet(
    x = x_train_scaled,
    y = y_train,
    family = "gaussian",
    alpha = 1,
    foldid = foldid,
    type.measure = "mae",
    standardize = FALSE
  )

  best_lambda <- cv_fit$lambda.min

  # Fit final outer-training model at selected lambda
  final_fit <- glmnet(
    x = x_train_scaled,
    y = y_train,
    family = "gaussian",
    alpha = 1,
    lambda = best_lambda,
    standardize = FALSE
  )

  # Predict untouched outer test samples
  pred <- as.numeric(
    predict(
      final_fit,
      newx = x_test_scaled,
      s = best_lambda
    )
  )

  # Performance
  pearson_r <- suppressWarnings(
    cor(y_test, pred, method = "pearson")
  )

  spearman_r <- suppressWarnings(
    cor(y_test, pred, method = "spearman")
  )

  r2 <- pearson_r^2
  mae <- mean(abs(pred - y_test))

  prediction_list[[rep]] <- data.frame(
    Sample = rownames(x_all)[test_idx],
    SampleIndex = test_idx,
    ActualAge = y_test,
    PredictedAge = pred,
    Repeat = rep,
    stringsAsFactors = FALSE
  )

  performance_list[[rep]] <- data.frame(
    Repeat = rep,
    N_train = length(train_idx),
    N_test = length(test_idx),
    N_candidate_genes = ncol(x_train_scaled),
    Lambda = best_lambda,
    Pearson_R = pearson_r,
    R2 = r2,
    Spearman_R = spearman_r,
    MAE = mae
  )

  # Clock genes = non-zero LASSO coefficients
  beta <- as.matrix(
    coef(final_fit, s = best_lambda)
  )

  beta <- beta[
    rownames(beta) != "(Intercept)",
    ,
    drop = FALSE
  ]

  clock_genes <- rownames(beta)[beta[, 1] != 0]

  clock_gene_list_bulk[[rep]] <- clock_genes

  model_info_list[[rep]] <- list(
    train_idx = train_idx,
    test_idx = test_idx,
    candidate_genes = colnames(x_train_scaled),
    center = scaled$center,
    scale = scaled$scale,
    lambda = best_lambda,
    clock_genes = clock_genes
  )

  cat("Lambda.min:", signif(best_lambda, 5), "\n")
  cat("Clock genes:", length(clock_genes), "\n")
  cat("Test R2:", round(r2, 3), "\n")
  cat("Test Spearman:", round(spearman_r, 3), "\n")
  cat("Test MAE:", round(mae, 3), "days\n")
}


# 5. Summarize predictive performance
performance_df <- bind_rows(performance_list)

performance_summary <- data.frame(
  Metric = c("R2", "Spearman_R", "MAE"),
  Median = c(
    median(performance_df$R2, na.rm = TRUE),
    median(performance_df$Spearman_R, na.rm = TRUE),
    median(performance_df$MAE, na.rm = TRUE)
  ),
  IQR = c(
    IQR(performance_df$R2, na.rm = TRUE),
    IQR(performance_df$Spearman_R, na.rm = TRUE),
    IQR(performance_df$MAE, na.rm = TRUE)
  )
)

cat("\n===== Performance summary =====\n")
print(performance_summary)


# 6. Stable bulk clock genes
all_clock_genes_bulk <- unlist(clock_gene_list_bulk)

clock_gene_frequency <- sort(
  table(all_clock_genes_bulk),
  decreasing = TRUE
)

clock_gene_frequency_df <- data.frame(
  Gene = names(clock_gene_frequency),
  Frequency = as.integer(clock_gene_frequency),
  stringsAsFactors = FALSE
)

stable_clock_genes_bulk <- clock_gene_frequency_df %>%
  filter(Frequency >= stable_min_repeats) %>%
  arrange(desc(Frequency), Gene)

cat(
  "\nStable bulk clock genes selected in >=",
  stable_min_repeats,
  "/",
  n_repeats,
  " repeats: ",
  nrow(stable_clock_genes_bulk),
  "\n",
  sep = ""
)


# 7. Consensus held-out prediction for visualization only
pred_all <- bind_rows(prediction_list)

# A sample may be held out in more than one repeat.
# Average its held-out predictions so every biological sample
# appears only once in the visualization.
pred_consensus <- pred_all %>%
  group_by(Sample, SampleIndex, ActualAge) %>%
  summarise(
    PredictedAge = mean(PredictedAge, na.rm = TRUE),
    N_test_repeats = n(),
    .groups = "drop"
  )

consensus_r <- cor(
  pred_consensus$ActualAge,
  pred_consensus$PredictedAge,
  method = "pearson"
)

consensus_r2 <- consensus_r^2

p_clock <- ggplot(
  pred_consensus,
  aes(x = ActualAge, y = PredictedAge)
) +
  geom_point(alpha = 0.75, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.6) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    hjust = 1.1,
    vjust = -0.6,
    label = paste0("R² = ", round(consensus_r2, 3))
  ) +
  labs(
    x = "Chronological developmental age (days)",
    y = "Predicted developmental age (days)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(colour = "black")
  )

print(p_clock)


# 8. Save outputs
write.csv(
  performance_df,
  file.path(out_dir, "repeated_test_performance.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  performance_summary,
  file.path(out_dir, "performance_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  pred_all,
  file.path(out_dir, "all_heldout_predictions.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  pred_consensus,
  file.path(out_dir, "consensus_heldout_predictions.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  stable_clock_genes_bulk,
  file.path(out_dir, "stable_bulk_clock_genes.csv"),
  row.names = FALSE,
  quote = FALSE
)

saveRDS(
  clock_gene_list_bulk,
  file.path(out_dir, "clock_gene_list_bulk.rds")
)

saveRDS(
  model_info_list,
  file.path(out_dir, "model_info_list.rds")
)

ggsave(
  file.path(out_dir, "bulk_developmental_clock.pdf"),
  p_clock,
  width = 4,
  height = 4,
  useDingbats = FALSE
)

ggsave(
  file.path(out_dir, "bulk_developmental_clock.png"),
  p_clock,
  width = 4,
  height = 4,
  dpi = 600
)

cat("\nFinished. Results written to:", out_dir, "\n")
