# Purpose: cross tissue maturation loto.
analysis_dir <- Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input")

# cross-tissue developmental maturation axis
# Minimal publishable revision:
# 1) leave-one-tissue-out (LOTO) split is preserved;
# 2) gene selection, subtype-centering, and scaling are fitted on training tissues only;
# 3) pseudocell construction and LOTO are repeated across independent seeds;
# 4) stable maturation genes are summarized across repeat × held-out tissue models, including zero coefficients.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(Matrix)
  library(glmnet)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

out_dir <- "result6_cross_tissue_maturation_axis"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stage_order <- c("e55d", "e90d", "0d", "30d", "90d", "180d")

stage_age_map <- c(
  "e55d" = -59,
  "e90d" = -24,
  "0d"   = 0,
  "30d"  = 30,
  "90d"  = 90,
  "180d" = 180
)

obj <- readRDS(file.path(analysis_dir, "immune/immune_anno_scVI.rds"))

obj$stage <- as.character(obj$stage)
obj$tissue <- as.character(obj$tissue)
obj$immune_subtype <- as.character(obj$immune_subtype)
obj$dev_age <- unname(stage_age_map[obj$stage])

stopifnot(!any(is.na(obj$dev_age)))

model_configs <- list(
  T_core = c(
    "CD4/activated T cell",
    "CD8 T cell",
    "Innate-like T cell"
  ),

  Resident_macrophage = c(
    "Resident macrophage"
  ),

  Shared_macrophage = c(
    "Resident macrophage",
    "Monocyte-derived macrophage",
    "IFN-response macrophage"
  )
)

get_assay_matrix <- function(obj, assay = "RNA", layer = "data") {
  DefaultAssay(obj) <- assay
  mat <- tryCatch(
    GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) GetAssayData(obj, assay = assay, slot = layer)
  )
  mat
}

# Function to create pseudocells
make_pseudocells <- function(obj,
                             subtypes_use,
                             assay = "RNA",
                             layer = "data",
                             cells_per_pseudocell = 30,
                             n_pseudocells_per_group = 30,
                             min_cells_per_group = 80,
                             seed = 123,
                             repeat_id = 1) {
  set.seed(seed)

  meta <- obj@meta.data %>%
    mutate(cell = rownames(.)) %>%
    filter(immune_subtype %in% subtypes_use)

  expr <- get_assay_matrix(obj, assay = assay, layer = layer)

  meta$group_id <- paste(
    meta$immune_subtype,
    meta$tissue,
    meta$stage,
    sep = "|||"
  )

  group_info <- meta %>%
    count(group_id, immune_subtype, tissue, stage, dev_age, name = "n_cells") %>%
    filter(n_cells >= min_cells_per_group)

  if (nrow(group_info) < 6) {
    stop("Too few eligible subtype × tissue × stage groups. Lower min_cells_per_group or change subtypes.")
  }

  expr_list <- list()
  meta_list <- list()
  idx <- 1

  for (gid in group_info$group_id) {
    cells <- meta$cell[meta$group_id == gid]
    group_meta <- group_info[group_info$group_id == gid, ]

    for (i in seq_len(n_pseudocells_per_group)) {
      sampled_cells <- sample(cells, cells_per_pseudocell, replace = TRUE)
      avg_expr <- Matrix::rowMeans(expr[, sampled_cells, drop = FALSE])

      pc_id <- paste0("R", repeat_id, "_PC", idx)
      expr_list[[idx]] <- avg_expr

      meta_list[[idx]] <- data.frame(
        pseudocell_id = pc_id,
        repeat_id = repeat_id,
        immune_subtype = group_meta$immune_subtype,
        tissue = group_meta$tissue,
        stage = group_meta$stage,
        dev_age = group_meta$dev_age,
        n_cells_group = group_meta$n_cells,
        stringsAsFactors = FALSE
      )

      idx <- idx + 1
    }
  }

  pc_expr <- do.call(cbind, expr_list)
  pc_meta <- bind_rows(meta_list)

  colnames(pc_expr) <- pc_meta$pseudocell_id
  rownames(pc_meta) <- pc_meta$pseudocell_id

  list(expr = pc_expr, meta = pc_meta)
}

# Select training genes from training tissues only
select_training_genes <- function(pc_expr_train, n_hvg = 3000) {
  genes <- rownames(pc_expr_train)

  exclude_genes <- unique(c(
    grep("^MT-", genes, value = TRUE),
    grep("^mt-", genes, value = TRUE),
    grep("^RPL", genes, value = TRUE),
    grep("^RPS", genes, value = TRUE),
    grep("^Rpl", genes, value = TRUE),
    grep("^Rps", genes, value = TRUE),
    c("MKI67", "TOP2A", "UBE2C", "STMN1", "PCNA", "HMGB2")
  ))

  genes_use <- setdiff(genes, exclude_genes)

  mat <- pc_expr_train[genes_use, , drop = FALSE]
  gene_mean <- Matrix::rowMeans(mat)
  gene_sd <- apply(mat, 1, sd)

  keep <- names(gene_mean)[gene_mean > 0.01 & gene_sd > quantile(gene_sd, 0.50, na.rm = TRUE)]
  gene_sd_keep <- gene_sd[keep]

  genes_use <- names(sort(gene_sd_keep, decreasing = TRUE))
  genes_use <- genes_use[seq_len(min(n_hvg, length(genes_use)))]

  genes_use
}

# Fit subtype-centering means in training tissues only, then apply to train and test.
fit_subtype_center <- function(expr_train, meta_train) {
  global_mean <- Matrix::rowMeans(expr_train)
  subtype_means <- lapply(unique(meta_train$immune_subtype), function(st) {
    cells <- rownames(meta_train)[meta_train$immune_subtype == st]
    if (length(cells) < 2) return(global_mean)
    Matrix::rowMeans(expr_train[, cells, drop = FALSE])
  })
  names(subtype_means) <- unique(meta_train$immune_subtype)
  list(global_mean = global_mean, subtype_means = subtype_means)
}

apply_subtype_center <- function(expr_mat, meta, center_fit) {
  expr_centered <- expr_mat
  for (st in unique(meta$immune_subtype)) {
    cells <- rownames(meta)[meta$immune_subtype == st]
    mu <- center_fit$subtype_means[[st]]
    if (is.null(mu)) mu <- center_fit$global_mean
    expr_centered[, cells] <- expr_mat[, cells, drop = FALSE] - mu
  }
  expr_centered
}

fit_scale <- function(x_train) {
  mu <- colMeans(x_train)
  sigma <- apply(x_train, 2, sd)
  sigma[is.na(sigma) | sigma == 0] <- 1
  list(center = mu, scale = sigma)
}

apply_scale <- function(x, scale_fit) {
  sweep(sweep(x, 2, scale_fit$center, FUN = "-"), 2, scale_fit$scale, FUN = "/")
}

# Run leave-one-tissue-out elastic net with train-only preprocessing
run_loto_elasticnet <- function(pc_expr,
                                pc_meta,
                                n_hvg = 3000,
                                alpha = 0.5,
                                nfolds = 5,
                                seed = 123,
                                repeat_id = 1,
                                subtype_center = TRUE) {
  set.seed(seed)

  y_all <- pc_meta$dev_age
  tissues <- unique(pc_meta$tissue)

  pred_list <- list()
  coef_list <- list()
  gene_list <- list()

  for (tt in tissues) {
    train_idx <- pc_meta$tissue != tt
    test_idx  <- pc_meta$tissue == tt

    if (sum(test_idx) < 10 || length(unique(y_all[test_idx])) < 2) next
    if (sum(train_idx) < 20 || length(unique(y_all[train_idx])) < 3) next

    meta_train <- pc_meta[train_idx, , drop = FALSE]
    meta_test <- pc_meta[test_idx, , drop = FALSE]

    genes_use <- select_training_genes(pc_expr[, rownames(meta_train), drop = FALSE], n_hvg = n_hvg)
    if (length(genes_use) < 100) next

    expr_train <- pc_expr[genes_use, rownames(meta_train), drop = FALSE]
    expr_test <- pc_expr[genes_use, rownames(meta_test), drop = FALSE]

    if (subtype_center) {
      center_fit <- fit_subtype_center(expr_train, meta_train)
      expr_train <- apply_subtype_center(expr_train, meta_train, center_fit)
      expr_test <- apply_subtype_center(expr_test, meta_test, center_fit)
    }

    x_train <- t(as.matrix(expr_train))
    x_test <- t(as.matrix(expr_test))

    scale_fit <- fit_scale(x_train)
    x_train <- apply_scale(x_train, scale_fit)
    x_test <- apply_scale(x_test, scale_fit)

    y_train <- y_all[train_idx]
    y_test <- y_all[test_idx]

    cv_nfolds <- min(nfolds, nrow(x_train))
    if (cv_nfolds < 3) next

    cvfit <- cv.glmnet(
      x = x_train,
      y = y_train,
      family = "gaussian",
      alpha = alpha,
      nfolds = cv_nfolds,
      standardize = FALSE
    )

    pred <- as.numeric(predict(cvfit, newx = x_test, s = "lambda.min"))
    model_id <- paste0("R", repeat_id, "__heldout__", tt)

    pred_list[[model_id]] <- data.frame(
      pseudocell_id = rownames(x_test),
      repeat_id = repeat_id,
      model_id = model_id,
      heldout_tissue = tt,
      actual_age = y_test,
      predicted_age = pred,
      stringsAsFactors = FALSE
    )

    coef_mat <- coef(cvfit, s = "lambda.min")

    coef_list[[model_id]] <- data.frame(
      gene = rownames(coef_mat),
      coefficient = as.numeric(coef_mat[, 1]),
      repeat_id = repeat_id,
      model_id = model_id,
      heldout_tissue = tt,
      stringsAsFactors = FALSE
    ) %>%
      filter(gene != "(Intercept)")

    gene_list[[model_id]] <- data.frame(
      model_id = model_id,
      repeat_id = repeat_id,
      heldout_tissue = tt,
      gene = genes_use,
      stringsAsFactors = FALSE
    )
  }

  pred_df <- bind_rows(pred_list) %>%
    left_join(
      pc_meta %>%
        mutate(pseudocell_id = rownames(.)) %>%
        select(pseudocell_id, repeat_id, immune_subtype, tissue, stage, dev_age),
      by = c("pseudocell_id", "repeat_id")
    )

  coef_df <- bind_rows(coef_list)
  gene_df <- bind_rows(gene_list)

  list(pred = pred_df, coef = coef_df, genes = gene_df)
}

safe_cor <- function(x, y, method = "spearman") {
  if (length(unique(x)) < 2 || length(unique(y)) < 2) return(NA_real_)
  suppressWarnings(cor(x, y, method = method, use = "complete.obs"))
}

# Calculate model performance metrics
calc_model_performance <- function(pred_df) {
  by_tissue_repeat <- pred_df %>%
    group_by(repeat_id, heldout_tissue) %>%
    summarise(
      n = n(),
      n_stages = n_distinct(stage),
      pearson = safe_cor(actual_age, predicted_age, method = "pearson"),
      spearman = safe_cor(actual_age, predicted_age, method = "spearman"),
      mae = mean(abs(actual_age - predicted_age)),
      .groups = "drop"
    )

  by_tissue <- by_tissue_repeat %>%
    group_by(heldout_tissue) %>%
    summarise(
      n_repeats = n_distinct(repeat_id),
      n = sum(n),
      median_pearson = median(pearson, na.rm = TRUE),
      q25_pearson = quantile(pearson, 0.25, na.rm = TRUE),
      q75_pearson = quantile(pearson, 0.75, na.rm = TRUE),
      median_spearman = median(spearman, na.rm = TRUE),
      q25_spearman = quantile(spearman, 0.25, na.rm = TRUE),
      q75_spearman = quantile(spearman, 0.75, na.rm = TRUE),
      median_mae = median(mae, na.rm = TRUE),
      q25_mae = quantile(mae, 0.25, na.rm = TRUE),
      q75_mae = quantile(mae, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  overall_repeat <- pred_df %>%
    group_by(repeat_id) %>%
    summarise(
      n = n(),
      n_tissues = n_distinct(heldout_tissue),
      pearson = safe_cor(actual_age, predicted_age, method = "pearson"),
      spearman = safe_cor(actual_age, predicted_age, method = "spearman"),
      mae = mean(abs(actual_age - predicted_age)),
      .groups = "drop"
    )

  overall <- overall_repeat %>%
    summarise(
      n_repeats = n_distinct(repeat_id),
      n = sum(n),
      n_tissues_median = median(n_tissues, na.rm = TRUE),
      median_pearson = median(pearson, na.rm = TRUE),
      q25_pearson = quantile(pearson, 0.25, na.rm = TRUE),
      q75_pearson = quantile(pearson, 0.75, na.rm = TRUE),
      median_spearman = median(spearman, na.rm = TRUE),
      q25_spearman = quantile(spearman, 0.25, na.rm = TRUE),
      q75_spearman = quantile(spearman, 0.75, na.rm = TRUE),
      median_mae = median(mae, na.rm = TRUE),
      q25_mae = quantile(mae, 0.25, na.rm = TRUE),
      q75_mae = quantile(mae, 0.75, na.rm = TRUE)
    )

  list(
    by_tissue_repeat = by_tissue_repeat,
    by_tissue = by_tissue,
    overall_repeat = overall_repeat,
    overall = overall
  )
}

# Extract stable maturation genes across repeat × held-out tissue models.
extract_stable_maturation_genes <- function(coef_df,
                                            gene_df,
                                            min_model_frac = 0.5,
                                            sign_consistency_cutoff = 0.75) {
  model_ids <- unique(gene_df$model_id)
  n_models_total <- length(model_ids)

  coef_complete <- gene_df %>%
    left_join(
      coef_df %>% select(model_id, gene, coefficient),
      by = c("model_id", "gene")
    ) %>%
    mutate(
      coefficient = ifelse(is.na(coefficient), 0, coefficient),
      selected = coefficient != 0,
      sign = case_when(
        coefficient > 0 ~ 1L,
        coefficient < 0 ~ -1L,
        TRUE ~ 0L
      )
    )

  stable_coef <- coef_complete %>%
    group_by(gene) %>%
    summarise(
      n_models_tested = n_distinct(model_id),
      n_models_selected = n_distinct(model_id[selected]),
      model_frac_tested = n_models_tested / n_models_total,
      selection_freq = n_models_selected / n_models_tested,
      mean_coef_all = mean(coefficient),
      median_coef_all = median(coefficient),
      mean_coef_nonzero = ifelse(any(selected), mean(coefficient[selected]), 0),
      median_coef_nonzero = ifelse(any(selected), median(coefficient[selected]), 0),
      positive_frac_nonzero = ifelse(any(selected), mean(coefficient[selected] > 0), NA_real_),
      negative_frac_nonzero = ifelse(any(selected), mean(coefficient[selected] < 0), NA_real_),
      sign_consistency = ifelse(
        any(selected),
        max(mean(coefficient[selected] > 0), mean(coefficient[selected] < 0)),
        0
      ),
      direction = case_when(
        n_models_selected == 0 ~ "not_selected",
        positive_frac_nonzero >= negative_frac_nonzero ~ "increase",
        TRUE ~ "decrease"
      ),
      .groups = "drop"
    ) %>%
    filter(
      n_models_tested >= ceiling(n_models_total * min_model_frac),
      selection_freq >= min_model_frac,
      sign_consistency >= sign_consistency_cutoff
    ) %>%
    arrange(desc(selection_freq), desc(abs(mean_coef_nonzero)))

  list(stable = stable_coef, complete = coef_complete)
}

# Run one maturation model with repeated pseudocell resampling
run_one_maturation_model <- function(obj,
                                     model_name,
                                     subtypes_use,
                                     min_cells_per_group = 80,
                                     cells_per_pseudocell = 30,
                                     n_pseudocells_per_group = 30,
                                     n_hvg = 3000,
                                     n_repeats = 30,
                                     seed = 123,
                                     assay = "RNA",
                                     layer = "data") {

  message("Running model: ", model_name)

  model_dir <- file.path(out_dir, model_name)
  dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)

  pred_all <- list()
  coef_all <- list()
  gene_all <- list()
  pc_first <- NULL

  for (rr in seq_len(n_repeats)) {
    message("  repeat ", rr, "/", n_repeats)
    seed_rr <- seed + rr - 1

    pc <- make_pseudocells(
      obj = obj,
      subtypes_use = subtypes_use,
      assay = assay,
      layer = layer,
      cells_per_pseudocell = cells_per_pseudocell,
      n_pseudocells_per_group = n_pseudocells_per_group,
      min_cells_per_group = min_cells_per_group,
      seed = seed_rr,
      repeat_id = rr
    )

    if (is.null(pc_first)) pc_first <- pc

    subtype_center <- length(unique(pc$meta$immune_subtype)) > 1

    loto <- run_loto_elasticnet(
      pc_expr = pc$expr,
      pc_meta = pc$meta,
      n_hvg = n_hvg,
      subtype_center = subtype_center,
      seed = seed_rr,
      repeat_id = rr
    )

    pred_all[[rr]] <- loto$pred
    coef_all[[rr]] <- loto$coef
    gene_all[[rr]] <- loto$genes
  }

  pred_df <- bind_rows(pred_all)
  coef_df <- bind_rows(coef_all)
  gene_df <- bind_rows(gene_all)

  if (nrow(pred_df) == 0) stop("No valid LOTO predictions were generated for model: ", model_name)

  perf <- calc_model_performance(pred_df)
  stable_obj <- extract_stable_maturation_genes(coef_df, gene_df)
  stable_coef <- stable_obj$stable
  coef_complete <- stable_obj$complete

  saveRDS(pc_first, file.path(model_dir, paste0(model_name, "_pseudocells_repeat1.rds")))

  write.table(
    pred_df,
    file.path(model_dir, paste0(model_name, "_LOTO_predictions.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    coef_df,
    file.path(model_dir, paste0(model_name, "_LOTO_coefficients_all_genes.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    coef_complete,
    file.path(model_dir, paste0(model_name, "_LOTO_coefficients_complete_train_gene_space.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    perf$by_tissue_repeat,
    file.path(model_dir, paste0(model_name, "_performance_by_tissue_repeat.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    perf$by_tissue,
    file.path(model_dir, paste0(model_name, "_performance_by_tissue_summary.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    perf$overall_repeat,
    file.path(model_dir, paste0(model_name, "_performance_overall_repeat.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    perf$overall,
    file.path(model_dir, paste0(model_name, "_performance_overall_summary.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    stable_coef,
    file.path(model_dir, paste0(model_name, "_stable_maturation_genes.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  list(
    model_name = model_name,
    pc = pc_first,
    pred = pred_df,
    coef = coef_df,
    coef_complete = coef_complete,
    performance = perf,
    stable_coef = stable_coef
  )
}

n_repeats_final <- 20

model_results <- list()

model_results$T_core <- run_one_maturation_model(
  obj = obj,
  model_name = "T_core",
  subtypes_use = model_configs$T_core,
  min_cells_per_group = 80,
  n_repeats = n_repeats_final
)

model_results$Resident_macrophage <- run_one_maturation_model(
  obj = obj,
  model_name = "Resident_macrophage",
  subtypes_use = model_configs$Resident_macrophage,
  min_cells_per_group = 80,
  n_repeats = n_repeats_final
)

model_results$Shared_macrophage <- run_one_maturation_model(
  obj = obj,
  model_name = "Shared_macrophage",
  subtypes_use = model_configs$Shared_macrophage,
  min_cells_per_group = 50,
  n_repeats = n_repeats_final
)

# plots
theme_pub <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black", size = base_size),
      axis.title = element_text(color = "black", size = base_size + 1),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold")
    )
}

stage_cols <- c(
  "e55d" = "#443A83",
  "e90d" = "#31688E",
  "0d"   = "#21908C",
  "30d"  = "#35B779",
  "90d"  = "#8FD744",
  "180d" = "#FDE725"
)

tissue_cols <- c(
  "adipose" = "#EEA236FF",
  "cerebrum" = "#357EBDFF",
  "duodenum" = "#5CB85CFF",
  "heart" = "#D43F3A99",
  "hypothalamus" = "#46B8DAFF",
  "liver" = "#20854E99",
  "muscle" = "#D43F3AFF",
  "skeletal muscle" = "#D43F3AFF"
)

# plot predicted vs actual ages
plot_predicted_vs_actual <- function(pred_df, model_name, output_pdf = NULL) {
  p <- ggplot(pred_df, aes(x = actual_age, y = predicted_age, color = heldout_tissue)) +
    geom_point(size = 0.35, alpha = 0.25) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3) +
    scale_color_manual(values = tissue_cols) +
    labs(
      x = "Actual developmental age",
      y = "Predicted developmental age",
      color = "Held-out tissue",
      title = paste0(model_name, ": repeated leave-one-tissue-out prediction")
    ) +
    theme_pub(base_size = 8)

  if (!is.null(output_pdf)) {
    ggsave(output_pdf, p, width = 4.5, height = 4)
  }

  p
}

plot_predicted_vs_actual(
  model_results$T_core$pred,
  model_name = "T-core",
  output_pdf = file.path(out_dir, "T_core_predicted_vs_actual.pdf")
)

plot_predicted_vs_actual(
  model_results$Resident_macrophage$pred,
  model_name = "Resident macrophage",
  output_pdf = file.path(out_dir, "Resident_macrophage_predicted_vs_actual.pdf")
)

plot_predicted_vs_actual(
  model_results$Shared_macrophage$pred,
  model_name = "Shared macrophage",
  output_pdf = file.path(out_dir, "Shared_macrophage_predicted_vs_actual.pdf")
)

# plot leave-one-tissue-out performance
plot_loto_performance <- function(perf_by_tissue, model_name, output_pdf = NULL) {
  df <- perf_by_tissue %>%
    arrange(desc(median_spearman)) %>%
    mutate(heldout_tissue = factor(heldout_tissue, levels = heldout_tissue))

  p <- ggplot(df, aes(x = median_spearman, y = heldout_tissue)) +
    geom_errorbarh(aes(xmin = q25_spearman, xmax = q75_spearman), height = 0.18, linewidth = 0.35) +
    geom_point(size = 1.4) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    labs(
      x = "Spearman correlation across resampling repeats",
      y = NULL,
      title = paste0(model_name, ": cross-tissue performance")
    ) +
    theme_pub(base_size = 8)

  if (!is.null(output_pdf)) {
    ggsave(output_pdf, p, width = 3.8, height = 3)
  }

  p
}

plot_loto_performance(
  model_results$T_core$performance$by_tissue,
  model_name = "T-core",
  output_pdf = file.path(out_dir, "T_core_LOTO_performance.pdf")
)

plot_loto_performance(
  model_results$Resident_macrophage$performance$by_tissue,
  model_name = "Resident macrophage",
  output_pdf = file.path(out_dir, "Resident_macrophage_LOTO_performance.pdf")
)

plot_loto_performance(
  model_results$Shared_macrophage$performance$by_tissue,
  model_name = "Shared macrophage",
  output_pdf = file.path(out_dir, "Shared_macrophage_LOTO_performance.pdf")
)

# plot model comparison
model_compare_df <- bind_rows(lapply(names(model_results), function(mn) {
  x <- model_results[[mn]]$performance$overall
  x$model <- mn
  x
}))

# plot maturation gene heatmap
plot_maturation_gene_heatmap <- function(model_result,
                                         model_name,
                                         top_n = 80,
                                         output_pdf = NULL) {
  pc <- model_result$pc
  stable_coef <- model_result$stable_coef

  pc_expr <- pc$expr
  pc_meta <- pc$meta

  if (nrow(stable_coef) == 0) {
    warning("No stable genes for heatmap: ", model_name)
    return(NULL)
  }

  genes_top <- stable_coef %>%
    arrange(desc(selection_freq), desc(abs(mean_coef_nonzero))) %>%
    slice_head(n = top_n) %>%
    pull(gene)

  genes_top <- intersect(genes_top, rownames(pc_expr))
  if (length(genes_top) < 2) {
    warning("Too few stable genes for heatmap: ", model_name)
    return(NULL)
  }

  pc_meta$stage <- factor(pc_meta$stage, levels = stage_order)
  pc_meta$group <- paste(pc_meta$immune_subtype, pc_meta$stage, sep = " | ")

  group_levels <- pc_meta %>%
    distinct(immune_subtype, stage, group) %>%
    arrange(immune_subtype, stage) %>%
    pull(group)

  avg_mat <- sapply(group_levels, function(g) {
    cells <- rownames(pc_meta)[pc_meta$group == g]
    Matrix::rowMeans(pc_expr[genes_top, cells, drop = FALSE])
  })

  avg_mat_z <- t(scale(t(avg_mat)))
  avg_mat_z[is.na(avg_mat_z)] <- 0
  avg_mat_z[avg_mat_z > 2.5] <- 2.5
  avg_mat_z[avg_mat_z < -2.5] <- -2.5

  gene_direction <- stable_coef %>%
    filter(gene %in% genes_top) %>%
    mutate(direction = factor(direction, levels = c("increase", "decrease")))

  row_split <- gene_direction$direction[
    match(rownames(avg_mat_z), gene_direction$gene)
  ]

  ht <- Heatmap(
    avg_mat_z,
    name = "z-score",
    col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    row_split = row_split,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 5),
    column_names_rot = 45,
    row_title_gp = gpar(fontsize = 8, fontface = "bold"),
    column_title = paste0(model_name, ": stable maturation-axis genes"),
    column_title_gp = gpar(fontsize = 9, fontface = "bold"),
    rect_gp = gpar(col = "white", lwd = 0.1)
  )

  if (!is.null(output_pdf)) {
    pdf(output_pdf, width = 7.5, height = 5)
    draw(ht)
    dev.off()
  }

  ht
}

plot_maturation_gene_heatmap(
  model_results$T_core,
  model_name = "T-core",
  output_pdf = file.path(out_dir, "T_core_maturation_gene_heatmap.pdf")
)

plot_maturation_gene_heatmap(
  model_results$Resident_macrophage,
  model_name = "Resident macrophage",
  output_pdf = file.path(out_dir, "Resident_macrophage_maturation_gene_heatmap.pdf")
)

plot_maturation_gene_heatmap(
  model_results$Shared_macrophage,
  model_name = "Shared macrophage",
  output_pdf = file.path(out_dir, "Shared_macrophage_maturation_gene_heatmap.pdf")
)

saveRDS(model_results, file.path(out_dir, "cross_tissue_maturation_axis_model_results.rds"))
