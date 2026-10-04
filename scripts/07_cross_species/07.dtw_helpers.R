# Purpose: dtw helpers.

library(dtw)
library(ggplot2)

# basic

scale_pseudobulk <- function(pb) {
  pb <- as.matrix(pb)
  out <- t(scale(t(pb)))
  out[is.na(out)] <- 0
  return(out)
}

compute_similarity_matrix <- function(x, y, method = "spearman") {
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  sim_mat <- matrix(NA_real_, ncol(x), ncol(y))
  rownames(sim_mat) <- colnames(x)
  colnames(sim_mat) <- colnames(y)
  
  if (method %in% c("cosine", "cosine_distance")) {
    normalize <- function(v) {
      s <- sqrt(sum(v^2))
      if (is.na(s) || s == 0) return(rep(0, length(v)))
      v / s
    }
    x <- apply(x, 2, normalize)
    y <- apply(y, 2, normalize)
  }
  
  for (i in seq_len(ncol(x))) {
    for (j in seq_len(ncol(y))) {
      xi <- x[, i]
      yj <- y[, j]
      
      if (method %in% c("spearman", "pearson")) {
        sim_mat[i, j] <- suppressWarnings(cor(xi, yj, method = method, use = "complete.obs"))
      } else if (method == "1-correlation") {
        sim_mat[i, j] <- 1 - suppressWarnings(cor(xi, yj, method = "spearman", use = "complete.obs"))
      } else if (method == "cosine") {
        sim_mat[i, j] <- sum(xi * yj)
      } else if (method == "cosine_distance") {
        sim_mat[i, j] <- 1 - sum(xi * yj)
      } else {
        stop("Unsupported similarity method: ", method)
      }
    }
  }
  
  sim_mat[is.na(sim_mat)] <- 0
  return(sim_mat)
}

aggregate_stage_pseudobulk <- function(df,
                                       sep = "//",
                                       agg_fun = c("sum", "mean", "median"),
                                       species_col = 1,
                                       sample_col = 2,
                                       stage_col = 3,
                                       verbose = TRUE) {
  agg_fun <- match.arg(agg_fun)
  df <- as.matrix(df)
  
  cn <- colnames(df)
  if (is.null(cn)) {
    stop("Input matrix must have column names.")
  }
  
  split_info <- strsplit(cn, split = sep, fixed = TRUE)
  split_len <- lengths(split_info)
  
  if (any(split_len < max(species_col, sample_col, stage_col))) {
    bad_cols <- cn[split_len < max(species_col, sample_col, stage_col)]
    stop(
      "Some column names do not contain enough fields after splitting by '", sep, "'.\n",
      "Examples: ", paste(head(bad_cols, 5), collapse = ", ")
    )
  }
  
  meta <- data.frame(
    colname = cn,
    species = vapply(split_info, `[`, character(1), species_col),
    sample  = vapply(split_info, `[`, character(1), sample_col),
    stage   = vapply(split_info, `[`, character(1), stage_col),
    stringsAsFactors = FALSE
  )
  
  meta$stage_group <- paste(meta$species, meta$stage, sep = sep)
  
  if (verbose) {
    message("Input columns: ", ncol(df))
    message("Unique species-stage groups: ", length(unique(meta$stage_group)))
    tb <- table(meta$stage_group)
    message("Samples per species-stage group: min=", min(tb), ", median=", stats::median(tb), ", max=", max(tb))
  }
  
  group_levels <- unique(meta$stage_group)
  
  out <- sapply(group_levels, function(g) {
    idx <- which(meta$stage_group == g)
    sub <- df[, idx, drop = FALSE]
    
    if (ncol(sub) == 1) {
      return(sub[, 1])
    }
    
    if (agg_fun == "sum") {
      return(rowSums(sub, na.rm = TRUE))
    } else if (agg_fun == "mean") {
      return(rowMeans(sub, na.rm = TRUE))
    } else if (agg_fun == "median") {
      return(apply(sub, 1, stats::median, na.rm = TRUE))
    }
  })
  
  out <- as.matrix(out)
  rownames(out) <- rownames(df)
  
  colnames(out) <- group_levels
  
  return(out)
}


# dynamic gene selection

select_dynamic_genes_one_species <- function(pb,
                                             stage_order,
                                             top_n = 2000,
                                             method = c("variance", "spearman")) {
  method <- match.arg(method)
  pb <- as.matrix(pb)
  
  use_stage <- intersect(stage_order, colnames(pb))
  pb <- pb[, use_stage, drop = FALSE]
  
  if (ncol(pb) < 2) {
    return(rownames(pb))
  }
  
  if (method == "variance") {
    score <- apply(pb, 1, var, na.rm = TRUE)
  } else if (method == "spearman") {
    stage_num <- seq_len(ncol(pb))
    score <- apply(pb, 1, function(v) {
      abs(suppressWarnings(cor(v, stage_num, method = "spearman", use = "complete.obs")))
    })
  }
  
  score[is.na(score)] <- 0
  score <- sort(score, decreasing = TRUE)
  top_n <- min(top_n, length(score))
  
  return(names(score)[seq_len(top_n)])
}

select_dynamic_genes_multi_species <- function(pseudobulk_list,
                                               stage_order_list,
                                               top_n = 2000,
                                               method = c("variance", "spearman"),
                                               combine = c("freq", "union", "intersect"),
                                               min_species = 2,
                                               anchor_genes = NULL) {
  method <- match.arg(method)
  combine <- match.arg(combine)
  
  gene_sets <- lapply(names(pseudobulk_list), function(org) {
    select_dynamic_genes_one_species(
      pb = pseudobulk_list[[org]],
      stage_order = stage_order_list[[org]],
      top_n = top_n,
      method = method
    )
  })
  names(gene_sets) <- names(pseudobulk_list)
  
  if (combine == "union") {
    selected <- Reduce(union, gene_sets)
  } else if (combine == "intersect") {
    selected <- Reduce(intersect, gene_sets)
  } else if (combine == "freq") {
    tab <- table(unlist(gene_sets))
    selected <- names(tab)[tab >= min_species]
  }
  
  # Add anchor genes directly to the selected genes.
  if (!is.null(anchor_genes)) {
    selected <- union(selected, anchor_genes)
  }
  
  return(selected)
}

# pairwise dtw

run_dtw_pair_to_plot <- function(pb1, pb2,
                                 stage_order1, stage_order2,
                                 org1_name, org2_name,
                                 feature_genes,
                                 similarity_method = "1-correlation") {
  pb1 <- as.matrix(pb1)
  pb2 <- as.matrix(pb2)
  
  common_genes <- intersect(rownames(pb1), rownames(pb2))
  hvg <- intersect(common_genes, feature_genes)
  
  if (length(hvg) < 10) {
    stop("Too few shared selected genes in ", org1_name, " vs ", org2_name,
         ". n = ", length(hvg))
  }
  
  stage_order1_use <- intersect(stage_order1, colnames(pb1))
  stage_order2_use <- intersect(stage_order2, colnames(pb2))
  
  pb1 <- pb1[hvg, stage_order1_use, drop = FALSE]
  pb2 <- pb2[hvg, stage_order2_use, drop = FALSE]
  
  pb1 <- scale_pseudobulk(pb1)
  pb2 <- scale_pseudobulk(pb2)
  
  sim_mat <- compute_similarity_matrix(pb1, pb2, method = similarity_method)
  
  if (similarity_method %in% c("spearman", "pearson")) {
    cdst <- acos(pmin(pmax(sim_mat, -1), 1))
  } else {
    cdst <- sim_mat
  }
  
  # Apply DTW without an additional window constraint.
  aln <- dtw::dtw(
    cdst,
    step.pattern = dtw::asymmetric,
    open.begin = TRUE,
    open.end = TRUE
  )
  
  df_dist <- as.data.frame(as.table(cdst))
  colnames(df_dist) <- c("Stage1", "Stage2", "Distance")
  df_dist$Stage1 <- factor(df_dist$Stage1, levels = rownames(cdst))
  df_dist$Stage2 <- factor(df_dist$Stage2, levels = colnames(cdst))
  
  df_path <- data.frame(
    Stage1 = factor(rownames(cdst)[aln$index1], levels = rownames(cdst)),
    Stage2 = factor(colnames(cdst)[aln$index2], levels = colnames(cdst)),
    Distance = NA_real_,
    index = seq_along(aln$index1)
  )
  
  plot_obj <- ggplot(df_dist, aes(x = Stage1, y = Stage2, fill = Distance)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradientn(
      colours = c("#0B2E6B", "#1F5E9C", "#4F93C6", "#8FB8D8", "#BFD4E6", "#DCE6F0", "#EEF1F5")
    )+
    geom_path(
      data = df_path,
      aes(x = Stage1, y = Stage2, group = 1),
      inherit.aes = FALSE,
      color = "cyan",
      linewidth = 1.2
    ) +
    theme_minimal(base_size = 12) +
    coord_fixed(ratio = 1) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(),
      panel.grid = element_blank(),
      legend.key.width = unit(8, "pt"),
      legend.title = element_blank()
    ) +
    labs(title = paste(org1_name, "vs", org2_name), x = org1_name, y = org2_name)
  
  alignment_df <- data.frame(
    stage1 = rownames(cdst)[aln$index1],
    stage2 = colnames(cdst)[aln$index2],
    stringsAsFactors = FALSE
  )
  
  return(list(
    plot = plot_obj,
    alignment = alignment_df,
    score = aln$distance,
    name = paste0(org1_name, "_vs_", org2_name),
    feature_genes_used = hvg
  ))
}

# main pipeline

run_dtw_all_pairs_to_list <- function(pseudobulk_list,
                                      stage_order_list,
                                      top_n_hvgs = 2000,
                                      anchor_genes = NULL,
                                      similarity_method = "1-correlation",
                                      dynamic_method = c("spearman", "variance"),
                                      gene_combine = c("freq", "union", "intersect"),
                                      min_species = 2) {
  dynamic_method <- match.arg(dynamic_method)
  gene_combine <- match.arg(gene_combine)
  
  org_names <- names(pseudobulk_list)
  
  common_genes <- Reduce(intersect, lapply(pseudobulk_list, rownames))
  pseudobulk_list <- lapply(pseudobulk_list, function(pb) {
    as.matrix(pb[common_genes, , drop = FALSE])
  })
  
  feature_genes <- select_dynamic_genes_multi_species(
    pseudobulk_list = pseudobulk_list,
    stage_order_list = stage_order_list,
    top_n = top_n_hvgs,
    method = dynamic_method,
    combine = gene_combine,
    min_species = min_species,
    anchor_genes = anchor_genes
  )
  
  plot_list <- list()
  alignment_list <- list()
  score_df <- data.frame()
  
  for (i in 1:(length(org_names) - 1)) {
    for (j in (i + 1):length(org_names)) {
      cat("Processing:", org_names[i], "vs", org_names[j], "\n")
      
      res <- run_dtw_pair_to_plot(
        pb1 = pseudobulk_list[[i]],
        pb2 = pseudobulk_list[[j]],
        stage_order1 = stage_order_list[[i]],
        stage_order2 = stage_order_list[[j]],
        org1_name = org_names[i],
        org2_name = org_names[j],
        feature_genes = feature_genes,
        similarity_method = similarity_method
      )
      
      plot_list[[res$name]] <- res$plot
      alignment_list[[res$name]] <- res$alignment
      score_df <- rbind(
        score_df,
        data.frame(
          pair = res$name,
          score = res$score,
          n_feature_genes = length(res$feature_genes_used),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  return(list(
    plots = plot_list,
    alignments = alignment_list,
    scores = score_df,
    feature_genes = feature_genes
  ))
}
