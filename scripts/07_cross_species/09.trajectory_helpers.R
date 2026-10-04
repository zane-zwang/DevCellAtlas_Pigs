# Purpose: trajectory helpers.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(edgeR)
  library(Mfuzz)
  library(Biobase)
  library(dtw)
  library(tibble)
})

# 0) helper

# 0.1 equal-number pseudotime binning
bin_pseudotime_equal_n <- function(pt, n_bins = 9) {
  stopifnot(length(pt) > 0)
  rk <- rank(pt, ties.method = "average", na.last = "keep")
  q <- ceiling(rk / (sum(!is.na(pt)) / n_bins))
  q[q < 1] <- 1
  q[q > n_bins] <- n_bins
  paste0("bin", q)
}

# 0.2 pseudobulk by (bin, sample), requiring >= min_cells
pseudobulk_sum <- function(counts, meta, min_cells = 50) {
  stopifnot(all(c("cell_id", "bin", "sample") %in% colnames(meta)))

  keep_grp <- meta %>%
    count(bin, sample, name = "n_cells") %>%
    filter(n_cells >= min_cells)

  meta2 <- meta %>%
    inner_join(keep_grp, by = c("bin", "sample")) %>%
    mutate(group = paste(bin, sample, sep = "|"))

  groups <- unique(meta2$group)

  pb <- sapply(groups, function(g) {
    cells <- meta2$cell_id[meta2$group == g]
    Matrix::rowSums(counts[, cells, drop = FALSE])
  })
  pb <- as.matrix(pb)
  colnames(pb) <- groups

  group_df <- data.frame(
    group = groups,
    bin = sub("\\|.*$", "", groups),
    sample = sub("^.*\\|", "", groups),
    stringsAsFactors = FALSE
  )

  list(pb = pb, group_df = group_df)
}

# 0.3 mean across samples per bin -> gene x bin
# keep optional n_reps output for QC
mean_across_reps <- function(pb_mat, group_df) {
  bins <- unique(group_df$bin)
  bins <- bins[order(readr::parse_number(bins))]

  out <- sapply(bins, function(b) {
    cols <- group_df$group[group_df$bin == b]
    if (length(cols) == 1) {
      pb_mat[, cols]
    } else {
      rowMeans(pb_mat[, cols, drop = FALSE])
    }
  })
  out <- as.matrix(out)
  colnames(out) <- bins

  nrep <- sapply(bins, function(b) sum(group_df$bin == b))

  list(
    mat = out,
    n_reps = nrep
  )
}

# 0.4 CPM
to_cpm <- function(count_mat, log = FALSE, prior.count = 1) {
  edgeR::cpm(count_mat, log = log, prior.count = prior.count)
}

# 0.5 dynamic-gene score on binned expression
# Combine variance and trajectory range for dynamic-gene selection.
# score combines:
#   - row variance on logCPM
#   - trajectory range on logCPM
#   - optional smoothness penalty can be added later
select_dynamic_genes <- function(mat_gene_bin,
                                 top_n = 5000,
                                 min_mean = 0.5,
                                 use_log = TRUE) {
  # mat_gene_bin: gene x bin, preferably CPM or pseudobulk mean counts
  x <- mat_gene_bin

  if (use_log) {
    x <- log2(x + 1)
  }

  mu <- rowMeans(x)
  va <- apply(x, 1, var)
  rg <- apply(x, 1, function(v) diff(range(v)))

  score <- scale(va)[, 1] + 0.5 * scale(rg)[, 1]
  score[mu < min_mean] <- -Inf
  score[is.na(score)] <- -Inf

  ord <- order(score, decreasing = TRUE)
  sel <- rownames(x)[ord][seq_len(min(top_n, sum(is.finite(score))))]

  data.frame(
    gene = rownames(x),
    mean_expr = mu,
    variance = va,
    range = rg,
    score = score,
    selected = rownames(x) %in% sel,
    stringsAsFactors = FALSE
  )
}


select_dynamic_genes_vmr_cor <- function(mat_gene_bin,
                                         pseudotime = NULL,
                                         min_mean = 0.5,
                                         vmr_quantile = 0.75,
                                         cor_cut = 0.2,
                                         use_log = TRUE) {
  # mat_gene_bin: gene x bin
  # pseudotime: numeric vector, length = ncol(mat_gene_bin)

  x <- mat_gene_bin

  if (use_log) {
    x <- log2(x + 1)
  }

  # 1. basic stats
  mu <- rowMeans(x)
  va <- apply(x, 1, var)

  # avoid division by zero
  vmr <- va / (mu + 1e-6)

  # 2. filter low expression
  keep_expr <- mu >= min_mean

  # 3. VMR-based selection
  vmr_thr <- quantile(vmr[keep_expr], vmr_quantile, na.rm = TRUE)
  keep_vmr <- vmr >= vmr_thr

  # 4. pseudotime correlation (optional but strongly recommended)
  if (!is.null(pseudotime)) {
    cor_val <- apply(x, 1, function(g) {
      suppressWarnings(cor(g, pseudotime, method = "spearman"))
    })
    cor_val[is.na(cor_val)] <- 0

    keep_cor <- abs(cor_val) >= cor_cut
  } else {
    cor_val <- rep(NA_real_, nrow(x))
    keep_cor <- rep(TRUE, nrow(x))
  }

  # Combine expression, variance-to-mean ratio and pseudotime criteria.
  selected <- keep_expr & keep_vmr & keep_cor

  # 6. output
  data.frame(
    gene = rownames(x),
    mean_expr = mu,
    variance = va,
    vmr = vmr,
    cor_pseudotime = cor_val,
    selected = selected,
    stringsAsFactors = FALSE
  )
}

# 0.6 consensus gene set
# instead of strict 3-way topN intersection
# modes:
#   "intersect3" = selected in all three species
#   "atleast2"   = selected in >=2 species
#   "union_rank" = union then rank by mean z-score

build_consensus_gene_set <- function(dynamic_stats_list,
                                     mode = c("atleast2", "intersect3", "union_rank"),
                                     top_union = 6000) {
  mode <- match.arg(mode)

  sel_list <- lapply(dynamic_stats_list, function(df) df$gene[df$selected])

  if (mode == "intersect3") {
    genes <- Reduce(intersect, sel_list)
    return(genes)
  }

  if (mode == "atleast2") {
    all_genes <- sort(unique(unlist(sel_list)))
    hit_n <- sapply(all_genes, function(g) sum(sapply(sel_list, function(v) g %in% v)))
    genes <- all_genes[hit_n >= 2]
    return(genes)
  }

  if (mode == "union_rank") {

    score_list <- lapply(dynamic_stats_list, function(df) {

      if ("score" %in% colnames(df)) {
        x <- df$score
      } else if ("vmr" %in% colnames(df)) {
        x <- df$vmr
      } else {
        stop("union_rank mode requires a numeric score column (e.g. score or vmr).")
      }
      names(x) <- df$gene
      x
    })

    all_genes <- sort(unique(unlist(sel_list)))
    score_mat <- sapply(score_list, function(sc) sc[all_genes])
    score_mat[is.na(score_mat)] <- -Inf

    tmp <- score_mat
    tmp[!is.finite(tmp)] <- NA_real_

    zmat <- apply(tmp, 2, function(v) {
      if (all(is.na(v))) return(rep(NA_real_, length(v)))
      as.numeric(scale(v))
    })
    if (is.vector(zmat)) zmat <- matrix(zmat, ncol = 1)

    rownames(zmat) <- all_genes
    colnames(zmat) <- names(score_list)

    mean_z <- rowMeans(zmat, na.rm = TRUE)
    hit_n <- rowSums(!is.na(score_mat) & is.finite(score_mat))

    out <- data.frame(
      gene = all_genes,
      mean_z = mean_z,
      hit_n = hit_n,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::arrange(desc(hit_n), desc(mean_z))

    return(out$gene[seq_len(min(top_union, nrow(out)))])
  }
}

# 0.7 ortholog-feature matrix
# here rownames are already unified orthogroup IDs,
# so mapping is trivial
# output rows = gene|species
build_feature_matrix_same_id <- function(expr_by_species, genes_use) {
  feats <- list()
  for (sp in names(expr_by_species)) {
    m <- expr_by_species[[sp]]
    g <- intersect(genes_use, rownames(m))
    m2 <- m[g, , drop = FALSE]
    rownames(m2) <- paste0(g, "|", sp)
    feats[[sp]] <- m2
  }
  do.call(rbind, feats)
}

# 0.8 Mfuzz
run_mfuzz <- function(mat_feature_bin,
                      k = 8,
                      m_est = NULL,
                      standardize = TRUE,
                      min_sd = 1e-8,
                      verbose = TRUE) {
  x <- as.matrix(mat_feature_bin)

  keep_finite <- apply(x, 1, function(v) all(is.finite(v)))
  if (verbose) {
    cat("Rows with all finite values:", sum(keep_finite), "/", nrow(x), "\n")
  }
  x <- x[keep_finite, , drop = FALSE]

  row_sd <- apply(x, 1, sd, na.rm = TRUE)
  keep_sd <- is.finite(row_sd) & row_sd > min_sd
  if (verbose) {
    cat("Rows with SD >", min_sd, ":", sum(keep_sd), "/", nrow(x), "\n")
    cat("Removed near-constant rows:", sum(!keep_sd), "\n")
  }
  x <- x[keep_sd, , drop = FALSE]

  stopifnot(nrow(x) > k)
  stopifnot(all(apply(x, 1, function(v) all(is.finite(v)))))

  eset <- Biobase::ExpressionSet(x)

  if (standardize) {
    eset <- Mfuzz::standardise(eset)
    expr_std <- Biobase::exprs(eset)

    keep_std <- apply(expr_std, 1, function(v) all(is.finite(v)))
    if (verbose) {
      cat("Rows finite after standardise:", sum(keep_std), "/", nrow(expr_std), "\n")
      cat("Removed after standardise:", sum(!keep_std), "\n")
    }

    if (!all(keep_std)) {
      expr_std <- expr_std[keep_std, , drop = FALSE]
      eset <- Biobase::ExpressionSet(expr_std)
    }
  }

  if (is.null(m_est)) {
    m_est <- Mfuzz::mestimate(eset)
  }

  cl <- Mfuzz::mfuzz(eset, c = k, m = m_est)
  membership <- cl$membership

  list(
    cl = cl,
    membership = membership,
    m = m_est,
    eset = eset
  )
}

# 0.9 center of mass
center_of_mass <- function(x) {
  w <- pmax(x, 0)
  if (sum(w) == 0) return(NA_real_)
  sum(seq_along(w) * w) / sum(w)
}

# 0.10 membership similarity
membership_similarity <- function(memb_x, memb_y) {
  sum(memb_x * memb_y)
}

# 0.11 direct rule-based classification
# replace p-value-based classify_orthogroup
# preserved:
#   same dominant cluster in all 3 species
#   all 3 pairwise similarities >= sim_cut
#   all max memberships >= memb_cutoff
# x_specific:
#   x differs from y and z, while y~z is similar and same cluster
# diverse:
#   all 3 pairs dissimilar
classify_orthogroup_direct <- function(cx, cy, cz,
                                       maxm_x, maxm_y, maxm_z,
                                       sim_xy, sim_xz, sim_yz,
                                       memb_cutoff = 0.5,
                                       sim_cut = 0.60,
                                       diff_cut = 0.30) {
  if (any(c(maxm_x, maxm_y, maxm_z) <= memb_cutoff)) {
    return(NA_character_)
  }

  # preserved
  if (cx == cy && cx == cz &&
      sim_xy >= sim_cut && sim_xz >= sim_cut && sim_yz >= sim_cut) {
    return("preserved")
  }

  if (sim_xy <= diff_cut && sim_xz <= diff_cut &&
      sim_yz >= sim_cut &&
      cy == cz && cx != cy) {
    return("human_specific")
  }

  if (sim_xy <= diff_cut && sim_yz <= diff_cut &&
      sim_xz >= sim_cut &&
      cx == cz && cy != cx) {
    return("mouse_specific")
  }

  if (sim_xz <= diff_cut && sim_yz <= diff_cut &&
      sim_xy >= sim_cut &&
      cx == cy && cz != cx) {
    return("pig_specific")
  }

  # diverse
  if (sim_xy <= diff_cut && sim_xz <= diff_cut && sim_yz <= diff_cut) {
    return("diverse")
  }

  "intermediate"
}

# 0.12 DTW
dtw_distance <- function(v1, v2, eps = 1e-8) {
  v1 <- as.numeric(v1)
  v2 <- as.numeric(v2)

  if (!all(is.finite(v1)) || !all(is.finite(v2))) {
    return(NA_real_)
  }

  # Exclude near-zero variance rows before scaling.
  sd1 <- sd(v1)
  sd2 <- sd(v2)

  if (is.na(sd1) || is.na(sd2)) return(NA_real_)

  if (sd1 < eps) {
    v1 <- rep(0, length(v1))
  } else {
    v1 <- (v1 - mean(v1)) / sd1
  }

  if (sd2 < eps) {
    v2 <- rep(0, length(v2))
  } else {
    v2 <- (v2 - mean(v2)) / sd2
  }

  out <- tryCatch(
    dtw::dtw(v1, v2, keep = FALSE)$distance,
    error = function(e) NA_real_
  )

  out
}

# 0.13 classify with membership + DTW
classify_orthogroup_joint <- function(
    cx, cy, cz,
    maxm_x, maxm_y, maxm_z,
    sim_xy, sim_xz, sim_yz,
    dtw_xy, dtw_xz, dtw_yz,
    memb_cutoff = 0.5,
    sim_cut = 0.60,
    diff_cut = 0.30,
    dtw_same_cut = 2.5,
    dtw_diff_cut = 4.5,
    require_same_cluster_for_preserved = TRUE
) {
  # 1) membership confidence
  if (any(c(maxm_x, maxm_y, maxm_z) <= memb_cutoff)) {
    return(NA_character_)
  }

  # 2) DTW availability
  if (any(!is.finite(c(dtw_xy, dtw_xz, dtw_yz)))) {
    return(NA_character_)
  }

  same_cluster_all <- (cx == cy && cx == cz)

  # preserved
  # membership high similarity + DTW small
  # optionally keep same-cluster constraint

  #  }
  #}

  n_sim_high <- sum(c(sim_xy, sim_xz, sim_yz) >= sim_cut)
  n_dtw_low  <- sum(c(dtw_xy, dtw_xz, dtw_yz) <= dtw_same_cut)

  if (n_sim_high >= 2 && n_dtw_low >= 2) {
    if (!require_same_cluster_for_preserved || same_cluster_all) {
      return("preserved")
    }
  }

  # human_specific
  # human differs from mouse/pig;
  # mouse and pig similar
  if (sim_xy <= diff_cut && sim_xz <= diff_cut &&
      sim_yz >= sim_cut &&
      dtw_xy >= dtw_diff_cut && dtw_xz >= dtw_diff_cut &&
      dtw_yz <= dtw_same_cut &&
      cy == cz && cx != cy) {
    return("human_specific")
  }

  # mouse_specific
  if (sim_xy <= diff_cut && sim_yz <= diff_cut &&
      sim_xz >= sim_cut &&
      dtw_xy >= dtw_diff_cut && dtw_yz >= dtw_diff_cut &&
      dtw_xz <= dtw_same_cut &&
      cx == cz && cy != cx) {
    return("mouse_specific")
  }

  # pig_specific
  if (sim_xz <= diff_cut && sim_yz <= diff_cut &&
      sim_xy >= sim_cut &&
      dtw_xz >= dtw_diff_cut && dtw_yz >= dtw_diff_cut &&
      dtw_xy <= dtw_same_cut &&
      cx == cy && cz != cx) {
    return("pig_specific")
  }

  # diverse
  # both membership and DTW support divergence
  #    dtw_xy >= dtw_diff_cut && dtw_xz >= dtw_diff_cut && dtw_yz >= dtw_diff_cut) {
  #}

  n_sim_low  <- sum(c(sim_xy, sim_xz, sim_yz) <= diff_cut)
  n_dtw_high <- sum(c(dtw_xy, dtw_xz, dtw_yz) >= dtw_diff_cut)

  if (n_sim_low >= 2 && n_dtw_high >= 2) {
    return("diverse")
  }
  
  "intermediate"
}

# 0.14 optional helper: summarize pairwise DTW to derive empirical cutoffs
summarise_dtw_distribution <- function(res_df) {
  all_dtw <- c(res_df$dtw_hm, res_df$dtw_hc, res_df$dtw_mc)
  all_dtw <- all_dtw[is.finite(all_dtw)]
  quantile(all_dtw, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
}


# A1. equal-width pseudotime binning
# Equal-width pseudotime bins.
bin_pseudotime_equal_width <- function(pt, n_bins = 9) {
  stopifnot(length(pt) > 0)
  out <- cut(pt, breaks = n_bins, include.lowest = TRUE, labels = FALSE)
  paste0("bin", out)
}

# Dynamic genes selected by the variance-to-mean ratio.
# per species, variance-based, no pseudotime-correlation filter
select_dynamic_genes_vmr <- function(mat_gene_bin,
                                               min_mean = 0.5,
                                               vmr_quantile = 0.75,
                                               use_log = FALSE) {
  x <- mat_gene_bin

  if (use_log) {
    x <- log2(x + 1)
  }

  mu <- rowMeans(x)
  va <- apply(x, 1, var)
  vmr <- va / (mu + 1e-6)

  keep_expr <- mu >= min_mean
  vmr_thr <- quantile(vmr[keep_expr], vmr_quantile, na.rm = TRUE)
  selected <- keep_expr & vmr >= vmr_thr

  data.frame(
    gene = rownames(x),
    mean_expr = mu,
    variance = va,
    vmr = vmr,
    selected = selected,
    stringsAsFactors = FALSE
  )
}

# Classify orthogroups by pairwise similarity thresholds.
classify_orthogroup_by_similarity_thresholds <- function(
    cx, cy, cz,
    maxm_x, maxm_y, maxm_z,
    sim_xy, sim_xz, sim_yz,
    memb_cutoff = 0.5,
    diff_pval = 0.05,
    same_pval = 0.5
) {
  # 1. all three ortholog-features must have confident membership
  if (any(c(maxm_x, maxm_y, maxm_z) <= memb_cutoff)) {
    return(NA_character_)
  }

  # 2. species-specific
  # x-specific: x differs from y and z, while y and z are similar and same cluster
  if (sim_xy < diff_pval &&
      sim_xz < diff_pval &&
      sim_yz > same_pval &&
      cy == cz) {
    return("human_specific")
  }

  if (sim_xy < diff_pval &&
      sim_xz > same_pval &&
      sim_yz < diff_pval &&
      cx == cz) {
    return("mouse_specific")
  }

  if (sim_xy > same_pval &&
      sim_xz < diff_pval &&
      sim_yz < diff_pval &&
      cx == cy) {
    return("pig_specific")
  }

  # 3. preserved / conserved
  if (sim_xy > same_pval &&
      sim_xz > same_pval &&
      sim_yz > same_pval &&
      cx == cy && cx == cz) {
    return("preserved")
  }

  # 4. diverged
  if (sim_xy < diff_pval &&
      sim_xz < diff_pval &&
      sim_yz < diff_pval) {
    return("diverse")
  }

  # 5. intermediate / not assigned
  return("intermediate")
}


# A3. similarity classify
estimate_similarity_cutoffs <- function(res_df,
                                        sim_cols = c("sim_xy", "sim_xz", "sim_yz"),
                                        low_q = 0.20,
                                        high_q = 0.80) {
  sim_all <- unlist(res_df[, sim_cols], use.names = FALSE)
  sim_all <- sim_all[is.finite(sim_all)]

  if (length(sim_all) == 0) {
    stop("No finite similarity values found.")
  }

  qs <- stats::quantile(sim_all, probs = c(low_q, 0.5, high_q), na.rm = TRUE)

  out <- list(
    diff_cut = as.numeric(qs[1]),
    mid_cut  = as.numeric(qs[2]),
    sim_cut  = as.numeric(qs[3]),
    low_q    = low_q,
    high_q   = high_q,
    sim_all  = sim_all
  )
  class(out) <- "traj_similarity_cutoffs"
  out
}

classify_orthogroup_direct_quantile <- function(
    cx, cy, cz,
    maxm_x, maxm_y, maxm_z,
    sim_xy, sim_xz, sim_yz,
    memb_cutoff = 0.5,
    sim_cut,
    diff_cut,
    require_same_cluster_for_preserved = TRUE
) {
  # 1. membership confidence
  if (any(c(maxm_x, maxm_y, maxm_z) <= memb_cutoff)) {
    return(NA_character_)
  }

  same_cluster_all <- (cx == cy && cx == cz)

  # 2. preserved
  if (sim_xy >= sim_cut &&
      sim_xz >= sim_cut &&
      sim_yz >= sim_cut) {
    if (!require_same_cluster_for_preserved || same_cluster_all) {
      return("preserved")
    }
  }

  # 3. human-specific
  if (sim_xy <= diff_cut &&
      sim_xz <= diff_cut &&
      sim_yz >= sim_cut &&
      cy == cz && cx != cy) {
    return("human_specific")
  }

  # 4. mouse-specific
  if (sim_xy <= diff_cut &&
      sim_yz <= diff_cut &&
      sim_xz >= sim_cut &&
      cx == cz && cy != cx) {
    return("mouse_specific")
  }

  # 5. pig-specific
  if (sim_xz <= diff_cut &&
      sim_yz <= diff_cut &&
      sim_xy >= sim_cut &&
      cx == cy && cz != cx) {
    return("pig_specific")
  }

  # 6. diverse
  if (sim_xy <= diff_cut &&
      sim_xz <= diff_cut &&
      sim_yz <= diff_cut) {
    return("diverse")
  }

  # 7. intermediate
  return("intermediate")
}
