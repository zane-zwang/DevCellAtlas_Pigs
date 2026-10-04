# Purpose: trajectory classification helpers.

run_trajectory_observed <- function(
    obj,
    species = c("human", "mouse", "pig"),
    pseudotime_col = "pseudotime_slingshot",
    sample_col = "sample",
    species_col = "species",
    n_bins = 9,
    min_cells = 50,
    k = 8,
    dynamic_gene_fun = select_dynamic_genes_vmr,
    dynamic_gene_args = list(min_mean = 0.5, vmr_quantile = 0.75, use_log = FALSE),
    consensus_mode = "intersect3"
) {
  get_meta <- function(obj) obj@meta.data
  get_raw_counts <- function(obj) {
    GetAssayData(object = obj, assay = "RNA", layer = "counts")
  }

  sp_list <- list()
  expr_bin_mean <- list()
  expr_cpm_all <- list()
  dynamic_stats <- list()
  qc_bin_reps <- list()

  for (sp in species) {
    sub <- subset(obj, cells = colnames(obj)[obj[[species_col, drop = TRUE]] == sp])
    sp_list[[sp]] <- list(
      counts = get_raw_counts(sub),
      meta = get_meta(sub)
    )
    sp_list[[sp]]$meta$cell_id <- rownames(sp_list[[sp]]$meta)
  }

  for (sp in names(sp_list)) {
    counts <- sp_list[[sp]]$counts
    meta <- sp_list[[sp]]$meta

    meta$bin <- bin_pseudotime_equal_width(meta[[pseudotime_col]], n_bins = n_bins)

    pb <- pseudobulk_sum(counts, meta, min_cells = min_cells)
    ma <- mean_across_reps(pb$pb, pb$group_df)

    expr_bin_mean[[sp]] <- ma$mat
    qc_bin_reps[[sp]] <- ma$n_reps
    expr_cpm_all[[sp]] <- to_cpm(expr_bin_mean[[sp]], log = FALSE)

    dynamic_stats[[sp]] <- do.call(
      dynamic_gene_fun,
      c(list(mat_gene_bin = expr_cpm_all[[sp]]), dynamic_gene_args)
    )
  }

  genes_use <- build_consensus_gene_set(
    dynamic_stats_list = dynamic_stats,
    mode = consensus_mode
  )

  expr_cpm_hvg <- lapply(expr_cpm_all, function(m) {
    g <- intersect(genes_use, rownames(m))
    m[g, , drop = FALSE]
  })

  genes_common3 <- Reduce(intersect, lapply(expr_cpm_hvg, rownames))
  expr_cpm_hvg <- lapply(expr_cpm_hvg, function(m) m[genes_common3, , drop = FALSE])

  mat_feat_bin <- build_feature_matrix_same_id(expr_cpm_hvg, genes_common3)

  mf <- run_mfuzz(mat_feat_bin, k = k, standardize = TRUE)
  membership <- mf$membership

  get_feat <- function(g, sp) paste0(g, "|", sp)

  res <- lapply(genes_common3, function(g) {
    fx <- get_feat(g, "human")
    fy <- get_feat(g, "mouse")
    fz <- get_feat(g, "pig")

    if (!all(c(fx, fy, fz) %in% rownames(membership))) return(NULL)

    mx <- membership[fx, ]
    my <- membership[fy, ]
    mz <- membership[fz, ]

    vh <- expr_cpm_hvg$human[g, ]
    vm <- expr_cpm_hvg$mouse[g, ]
    vc <- expr_cpm_hvg$pig[g, ]

    data.frame(
      gene = g,
      cx = which.max(mx),
      cy = which.max(my),
      cz = which.max(mz),
      maxm_x = max(mx),
      maxm_y = max(my),
      maxm_z = max(mz),
      sim_xy = membership_similarity(mx, my),
      sim_xz = membership_similarity(mx, mz),
      sim_yz = membership_similarity(my, mz),
      dtw_xy = dtw_distance(vh, vm),
      dtw_xz = dtw_distance(vh, vc),
      dtw_yz = dtw_distance(vm, vc),
      stringsAsFactors = FALSE
    )
  }) %>% dplyr::bind_rows()

  list(
    res = res,
    expr_cpm_hvg = expr_cpm_hvg,
    membership = membership,
    dynamic_stats = dynamic_stats,
    qc_bin_reps = qc_bin_reps,
    params = list(
      n_bins = n_bins,
      min_cells = min_cells,
      k = k,
      consensus_mode = consensus_mode
    )
  )
}

estimate_empirical_thresholds <- function(res_df,
                                          sim_cols = c("sim_xy", "sim_xz", "sim_yz"),
                                          dtw_cols = c("dtw_xy", "dtw_xz", "dtw_yz"),
                                          sim_low_q = 0.20,
                                          sim_high_q = 0.80,
                                          dtw_low_q = 0.20,
                                          dtw_high_q = 0.80) {
  sim_all <- unlist(res_df[, sim_cols], use.names = FALSE)
  sim_all <- sim_all[is.finite(sim_all)]

  dtw_all <- unlist(res_df[, dtw_cols], use.names = FALSE)
  dtw_all <- dtw_all[is.finite(dtw_all)]

  list(
    sim_diff_cut = as.numeric(quantile(sim_all, sim_low_q, na.rm = TRUE)),
    sim_same_cut = as.numeric(quantile(sim_all, sim_high_q, na.rm = TRUE)),
    dtw_same_cut = as.numeric(quantile(dtw_all, dtw_low_q, na.rm = TRUE)),
    dtw_diff_cut = as.numeric(quantile(dtw_all, dtw_high_q, na.rm = TRUE)),
    sim_low_q = sim_low_q,
    sim_high_q = sim_high_q,
    dtw_low_q = dtw_low_q,
    dtw_high_q = dtw_high_q
  )
}

classify_orthogroup_direct_final <- function(
    cx, cy, cz,
    maxm_x, maxm_y, maxm_z,
    sim_xy, sim_xz, sim_yz,
    dtw_xy, dtw_xz, dtw_yz,
    memb_cutoff = 0.5,
    sim_same_cut,
    sim_diff_cut,
    dtw_same_cut,
    dtw_diff_cut,
    require_same_cluster_for_preserved = TRUE
) {
  if (any(c(maxm_x, maxm_y, maxm_z) <= memb_cutoff)) {
    return(NA_character_)
  }

  same_cluster_all <- (cx == cy && cy == cz)

  sim_high <- c(sim_xy, sim_xz, sim_yz) >= sim_same_cut
  sim_low  <- c(sim_xy, sim_xz, sim_yz) <= sim_diff_cut

  dtw_low  <- c(dtw_xy, dtw_xz, dtw_yz) <= dtw_same_cut
  dtw_high <- c(dtw_xy, dtw_xz, dtw_yz) >= dtw_diff_cut

  # preserved
  if (all(sim_high) && all(dtw_low)) {
    if (!require_same_cluster_for_preserved || same_cluster_all) {
      return("preserved")
    }
  }

  # human-specific
  if (sim_xy <= sim_diff_cut && sim_xz <= sim_diff_cut &&
      sim_yz >= sim_same_cut &&
      dtw_xy >= dtw_diff_cut && dtw_xz >= dtw_diff_cut &&
      dtw_yz <= dtw_same_cut &&
      cy == cz && cx != cy) {
    return("human_specific")
  }

  # mouse-specific
  if (sim_xy <= sim_diff_cut && sim_yz <= sim_diff_cut &&
      sim_xz >= sim_same_cut &&
      dtw_xy >= dtw_diff_cut && dtw_yz >= dtw_diff_cut &&
      dtw_xz <= dtw_same_cut &&
      cx == cz && cy != cx) {
    return("mouse_specific")
  }

  # pig-specific
  if (sim_xz <= sim_diff_cut && sim_yz <= sim_diff_cut &&
      sim_xy >= sim_same_cut &&
      dtw_xz >= dtw_diff_cut && dtw_yz >= dtw_diff_cut &&
      dtw_xy <= dtw_same_cut &&
      cx == cy && cz != cx) {
    return("pig_specific")
  }

  # diverse
  if (all(sim_low) && all(dtw_high)) {
    return("diverse")
  }

  return("intermediate")
}


classify_orthogroup_direct_final_wide <- function(
    cx, cy, cz,
    maxm_x, maxm_y, maxm_z,
    sim_xy, sim_xz, sim_yz,
    dtw_xy, dtw_xz, dtw_yz,
    memb_cutoff = 0.5,
    sim_same_cut,
    sim_diff_cut,
    dtw_same_cut,
    dtw_diff_cut,
    require_same_cluster_for_preserved = TRUE
) {
  # 1. membership confidence
  if (any(c(maxm_x, maxm_y, maxm_z) <= memb_cutoff)) {
    return(NA_character_)
  }

  same_cluster_all <- (cx == cy && cy == cz)

  # pairwise evidence counts
  sim_high <- c(sim_xy, sim_xz, sim_yz) >= sim_same_cut
  sim_low  <- c(sim_xy, sim_xz, sim_yz) <= sim_diff_cut

  dtw_low  <- c(dtw_xy, dtw_xz, dtw_yz) <= dtw_same_cut
  dtw_high <- c(dtw_xy, dtw_xz, dtw_yz) >= dtw_diff_cut

  n_sim_high <- sum(sim_high, na.rm = TRUE)
  n_sim_low  <- sum(sim_low,  na.rm = TRUE)
  n_dtw_low  <- sum(dtw_low,  na.rm = TRUE)
  n_dtw_high <- sum(dtw_high, na.rm = TRUE)

  # 2. preserved
  # balanced rule:
  # at least 2/3 high-sim pairs + at least 2/3 low-DTW pairs
  # plus same dominant cluster across all 3 species
  if (n_sim_high >= 2 && n_dtw_low >= 2) {
    if (!require_same_cluster_for_preserved || same_cluster_all) {
      return("preserved")
    }
  }

  # 3. species-specific
  # require one strongly supported similar pair,
  # and the focal species should differ from the other two
  # by majority evidence rather than all-or-none evidence

  # human-specific: mouse ~ pig, human differs
  yz_support <- (sim_yz >= sim_same_cut) + (dtw_yz <= dtw_same_cut)
  x_diff_y   <- (sim_xy <= sim_diff_cut) + (dtw_xy >= dtw_diff_cut)
  x_diff_z   <- (sim_xz <= sim_diff_cut) + (dtw_xz >= dtw_diff_cut)

  if (cy == cz && cx != cy &&
      yz_support >= 1 &&
      (x_diff_y + x_diff_z) >= 3) {
    return("human_specific")
  }

  # mouse-specific: human ~ pig, mouse differs
  xz_support <- (sim_xz >= sim_same_cut) + (dtw_xz <= dtw_same_cut)
  y_diff_x   <- (sim_xy <= sim_diff_cut) + (dtw_xy >= dtw_diff_cut)
  y_diff_z   <- (sim_yz <= sim_diff_cut) + (dtw_yz >= dtw_diff_cut)

  if (cx == cz && cy != cx &&
      xz_support >= 1 &&
      (y_diff_x + y_diff_z) >= 3) {
    return("mouse_specific")
  }

  # pig-specific: human ~ mouse, pig differs
  xy_support <- (sim_xy >= sim_same_cut) + (dtw_xy <= dtw_same_cut)
  z_diff_x   <- (sim_xz <= sim_diff_cut) + (dtw_xz >= dtw_diff_cut)
  z_diff_y   <- (sim_yz <= sim_diff_cut) + (dtw_yz >= dtw_diff_cut)

  if (cx == cy && cz != cx &&
      xy_support >= 1 &&
      (z_diff_x + z_diff_y) >= 3) {
    return("pig_specific")
  }

  # 4. diverse
  # balanced rule:
  # at least 2/3 low-sim pairs + at least 2/3 high-DTW pairs
  if (n_sim_low >= 2 && n_dtw_high >= 2) {
    return("diverse")
  }

  # 5. intermediate
  return("intermediate")
}


bootstrap_resample_obj <- function(obj,
                                   species_col = "species",
                                   sample_col = "sample",
                                   within_group_replace = TRUE) {
  meta <- obj@meta.data
  meta$cell_id <- rownames(meta)

  sampled_cells <- meta %>%
    dplyr::group_by(.data[[species_col]], .data[[sample_col]]) %>%
    dplyr::group_modify(~ {
      n <- nrow(.x)
      idx <- sample(seq_len(n), size = n, replace = within_group_replace)
      .x[idx, , drop = FALSE]
    }) %>%
    dplyr::ungroup() %>%
    dplyr::pull(cell_id)

  obj[, sampled_cells]
}

run_trajectory_bootstrap <- function(
    obj,
    n_boot = 100,
    seed = 123,
    ...
) {
  set.seed(seed)

  boot_list <- vector("list", n_boot)

  for (i in seq_len(n_boot)) {
    message("Bootstrap ", i, "/", n_boot)
    obj_b <- bootstrap_resample_obj(obj)

    obs_b <- tryCatch(
      run_trajectory_observed(obj = obj_b, ...),
      error = function(e) NULL
    )

    if (is.null(obs_b)) {
      boot_list[[i]] <- NULL
      next
    }

    tmp <- obs_b$res
    tmp$boot_id <- i
    boot_list[[i]] <- tmp
  }

  dplyr::bind_rows(boot_list)
}

summarise_bootstrap_stability <- function(boot_res_df) {
  # dominant cluster stability per species
  cluster_stab <- boot_res_df %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      stab_x = max(table(cx)) / dplyr::n(),
      stab_y = max(table(cy)) / dplyr::n(),
      stab_z = max(table(cz)) / dplyr::n(),
      .groups = "drop"
    )

  # class recurrence
  class_stab <- boot_res_df %>%
    dplyr::group_by(gene, class_boot) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup()

  class_major <- class_stab %>%
    dplyr::group_by(gene) %>%
    dplyr::slice_max(order_by = prop, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::rename(class_final = class_boot,
                  class_confidence = prop)

  # bootstrap means of continuous stats
  cont_sum <- boot_res_df %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      sim_xy_boot_mean = mean(sim_xy, na.rm = TRUE),
      sim_xz_boot_mean = mean(sim_xz, na.rm = TRUE),
      sim_yz_boot_mean = mean(sim_yz, na.rm = TRUE),
      dtw_xy_boot_mean = mean(dtw_xy, na.rm = TRUE),
      dtw_xz_boot_mean = mean(dtw_xz, na.rm = TRUE),
      dtw_yz_boot_mean = mean(dtw_yz, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    cluster_stab = cluster_stab,
    class_stab = class_stab,
    class_major = class_major,
    cont_sum = cont_sum
  )
}


finalise_orthogroup_classification <- function(
    observed_res,
    boot_summary,
    min_class_confidence = 0.70,
    min_cluster_stability = 0.60
) {
  out <- observed_res %>%
    dplyr::left_join(boot_summary$cluster_stab, by = "gene") %>%
    dplyr::left_join(boot_summary$class_major, by = "gene") %>%
    dplyr::left_join(boot_summary$cont_sum, by = "gene")

  out <- out %>%
    dplyr::mutate(
      cluster_pattern = dplyr::case_when(
        cx == cy & cy == cz ~ "same_all",
        cx == cy & cz != cx ~ "xy_same",
        cx == cz & cy != cx ~ "xz_same",
        cy == cz & cx != cy ~ "yz_same",
        TRUE ~ "all_diff"
      ),
      min_stability = pmin(stab_x, stab_y, stab_z, na.rm = TRUE),
      is_high_confidence = !is.na(class_final) &
        class_confidence >= min_class_confidence &
        min_stability >= min_cluster_stability
    )

  out
}
