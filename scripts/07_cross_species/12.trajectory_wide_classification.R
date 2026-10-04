# Purpose: trajectory wide classification.
publication_root <- normalizePath(Sys.getenv("DEVPIGATLAS_ROOT", "."), mustWork = TRUE)
source(file.path(publication_root, "scripts/07_cross_species/09.trajectory_helpers.R"))
source(file.path(publication_root, "scripts/07_cross_species/10.trajectory_classification_helpers.R"))
trajectory_output_dir <- Sys.getenv("OPC_OL_OUTPUT_DIR", file.path(publication_root, "results/cross_species/olig"))
trajectory_seurat_file <- Sys.getenv("OPC_OL_SLINGSHOT_SEURAT_RDS", file.path(trajectory_output_dir, "olig_slingshot_seurat.rds"))
dir.create(trajectory_output_dir, recursive = TRUE, showWarnings = FALSE)
obj <- readRDS(trajectory_seurat_file)
if (!"pseudotime_slingshot" %in% colnames(obj@meta.data)) {
  stop("The Slingshot Seurat input lacks pseudotime_slingshot.")
}

# -------- step 1 --------
obs <- run_trajectory_observed(
  obj = obj,
  n_bins = 9,
  min_cells = 30,
  k = 6,
  consensus_mode = "intersect3"
)

res_obs <- obs$res

# -------- step 2 --------
thr <- estimate_empirical_thresholds(
  res_obs,
  sim_low_q = 0.25,
  sim_high_q = 0.75,
  dtw_low_q = 0.25,
  dtw_high_q = 0.75
)

print(thr)

res_obs$class_observed <- mapply(
  classify_orthogroup_direct_final_wide,
  cx = res_obs$cx, cy = res_obs$cy, cz = res_obs$cz,
  maxm_x = res_obs$maxm_x, maxm_y = res_obs$maxm_y, maxm_z = res_obs$maxm_z,
  sim_xy = res_obs$sim_xy, sim_xz = res_obs$sim_xz, sim_yz = res_obs$sim_yz,
  dtw_xy = res_obs$dtw_xy, dtw_xz = res_obs$dtw_xz, dtw_yz = res_obs$dtw_yz,
  MoreArgs = list(
    memb_cutoff = 0.5,
    sim_same_cut = thr$sim_same_cut,
    sim_diff_cut = thr$sim_diff_cut,
    dtw_same_cut = thr$dtw_same_cut,
    dtw_diff_cut = thr$dtw_diff_cut,
    require_same_cluster_for_preserved = TRUE
  )
)

table(res_obs$class_observed, useNA = "ifany")


obs$res <- res_obs
saveRDS(obs, file.path(trajectory_output_dir, "traj_analysis_bundle.rds"))

# Trajectory classification table.
suppressPackageStartupMessages({
  library(tidyverse)
})

traj_obj <- readRDS(file.path(trajectory_output_dir, "traj_analysis_bundle.rds"))

res <- traj_obj$res
expr_cpm_hvg <- traj_obj$expr_cpm_hvg
membership <- traj_obj$membership

species_levels <- c("human", "mouse", "pig")

# 1. row-wise scaling
scale_by_row <- function(mat) {
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  z
}

center_of_mass <- function(x) {
  x <- as.numeric(x)
  x2 <- x - min(x, na.rm = TRUE)
  if (sum(x2, na.rm = TRUE) == 0) return(NA_real_)
  sum(seq_along(x2) * x2, na.rm = TRUE) / sum(x2, na.rm = TRUE)
}

peak_bin <- function(x) {
  x <- as.numeric(x)
  if (all(!is.finite(x))) return(NA_integer_)
  which.max(x)
}

stage_group_from_bin <- function(bin, n_bins = 9) {
  case_when(
    is.na(bin) ~ NA_character_,
    bin <= floor(n_bins / 3) ~ "Early",
    bin <= ceiling(2 * n_bins / 3) ~ "Mid",
    TRUE ~ "Late"
  )
}

expr_scaled <- lapply(expr_cpm_hvg, scale_by_row)

genes_common <- Reduce(intersect, lapply(expr_scaled, rownames))
expr_scaled <- lapply(expr_scaled, function(m) m[genes_common, , drop = FALSE])

# 2. expression trajectory features
traj_features <- bind_rows(lapply(species_levels, function(sp) {
  mat <- expr_scaled[[sp]]
  tibble(
    gene = rownames(mat),
    species = sp,
    com = apply(mat, 1, center_of_mass),
    peak_bin = apply(mat, 1, peak_bin),
    peak_stage_group = stage_group_from_bin(peak_bin, n_bins = 9)
  )
})) %>%
  pivot_wider(
    names_from = species,
    values_from = c(com, peak_bin, peak_stage_group),
    names_glue = "{species}_{.value}"
  )

# 3. membership features
feature_info <- as.data.frame(membership) %>%
  rownames_to_column("feature") %>%
  mutate(
    gene = sub("\\|.*$", "", feature),
    species = sub("^.*\\|", "", feature),
    dominant_cluster = max.col(across(where(is.numeric)), ties.method = "first"),
    max_membership = apply(select(., where(is.numeric)), 1, max)
  ) %>%
  select(gene, species, dominant_cluster, max_membership)

membership_wide <- feature_info %>%
  pivot_wider(
    names_from = species,
    values_from = c(dominant_cluster, max_membership),
    names_glue = "{species}_{.value}"
  )

# 4. normalize column names from res
final_class_table <- res %>%
  rename(
    human_cluster = cx,
    mouse_cluster = cy,
    pig_cluster = cz,
    human_membership = maxm_x,
    mouse_membership = maxm_y,
    pig_membership = maxm_z,
    sim_human_mouse = sim_xy,
    sim_human_pig = sim_xz,
    sim_mouse_pig = sim_yz,
    dtw_human_mouse = dtw_xy,
    dtw_human_pig = dtw_xz,
    dtw_mouse_pig = dtw_yz
  ) %>%
  mutate(
    class_final = ifelse("class_final" %in% colnames(.), class_final, class_observed),
    class_3way = case_when(
      class_final == "preserved" ~ "preserved",
      class_final %in% c("human_specific", "mouse_specific", "pig_specific", "diverse") ~ "diverged",
      class_final == "intermediate" ~ "intermediate",
      TRUE ~ NA_character_
    ),
    max_dtw = pmax(dtw_human_mouse, dtw_human_pig, dtw_mouse_pig, na.rm = TRUE),
    min_dtw = pmin(dtw_human_mouse, dtw_human_pig, dtw_mouse_pig, na.rm = TRUE)
  ) %>%
  left_join(traj_features, by = "gene") %>%
  mutate(
    delta_com_human_mouse = human_com - mouse_com,
    delta_com_human_pig   = human_com - pig_com,
    delta_com_mouse_pig   = mouse_com - pig_com,
    
    human_shift_vs_mouse = case_when(
      delta_com_human_mouse > 1 ~ "human_later",
      delta_com_human_mouse < -1 ~ "human_earlier",
      TRUE ~ "similar_timing"
    ),
    human_shift_vs_pig = case_when(
      delta_com_human_pig > 1 ~ "human_later",
      delta_com_human_pig < -1 ~ "human_earlier",
      TRUE ~ "similar_timing"
    ),
    mouse_shift_vs_pig = case_when(
      delta_com_mouse_pig > 1 ~ "mouse_later",
      delta_com_mouse_pig < -1 ~ "mouse_earlier",
      TRUE ~ "similar_timing"
    )
  )

write.csv(
  final_class_table,
  file.path(trajectory_output_dir, "trajectory_final_classification_table.csv"),
  row.names = FALSE
)

saveRDS(
  final_class_table,
  file.path(trajectory_output_dir, "trajectory_final_classification_table.rds")
)
