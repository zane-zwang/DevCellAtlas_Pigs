# Purpose: trajectory conservation.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(Matrix)
  library(RColorBrewer)
  library(viridis)
  library(patchwork)
})

# 0. input
traj_obj <- readRDS('traj_analysis_bundle.rds')

res <- traj_obj$res
expr_cpm_hvg <- traj_obj$expr_cpm_hvg
membership <- traj_obj$membership

outdir <- "trajectory_exploration_output"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

species_levels <- c("mouse", "human", "pig")
species_cols <- c(
  mouse   = "#edab1c",
  human   = "#7DAEE0",
  pig = "#725f97"
)

class_cols <- c(
  preserved        = "#299D8F",
  intermediate     = "#E9C46A",
  diverse          = "#D87659",
  human_specific   = "#7DAEE0",
  mouse_specific   = "#edab1c",
  pig_specific = "#725f97"
)

# 1. helper functions

# z-score by row
scale_by_row <- function(mat) {
  t(scale(t(mat)))
}

# center of mass on expression vector
center_of_mass <- function(x) {
  x <- as.numeric(x)
  # Shift negative values before calculating centers of mass.
  x2 <- x - min(x, na.rm = TRUE)
  if (sum(x2, na.rm = TRUE) == 0) {
    return(NA_real_)
  }
  sum(seq_along(x2) * x2, na.rm = TRUE) / sum(x2, na.rm = TRUE)
}

# Parse gene and species from feature IDs.
split_feature_name <- function(x) {
  tibble(
    feature = x,
    gene = sub("\\|.*$", "", x),
    species = sub("^.*\\|", "", x)
  )
}

# Order clusters by mean center of mass.
get_cluster_order <- function(expr_scaled_list, membership_df) {
  # expr_scaled_list: list(species -> gene x bin)
  # membership_df: res, must contain cx/cy/cz or equivalent cluster columns
  
  long_com <- bind_rows(lapply(names(expr_scaled_list), function(sp) {
    mat <- expr_scaled_list[[sp]]
    tibble(
      gene = rownames(mat),
      species = sp,
      com = apply(mat, 1, center_of_mass)
    )
  }))
  
  cluster_df <- res %>%
    transmute(
      gene,
      human_cluster   = cx,
      mouse_cluster   = cy,
      pig_cluster = cz
    ) %>%
    pivot_longer(
      cols = c(human_cluster, mouse_cluster, pig_cluster),
      names_to = "cluster_species",
      values_to = "cluster"
    ) %>%
    mutate(
      species = case_when(
        cluster_species == "human_cluster"   ~ "human",
        cluster_species == "mouse_cluster"   ~ "mouse",
        cluster_species == "pig_cluster" ~ "pig"
      )
    ) %>%
    select(gene, species, cluster)
  
  out <- long_com %>%
    inner_join(cluster_df, by = c("gene", "species")) %>%
    group_by(cluster) %>%
    summarise(mean_com = mean(com, na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_com) %>%
    mutate(cluster = as.integer(cluster))
  
  out
}

# Group ordered clusters into early, mid and late phases.
assign_cluster_stage_group <- function(cluster_order_df) {
  k <- nrow(cluster_order_df)
  if (k < 3) stop("Need at least 3 clusters to define early/mid/late.")
  
  idx <- seq_len(k)
  stage_group <- case_when(
    idx <= floor(k/3) ~ "Early",
    idx <= ceiling(2*k/3) ~ "Mid",
    TRUE ~ "Late"
  )
  
  cluster_order_df %>%
    mutate(stage_group = stage_group)
}

build_long_expression_df <- function(expr_list, genes_use = NULL, scale_rows = TRUE) {
  mats <- expr_list
  
  if (!is.null(genes_use)) {
    mats <- lapply(mats, function(m) {
      g <- intersect(genes_use, rownames(m))
      m[g, , drop = FALSE]
    })
  }
  
  if (scale_rows) {
    mats <- lapply(mats, scale_by_row)
  }
  
  bind_rows(lapply(names(mats), function(sp) {
    mat <- mats[[sp]]
    as.data.frame(mat) %>%
      rownames_to_column("gene") %>%
      pivot_longer(
        cols = -gene,
        names_to = "bin",
        values_to = "value"
      ) %>%
      mutate(
        species = sp,
        bin_num = readr::parse_number(bin)
      )
  }))
}

get_cluster_order_from_features <- function(expr_scaled_list, feature_info, membership_cutoff = 0.5) {
  long_com <- bind_rows(lapply(names(expr_scaled_list), function(sp) {
    mat <- expr_scaled_list[[sp]]
    tibble(
      gene = rownames(mat),
      species = sp,
      com = apply(mat, 1, center_of_mass)
    )
  }))
  
  out <- long_com %>%
    inner_join(
      feature_info %>%
        select(gene, species, dom_cluster, max_membership),
      by = c("gene", "species")
    ) %>%
    filter(max_membership > membership_cutoff) %>%
    group_by(cluster = dom_cluster) %>%
    summarise(mean_com = mean(com, na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_com)
  
  out
}

# 2. build scaled expression matrices
expr_scaled <- lapply(expr_cpm_hvg, scale_by_row)

# Restrict to genes shared by all three species.
genes_common3 <- Reduce(intersect, lapply(expr_scaled, rownames))
expr_scaled <- lapply(expr_scaled, function(m) m[genes_common3, , drop = FALSE])

# 3. cluster order and stage grouping
cluster_order_df <- get_cluster_order(expr_scaled, res)
cluster_order_df <- assign_cluster_stage_group(cluster_order_df)

write.csv(cluster_order_df,
          file.path(outdir, "mfuzz_cluster_order_and_stage_group.csv"),
          row.names = FALSE)

# Assign gene-level stage groups using the human cluster.
res3 <- res %>%
  left_join(
    cluster_order_df %>%
      rename(cx = cluster,
             human_stage_group = stage_group,
             human_cluster_mean_com = mean_com),
    by = "cx"
  )

saveRDS(res3, file.path(outdir, "traj_classification_summary_with_stage_group.rds"))

# 4. plot 1: trajectory clusters panel

# Identify each gene's dominant cluster in each species.
gene_cluster_long <- res3 %>%
  transmute(
    gene,
    class_observed,
    human_cluster = cx,
    mouse_cluster = cy,
    pig_cluster = cz
  ) %>%
  pivot_longer(
    cols = c(human_cluster, mouse_cluster, pig_cluster),
    names_to = "cluster_species",
    values_to = "cluster"
  ) %>%
  mutate(
    species = case_when(
      cluster_species == "human_cluster"   ~ "human",
      cluster_species == "mouse_cluster"   ~ "mouse",
      cluster_species == "pig_cluster" ~ "pig"
    )
  ) %>%
  select(gene, species, cluster, class_observed)

expr_long <- build_long_expression_df(expr_cpm_hvg, genes_use = genes_common3, scale_rows = TRUE)

plot_df_cluster <- expr_long %>%
  inner_join(gene_cluster_long, by = c("gene", "species")) %>%
  left_join(cluster_order_df, by = "cluster") %>%
  mutate(
    cluster_f = factor(cluster, levels = cluster_order_df$cluster),
    species = factor(species, levels = species_levels)
  )


# Show all genes as a light-gray background.
p_cluster_all <- ggplot(plot_df_cluster, aes(x = bin_num, y = value, group = interaction(gene, species))) +
  geom_line(color = "grey75", alpha = 0.25, linewidth = 0.25) +
  stat_summary(
    aes(color = species, group = species),
    fun = mean,
    geom = "line",
    linewidth = 0.9
  ) +
  facet_wrap(~ cluster_f, scales = "free_y", nrow = 1) +
  scale_color_manual(values = species_cols) +
  labs(
    x = "Pseudotime bins",
    y = "Scaled expression",
    color = "Species",
    title = "Trajectory clusters across species"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    axis.text.x = element_text(angle = 0)
  )

ggsave(
  file.path(outdir, "trajectory_clusters_all_genes.pdf"),
  p_cluster_all, width = 18, height = 3
)

# plot preserved only
plot_df_preserved <- plot_df_cluster %>%
  filter(class_observed == "preserved")

if (nrow(plot_df_preserved) > 0) {
  p_cluster_pres <- ggplot(plot_df_preserved, aes(x = bin_num, y = value, group = interaction(gene, species))) +
    geom_line(color = "grey75", alpha = 0.25, linewidth = 0.25) +
    stat_summary(
      aes(color = species, group = species),
      fun = mean,
      geom = "line",
      linewidth = 1
    ) +
    facet_wrap(~ cluster_f, scales = "free_y", nrow = 1) +
    scale_color_manual(values = species_cols) +
    labs(
      x = "Pseudotime bins",
      y = "Scaled expression",
      color = "Species",
      title = "Preserved trajectories across ordered clusters"
    ) +
    theme_classic(base_size = 12)
  
  ggsave(
    file.path(outdir, "trajectory_clusters_preserved_only.pdf"),
    p_cluster_pres, width = 18, height = 3
  )
}

# Cluster-level trajectory summaries.
feature_info <- as.data.frame(membership) %>%
  rownames_to_column("feature") %>%
  mutate(
    dom_cluster = max.col(across(-feature), ties.method = "first"),
    max_membership = apply(select(., -feature), 1, max)
  ) %>%
  mutate(
    gene = sub("\\|.*$", "", feature),
    species = sub("^.*\\|", "", feature)
  )

cluster_order_df <- get_cluster_order_from_features(expr_scaled, feature_info, membership_cutoff = 0.5)
cluster_order_df <- assign_cluster_stage_group(cluster_order_df)

plot_df_cluster <- expr_long %>%
  inner_join(
    feature_info %>% select(gene, species, dom_cluster, max_membership),
    by = c("gene", "species")
  ) %>%
  filter(max_membership > 0.5) %>%
  rename(cluster = dom_cluster) %>%
  left_join(cluster_order_df, by = "cluster") %>%
  mutate(
    cluster_f = factor(cluster, levels = cluster_order_df$cluster),
    species = factor(species, levels = species_levels)
  )

p_cluster_all <- ggplot(plot_df_cluster, aes(x = bin_num, y = value, group = interaction(gene, species))) +
  geom_line(color = "grey75", alpha = 0.18, linewidth = 0.2) +
  stat_summary(
    aes(color = species, group = species),
    fun = mean,
    geom = "line",
    linewidth = 1
  ) +
  facet_wrap(~ cluster_f, scales = "free_y", nrow = 1) +
  scale_color_manual(values = species_cols) +
  labs(
    x = "Pseudotime bins",
    y = "Scaled expression",
    color = "Species",
    title = "Trajectory clusters across species (confident members only)"
  ) +
  theme_classic(base_size = 12)

ggsave(
  file.path(outdir, "trajectory_clusters_all_genes_confident_members.pdf"),
  p_cluster_all, width = 18, height = 3
)

# 5. plot 2: center of mass violin
com_long <- bind_rows(lapply(names(expr_scaled), function(sp) {
  mat <- expr_scaled[[sp]]
  tibble(
    gene = rownames(mat),
    species = sp,
    com = apply(mat, 1, center_of_mass)
  )
})) %>%
  inner_join(gene_cluster_long, by = c("gene", "species")) %>%
  left_join(cluster_order_df, by = "cluster") %>%
  mutate(
    cluster_f = factor(cluster, levels = cluster_order_df$cluster),
    species = factor(species, levels = species_levels)
  )

p_com <- ggplot(com_long, aes(x = cluster_f, y = com, fill = species)) +
  geom_violin(position = position_dodge(width = 0.8), scale = "width", trim = TRUE, width =  0.5) +
  scale_fill_manual(values = species_cols) +
  labs(
    x = "Trajectory clusters",
    y = "Center of mass",
    fill = "Species",
    title = "Center of mass distribution across ordered clusters"
  ) +
  theme_classic(base_size = 12)

ggsave(
  file.path(outdir, "trajectory_cluster_center_of_mass_violin.pdf"),
  p_com, width = 7, height = 5
)

# 6. plot 3: Early/Mid/Late prop

# class: Preserved / Intermediate / Diverged
res4 <- res3 %>%
  mutate(class_3way = case_when(
    class_observed == "preserved" ~ "Preserved",
    class_observed %in% c("diverse", "human_specific", "mouse_specific", "pig_specific") ~ "Diverged",
    class_observed == "intermediate" ~ "Intermediate",
    TRUE ~ NA_character_
  )) %>%
  left_join(
    cluster_order_df %>%
      rename(cx = cluster, stage_group = stage_group),
    by = "cx"
  )

prop_df <- res4 %>%
  filter(!is.na(class_3way), !is.na(stage_group)) %>%
  count(stage_group, class_3way, name = "n") %>%
  group_by(stage_group) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    stage_group = factor(stage_group, levels = c("Early", "Mid", "Late")),
    class_3way = factor(class_3way, levels = c("Preserved", "Intermediate", "Diverged"))
  )

class3_cols <- c(
  Preserved        = "#299D8F",
  Intermediate     = "#E9C46A",
  Diverged          = "#D87659"
)

p_prop <- ggplot(prop_df, aes(x = stage_group, y = prop, fill = class_3way)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = n),
            position = position_stack(vjust = 0.5),
            color = "white",
            size = 4) +
  scale_fill_manual(values = class3_cols) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = NULL,
    y = "Proportion",
    fill = NULL,
    title = "Composition of trajectory classes across Early / Mid / Late groups"
  ) +
  theme_classic(base_size = 12)

ggsave(
  file.path(outdir, "trajectory_class_composition_early_mid_late.pdf"),
  p_prop, width = 5.2, height = 4.2
)

# Export trajectory classification summaries.
write.csv(res4, file.path(outdir, "trajectory_classification_exploration_table.csv"), row.names = FALSE)

class_count <- res4 %>%
  count(class_observed, sort = TRUE)
write.csv(class_count, file.path(outdir, "trajectory_class_counts.csv"), row.names = FALSE)

stage_class_count <- res4 %>%
  filter(!is.na(stage_group)) %>%
  count(stage_group, class_observed)
write.csv(stage_class_count,
          file.path(outdir, "trajectory_class_counts_by_stage_group.csv"),
          row.names = FALSE)

write.csv(plot_df_cluster, file.path(outdir, "trajectory_cluster.csv"), row.names=F)

message("Done. Outputs saved in: ", outdir)
