# Purpose: similarity.
publication_root <- normalizePath(Sys.getenv("DEVPIGATLAS_ROOT", "."), mustWork = TRUE)

library(Seurat)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(Matrix)
library(edgeR)
library(ComplexHeatmap)
library(circlize)
source(file.path(publication_root, "scripts/07_cross_species/14.similarity_helpers.R"))

# paths
deg_dir <- "./figE_all_results"

human_deg <- readRDS(file.path(deg_dir, "human_deg.rds"))
mouse_deg <- readRDS(file.path(deg_dir, "mouse_deg.rds"))
pig_deg   <- readRDS(file.path(deg_dir, "pig_deg.rds"))

human_spec_deg <- readRDS(file.path(deg_dir, "human_spec_deg.rds"))
mouse_spec_deg <- readRDS(file.path(deg_dir, "mouse_spec_deg.rds"))
pig_spec_deg   <- readRDS(file.path(deg_dir, "pig_spec_deg.rds"))
conserved_deg  <- readRDS(file.path(deg_dir, "conserved_deg.rds"))

obj <- readRDS(file.path(deg_dir, "downsample_1000_RDS.rds"))

traj_dir <- "."
traj_bundle <- readRDS(file.path(traj_dir, "traj_analysis_bundle.rds"))
traj_table  <- readRDS(file.path(traj_dir, "trajectory_final_classification_table.rds"))

cells_order <- c(
  "astrocyte",
  "endothelial_cell",
  "excitatory_neuron",
  "inhibitory_interneuron",
  "microglia",
  "mural",
  "oligodendrocyte",
  "oligodendrocyte_progenitor_cell"
)

obj$celltype <- factor(obj$celltype, levels = cells_order)
obj$species  <- factor(obj$species, levels = c("human", "mouse", "pig"))

# ----------------------------- Cell-type identity conservation -----------------------------
# 1. construct species × celltype pseudobulk

pb <- make_sample_celltype_pseudobulk(
  obj = obj,
  species_col = "species",
  sample_col = "sample",
  celltype_col = "celltype",
  min_cells = 20
)

pb_mat <- pb$logcpm
pb_info <- pb$info


# 2. calculate the correlation of human-pig and human-mouse on sample-level
cells_order <- c(
  "astrocyte",
  "endothelial_cell",
  "excitatory_neuron",
  "inhibitory_interneuron",
  "microglia",
  "mural",
  "oligodendrocyte",
  "oligodendrocyte_progenitor_cell"
)

marker_genes_all <- unique(c(
  human_deg$gene,
  mouse_deg$gene,
  pig_deg$gene
))

celltype_cor_pairwise <- calc_pairwise_celltype_cor(
  pb_mat = pb_mat,
  pb_info = pb_info,
  celltypes = cells_order,
  gene_use = marker_genes_all,
  method = "spearman"
)

celltype_cor_sample <- celltype_cor_pairwise %>%
  group_by(celltype, comparison, human_sample) %>%
  summarise(
    mean_cor = mean(cor, na.rm = TRUE),
    median_cor = median(cor, na.rm = TRUE),
    n_pairs = n(),
    .groups = "drop"
  )


test_celltype_identity <- function(df, value_col = "median_cor") {
  wide <- df %>%
    select(celltype, human_sample, comparison, all_of(value_col)) %>%
    rename(value = all_of(value_col)) %>%
    pivot_wider(names_from = comparison, values_from = value)
  
  wide %>%
    group_by(celltype) %>%
    summarise(
      n_human_sample = sum(!is.na(human_pig) & !is.na(human_mouse)),
      median_cor_HP = median(human_pig, na.rm = TRUE),
      median_cor_HM = median(human_mouse, na.rm = TRUE),
      delta_HP_minus_HM = median(human_pig - human_mouse, na.rm = TRUE),
      prop_pig_closer = mean((human_pig - human_mouse) > 0, na.rm = TRUE),
      p_value = ifelse(
        n_human_sample >= 3,
        wilcox.test(human_pig, human_mouse, paired = TRUE)$p.value,
        NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"))
}

identity_stats <- test_celltype_identity(celltype_cor_sample)

identity_stats

p1 <- celltype_cor_sample %>%
  mutate(
    celltype = factor(celltype, levels = cells_order),
    comparison = factor(comparison, levels = c("human_mouse", "human_pig"))
  ) %>%
  ggplot(aes(x = comparison, y = median_cor, fill = comparison)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.5) +
  facet_wrap(~ celltype, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c('#edab1c','#725f97')) +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = "Median sample-level correlation",
    title = "Human-centered cell-type identity conservation"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )
p1

celltype_cor_sample %>%
  mutate(
    comparison = factor(
      comparison,
      levels = c("human_mouse", "human_pig"),
      labels = c("Human–mouse", "Human–pig")
    )
  ) %>%
  ggplot(aes(x = comparison, y = median_cor, group = human_sample)) +
  geom_line(alpha = 0.35, linewidth = 0.4) +
  geom_point(size = 1.6, alpha = 0.85, aes(color = comparison)) +
  facet_wrap(~ celltype, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c('#edab1c','#725f97')) +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = "Median sample-level correlation",
    title = "Human-centered cell-type identity conservation"
  ) +
  theme(
    aspect.ratio = 0.75,
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

human_spec_list <- human_spec_deg %>%
  group_by(cluster) %>%
  summarise(
    genes = list(unique(gene)),
    n_genes = n_distinct(gene),
    .groups = "drop"
  )

human_spec_list

human_spec_cor_pairwise <- calc_human_specific_program_cor(
  pb_mat = pb_mat,
  pb_info = pb_info,
  human_spec_list = human_spec_list,
  method = "spearman",
  min_genes = 10
)

human_spec_cor_sample <- human_spec_cor_pairwise %>%
  group_by(celltype, comparison, human_sample, n_genes) %>%
  summarise(
    mean_cor = mean(cor, na.rm = TRUE),
    median_cor = median(cor, na.rm = TRUE),
    n_pairs = n(),
    .groups = "drop"
  )

human_spec_stats <- test_celltype_identity(human_spec_cor_sample)

human_spec_stats

p2 <- human_spec_cor_sample %>%
  mutate(
    celltype = factor(celltype, levels = cells_order),
    comparison = factor(comparison, levels = c("human_mouse", "human_pig"))
  ) %>%
  ggplot(aes(x = comparison, y = median_cor, fill = comparison)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.5) +
  facet_wrap(~ celltype, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c('#edab1c','#725f97')) +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = "Median sample-level correlation",
    title = "Retention of human-specific cell-type programs"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )
p2

celltype_cols = c("astrocyte"="#FED43999", "endothelial_cell"="#709AE199", "excitatory_neuron"="#8A919799", 
                  "inhibitory_interneuron"="#D2AF8199", "microglia"="#FD744699", "mural"="#D5E4A299",
                  "oligodendrocyte"="#197EC099", "oligodendrocyte_progenitor_cell"="#F05C3B99", 
                  "choroid_plexus"="#46732E99", "fibroblast"="#71D0F599", 
                  "intermediate_progenitor_cell"="#37033599","mixed_neuron"="#07514999", "radial_glia"="#C8081399")

human_spec_stats %>%
  mutate(
    celltype = reorder(celltype, delta_HP_minus_HM),
    sig_label = case_when(
      p_adj < 0.05 & delta_HP_minus_HM >= 0.05 ~ "FDR < 0.05, meaningful effect",
      p_adj < 0.05 ~ "FDR < 0.05, small effect",
      TRUE ~ "Not significant"
    )
  ) %>%
  ggplot(aes(x = delta_HP_minus_HM, y = celltype)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(aes(size = n_human_sample, shape = sig_label, color = celltype), alpha = 0.9) +
  scale_color_manual(values = celltype_cols) +
  theme_classic(base_size = 12) +
  labs(
    x = "Δ median correlation\nHuman–pig minus human–mouse",
    y = NULL,
    size = "Human samples",
    shape = NULL,
    title = "Retention of human-specific cell-type programs"
  ) +
  theme(
    aspect.ratio = 1
  )


traj_long_sim <- traj_table %>%
  select(
    gene,
    class_final,
    sim_human_mouse,
    sim_human_pig
  ) %>%
  pivot_longer(
    cols = c(sim_human_mouse, sim_human_pig),
    names_to = "comparison",
    values_to = "similarity"
  ) %>%
  mutate(
    comparison = recode(
      comparison,
      sim_human_mouse = "human_mouse",
      sim_human_pig = "human_pig"
    )
  )

traj_long_dtw <- traj_table %>%
  select(
    gene,
    class_final,
    dtw_human_mouse,
    dtw_human_pig
  ) %>%
  pivot_longer(
    cols = c(dtw_human_mouse, dtw_human_pig),
    names_to = "comparison",
    values_to = "dtw_distance"
  ) %>%
  mutate(
    comparison = recode(
      comparison,
      dtw_human_mouse = "human_mouse",
      dtw_human_pig = "human_pig"
    )
  )


traj_stats <- tibble(
  metric = c("trajectory_similarity", "DTW_distance"),
  median_human_pig = c(
    median(traj_table$sim_human_pig, na.rm = TRUE),
    median(traj_table$dtw_human_pig, na.rm = TRUE)
  ),
  median_human_mouse = c(
    median(traj_table$sim_human_mouse, na.rm = TRUE),
    median(traj_table$dtw_human_mouse, na.rm = TRUE)
  ),
  delta_pig_minus_mouse = c(
    median(traj_table$sim_human_pig - traj_table$sim_human_mouse, na.rm = TRUE),
    median(traj_table$dtw_human_mouse - traj_table$dtw_human_pig, na.rm = TRUE)
  ),
  p_value = c(
    wilcox.test(
      traj_table$sim_human_pig,
      traj_table$sim_human_mouse,
      paired = TRUE
    )$p.value,
    wilcox.test(
      traj_table$dtw_human_pig,
      traj_table$dtw_human_mouse,
      paired = TRUE
    )$p.value
  )
) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

traj_stats

traj_delta <- traj_table %>%
  mutate(
    class_final = as.character(class_final),
    
    delta_similarity = sim_human_pig - sim_human_mouse,
    delta_DTW = dtw_human_mouse - dtw_human_pig,
    
    pig_closer_similarity = delta_similarity > 0,
    pig_closer_DTW = delta_DTW > 0
  )

traj_class_stats <- traj_delta %>%
  group_by(class_final) %>%
  summarise(
    n_gene = n(),
    
    median_delta_similarity = median(delta_similarity, na.rm = TRUE),
    prop_pig_closer_similarity = mean(pig_closer_similarity, na.rm = TRUE),
    p_similarity = ifelse(
      sum(!is.na(delta_similarity)) >= 10,
      wilcox.test(delta_similarity, mu = 0)$p.value,
      NA_real_
    ),
    
    median_delta_DTW = median(delta_DTW, na.rm = TRUE),
    prop_pig_closer_DTW = mean(pig_closer_DTW, na.rm = TRUE),
    p_DTW = ifelse(
      sum(!is.na(delta_DTW)) >= 10,
      wilcox.test(delta_DTW, mu = 0)$p.value,
      NA_real_
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_similarity = p.adjust(p_similarity, method = "BH"),
    p_adj_DTW = p.adjust(p_DTW, method = "BH")
  )

traj_class_stats

p3_dtw <- traj_long_dtw %>%
  mutate(comparison = factor(comparison, levels = c("human_mouse", "human_pig"))) %>%
  ggplot(aes(x = comparison, y = dtw_distance, fill = comparison)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 0.5, alpha = 0.25) +
  scale_fill_manual(values = c('#edab1c','#725f97')) +
  theme_classic(base_size = 13) +
  labs(
    x = NULL,
    y = "DTW distance",
    title = "OPC/OLG developmental trajectory distance"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none",
    aspect.ratio = 0.75
  )

p3_dtw

class_order <- c(
  "preserved",
  "intermediate",
  "diverse",
  "human_specific",
  "pig_specific",
  "mouse_specific"
)

class_cols <- c(
  "preserved"        = "#299D8F",
  "intermediate"     = "#E9C46A",
  "diverse"          = "#D87659",
  "human_specific"   = "#7DAEE0",
  "mouse_specific"   = "#edab1c",
  "pig_specific" = "#725f97"
)

p3_delta_dtw_by_class <- traj_delta %>%
  mutate(class_final = factor(class_final, levels = class_order)) %>%
  filter(!is.na(class_final)) %>%
  ggplot(aes(x = class_final, y = delta_DTW)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_boxplot(aes(fill = class_final), width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.45, alpha = 0.25) +
  scale_fill_manual(values = class_cols) +
  theme_classic(base_size = 13) +
  labs(
    x = NULL,
    y = "Δ DTW proximity\nhuman-mouse distance − human-pig distance",
    title = "Human-pig trajectory proximity across gene classes"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    aspect.ratio = 0.75
  )

p3_delta_dtw_by_class

traj_delta <- traj_table %>%
  mutate(
    class_final = as.character(class_final),
    phase = factor(human_peak_stage_group, levels = c("Early", "Mid", "Late")),
    
    delta_similarity = sim_human_pig - sim_human_mouse,
    delta_DTW = dtw_human_mouse - dtw_human_pig
  )

traj_phase_stats <- traj_delta %>%
  group_by(phase) %>%
  summarise(
    n_gene = n(),
    
    median_delta_similarity = median(delta_similarity, na.rm = TRUE),
    prop_pig_closer_similarity = mean(delta_similarity > 0, na.rm = TRUE),
    p_similarity = ifelse(
      sum(!is.na(delta_similarity)) >= 10,
      wilcox.test(delta_similarity, mu = 0)$p.value,
      NA_real_
    ),
    
    median_delta_DTW = median(delta_DTW, na.rm = TRUE),
    prop_pig_closer_DTW = mean(delta_DTW > 0, na.rm = TRUE),
    p_DTW = ifelse(
      sum(!is.na(delta_DTW)) >= 10,
      wilcox.test(delta_DTW, mu = 0)$p.value,
      NA_real_
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_similarity = p.adjust(p_similarity, method = "BH"),
    p_adj_DTW = p.adjust(p_DTW, method = "BH")
  )

traj_phase_stats

p3_delta_dtw_by_phase <- traj_delta %>%
  filter(!is.na(delta_DTW), !is.na(phase)) %>%
  ggplot(aes(x = phase, y = delta_DTW)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_boxplot(aes(fill = phase), width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.45, alpha = 0.25) +
  scale_fill_manual(values = c("#7f4ea860","#9fd4ed","#fedeab")) +
  theme_classic(base_size = 13) +
  labs(
    x = NULL,
    y = "Δ DTW proximity\nhuman-mouse distance − human-pig distance",
    title = "Human-pig trajectory proximity across developmental phases"
  ) +
  theme(
    aspect.ratio = 0.9
  )

p3_delta_dtw_by_phase


traj_table <- traj_table %>%
  mutate(
    phase = factor(human_peak_stage_group, levels = c("Early", "Mid", "Late"))
  )

traj_long_sim <- traj_table %>%
  select(
    gene,
    phase,
    sim_human_mouse,
    sim_human_pig
  ) %>%
  pivot_longer(
    cols = c(sim_human_mouse, sim_human_pig),
    names_to = "comparison",
    values_to = "similarity"
  ) %>%
  mutate(
    comparison = recode(
      comparison,
      sim_human_mouse = "human_mouse",
      sim_human_pig = "human_pig"
    )
  )

traj_long_dtw <- traj_table %>%
  select(
    gene,
    phase,
    dtw_human_mouse,
    dtw_human_pig
  ) %>%
  pivot_longer(
    cols = c(dtw_human_mouse, dtw_human_pig),
    names_to = "comparison",
    values_to = "dtw_distance"
  ) %>%
  mutate(
    comparison = recode(
      comparison,
      dtw_human_mouse = "human_mouse",
      dtw_human_pig = "human_pig"
    )
  )

ggplot(traj_long_sim, aes(x=comparison, y=similarity, fill=comparison)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 0.5, alpha = 0.25) +
  facet_wrap(~ phase, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c('#edab1c','#725f97')) +
  theme_classic(base_size = 13) +
  labs(
    x = NULL,
    y = "Trajectory similarity",
    title = "Human-pig trajectory proximity across developmental phases"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none",
    aspect.ratio = 0.75
  )

ggplot(traj_long_dtw, aes(x=comparison, y=dtw_distance, fill=comparison)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 0.5, alpha = 0.25) +
  facet_wrap(~ phase, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c('#edab1c','#725f97')) +
  theme_classic(base_size = 13) +
  labs(
    x = NULL,
    y = "DTW distance",
    title = "Human-pig trajectory proximity across developmental phases"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none",
    aspect.ratio = 0.75
  )

test_traj_group <- function(df, value_col = "dtw_distance") {
  wide <- df %>%
    select(phase, gene, comparison, all_of(value_col)) %>%
    rename(value = all_of(value_col)) %>%
    pivot_wider(names_from = comparison, values_from = value) %>%
    mutate(
      delta_pig_closer = if (value_col == "dtw_distance") {
        human_mouse - human_pig
      } else {
        human_pig - human_mouse
      }
    )
  
  wide %>%
    group_by(phase) %>%
    summarise(
      n_gene = sum(!is.na(human_pig) & !is.na(human_mouse)),
      median_HP = median(human_pig, na.rm = TRUE),
      median_HM = median(human_mouse, na.rm = TRUE),
      median_delta_pig_closer = median(delta_pig_closer, na.rm = TRUE),
      prop_pig_closer = mean(delta_pig_closer > 0, na.rm = TRUE),
      p_value = ifelse(
        n_gene >= 3,
        wilcox.test(human_pig, human_mouse, paired = TRUE)$p.value,
        NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"))
}

traj_group_stats_dtw <- test_traj_group(traj_long_dtw, value_col = "dtw_distance")
traj_group_stats_sim <- test_traj_group(traj_long_sim, value_col = "similarity")

traj_table <- traj_table %>%
  mutate(class_final = factor(class_final, 
                              levels = c("preserved","intermediate","diverse",
                                         "human_specific","pig_specific","mouse_specific")))

traj_long_sim <- traj_table %>%
  select(
    gene,
    class_final,
    sim_human_mouse,
    sim_human_pig
  ) %>%
  pivot_longer(
    cols = c(sim_human_mouse, sim_human_pig),
    names_to = "comparison",
    values_to = "similarity"
  ) %>%
  mutate(
    comparison = recode(
      comparison,
      sim_human_mouse = "human_mouse",
      sim_human_pig = "human_pig"
    )
  ) %>%
  filter(!is.na(class_final))

traj_long_dtw <- traj_table %>%
  select(
    gene,
    class_final,
    dtw_human_mouse,
    dtw_human_pig
  ) %>%
  pivot_longer(
    cols = c(dtw_human_mouse, dtw_human_pig),
    names_to = "comparison",
    values_to = "dtw_distance"
  ) %>%
  mutate(
    comparison = recode(
      comparison,
      dtw_human_mouse = "human_mouse",
      dtw_human_pig = "human_pig"
    )
  ) %>%
  filter(!is.na(class_final))

ggplot(traj_long_sim, aes(x=comparison, y=similarity, fill=comparison)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 0.5, alpha = 0.25) +
  facet_wrap(~ class_final, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = c('#edab1c','#725f97')) +
  theme_classic(base_size = 13) +
  labs(
    x = NULL,
    y = "Trajectory similarity",
    title = "Human-pig trajectory proximity across gene classes"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none",
    aspect.ratio = 0.75
  )

ggplot(traj_long_dtw, aes(x=comparison, y=dtw_distance, fill=comparison)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 0.5, alpha = 0.25) +
  facet_wrap(~ class_final, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = c('#edab1c','#725f97')) +
  theme_classic(base_size = 13) +
  labs(
    x = NULL,
    y = "DTW distance",
    title = "Human-pig trajectory proximity across gene classes"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none",
    aspect.ratio = 0.75
  )

test_traj_class <- function(df, value_col = "dtw_distance") {
  wide <- df %>%
    select(class_final, gene, comparison, all_of(value_col)) %>%
    rename(value = all_of(value_col)) %>%
    pivot_wider(names_from = comparison, values_from = value) %>%
    mutate(
      delta_pig_closer = if (value_col == "dtw_distance") {
        human_mouse - human_pig
      } else {
        human_pig - human_mouse
      }
    )
  
  wide %>%
    group_by(class_final) %>%
    summarise(
      n_gene = sum(!is.na(human_pig) & !is.na(human_mouse)),
      median_HP = median(human_pig, na.rm = TRUE),
      median_HM = median(human_mouse, na.rm = TRUE),
      median_delta_pig_closer = median(delta_pig_closer, na.rm = TRUE),
      prop_pig_closer = mean(delta_pig_closer > 0, na.rm = TRUE),
      p_value = ifelse(
        n_gene >= 3,
        wilcox.test(human_pig, human_mouse, paired = TRUE)$p.value,
        NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"))
}

traj_class_stats_dtw <- test_traj_class(traj_long_dtw, value_col = "dtw_distance")
traj_class_stats_sim <- test_traj_class(traj_long_sim, value_col = "similarity")


source(file.path(publication_root, "scripts/07_cross_species/15.robustness_inference.R"))
