# Purpose: subtype program enrichment.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(RColorBrewer)
})

# Input / output
res_dir <- "result4_tissue_directed_programs"
fig_dir <- file.path(res_dir, "figures_subtype_enrichment")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

min_cells <- 50
z_cap <- 2.5

# Palettes
lineage_cols <- c(
  "Macrophage" = "#B35806",
  "T_NK"       = "#E69253",
  "T-NK"       = "#E69253",
  "B_Plasma"  = "#41B6C4",
  "B-Plasma"  = "#41B6C4",
  "DC"         = "#756BB1",
  "Mast"       = "#CC79A7"
)

tissue_cols <- c(
  "adipose"      = "#EEA236FF",
  "cerebrum"     = "#357EBDFF",
  "duodenum"     = "#5CB85CFF",
  "heart"        = "#D43F3A99",
  "hypothalamus" = "#46B8DAFF",
  "liver"        = "#20854E99",
  "muscle"       = "#D43F3AFF",
  "skeletal muscle" = "#D43F3AFF"
)

tissue_class_cols <- c(
  "tissue_restricted" = "#B2182B",
  "tissue_shared"     = "#2166AC",
  "intermediate"      = "#999999"
)

source_cols <- c(
  "restricted_subtype_driven" = "#B2182B",
  "shared_subtype_associated" = "#2166AC",
  "mixed_or_intermediate"     = "#4D4D4D"
)

heat_col <- colorRamp2(
  c(-z_cap, 0, z_cap),
  c("#2166AC", "white", "#B2182B")
)

# Helper functions
clip_z <- function(x, cap = 2.5) {
  x[is.na(x)] <- 0
  x[x > cap] <- cap
  x[x < -cap] <- -cap
  x
}

cluster_order_fun <- function(x) {
  x[order(as.integer(sub("^C", "", x)))]
}

theme_pub <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black", size = base_size),
      axis.title = element_text(color = "black", size = base_size + 1),
      axis.line = element_line(linewidth = 0.3, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      panel.border = element_blank()
    )
}

# Load results
cluster_score_df <- read.table(
  file.path(res_dir, "tissue_gene_cluster_scores_by_subtype_tissue.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

subtype_class <- read.table(
  file.path(res_dir, "subtype_tissue_class.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

program_source_file <- file.path(res_dir, "tissue_gene_cluster_program_source_classification.tsv")
program_source <- NULL
if (file.exists(program_source_file)) {
  program_source <- read.table(
    program_source_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
  )
}

gene_cluster_file <- file.path(res_dir, "lineage_level_tissue_gene_clusters.tsv")
gene_cluster_df <- read.table(
  gene_cluster_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

cluster_size_df <- gene_cluster_df %>%
  count(gene_cluster, name = "n_genes")

cluster_levels <- cluster_order_fun(unique(cluster_score_df$gene_cluster))

# Optional: manually edit these labels after checking GO enrichment
module_labels <- setNames(cluster_levels, cluster_levels)

# Example:
# module_labels["C1"]  <- "C1 ER stress"
# module_labels["C2"]  <- "C2 B-cell activation"
# module_labels["C3"]  <- "C3 Antigen presentation"
# module_labels["C12"] <- "C12 TCR signaling"

# 1. Aggregate C1-Ck activity by immune_subtype
#    Weighted mean across tissues; weights = n_cells.
subtype_score <- cluster_score_df %>%
  filter(n_cells >= min_cells) %>%
  group_by(gene_cluster, immune_lineage, immune_subtype) %>%
  summarise(
    weighted_mean_score = weighted.mean(mean_score, w = n_cells, na.rm = TRUE),
    total_cells = sum(n_cells, na.rm = TRUE),
    n_tissues_observed = n_distinct(tissue),
    tissue_with_max_score = tissue[which.max(mean_score)],
    max_tissue_score = max(mean_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(subtype_class, by = "immune_subtype") %>%
  mutate(
    gene_cluster = factor(gene_cluster, levels = cluster_levels)
  )

# z-score within each gene cluster across subtypes
subtype_score <- subtype_score %>%
  group_by(gene_cluster) %>%
  mutate(
    subtype_enrichment_z = as.numeric(scale(weighted_mean_score)),
    subtype_enrichment_z = clip_z(subtype_enrichment_z, z_cap)
  ) %>%
  ungroup()

write.table(
  subtype_score,
  file.path(fig_dir, "C1_Ck_subtype_enrichment_scores.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Define subtype order
preferred_subtype_order <- c(
  "CNS-associated macrophage",
  "Cycling macrophage",
  "IFN-response macrophage",
  "Kupffer-like macrophage",
  "Macrophage-like mixed",
  "Microglia",
  "Monocyte-derived macrophage",
  "Resident macrophage",
  "CD4/activated T cell",
  "CD8 T cell",
  "Cycling T cell",
  "Duodenum-enriched activated T cell",
  "Innate-like T cell",
  "NK cell",
  "B cell",
  "Plasma cell",
  "pDC-like/pre-DC",
  "cDC1-like DC",
  "cDC2/activated DC",
  "Mast cell"
)

subtype_present <- unique(subtype_score$immune_subtype)
subtype_order <- c(
  preferred_subtype_order[preferred_subtype_order %in% subtype_present],
  setdiff(subtype_present, preferred_subtype_order)
)

# 2. Main figure: C1-Ck x subtype heatmap
subtype_mat <- subtype_score %>%
  select(gene_cluster, immune_subtype, subtype_enrichment_z) %>%
  pivot_wider(
    names_from = immune_subtype,
    values_from = subtype_enrichment_z,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(subtype_mat) <- as.character(subtype_mat$gene_cluster)
subtype_mat$gene_cluster <- NULL

subtype_mat <- as.matrix(subtype_mat)

subtype_mat <- subtype_mat[
  cluster_levels[cluster_levels %in% rownames(subtype_mat)],
  subtype_order[subtype_order %in% colnames(subtype_mat)],
  drop = FALSE
]

subtype_anno_df <- subtype_score %>%
  distinct(immune_subtype, immune_lineage, tissue_class) %>%
  right_join(
    tibble(immune_subtype = colnames(subtype_mat)),
    by = "immune_subtype"
  ) %>%
  mutate(
    immune_lineage = ifelse(is.na(immune_lineage), "Unknown", immune_lineage),
    tissue_class = ifelse(is.na(tissue_class), "intermediate", tissue_class)
  )

rownames(subtype_anno_df) <- subtype_anno_df$immune_subtype

top_anno <- HeatmapAnnotation(
  lineage = subtype_anno_df$immune_lineage,
  subtype_class = subtype_anno_df$tissue_class,
  col = list(
    lineage = lineage_cols,
    subtype_class = tissue_class_cols
  ),
  annotation_name_gp = gpar(fontsize = 7),
  annotation_legend_param = list(
    lineage = list(title = "Lineage"),
    subtype_class = list(title = "Subtype class")
  )
)

row_anno_list <- list(
  n_genes = anno_barplot(
    cluster_size_df$n_genes[
      match(rownames(subtype_mat), cluster_size_df$gene_cluster)
    ],
    gp = gpar(fill = "grey55", col = NA),
    border = FALSE,
    axis_param = list(gp = gpar(fontsize = 6))
  )
)

if (!is.null(program_source)) {
  source_vec <- program_source$program_source[
    match(rownames(subtype_mat), program_source$gene_cluster)
  ]
  source_vec[is.na(source_vec)] <- "mixed_or_intermediate"

  row_anno_list$source <- source_vec

  left_anno <- rowAnnotation(
    n_genes = row_anno_list$n_genes,
    source = row_anno_list$source,
    col = list(source = source_cols),
    annotation_name_gp = gpar(fontsize = 7),
    annotation_legend_param = list(
      source = list(title = "Program source")
    )
  )
} else {
  left_anno <- rowAnnotation(
    n_genes = row_anno_list$n_genes,
    annotation_name_gp = gpar(fontsize = 7)
  )
}

row_labels <- module_labels[rownames(subtype_mat)]
row_labels[is.na(row_labels)] <- rownames(subtype_mat)

ht_subtype <- Heatmap(
  subtype_mat,
  name = "Subtype\nenrichment\nz-score",
  col = heat_col,
  top_annotation = top_anno,
  left_annotation = left_anno,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = subtype_anno_df$immune_lineage,
  row_labels = row_labels,
  row_names_gp = gpar(fontsize = 8, fontface = "bold"),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 45,
  column_title = "C1-Ck tissue-directed programs across immune subtypes",
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  rect_gp = gpar(col = "white", lwd = 0.25),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6)
  )
)

pdf(
  file.path(fig_dir, "Fig_C1_Ck_subtype_enrichment_heatmap.pdf"),
  width = 9.2,
  height = 4.2,
  useDingbats = FALSE
)
draw(ht_subtype, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

png(
  file.path(fig_dir, "Fig_C1_Ck_subtype_enrichment_heatmap.png"),
  width = 9.2,
  height = 4.2,
  units = "in",
  res = 600
)
draw(ht_subtype, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# Top subtype enrichments per cluster.
#    This is more compact and often better for main figure panels.
top_n_subtype_per_cluster <- 6

dot_df <- subtype_score %>%
  group_by(gene_cluster) %>%
  arrange(desc(subtype_enrichment_z), desc(total_cells)) %>%
  slice_head(n = top_n_subtype_per_cluster) %>%
  ungroup() %>%
  mutate(
    gene_cluster = factor(as.character(gene_cluster), levels = cluster_levels),
    immune_subtype = factor(immune_subtype, levels = rev(subtype_order)),
    immune_lineage = factor(
      immune_lineage,
      levels = c("Macrophage", "T_NK", "B_Plasma", "DC", "Mast")
    )
  )

p_dot <- ggplot(
  dot_df,
  aes(
    x = gene_cluster,
    y = immune_subtype,
    size = total_cells,
    color = subtype_enrichment_z
  )
) +
  geom_point(alpha = 0.95) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-z_cap, z_cap),
    oob = scales::squish,
    name = "Enrichment\nz-score"
  ) +
  scale_size_continuous(
    range = c(1.2, 4.8),
    breaks = scales::pretty_breaks(n = 4),
    name = "Cells"
  ) +
  labs(
    x = "Tissue-directed gene program",
    y = NULL,
    title = "Subtype enrichment of tissue-directed programs",
    subtitle = paste0("Top ", top_n_subtype_per_cluster, " subtype enrichments per C module")
  ) +
  theme_pub(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_line(linewidth = 0.25, color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave(
  file.path(fig_dir, "Fig_C1_Ck_top_subtype_enrichment_dotplot.pdf"),
  p_dot,
  width = 7.2,
  height = 5.6,
  useDingbats = FALSE
)

ggsave(
  file.path(fig_dir, "Fig_C1_Ck_top_subtype_enrichment_dotplot.png"),
  p_dot,
  width = 7.2,
  height = 5.6,
  dpi = 600
)

# 4. Extended figure: C1-Ck x subtype-tissue heatmap
#    This tells whether enrichment is subtype-wide or tissue-specific.
subtype_tissue_score <- cluster_score_df %>%
  filter(n_cells >= min_cells) %>%
  mutate(
    group = paste(immune_subtype, tissue, sep = " | "),
    gene_cluster = factor(gene_cluster, levels = cluster_levels)
  ) %>%
  group_by(gene_cluster) %>%
  mutate(
    subtype_tissue_enrichment_z = as.numeric(scale(mean_score)),
    subtype_tissue_enrichment_z = clip_z(subtype_tissue_enrichment_z, z_cap)
  ) %>%
  ungroup() %>%
  left_join(
    subtype_class %>% select(immune_subtype, tissue_class),
    by = "immune_subtype"
  )

st_mat <- subtype_tissue_score %>%
  select(gene_cluster, group, subtype_tissue_enrichment_z) %>%
  pivot_wider(
    names_from = group,
    values_from = subtype_tissue_enrichment_z,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(st_mat) <- as.character(st_mat$gene_cluster)
st_mat$gene_cluster <- NULL
st_mat <- as.matrix(st_mat)

st_mat <- st_mat[
  cluster_levels[cluster_levels %in% rownames(st_mat)],
  ,
  drop = FALSE
]

group_anno_df <- data.frame(group = colnames(st_mat)) %>%
  separate(
    group,
    into = c("immune_subtype", "tissue"),
    sep = " \\| ",
    remove = FALSE
  ) %>%
  left_join(
    subtype_score %>% distinct(immune_subtype, immune_lineage),
    by = "immune_subtype"
  ) %>%
  left_join(
    subtype_class %>% select(immune_subtype, tissue_class),
    by = "immune_subtype"
  ) %>%
  mutate(
    immune_lineage = ifelse(is.na(immune_lineage), "Unknown", immune_lineage),
    tissue_class = ifelse(is.na(tissue_class), "intermediate", tissue_class)
  )

# order columns by lineage -> subtype -> tissue
group_order <- group_anno_df %>%
  mutate(
    immune_lineage = factor(
      immune_lineage,
      levels = c("Macrophage", "T_NK", "B_Plasma", "DC", "Mast", "Unknown")
    ),
    immune_subtype = factor(immune_subtype, levels = subtype_order),
    tissue = factor(tissue, levels = names(tissue_cols))
  ) %>%
  arrange(immune_lineage, immune_subtype, tissue) %>%
  pull(group)

group_order <- group_order[group_order %in% colnames(st_mat)]
st_mat <- st_mat[, group_order, drop = FALSE]

group_anno_df <- group_anno_df[match(colnames(st_mat), group_anno_df$group), ]
rownames(group_anno_df) <- group_anno_df$group

top_anno_st <- HeatmapAnnotation(
  lineage = group_anno_df$immune_lineage,
  tissue = group_anno_df$tissue,
  subtype_class = group_anno_df$tissue_class,
  col = list(
    lineage = lineage_cols,
    tissue = tissue_cols,
    subtype_class = tissue_class_cols
  ),
  annotation_name_gp = gpar(fontsize = 7),
  annotation_legend_param = list(
    lineage = list(title = "Lineage"),
    tissue = list(title = "Tissue"),
    subtype_class = list(title = "Subtype class")
  )
)

ht_subtype_tissue <- Heatmap(
  st_mat,
  name = "Subtype-tissue\nenrichment\nz-score",
  col = heat_col,
  top_annotation = top_anno_st,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = group_anno_df$immune_lineage,
  row_labels = row_labels[rownames(st_mat)],
  row_names_gp = gpar(fontsize = 8, fontface = "bold"),
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 4.6),
  column_names_rot = 45,
  column_title = "C1-Ck program activity across subtype-tissue combinations",
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  rect_gp = gpar(col = "white", lwd = 0.15),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6)
  )
)

pdf(
  file.path(fig_dir, "ExtFig_C1_Ck_subtype_tissue_enrichment_heatmap.pdf"),
  width = 13.5,
  height = 4.8,
  useDingbats = FALSE
)
draw(ht_subtype_tissue, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

png(
  file.path(fig_dir, "ExtFig_C1_Ck_subtype_tissue_enrichment_heatmap.png"),
  width = 13.5,
  height = 4.8,
  units = "in",
  res = 600
)
draw(ht_subtype_tissue, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# 5. Optional: identify representative subtype drivers for text
driver_summary <- subtype_score %>%
  group_by(gene_cluster) %>%
  arrange(desc(subtype_enrichment_z), desc(total_cells)) %>%
  summarise(
    top_subtypes = paste(head(immune_subtype, 5), collapse = "; "),
    top_lineages = paste(unique(head(immune_lineage, 5)), collapse = "; "),
    top_tissues = paste(unique(head(tissue_with_max_score, 5)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(as.integer(sub("^C", "", gene_cluster)))

write.table(
  driver_summary,
  file.path(fig_dir, "C1_Ck_top_subtype_driver_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Finished: C1-Ck subtype enrichment figures written to: ", fig_dir)
