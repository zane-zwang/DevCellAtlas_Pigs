# Purpose: tissue program activity.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(patchwork)
  library(RColorBrewer)
})

res_dir <- "result4_tissue_directed_programs"
fig_dir <- file.path(res_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

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

heat_col <- colorRamp2(
  c(-2, 0, 2),
  c("#321653", "white", "#823e1b")
)

source_cols <- c(
  "restricted_subtype_driven" = "#B2182B",
  "shared_subtype_associated" = "#2166AC",
  "mixed_or_intermediate" = "#4D4D4D"
)

tissue_class_cols <- c(
  "tissue_restricted" = "#B2182B",
  "tissue_shared" = "#2166AC",
  "intermediate" = "#999999"
)

lineage_cols <- c(
  "Macrophage" = "#B35806",
  "T_NK" = "#E69253",
  "T-NK" = "#E69253",
  "B_Plasma" = "#41B6C4",
  "B-Plasma" = "#41B6C4",
  "DC" = "#756BB1",
  "Mast" = "#CC79A7"
)

tissue_cols <- c(
  'adipose'='#EEA236FF', 
  'cerebrum'='#357EBDFF', 
  'duodenum'='#5CB85CFF', 
  'heart'='#D43F3A99',
  'hypothalamus'='#46B8DAFF', 
  'liver'='#20854E99' , 
  'muscle'='#D43F3AFF'
)

# ---- import results ----
obj <- readRDS(file.path(res_dir, "immune_with_tissue_gene_cluster_scores.rds"))

gene_cluster_df <- read.table(
  file.path(res_dir, "lineage_level_tissue_gene_clusters.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

cluster_score_df <- read.table(
  file.path(res_dir, "tissue_gene_cluster_scores_by_subtype_tissue.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

program_source <- read.table(
  file.path(res_dir, "tissue_gene_cluster_program_source_classification.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

program_drivers <- read.table(
  file.path(res_dir, "top_subtype_drivers_of_tissue_gene_cluster_programs.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

within_validation <- read.table(
  file.path(res_dir, "within_shared_subtype_cluster_validation.tsv"),
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

# ---- workflow schematic ----
workflow_df <- data.frame(
  step = factor(
    c(
      "Immune lineages",
      "Tissue-associated\nmarkers",
      "Gene clusters\nC1–Ck",
      "Subtype mapping",
      "Shared-subtype\nvalidation"
    ),
    levels = c(
      "Immune lineages",
      "Tissue-associated\nmarkers",
      "Gene clusters\nC1–Ck",
      "Subtype mapping",
      "Shared-subtype\nvalidation"
    )
  ),
  x = 1:5,
  y = 1
)

p_workflow <- ggplot(workflow_df, aes(x = x, y = y)) +
  geom_segment(
    data = data.frame(x = 1:4, xend = 2:5, y = 1, yend = 1),
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.12, "inches")),
    linewidth = 0.4,
    color = "grey35"
  ) +
  geom_label(
    aes(label = step),
    size = 2.8,
    label.size = 0.25,
    label.r = unit(0.08, "lines"),
    fill = "white"
  ) +
  xlim(0.5, 5.5) +
  ylim(0.7, 1.3) +
  theme_void() +
  ggtitle("Tissue-directed program discovery")

ggsave(
  file.path(fig_dir, "Fig4A_workflow.pdf"),
  p_workflow,
  width = 7.2,
  height = 1.2
)

# ---- lineage-level tissue gene cluster heatmap ----
DefaultAssay(obj) <- "RNA"

# Restrict to genes assigned to gene clusters.
genes_use <- intersect(gene_cluster_df$gene, rownames(obj))

obj$lineage_tissue <- paste(obj$immune_lineage, obj$tissue, sep = "|||")

avg_expr <- AverageExpression(
  obj,
  assays = "RNA",
  features = genes_use,
  group.by = "lineage_tissue",
  slot = "data"
)$RNA

avg_mat <- as.matrix(avg_expr)

# Match gene order to the cluster annotation.
gene_cluster_df2 <- gene_cluster_df %>%
  filter(gene %in% rownames(avg_mat)) %>%
  arrange(gene_cluster, gene)

avg_mat <- avg_mat[gene_cluster_df2$gene, , drop = FALSE]

# row z-score
avg_mat_z <- t(scale(t(avg_mat)))
avg_mat_z[is.na(avg_mat_z)] <- 0
avg_mat_z[avg_mat_z > 2.5] <- 2.5
avg_mat_z[avg_mat_z < -2.5] <- -2.5

# column annotation
col_anno_df <- data.frame(
  lineage_tissue = colnames(avg_mat_z)
) %>%
  tidyr::separate(
    lineage_tissue,
    into = c("immune_lineage", "tissue"),
    sep = "\\|\\|\\|",
    remove = FALSE
  )

rownames(col_anno_df) <- col_anno_df$lineage_tissue

top_anno <- HeatmapAnnotation(
  lineage = col_anno_df$immune_lineage,
  tissue = col_anno_df$tissue,
  col = list(
    lineage = lineage_cols,
    tissue = tissue_cols
  ),
  annotation_name_gp = gpar(fontsize = 7),
  annotation_legend_param = list(
    lineage = list(title = "Lineage"),
    tissue = list(title = "Tissue")
  )
)

# Order genes within each cluster by expression pattern.
new_gene_order <- c()
for (clust in unique(gene_cluster_df2$gene_cluster)) {
  genes_in_clust <- gene_cluster_df2$gene[gene_cluster_df2$gene_cluster == clust]
  mat_sub <- avg_mat_z[genes_in_clust, , drop = FALSE]
  
  if (nrow(mat_sub) > 1) {
    dist_sub <- dist(mat_sub, method = "euclidean")
    hc_sub <- hclust(dist_sub, method = "ward.D2")
    genes_sorted <- genes_in_clust[hc_sub$order]
  } else {
    genes_sorted <- genes_in_clust
  }
  
  new_gene_order <- c(new_gene_order, genes_sorted)
}

# Apply the same ordering to the matrix and metadata.
avg_mat_z <- avg_mat_z[new_gene_order, ]
gene_cluster_df2 <- gene_cluster_df2 %>%
  mutate(gene = factor(gene, levels = new_gene_order))


row_split <- gene_cluster_df2$gene_cluster
col_split <- col_anno_df$immune_lineage

heat_col <- colorRamp2(
  c(-2, 0, 2),
  c("#321653", "white", "#823e1b")
)

ht_gene_cluster <- Heatmap(
  avg_mat_z,
  name = "Expression\nz-score",
  col = heat_col,
  top_annotation = top_anno,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = row_split,
  column_split = col_split,
  border = TRUE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 45,
  row_title_gp = gpar(fontsize = 8, fontface = "bold"),
  column_title = "Lineage-level tissue-associated gene clusters",
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  rect_gp = gpar(col = NA, lwd = 0),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6)
  )
)

pdf(file.path(fig_dir, "Fig4B_lineage_level_gene_cluster_heatmap_2.pdf"), width = 7.8, height = 6.5)
draw(ht_gene_cluster)
dev.off()

plot_df <- cluster_score_df %>%
  filter(n_cells >= 50) %>%
  mutate(
    group = paste(immune_subtype, tissue, sep = " | ")
  )

score_mat <- plot_df %>%
  select(gene_cluster, group, mean_score) %>%
  pivot_wider(
    names_from = group,
    values_from = mean_score,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(score_mat) <- score_mat$gene_cluster
score_mat$gene_cluster <- NULL

score_mat <- as.matrix(score_mat)

score_mat_z <- t(scale(t(score_mat)))
score_mat_z[is.na(score_mat_z)] <- 0
score_mat_z[score_mat_z > 2.5] <- 2.5
score_mat_z[score_mat_z < -2.5] <- -2.5

# column annotation: subtype, tissue, tissue_class
group_anno_df <- data.frame(group = colnames(score_mat_z)) %>%
  tidyr::separate(
    group,
    into = c("immune_subtype", "tissue"),
    sep = " \\| ",
    remove = FALSE
  ) %>%
  left_join(subtype_class, by = "immune_subtype")

rownames(group_anno_df) <- group_anno_df$group

top_anno2 <- HeatmapAnnotation(
  tissue = group_anno_df$tissue,
  subtype_class = group_anno_df$tissue_class,
  col = list(
    tissue = tissue_cols,
    subtype_class = tissue_class_cols
  ),
  annotation_name_gp = gpar(fontsize = 7)
)

row_anno_source <- program_source %>%
  select(gene_cluster, program_source) %>%
  distinct()

row_source <- row_anno_source$program_source[
  match(rownames(score_mat_z), row_anno_source$gene_cluster)
]

left_anno <- rowAnnotation(
  source = row_source,
  col = list(source = source_cols),
  annotation_name_gp = gpar(fontsize = 7)
)

ht_cluster_score <- Heatmap(
  score_mat_z,
  name = "Score\nz-score",
  col = heat_col,
  top_annotation = top_anno2,
  left_annotation = left_anno,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8, fontface = "bold"),
  column_names_gp = gpar(fontsize = 4.8),
  column_names_rot = 45,
  column_title = "Tissue gene cluster activity across immune subtypes",
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  rect_gp = gpar(col = "white", lwd = 0.1),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6)
  )
)

pdf(file.path(fig_dir, "Fig4C_cluster_scores_by_subtype_tissue_heatmap.pdf"), width = 9.5, height = 4.8)
draw(ht_cluster_score)
dev.off()

ht_cluster_score
