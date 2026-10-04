# Purpose: csi modules activity targets.
publication_root <- normalizePath(Sys.getenv("DEVPIGATLAS_ROOT", "."), mustWork = TRUE)
scenic_output_dir <- Sys.getenv("SCENIC_OUTPUT_DIR", file.path(publication_root, "results/scenic"))
scenic_metadata_file <- Sys.getenv("SCENIC_METADATA_CSV", "path/to/input/all_obs.csv")
dir.create(scenic_output_dir, recursive = TRUE, showWarnings = FALSE)

library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(patchwork)
library(ggplot2)
library(stringr)
library(circlize)
library(doParallel)
library(tidyverse)
library(igraph)
library(ggraph)
library(ggnetwork)
library(reshape2)
library(scico)
library(ggsci)

meta <- read.csv(scenic_metadata_file, row.names = 1)
loom <- open_loom(file.path(scenic_output_dir, 'out_SCENIC.loom'))
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")

regulons <- SCENIC::regulonsToGeneLists(regulons_incidMat)

regulonAUC <- SCopeLoomR::get_regulons_AUC(loom,column.attr.name='RegulonsAUC')

close_loom(loom)


# --------- auc heatmap ------------
cellinfo <- meta
cellTypes <-  as.data.frame(subset(cellinfo,select = 'celllineage'))

selectedResolution <- "celllineage"
sub_regulonAUC <- regulonAUC

sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),]
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,selectedResolution])

regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells)
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

saveRDS(regulonActivity_byGroup, file.path(scenic_output_dir, 'regulonActivity_byGroup.rds'))
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),center = T, scale=T))

# ---------- csi ----------------
auc <- getAUC(regulonAUC)
auc <- auc[complete.cases(auc), ]

cor_mat <- cor(t(auc), method = "pearson")

calc_csi <- function(cor_mat, eps = 0.05) {
  n <- ncol(cor_mat)
  csi <- matrix(0, n, n)
  colnames(csi) <- rownames(csi) <- colnames(cor_mat)
  
  for (i in seq_len(n)) {
    for (j in seq_len(i)) {
      if (i == j) {
        csi[i, j] <- 1
      } else {
        thr <- cor_mat[i, j] - eps
        nodeA <- which(cor_mat[i, ] >= thr)
        nodeB <- which(cor_mat[, j] >= thr)
        csi_ij <- 1 - (length(unique(c(nodeA, nodeB))) / n)
        csi[i, j] <- csi_ij
        csi[j, i] <- csi_ij
      }
    }
  }
  return(csi)
}

CSI <- calc_csi(cor_mat, eps = 0.05)

# plot csi df
library(dynamicTreeCut)
D <- as.dist(1 - CSI)
hc <- hclust(D, method = "ward.D2")

modules <- cutreeDynamic(
  dendro = hc,
  distM  = as.matrix(D),
  deepSplit = 2,
  pamRespectsDendro = TRUE,
  minClusterSize = 15
)

module_df <- data.frame(
  regulon = hc$labels,
  module  = paste0("M", modules),
  stringsAsFactors = FALSE
)

# module activity per cell (mean AUC of regulons in module)
module_list <- split(module_df$regulon, module_df$module)

module_auc_cell <- sapply(module_list, function(regs){
  colMeans(auc[regs, , drop = FALSE])
})

# summarize by cell type
module_auc_ct <- as.data.frame(module_auc_cell) %>%
  mutate(celltype = meta[rownames(.), "celltype"]) %>%
  group_by(celltype) %>%
  summarise(across(starts_with("M"), mean, na.rm = TRUE))


hub_by_module <- lapply(module_list, function(regs){
  sub <- CSI[regs, regs, drop = FALSE]
  score <- rowMeans(sub, na.rm = TRUE)
  names(sort(score, decreasing = TRUE))[1:3]
})

ord <- order(module_df$module)
reg_ord <- module_df$regulon[ord]
CSI_ord <- CSI[reg_ord, reg_ord]

mod_vec <- module_df$module[match(reg_ord, module_df$regulon)]
mod_fac <- factor(mod_vec, levels = unique(mod_vec))

col_fun <- colorRamp2(
  c(0, 0.2, 0.45, 0.7, 1),
  c("#fbfaf7", "#eee7d3", "#d0b36f", "#8c6d31", "#2f2415")
)

mod_levels <- levels(mod_fac)
mod_colors <-pal_jama("default")(7)
names(mod_colors) <- mod_levels

ha_row <- rowAnnotation(
  Module = mod_fac,
  col = list(Module = mod_colors),
  show_annotation_name = FALSE
)

ht <- Heatmap(
  CSI_ord,
  name = "CSI",
  col = col_fun,
  cluster_rows = F,
  cluster_columns = F,
  show_row_names = FALSE,
  show_column_names = FALSE,
  width  = unit(140, "mm"),
  height = unit(140, "mm"),
  left_annotation = ha_row,
  heatmap_legend_param = list(
    title = "Connection Specificity Index (CSI)",
    at = c(0, 0.5, 1),
    labels = c("low", "mid", "high")
  )
)

draw(ht)

rle_mod <- rle(as.character(mod_fac))
ends <- cumsum(rle_mod$lengths)
starts <- c(1, head(ends + 1, -1))

decorate_heatmap_body("CSI", {
  for (k in seq_along(starts)) {
    s <- starts[k]; e <- ends[k]
    grid.rect(
      x = unit((s - 0.5) / nrow(CSI_ord), "npc"),
      y = unit(1 - (s - 0.5) / nrow(CSI_ord), "npc"),
      width  = unit((e - s + 1) / nrow(CSI_ord), "npc"),
      height = unit((e - s + 1) / nrow(CSI_ord), "npc"),
      just = c("left", "top"),
      gp = gpar(col = "#E69F00", lwd = 3, fill = NA)
    )
  }
})

# ------------- regulon cross stages ---------------
meta$celllineage_stage <- paste0(meta$celllineage, "_", meta$stage)
csi_helper_file <- file.path(publication_root, "scripts/02_developmental_dynamics/helpers/03.calc_csi_module_activity.R")
source(csi_helper_file)

clusters_df <- data.frame("regulon" = module_df$regulon,
                          "csi_cluster" = module_df$module)

lineage_stage_activity <- calc_csi_module_activity(clusters_df,
                                                   regulonAUC,
                                                   meta, "celllineage_stage")

celllineage <- unique(meta$celllineage)
stage <- unique(meta$stage)
lineage_stage_order <- as.vector(sapply(celllineage, function(c) paste(c, stage, sep = "_")))

lineage_stage_activity <- lineage_stage_activity[, intersect(lineage_stage_order, colnames(lineage_stage_activity)), drop = FALSE]

lineage_stage_activity_scaled <- t(scale(t(lineage_stage_activity), center = TRUE, scale = TRUE))
saveRDS(lineage_stage_activity, file.path(scenic_output_dir, 'lineage_stage_module_activity.rds'))
saveRDS(lineage_stage_activity_scaled, file.path(scenic_output_dir, 'lineage_stage_module_activity_scaled.rds'))

lineage_activity <- calc_csi_module_activity(clusters_df,
                                             regulonAUC,
                                             meta, "celllineage")

lineage_order <- unique(meta$celllineage)

lineage_activity <- lineage_activity[, intersect(lineage_order, colnames(lineage_activity)), drop = FALSE]

lineage_activity_scaled <- t(scale(t(lineage_activity), center = TRUE, scale = TRUE))
saveRDS(lineage_activity, file.path(scenic_output_dir, 'lineage_module_activity.rds'))
saveRDS(lineage_activity_scaled, file.path(scenic_output_dir, 'lineage_module_activity_scaled.rds'))

# ------------ regulon modules downstream --------------

common_regs <- Reduce(intersect, list(rownames(CSI), rownames(auc), module_df$regulon))
CSI_use <- CSI[common_regs, common_regs]
auc_use <- auc[common_regs, , drop = FALSE]
module_df2 <- module_df %>% filter(regulon %in% common_regs)
module_list <- split(module_df2$regulon, module_df2$module)

# module activity
module_auc_cell <- sapply(module_list, function(regs){
  regs2 <- intersect(regs, rownames(auc_use))
  if (length(regs2) == 0) return(rep(NA_real_, ncol(auc_use)))
  colMeans(auc_use[regs2, , drop = FALSE], na.rm = TRUE)
})

module_auc_cell <- as.data.frame(module_auc_cell)
rownames(module_auc_cell) <- colnames(auc_use)

module_mat <- t(as.matrix(module_auc_cell))  # module × cell
cell_annot <- meta[colnames(module_mat), "celltype"]

saveRDS(module_mat, file.path(scenic_output_dir, 'module_mat.rds'))
saveRDS(cell_annot, file.path(scenic_output_dir, 'cell_annot.rds'))

# module cross stages
df_long <- module_auc_cell %>%
  tibble::rownames_to_column("cell") %>%
  left_join(meta %>% tibble::rownames_to_column("cell"), by = "cell") %>%
  pivot_longer(cols = all_of(names(module_list)),
               names_to = "module", values_to = "activity")

module_cl_stage <- df_long %>%
  group_by(celltype, celllineage, stage, module) %>%
  summarise(
    mean_activity   = mean(activity, na.rm = TRUE),
    median_activity = median(activity, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )


stages <- c("e55d", "e90d", "0d", "30d", "90d", "180d")
module_cl_stage$stage <- factor(module_cl_stage$stage, levels = stages)

stage_order <- c("e55d", "e90d", "0d", "30d", "90d", "180d")

plot_df <- module_cl_stage %>%
  mutate(stage = factor(stage, levels = stage_order)) %>%
  group_by(celllineage, module, stage) %>%
  summarise(
    mean_activity = mean(mean_activity, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(plot_df, aes(x = stage, y = mean_activity,
                    group = module, colour = module)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = mod_colors) +
  facet_wrap(~ celllineage, scales = "free_y") +
  theme_bw() +
  theme(
    aspect.ratio = 0.6,
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "stage", y = "mean activity")

sub1 <- subset(plot_df, subset=celllineage=='Stromal')
sub2 <- subset(sub1, subset=module %in% c('M3', 'M4'))

p1 <- ggplot(sub2, aes(x = stage, y = mean_activity,
                    group = module, colour = module)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = mod_colors) +
  theme_bw() +
  theme(
    aspect.ratio = 0.5,
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "stage", y = "Average activity score")

sub1 <- subset(plot_df, subset=celllineage %in% c('Epithelial','Muscle','Stromal'))
sub2 <- subset(sub1, subset=module %in% c('M4'))

lineage_colors <- c('Epithelial'='#5f2d9f','Muscle'='#9f2d2d','Stromal'='#9f5f2d','Endothelial'='#2d9f5f')

p2 <- ggplot(sub2, aes(x = stage, y = mean_activity,
                 group = celllineage, colour = celllineage)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = lineage_colors) +
  theme_bw() +
  theme(
    aspect.ratio = 0.5,
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "stage", y = "Average activity score")

p1 / p2

# target gene enrichment
source(file.path(publication_root, "scripts/02_developmental_dynamics/02.perform_go_kegg_pig.R"))
library(clusterProfiler)
library(org.Ss.eg.db)

# --------- target ------------
reg <- read.csv(file.path(scenic_output_dir, 'reg.csv'))
row1 <- as.character(reg[1, ])
row2 <- as.character(reg[2, ])

new_names <- mapply(function(x, y) {
  if(nchar(x) > 0 & nchar(y) > 0) {
    paste0(x, "_", y)
  } else if(nchar(x) > 0) {
    x
  } else {
    y
  }
}, row1, row2, USE.NAMES = FALSE)

colnames(reg) <- new_names

reg <- reg[-c(1,2), ]


parse_target_info <- function(tf, s) {
  pattern <- "\\('([^']+)',\\s*([0-9\\.eE+-]+)\\)"
  matches <- str_match_all(s, pattern)[[1]]
  if(nrow(matches) == 0) return(data.frame(TF = character(), Gene = character(), Score = numeric()))
  
  data.frame(
    TF = rep(tf, nrow(matches)),
    Gene = matches[,2],
    Score = as.numeric(matches[,3]),
    stringsAsFactors = FALSE
  )
}

TF_Genes_score <- do.call(rbind, lapply(seq_len(nrow(reg)), function(i) {
  parse_target_info(reg$TF[i], reg$TargetGenes[i])
}))
TF_Genes_score = unique(TF_Genes_score)

write.csv(TF_Genes_score, file.path(scenic_output_dir, "All_TF_targets_Genes_score.csv"), row.names = F)


TF_GeneNum = TF_Genes_score %>%
  ungroup() %>%
  group_by(TF) %>%
  summarise(Freq = n(), .groups = "drop")
TF_GeneNum <- as.data.frame(TF_GeneNum)

TF_GeneNum$TF_genenum = paste0(TF_GeneNum$TF, "(", TF_GeneNum$Freq, ")")
write.csv(TF_GeneNum, file.path(scenic_output_dir, "TF_GeneNum_score.csv"), row.names = F)

K <- 50

TF_Genes_filt <- TF_Genes_score %>%
  group_by(TF) %>%
  slice_max(order_by = Score, n = K, with_ties = FALSE) %>%
  ungroup()

write.csv(TF_Genes_filt, file.path(scenic_output_dir, 'TF_Genes_filt.csv'), quote = F)

regulon2targets <- split(TF_Genes_filt$Gene, TF_Genes_filt$TF)
regulon2targets <- lapply(regulon2targets, unique)


hub_table <- module_df2 %>%
  group_by(module) %>%
  group_modify(~{
    regs <- .x$regulon
    sub <- CSI_use[regs, regs, drop = FALSE]
    hub_score <- rowMeans(sub, na.rm = TRUE)
    tibble(regulon = names(hub_score), hub_score = hub_score)
  }) %>%
  ungroup()

hub_TFs <- hub_table %>% 
  group_by(module) %>% 
  summarise(hubs = list(regulon), .groups = "drop")

strip_regulon <- function(x) sub("\\(.*\\)$", "", x)
module_targets <- setNames(
  lapply(hub_TFs$hubs, function(regs){
    tf <- strip_regulon(regs)
    tf <- intersect(tf, names(regulon2targets))
    unique(unlist(regulon2targets[tf]))
  }),
  hub_TFs$module
)

module_targets_df <- stack(module_targets)
colnames(module_targets_df) <- c("gene", "module")

write.csv(module_targets_df, file.path(scenic_output_dir, 'module_targets_df.csv'), quote=F)

# perform enrichment (BP)
BuGn_cols <- brewer.pal(9, "Greens")

bp_res_list <- list()
for (i in unique(module_targets_df$module)) {
  message("Running module: ", i)
  
  bp_res <- perform_enrichment(
    module_targets_df[module_targets_df$module == i, ],
    method = "GO",
    universe_genes = unique(module_targets_df$gene),
    ont = "BP",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    final_padj_cutoff = 0.1,
    simplify_go = TRUE
  )
  
  if (!is.null(bp_res) && nrow(bp_res) > 0) {
    bp_res$module <- i
  }
  
  bp_res_list[[i]] <- bp_res
}

bp_all <- bind_rows(bp_res_list)
bp_top5 <- bp_all %>%
  filter(p.adjust < 0.05) %>%
  group_by(module) %>%
  arrange(p.adjust, desc(Count)) %>%
  slice_head(n = 5) %>%
  ungroup()

bp_top5 <- bp_top5 %>%
  mutate(
    GeneRatio_num = as.numeric(sub("/.*", "", GeneRatio)) /
      as.numeric(sub(".*/", "", GeneRatio)),
    logP = -log10(p.adjust)
  )

module_order <- paste0("M", 1:7)

bp_top5_plot <- bp_top5 %>%
  mutate(module = factor(module, levels = module_order)) %>%
  arrange(module, desc(logP)) %>%
  mutate(
    pathway_module = paste(module, Description, sep = "___"),
    pathway_module = factor(pathway_module, levels = rev(unique(pathway_module)))
  )


ggplot(bp_top5_plot,
       aes(x = module,
           y = pathway_module)) +
  geom_point(aes(size = GeneRatio_num, fill = logP), shape = 21, colour = "grey30") +
  scale_y_discrete(
    labels = function(x) sub("^.*___", "", x)
  ) +
  scale_fill_gradientn(
    colors = BuGn_cols,
    name = "-log10(adj. P)"
  ) +
  scale_size(range = c(2, 8)) +
  labs(
    x = "CSI module",
    y = "KEGG pathway",
    fill = "-log10(adj. P)",
    size = "Gene Ratio"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_blank(),
    aspect.ratio = 3
  )

# perform enrichment (KEGG)
kegg_res_list <- list()

for (i in unique(module_targets_df$module)) {
  message("Running module: ", i)
  
  kegg_res <- perform_enrichment(
    module_targets_df[module_targets_df$module == i, ],
    method = "KEGG",
    universe_genes = unique(module_targets_df$gene),
    organism = "ssc",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    final_padj_cutoff = 0.1
  )
  
  if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
    kegg_res$module <- i
  }
  
  kegg_res_list[[i]] <- kegg_res
}

kegg_all <- bind_rows(kegg_res_list)
kegg_top5 <- kegg_all %>%
  filter(p.adjust < 0.05) %>%
  group_by(module) %>%
  arrange(p.adjust, desc(Count)) %>%
  slice_head(n = 5) %>%
  ungroup()

kegg_top5 <- kegg_top5 %>%
  mutate(
    GeneRatio_num = as.numeric(sub("/.*", "", GeneRatio)) /
      as.numeric(sub(".*/", "", GeneRatio)),
    logP = -log10(p.adjust)
  )


module_order <- paste0("M", 1:7)

kegg_top5_plot <- kegg_top5 %>%
  mutate(module = factor(module, levels = module_order)) %>%
  arrange(module, desc(logP)) %>%
  mutate(
    pathway_module = paste(module, Description, sep = "___"),
    pathway_module = factor(pathway_module, levels = rev(unique(pathway_module)))
  )


ggplot(kegg_top5_plot,
       aes(x = module,
           y = pathway_module)) +
  geom_point(aes(size = GeneRatio_num, fill = logP), shape = 21, colour = "grey30") +
  scale_y_discrete(
    labels = function(x) sub("^.*___", "", x)
  ) +
  scale_fill_gradientn(
    colors = BuGn_cols,
    name = "-log10(adj. P)"
  ) +
  scale_size(range = c(2, 8)) +
  labs(
    x = "CSI module",
    y = "KEGG pathway",
    fill = "-log10(adj. P)",
    size = "Gene Ratio"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_blank(),
    aspect.ratio = 3
  )
