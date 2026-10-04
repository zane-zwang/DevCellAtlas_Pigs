# Purpose: immune composition.

library(ggplot2)
library(dplyr)


meta <- readRDS('immune_meta.rds')

colors_tissue <- c('adipose'='#EEA236FF', 'cerebrum'='#357EBDFF', 'duodenum'='#5CB85CFF', 'heart'='#D43F3A99',
                   'hypothalamus'='#46B8DAFF', 'liver'='#20854E99' , 'muscle'='#D43F3AFF')

colors_stage <- c("e55d"="#443a83","e90d"="#31688e", "0d"="#21908c","30d"="#35b779","90d"="#8fd744","180d"="#fde725")

meta$immune_subtype <- factor(meta$immune_subtype, levels = cell_sort)

ggplot(meta,aes(x=immune_subtype,fill=stage))+
  geom_bar(position="fill")+
  scale_fill_manual(values = colors_stage)+
  theme_classic()+
  xlab("")+ylab("Abundance")+
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.title = element_blank(),
  )


library(tidyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(colorRamp2)
library(circlize)
library(grid)

# calculate Ro/e

mat <- table(meta$tissue, meta$immune_subtype)
mat <- table(meta$stage, meta$immune_subtype)
roe <- t(t(mat / rowSums(mat)) / colSums(mat) * sum(mat))


roe_plot <- as.matrix(roe)
cell_sort_use <- intersect(cell_sort, colnames(roe_plot))
roe_plot <- roe_plot[, cell_sort_use, drop = FALSE]

cap <- 4

col_fun <- colorRamp2(
  c(0, 0.5, 1, 2, cap),
  c("white", "#FEE5D9", "#FC9272", "#DE2D26", "#67000D")
)

Heatmap(
  roe_plot,
  name = "R_o/e",
  show_heatmap_legend = TRUE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_side = "right",
  show_column_names = TRUE,
  show_row_names = TRUE,
  col = col_fun,
  rect_gp = gpar(col = "dimgray", lwd = 0.5),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(
    title = "R_o/e",
    at = c(0, 0.5, 1, 2, 4),
    labels = c("0", "0.5", "1", "2", "≥4")
  )
)
