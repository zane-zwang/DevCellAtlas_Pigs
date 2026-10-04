# Purpose: neighborhood abundance.

library(miloR)
library(Seurat)
library(ggplot2)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(scuttle)


da_results <- readRDS('da_results_2.rds')
milo_obj <- readRDS('milo_obj_buildNhoodGraph.rds')

da_results %>%
  arrange(SpatialFDR) %>%
  head() 

ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1)

plotNhoodGraphDA(milo_obj, da_results, layout="X_UMAP",alpha=0.1) + scale_fill_gradient2(name = "logFC", low = "#913a35", high = "#6a61ae")

da_results <- annotateNhoods(milo_obj, da_results, coldata_col = "immune_major")


da_results <- readRDS('da_results_3.rds')
plotDAbeeswarm(da_results, group.by = "immune_subtype")

library(ggbeeswarm)

da.res <- da_results

cell_sort <- c(
  "CNS-associated macrophage", "Cycling macrophage", "IFN-response macrophage",
  "Kupffer-like macrophage", "Macrophage-like mixed", "Microglia",
  "Monocyte-derived macrophage", "Resident macrophage",
  "CD4/activated T cell", "CD8 T cell", "Cycling T cell",
  "Duodenum-enriched activated T cell", "Innate-like T cell", "NK cell",
  "B cell", "Plasma cell",
  "pDC-like/pre-DC", "cDC1-like DC", "cDC2/activated DC", "Mast cell"
)

group.by <- "immune_subtype"
alpha <- 0.1
subset.nhoods <- NULL

if (!is.null(group.by)) {
  if (!group.by %in% colnames(da.res)) {
    stop(
      group.by,
      " is not a column in da.res. Have you forgotten to run annotateNhoods(x, da.res, ",
      group.by,
      ")?"
    )
  }
  da.res <- da.res %>%
    mutate(group_by = .data[[group.by]])
} else {
  da.res <- da.res %>%
    mutate(group_by = "g1")
}

if (!is.null(subset.nhoods)) {
  da.res <- da.res[subset.nhoods, ]
}

plot_df <- da.res %>%
  mutate(
    is_signif = SpatialFDR < alpha,
    logFC_color = ifelse(is_signif, logFC, NA_real_),
    immune_subtype = factor(immune_subtype, levels = cell_sort),
    group_by = factor(group_by, levels = cell_sort)
  ) %>%
  arrange(group_by, Nhood)


beeswarm_pos <- ggplot_build(
  plot_df %>%
    ggplot(aes(group_by, logFC)) +
    geom_quasirandom()
)

plot_df$pos_x <- beeswarm_pos$data[[1]]$x
plot_df$pos_y <- beeswarm_pos$data[[1]]$y

n_groups <- length(levels(plot_df$group_by))

ggplot(plot_df, aes(pos_x, pos_y, color = logFC_color)) +
  geom_point(size = 1.2) +
  scale_color_gradient2() +
  guides(color = "none") +
  xlab(group.by) +
  ylab("Log Fold Change") +
  scale_x_continuous(
    breaks = seq_len(n_groups),
    labels = levels(plot_df$group_by)
  ) +
  theme_bw(base_size = 20) +
  theme(
    strip.text.y = element_text(angle = 0),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )
