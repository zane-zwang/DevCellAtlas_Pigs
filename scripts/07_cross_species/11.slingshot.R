# Purpose: slingshot.
analysis_dir <- Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input")
publication_root <- normalizePath(Sys.getenv("DEVPIGATLAS_ROOT", "."), mustWork = TRUE)
trajectory_output_dir <- Sys.getenv("OPC_OL_OUTPUT_DIR", file.path(publication_root, "results/cross_species/olig"))
trajectory_seurat_file <- Sys.getenv("OPC_OL_SLINGSHOT_SEURAT_RDS", file.path(trajectory_output_dir, "olig_slingshot_seurat.rds"))
dir.create(trajectory_output_dir, recursive = TRUE, showWarnings = FALSE)

library(slingshot)
library(Seurat)
library(cowplot)
library(ggplot2)
library(Matrix)
library(tradeSeq)
library(RColorBrewer)
library(scales)
library(dplyr)

obj <- readRDS(file.path(analysis_dir, "cross_species/olig/adata_CAMEX_umap_cluster_flt.rds"))
sce <- as.SingleCellExperiment(obj, assay = "RNA")

sce_slingshot <- slingshot(sce,
                     reducedDim = 'UMAP',
                     clusterLabels = sce$celltype,
                     start.clus = 'OPC-1',
                     approx_points = 150)

saveRDS(sce_slingshot, file.path(trajectory_output_dir, 'olig_slingshot.rds'))

SlingshotDataSet(sce_slingshot)

cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

cell_colors <- cell_pal(sce_slingshot$celltype, brewer_pal("qual", "Set2"))

celltype_label <- obj@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(celltype = obj@meta.data$celltype) %>%
  group_by(celltype) %>%
  summarise(UMAP1 = median(UMAP_1),
            UMAP2 = median(UMAP_2))

pdf(file.path(trajectory_output_dir, '3_species_olig_slingshot_cluster.pdf'))
plot(reducedDims(sce_slingshot)$UMAP, col = cell_colors, pch=16, asp = 1, cex = 0.2)
lines(SlingshotDataSet(sce_slingshot), lwd=2, col='black')
dev.off()

pdf(file.path(trajectory_output_dir, "3_species_olig_slingshot_pseudotime.pdf"), width = 6, height = 6)
par(pty = "s", mar = c(4, 4, 1, 1))
plot(reducedDims(sce_slingshot)$UMAP, pch=16, asp = 1, cex = 0.2, col = hcl.colors(100, alpha = 1)[cut(sce_slingshot$slingPseudotime_1, breaks = 100)])
lines(SlingshotDataSet(sce_slingshot), lwd=2, col='black')
dev.off()


pseudotime <- as.data.frame(slingPseudotime(sce_slingshot))
if (!"pseudotime_slingshot" %in% colnames(pseudotime)) {
  if (ncol(pseudotime) != 1L) stop("Specify the intended Slingshot lineage before trajectory classification.")
  colnames(pseudotime) <- "pseudotime_slingshot"
}
obj <- AddMetaData(obj, metadata = pseudotime)
saveRDS(obj, trajectory_seurat_file)

FeaturePlot(obj, features='pseudotime_slingshot') + scale_color_viridis_c() + theme(aspect.ratio = 1)
ggsave(file.path(trajectory_output_dir, '3_species_olig_slingshot_pseudotime.pdf'))
ggsave(file.path(trajectory_output_dir, '3_species_olig_slingshot_pseudotime.png'), dpi=300)

FeaturePlot(obj, features='pseudotime_slingshot', split.by='species') + theme(aspect.ratio = 1) &
  scale_color_viridis_c()
ggsave(file.path(trajectory_output_dir, '3_species_olig_slingshot_pseudotime_split.pdf'), width=21, height=7)
