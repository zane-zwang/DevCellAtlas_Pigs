# Purpose: immune differential abundance.
publication_root <- normalizePath(Sys.getenv("DEVPIGATLAS_ROOT", "."), mustWork = TRUE)
analysis_dir <- Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input")

# devtools::install_local(file.path(publication_root, "resources/software/miloR-devel.zip"))
library(miloR)
library(Seurat)
library(ggplot2)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(scuttle)

setwd(file.path(analysis_dir, "immune/milo"))
obj <- readRDS(file.path(analysis_dir, "immune/immune_anno_scVI.rds"))
obj@meta.data <- obj@meta.data %>%  
  mutate(group = case_when(  
    stage == "e55d" ~ "stage1",  
    stage == "e90d" ~ "stage1",  
    stage == "0d" ~ "stage1",  
    stage == "30d" ~ "stage2",  
    stage == "90d" ~ "stage2",  
    stage == "180d" ~ "stage2",  
    #TRUE ~ as.numeric(as.character(group))  
  ))

obj <- obj %>% RunPCA(npcs = 50, verbose = F)
sce1 = as.SingleCellExperiment(obj)


milo_obj = Milo(sce1)
milo_obj = buildGraph(milo_obj, k = 30, d = 30,reduced.dim = "X_SCVI")
milo_obj = makeNhoods(milo_obj, prop = 0.1, k = 30, d=30, refined = TRUE)
milo_obj = countCells(milo_obj, meta.data = as.data.frame(colData(milo_obj)), sample="sample")

sample_order = colnames(nhoodCounts(milo_obj))

exp_design <- data.frame(colData(milo_obj))[,c("tissue", "group", "sample", "stage")]

exp_design = distinct(exp_design)
rownames(exp_design) = exp_design$sample
exp_design = exp_design[sample_order,]
exp_design
saveRDS(exp_design,'exp_design.rds')

milo_obj = calcNhoodDistance(milo_obj, d=30, reduced.dim = "X_SCVI")
saveRDS(milo_obj,'milo_obj.rds')


da_results = testNhoods(milo_obj, design = ~ group, design.df = exp_design, reduced.dim='X_UMAP')
saveRDS(da_results,'da_results_2.rds')


da_results %>%
  arrange(SpatialFDR) %>%
  head() 

milo_obj = buildNhoodGraph(milo_obj)
colData(milo_obj)
saveRDS(milo_obj,'milo_obj_buildNhoodGraph.rds')

umap_pl1 <- plotReducedDim(milo_obj, dimred = "X_UMAP", colour_by="sample", text_by = "sample", text_size = 3, point_size=0.5) +  guides(fill="none")
umap_pl2 <- plotReducedDim(milo_obj, dimred = "X_UMAP", colour_by="immune_subtype", text_by = "immune_subtype", text_size = 3, point_size=0.5) +  guides(fill="none")
umap_pl3 <- plotReducedDim(milo_obj, dimred = "X_UMAP", colour_by="group", text_by = "group", text_size = 3, point_size=0.5) +  guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(milo_obj, da_results, layout="X_UMAP",alpha=0.1) + scale_fill_gradient2(name = "logFC", low = "#913a35", high = "#6a61ae") 

nh_plt_comb = umap_pl1+umap_pl2+umap_pl3 + nh_graph_pl + plot_layout(guides="collect")
ggsave(plot=nh_plt_comb, filename='nh_plt_comb.png', width=14.5, height=10, dpi=300)
ggsave(plot=nh_plt_comb, filename='nh_plt_comb.pdf', width=14.5, height=10, dpi=300)

da_results <- annotateNhoods(milo_obj, da_results, coldata_col = "immune_subtype")

#da_results$celltype <- ifelse(da_results$immune_subtype_fraction < 0.7, "Mixed", da_results$immune_subtype)

plotDAbeeswarm(da_results, group.by = "immune_subtype")
ggsave('swarm.pdf', height=9, width=9, dpi=300)


colors <- c('#1f77b4', '#aec7e8', '#AA336A', '#ff9896', '#ff7f0e', '#ffbb78', '#2ca02c', 
'#90EE90', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', 
'#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5')

Idents(object = obj) <- "immune_subtype"

DimPlot(obj, reduction = "X_umap", cols=colors, alpha = 0.75, pt.size = 0.2, shuffle = TRUE, label = TRUE) + 
theme(text=element_text(size=20, color="black"), aspect.ratio=1)
ggsave("immune_celltypes.png", width = 14, height = 7)
