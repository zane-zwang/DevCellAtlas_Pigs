# Purpose: identify marker signatures.
analysis_dir <- Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input")

library(Seurat)
library(dplyr)
library(purrr)
library(ComplexHeatmap)
library(circlize)
library(ggsci)

# ------------- functions -----------------
## add_average_expression
add_average_expression <- function(obj,
                                   layer = "data",
                                   colname = "average_expression") {
  mat <- GetAssayData(obj, layer = layer)
  avg <- Matrix::colMeans(mat)
  obj <- AddMetaData(obj, metadata = avg, col.name = colname)
  return(obj)
}

## downsample
downsample_by_group <- function(obj,
                                group_var,
                                max_cells = 1000) {
  meta <- obj@meta.data
  if (!group_var %in% colnames(meta)) {
    stop("group_var not found in meta.data: ", group_var)
  }
  
  meta$cell_id <- rownames(meta)
  
  ds_meta <- meta %>%
    group_by(.data[[group_var]]) %>%
    group_modify(~ {
      n <- nrow(.x)
      if (n <= max_cells) {
        .x
      } else {
        dplyr::slice_sample(.x, n = max_cells)
      }
    }) %>%
    ungroup()
  
  obj_ds <- subset(obj, cells = ds_meta$cell_id)
  return(obj_ds)
}

## run deg
run_deg_by_species <- function(obj,
                               species_var = "species",
                               species_levels = NULL,
                               ident_var = "celltype",
                               only.pos = TRUE,
                               min.pct = 0.1,
                               logfc.threshold = 0.5,
                               padj_cutoff = 0.05,
                               ...) {
  meta <- obj@meta.data
  if (!species_var %in% colnames(meta)) {
    stop("species_var not found in meta.data: ", species_var)
  }
  if (!ident_var %in% colnames(meta)) {
    stop("ident_var not found in meta.data: ", ident_var)
  }
  
  Idents(obj) <- obj[[ident_var]][, 1]
  
  if (is.null(species_levels)) {
    species_levels <- unique(meta[[species_var]])
  }
  
  deg_list <- list()
  for (sp in species_levels) {
    message("Running DEG for species: ", sp)
    cells_use <- rownames(meta)[meta[[species_var]] == sp]
    obj_sp <- subset(obj, cells = cells_use)
    
    deg <- FindAllMarkers(
      obj_sp,
      only.pos        = only.pos,
      min.pct         = min.pct,
      logfc.threshold = logfc.threshold,
      ...
    )
    
    if (!"p_val_adj" %in% colnames(deg)) {
      warning("p_val_adj not found in DEG result for ", sp, ", returning unfiltered.")
    } else {
      deg <- deg[deg$p_val_adj < padj_cutoff, ]
    }
    
    deg_list[[sp]] <- deg
  }
  
  return(deg_list)
}

get_conserved_deg <- function(deg_list) {
  if (length(deg_list) < 2) {
    stop("Need at least 2 species to compute conserved DEG.")
  }
  conserved <- reduce(
    deg_list,
    ~ dplyr::inner_join(.x, .y, by = c("cluster", "gene"))
  )
  conserved <- conserved %>% dplyr::distinct(gene, .keep_all = TRUE)
  return(conserved)
}


# ------------------ get species specific genes --------------------
data <- readRDS(file.path(analysis_dir, "cross_species/camex/adata_CAMEX_umap.rds"))

cells_order <- c('astrocyte', 'endothelial_cell', 'excitatory_neuron', 'inhibitory_interneuron', 'microglia', 'mural', 'oligodendrocyte', 'oligodendrocyte_progenitor_cell')
cells_keep <- cells_order
keep_cells <- data$celltype %in% cells_keep
data <- subset(data, cells = colnames(data)[keep_cells])

data$celltype_FIG1_species <- paste0(data$species, "_", data$celltype)
data <- data %>% NormalizeData() 

data <- add_average_expression(data, layer = "data", colname = "average_expression")

downsample_1000_RDS <- downsample_by_group(
  data,
  group_var = "celltype_FIG1_species",
  max_cells = 1000
)


downsample_1000_RDS$celltype <- factor(
  downsample_1000_RDS$celltype,
  levels = cells_order
)

downsample_1000_RDS <- ScaleData(
  downsample_1000_RDS,
  features = rownames(downsample_1000_RDS)
)

deg_list <- run_deg_by_species(
  downsample_1000_RDS,
  species_var      = "species",
  species_levels   = c("human", "mouse", "pig"),
  ident_var        = "celltype",
  only.pos         = TRUE,
  min.pct          = 0.1,
  logfc.threshold  = 0.5,
  padj_cutoff      = 0.05
)

human_deg   <- deg_list[["human"]]
mouse_deg   <- deg_list[["mouse"]]
pig_deg <- deg_list[["pig"]]

conserved_deg <- get_conserved_deg(deg_list)


# --------------------- functions ---------------------
## order cell
order_cells_by_type_and_avg <- function(obj,
                                        celltype_var   = "celltype",
                                        avg_expr_col   = "average_expression",
                                        cells_order    = NULL) {
  meta <- obj@meta.data
  
  if (!celltype_var %in% colnames(meta)) {
    stop("celltype_var not found: ", celltype_var)
  }
  if (!avg_expr_col %in% colnames(meta)) {
    stop("avg_expr_col not found: ", avg_expr_col)
  }
  
  meta$cell_id <- rownames(meta)
  
  if (!is.null(cells_order)) {
    meta[[celltype_var]] <- factor(meta[[celltype_var]], levels = cells_order)
  }
  
  meta <- meta %>%
    arrange(.data[[celltype_var]]) %>%
    group_by(.data[[celltype_var]]) %>%
    arrange(desc(.data[[avg_expr_col]]), .by_group = TRUE) %>%
    ungroup()
  
  return(meta$cell_id)
}

get_scaled_matrix <- function(obj,
                              genes,
                              ordered_cells = NULL) {
  mat <- GetAssayData(obj, layer = "scale.data")
  genes <- intersect(genes, rownames(mat))
  if (length(genes) == 0) {
    stop("No genes found in object for the given gene list.")
  }
  mat <- mat[genes, , drop = FALSE]
  
  if (!is.null(ordered_cells)) {
    ordered_cells <- intersect(ordered_cells, colnames(mat))
    mat <- mat[, ordered_cells, drop = FALSE]
  }
  return(mat)
}

## get_species_specific_deg
get_species_specific_deg <- function(target_deg,
                                     other_deg_list) {
  other_genes <- unique(unlist(lapply(other_deg_list, function(x) x$gene)))
  
  spec <- target_deg[!(target_deg$gene %in% other_genes), ]
  spec <- spec %>% dplyr::distinct(gene, .keep_all = TRUE)
  return(spec)
}

human_obj   <- subset(downsample_1000_RDS, species == "human")
mouse_obj   <- subset(downsample_1000_RDS, species == "mouse")
pig_obj <- subset(downsample_1000_RDS, species == "pig")

human_cells_order   <- order_cells_by_type_and_avg(human_obj,   "celltype", "average_expression", cells_order)
mouse_cells_order   <- order_cells_by_type_and_avg(mouse_obj,   "celltype", "average_expression", cells_order)
pig_cells_order <- order_cells_by_type_and_avg(pig_obj, "celltype", "average_expression", cells_order)

## conserved genes
genes_conserved <- conserved_deg$gene

human_conserved_mat   <- get_scaled_matrix(human_obj,   genes_conserved, human_cells_order)
mouse_conserved_mat   <- get_scaled_matrix(mouse_obj,   genes_conserved, mouse_cells_order)
pig_conserved_mat     <- get_scaled_matrix(pig_obj, genes_conserved, pig_cells_order)
 
## specfic genes
human_spec_deg <- get_species_specific_deg(
  human_deg,
  list(mouse_deg, pig_deg)
)

mouse_spec_deg <- get_species_specific_deg(
  mouse_deg,
  list(human_deg, pig_deg)
)

pig_spec_deg <- get_species_specific_deg(
  pig_deg,
  list(human_deg, mouse_deg)
)

human_spec_genes <- human_spec_deg$gene

human_spec_mat_human   <- get_scaled_matrix(human_obj,   human_spec_genes, human_cells_order)
human_spec_mat_mouse   <- get_scaled_matrix(mouse_obj,   human_spec_genes, mouse_cells_order)
human_spec_mat_pig <- get_scaled_matrix(pig_obj, human_spec_genes, pig_cells_order)

mouse_spec_genes <- mouse_spec_deg$gene

mouse_spec_mat_mouse <- get_scaled_matrix(
  mouse_obj,
  genes         = mouse_spec_genes,
  ordered_cells = mouse_cells_order
)

mouse_spec_mat_human <- get_scaled_matrix(
  human_obj,
  genes         = mouse_spec_genes,
  ordered_cells = human_cells_order
)

mouse_spec_mat_pig <- get_scaled_matrix(
  pig_obj,
  genes         = mouse_spec_genes,
  ordered_cells = pig_cells_order
)

pig_spec_genes <- pig_spec_deg$gene

pig_spec_mat_pig <- get_scaled_matrix(
  pig_obj,
  genes         = pig_spec_genes,
  ordered_cells = pig_cells_order
)

pig_spec_mat_human <- get_scaled_matrix(
  human_obj,
  genes         = pig_spec_genes,
  ordered_cells = human_cells_order
)

pig_spec_mat_mouse <- get_scaled_matrix(
  mouse_obj,
  genes         = pig_spec_genes,
  ordered_cells = mouse_cells_order
)


n_human_specific <- length(unique(human_spec_deg$gene))
n_mouse_specific <- length(unique(mouse_spec_deg$gene))
n_pig_specific   <- length(unique(pig_spec_deg$gene))
n_conserved <- length(unique(conserved_deg$gene))

n_gene <- data.frame(
  Category = c("Human_specific", "Mouse_specific", "Pig_specific", "Conserved"),
  Count    = c(n_human_specific, n_mouse_specific, n_pig_specific, n_conserved)
)

# ----------------------- plot heatmap --------------------------

celltype_colors <- pal_simpsons("springfield", alpha = 0.6)(length(cells_order))
names(celltype_colors) <- cells_order

# human
human_col_annot_df <- data.frame(
  celltype = human_obj$celltype[human_cells_order]
)
rownames(human_col_annot_df) <- human_cells_order

ha_human <- HeatmapAnnotation(
  celltype = human_col_annot_df$celltype,
  col = list(celltype = celltype_colors),
  show_annotation_name = FALSE,
  show_legend = TRUE, which = "row"
)

mouse_col_annot_df <- data.frame(
  celltype = mouse_obj$celltype[mouse_cells_order]
)
rownames(mouse_col_annot_df) <- mouse_cells_order

ha_mouse <- HeatmapAnnotation(
  celltype = mouse_col_annot_df$celltype,
  col = list(celltype = celltype_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE, which = "row"
)

pig_col_annot_df <- data.frame(
  celltype = pig_obj$celltype[pig_cells_order]
)
rownames(pig_col_annot_df) <- pig_cells_order

ha_pig <- HeatmapAnnotation(
  celltype = pig_col_annot_df$celltype,
  col = list(celltype = celltype_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE, which = "row"
)


## plot
col_fun_human <- colorRamp2(c(0, 2, 4), c("white", "#8EC5F4", "#066ed5"))
col_fun_pig  <- colorRamp2(c(0, 2, 4), c("white", "#B8A6E0", "#5922c7"))
col_fun_mouse  <- colorRamp2(c(0, 2, 4), c("white", "#F0C46A", "#d49306"))
col_fun_con <- colorRamp2(c(0, 2, 4), c("white", "#9E9E9E", "#050505"))

# conserved
human_conserverd_Heatmap <- Heatmap(t(human_conserved_mat), show_row_names = FALSE,
        col = col_fun_con, show_column_names = FALSE, cluster_rows = FALSE,
        cluster_columns = FALSE, show_heatmap_legend = TRUE,
        left_annotation = ha_human)

mouse_conserverd_Heatmap <- Heatmap(t(mouse_conserved_mat), show_row_names = FALSE,
        col = col_fun_con, show_column_names = FALSE, cluster_rows = FALSE,
        cluster_columns = FALSE, show_heatmap_legend = TRUE,
        left_annotation = ha_mouse)

pig_conserverd_Heatmap <- Heatmap(t(pig_conserved_mat), show_row_names = FALSE,
        col = col_fun_con, show_column_names = FALSE, cluster_rows = FALSE,
        cluster_columns = FALSE, show_heatmap_legend = TRUE,
        left_annotation = ha_pig)

merged_conserverd_Heatmap <- human_conserverd_Heatmap %v% mouse_conserverd_Heatmap %v% pig_conserverd_Heatmap

pdf("figE_all_results/conserved_gene.pdf",width=8,height=15)
draw(merged_conserverd_Heatmap)
dev.off()

# human specific
human_specific_Heatmap_human <- Heatmap(t(human_spec_mat_human), show_row_names = FALSE,
        col = col_fun_human, show_column_names = FALSE, cluster_rows = FALSE,
        cluster_columns = FALSE, show_heatmap_legend = TRUE,
        left_annotation = ha_human)

human_specific_Heatmap_mouse <- Heatmap(t(human_spec_mat_mouse), show_row_names = FALSE,
        col = col_fun_human, show_column_names = FALSE, cluster_rows = FALSE,
        cluster_columns = FALSE, show_heatmap_legend = TRUE,
        left_annotation = ha_mouse)

human_specific_Heatmap_pig <- Heatmap(t(human_spec_mat_pig), show_row_names = FALSE,
        col = col_fun_human, show_column_names = FALSE, cluster_rows = FALSE,
        cluster_columns = FALSE, show_heatmap_legend = TRUE,
        left_annotation = ha_pig)

merged_human_specific_Heatmap <- human_specific_Heatmap_human %v% human_specific_Heatmap_mouse %v% human_specific_Heatmap_pig

pdf("figE_all_results/human_specific_gene.pdf",width=8,height=15)
draw(merged_human_specific_Heatmap)
dev.off()

# mouse specific
mouse_specific_Heatmap_human <- Heatmap(t(mouse_spec_mat_human), show_row_names = FALSE,
        col = col_fun_mouse, show_column_names = FALSE, cluster_rows = FALSE,
        cluster_columns = FALSE, show_heatmap_legend = TRUE,
        left_annotation = ha_human)

mouse_specific_Heatmap_mouse <- Heatmap(t(mouse_spec_mat_mouse), show_row_names = FALSE,
        col = col_fun_mouse, show_column_names = FALSE, cluster_rows = FALSE,
        cluster_columns = FALSE, show_heatmap_legend = TRUE,
        left_annotation = ha_mouse)

mouse_specific_Heatmap_pig <- Heatmap(t(mouse_spec_mat_pig), show_row_names = FALSE,
        col = col_fun_mouse, show_column_names = FALSE, cluster_rows = FALSE,
        cluster_columns = FALSE, show_heatmap_legend = TRUE,
        left_annotation = ha_pig)

merged_mouse_specific_Heatmap <- mouse_specific_Heatmap_human %v% mouse_specific_Heatmap_mouse %v% mouse_specific_Heatmap_pig

pdf("figE_all_results/mouse_specific_gene.pdf",width=8,height=15)
draw(merged_mouse_specific_Heatmap)
dev.off()

# pig specific
pig_specific_Heatmap_human <- Heatmap(t(pig_spec_mat_human), show_row_names = FALSE,
        col = col_fun_pig, show_column_names = FALSE, cluster_rows = FALSE,
        cluster_columns = FALSE, show_heatmap_legend = TRUE,
        left_annotation = ha_human)

pig_specific_Heatmap_mouse <- Heatmap(t(pig_spec_mat_mouse), show_row_names = FALSE,
        col = col_fun_pig, show_column_names = FALSE, cluster_rows = FALSE,
        cluster_columns = FALSE, show_heatmap_legend = TRUE,
        left_annotation = ha_mouse)

pig_specific_Heatmap_pig <- Heatmap(t(pig_spec_mat_pig), show_row_names = FALSE,
        col = col_fun_pig, show_column_names = FALSE, cluster_rows = FALSE,
        cluster_columns = FALSE, show_heatmap_legend = TRUE,
        left_annotation = ha_pig)

merged_pig_specific_Heatmap <- pig_specific_Heatmap_human %v% pig_specific_Heatmap_mouse %v% pig_specific_Heatmap_pig

pdf("figE_all_results/pig_specific_gene.pdf",width=8,height=15)
draw(merged_pig_specific_Heatmap)
dev.off()


################ save results
save_all_results <- function(dir = "all_results") {
  if (!dir.exists(dir)) dir.create(dir)
  
  # Helper
  save_obj <- function(x, name) {
    saveRDS(x, file = file.path(dir, paste0(name, ".rds")))
  }
  
  message("Saving core DEG results...")
  save_obj(human_deg,  "human_deg")
  save_obj(mouse_deg,  "mouse_deg")
  save_obj(pig_deg,    "pig_deg")
  
  message("Saving species-specific DEGs...")
  save_obj(human_spec_deg, "human_spec_deg")
  save_obj(mouse_spec_deg, "mouse_spec_deg")
  save_obj(pig_spec_deg,   "pig_spec_deg")
  
  message("Saving conserved DEGs...")
  save_obj(conserved_deg, "conserved_deg")
  save_obj(genes_conserved, "genes_conserved")
  
  message("Saving downsampled Seurat object...")
  save_obj(downsample_1000_RDS, "downsample_1000_RDS")
  
  message("Saving conserved matrices...")
  save_obj(human_conserved_mat, "human_conserved_mat")
  save_obj(mouse_conserved_mat, "mouse_conserved_mat")
  save_obj(pig_conserved_mat,   "pig_conserved_mat")
  
  message("Saving Human-specific matrices...")
  save_obj(human_spec_mat_human, "human_spec_mat_human")
  save_obj(human_spec_mat_mouse, "human_spec_mat_mouse")
  save_obj(human_spec_mat_pig,   "human_spec_mat_pig")
  
  message("Saving Mouse-specific matrices...")
  save_obj(mouse_spec_mat_mouse, "mouse_spec_mat_mouse")
  save_obj(mouse_spec_mat_human, "mouse_spec_mat_human")
  save_obj(mouse_spec_mat_pig,   "mouse_spec_mat_pig")
  
  message("Saving Pig-specific matrices...")
  save_obj(pig_spec_mat_pig,    "pig_spec_mat_pig")
  save_obj(pig_spec_mat_human,  "pig_spec_mat_human")
  save_obj(pig_spec_mat_mouse,  "pig_spec_mat_mouse")

  message("Saving counts...")
  save_obj(n_gene, "n_gene_conserved_specific")
  
  message("Saving cell orders...")
  save_obj(human_cells_order, "human_cells_order")
  save_obj(mouse_cells_order, "mouse_cells_order")
  save_obj(pig_cells_order,   "pig_cells_order")
  
  message("Saving species subset Seurat objects (optional)...")
  save_obj(human_obj, "human_obj")
  save_obj(mouse_obj, "mouse_obj")
  save_obj(pig_obj,   "pig_obj")

  message("== All important results saved successfully ==")
}

save_all_results("figE_all_results")
