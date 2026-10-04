# Purpose: pseudobulk tai.
publication_root <- normalizePath(Sys.getenv("DEVPIGATLAS_ROOT", "."), mustWork = TRUE)

library(ggplot2)
library(myTAI)
library(cowplot)
library(DESeq2)
library(patchwork)
source(file.path(publication_root, "figures/Figure6/02.plot_signature_helper.R"))
source(file.path(publication_root, "figures/Figure6/03.plot_contribution_helper.R"))
source(file.path(publication_root, "figures/Figure6/04.plot_relative_expression_helper.R"))
custom.myTAI.cols <- function(n) {
  cols <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
            "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
            "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
            "#8A7C64", "#599861")
  return(cols[1:n])
}

### tissue tai
sc.data <- readr::read_tsv('gene18668_ages.tsv')
sc.PhyloMap <- dplyr::select(sc.data, Phylostratum = rank, GeneID = 'SYMBOL')
tissues <- c('adipose','cerebrum','duodenum','heart','hypothalamus','liver','muscle')

read_exp_csv <- function(name){
  df <- read.csv(name, header = T, row.names = 1)
  df$GeneID <- rownames(df)
  return(df)
}

folder_path <- 'pseudobulk_exp_for_tai/01_tissues/'
file_names <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)
ts_pseudo_exp_list <- lapply(file_names, read_exp_csv)
names(ts_pseudo_exp_list) <- tissues

folder_path <- 'pseudobulk_exp_for_tai/02_tissues_cells/'
ts_cell_pseudo_exp_list <- list()
for(tissue in tissues){
  pattern <- paste0("^", tissue, ".*\\.csv$")
  file_names <- list.files(path = folder_path, pattern = pattern, full.names = TRUE)
  ts_cell_pseudo_exp_list[[tissue]] <- lapply(file_names, read_exp_csv)
  cell_names <- gsub("_stage_pseudo_aveExp\\.csv$", "", basename(file_names))
  names(ts_cell_pseudo_exp_list[[tissue]]) <- cell_names
}

folder_path <- 'animal_age_data/DevCortex_pseudo_exp/04_celltypes/'
file_names <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)
ct_pseudo_exp_list <- lapply(file_names, read_exp_csv)
names(ct_pseudo_exp_list) <- gsub("_stage_pseudo_aveExp\\.csv$", "", basename(file_names))

tai_process_log <- function(df){
  target_order <- c('Phylostratum', 'GeneID', 'e55d', 'e90d', 'X0d', 'X30d', 'X90d', 'X180d')
  actual_order <- intersect(target_order, colnames(df))
  df <- df[, actual_order]
  df <- na.omit(df)
  df_log2 <- tf(df, log2, pseudocount = 1)
  return(df_log2)
}

pseudo_ts_plot_list <- list()
pseudo_ts_tai_list <- list()
pseudo_ts_flat_list <- list()
for(t in tissues){
  message("Processing: ", t)
  sc.PES <- merge(sc.PhyloMap, ts_pseudo_exp_list[[t]], by = 'GeneID')
  sc.PES <- tai_process_log(sc.PES)
  pseudo_ts_plot_list[[t]] <- PlotSignature_z(sc.PES)+
    xlab('')+labs(title = t)
  pseudo_ts_tai_list[[t]] <- TAI(sc.PES)
  pseudo_ts_flat_list[[t]] <- FlatLineTest(ExpressionSet = sc.PES)
}

pseudo_ts_plot_list[[7]]
saveRDS(pseudo_ts_plot_list, 'tai_ensembl/pseudo_exp_raw/tissue/pseudo_ts_plot_list_log.rds')
saveRDS(pseudo_ts_tai_list, 'tai_ensembl/pseudo_exp_raw/tissue/pseudo_ts_tai_list_log.rds')
saveRDS(pseudo_ts_flat_list, 'tai_ensembl/pseudo_exp_raw/tissue/pseudo_ts_flat_list_log.rds')

plots <- readRDS('tai_ensembl/pseudo_exp_cus/tissue/pseudo_ts_plot_list_log.rds')


wrap_plots(pseudo_ts_plot_list, ncol = 3)

sc.PES <- merge(sc.PhyloMap, ts_pseudo_exp_list[[6]], by = 'GeneID')
sc.PES <- tai_process_log(sc.PES)
PlotSignature(sc.PES, measure = "TAI")
PlotContribution_z(sc.PES, legendName = "PS")
PlotRE_z(sc.PES, Groups = list(1:3))

############################## tissue cell tai

pseudoExp_cell_plot_list <- list()
pseudoExp_cell_tai_list <- list()
pseudoExp_cell_flat_list <- list()
for(t in tissues){
  message("Processing tissue: ", t)
  cells <- names(ts_cell_pseudo_exp_list[[t]])
  pseudoExp_cell_plot_list[[t]] <- list()
  pseudoExp_cell_tai_list[[t]] <- list()
  pseudoExp_cell_flat_list[[t]] <- list()
  
  for(c in cells){
    message("  Processing cell type: ", c)
    
    tryCatch({
      sc.PES <- merge(sc.PhyloMap, ts_cell_pseudo_exp_list[[t]][[c]], by = 'GeneID')
      sc.PES <- tai_process_log(sc.PES)
      
      pseudoExp_cell_plot_list[[t]][[c]] <- PlotSignature(sc.PES) +
        xlab('') + labs(title = c)
      pseudoExp_cell_tai_list[[t]][[c]] <- TAI(sc.PES)
      pseudoExp_cell_flat_list[[t]][[c]] <- FlatLineTest(ExpressionSet = sc.PES)
    }, error = function(e){
      warning(sprintf("    Skipped cell type '%s' in tissue '%s' due to error: %s", c, t, e$message))
    })
  }
}
saveRDS(pseudoExp_cell_plot_list, 'tai_ensembl/pseudo_exp_raw/cell/pseudoExp_cell_plot_list_log.rds')
saveRDS(pseudoExp_cell_tai_list, 'tai_ensembl/pseudo_exp_raw/cell/pseudoExp_cell_tai_list_log.rds')
saveRDS(pseudoExp_cell_flat_list, 'tai_ensembl/pseudo_exp_raw/cell/pseudoExp_cell_flat_list_log.rds')

pseudo_ts_plot_list$muscle


# celltype
pseudo_ct_plot_list <- list()
pseudo_ct_tai_list <- list()
pseudo_ct_flat_list <- list()

cts <- names(ct_pseudo_exp_list)
for(t in cts){
  message("Processing: ", t)
  
  tryCatch({
    sc.PES <- merge(sc.PhyloMap, ct_pseudo_exp_list[[t]], by = 'GeneID')
    sc.PES <- tai_process_log(sc.PES)
    
    pseudo_ct_plot_list[[t]] <- PlotSignature(sc.PES) +
      xlab('') + labs(title = t)
    pseudo_ct_tai_list[[t]] <- TAI(sc.PES)
    pseudo_ct_flat_list[[t]] <- FlatLineTest(ExpressionSet = sc.PES)
  }, error = function(e){
    warning(sprintf("  Skipped cell type '%s' due to error: %s", t, e$message))
  })
}

pseudo_ct_plot_list[[1]]
saveRDS(pseudo_ct_plot_list, 'tai_ensembl/pseudo_exp_raw/ct/pseudo_ct_plot_list_log.rds')
saveRDS(pseudo_ct_tai_list, 'tai_ensembl/pseudo_exp_raw/ct/pseudo_ct_tai_list_log.rds')
saveRDS(pseudo_ct_flat_list, 'tai_ensembl/pseudo_exp_raw/ct/pseudo_ct_flat_list_log.rds')
