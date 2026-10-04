# Purpose: filter cse associations.

library(data.table)
analysis_dir <- Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/analysis_data")
input_dir <- Sys.getenv("QTL_CSE_ASSOCIATION_DIR", file.path(analysis_dir, "qtl/mapping/cse_omiga/01_acat"))

groups <- c("bin1", "bin2", "bin3", "bin4")
cells  <- c("adipocyte", "fibro_adipogenic_progenitor_cell", "lymphatic_endothelial_cell", "muscle_stem_cell", "peripheral_glial", "tenocyte", "type_ii_x_myonuclei", "capillary_endothelial_cell", "macrophage", "pericyte", "t_cell", "type_ii_a_b_myonuclei", "type_i_myonuclei")
output_dir <- Sys.getenv("QTL_CSE_FILTER_OUTPUT_DIR", file.path(analysis_dir, "qtl/cse/02_filter"))

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

all_significant_list <- list()
for (grp in groups) {
  for (cl in cells) {
    
    file_name <- file.path(input_dir, grp, cl, paste0("muscle_", grp, "_", cl, ".cis_qtl.txt.gz"))
    
    if (!file.exists(file_name)) {
      warning("File not found: ", file_name)
      next
    }
    
    message("Processing: ", file_name)
    df <- fread(file_name)
    
    df$qval_g1 <- p.adjust(df$pval_g1, method = "BH")

    df$is_ieGene <- (df$pval_g1_acat < 0.05) & (df$qval_g1 < 0.05)
    df$group <- grp
    df$cell  <- cl
    
    all_significant_list[[paste0(grp,'_',cl)]] <- df[is_ieGene == TRUE]

    out_file <- file.path(output_dir, paste0(grp, "_", cl, "_significant_cseQTLs.tsv"))
    fwrite(df[is_ieGene == TRUE], out_file, sep = "\t")
  }
}

all_res <- rbindlist(all_significant_list, fill = TRUE)
fwrite(all_res, file.path(output_dir, "all_significant_cseQTLs.tsv"), sep = "\t")
