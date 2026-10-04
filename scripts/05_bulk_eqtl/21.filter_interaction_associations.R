# Purpose: filter interaction associations.

library(data.table)
analysis_dir <- Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/analysis_data")
input_dir <- Sys.getenv("QTL_INTERACTION_ASSOCIATION_DIR", file.path(analysis_dir, "qtl/mapping/bisque_omiga/01_acat"))

groups <- c("bin1", "bin2", "bin3", "bin4")
cells  <- c("adip", "cap", "fap", "lec", "mac", "msc", "peri", "glial", "tc", "teno", "mf1", "mf2ab", "mf2x")
output_dir <- Sys.getenv("QTL_INTERACTION_FILTER_OUTPUT_DIR", file.path(analysis_dir, "qtl/interaction/02_filter_2"))

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
    
    df$is_ieGene <- df$pval_g2 < df$pval_g2_acat & df$qval_g2 < 0.05
    df$group <- grp
    df$cell  <- cl
    
    all_significant_list[[paste0(grp,'_',cl)]] <- df[is_ieGene == TRUE]

    out_file <- file.path(output_dir, paste0(grp, "_", cl, "_significant_ieQTLs.tsv"))
    fwrite(df[is_ieGene == TRUE], out_file, sep = "\t")
  }
}

all_res <- rbindlist(all_significant_list, fill = TRUE)
fwrite(all_res, file.path(output_dir, "all_significant_ieQTLs.tsv"), sep = "\t")
