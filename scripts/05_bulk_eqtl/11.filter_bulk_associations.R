# Purpose: filter bulk associations.

library(data.table)
analysis_dir <- Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/analysis_data")
input_dir <- Sys.getenv("QTL_BULK_ASSOCIATION_DIR", file.path(analysis_dir, "qtl/mapping/bulk_omiga/01_acat"))

groups <- c("bin1", "bin2", "bin3", "bin4")
output_dir <- Sys.getenv("QTL_BULK_FILTER_OUTPUT_DIR", file.path(input_dir, "beQTL/02_filter"))

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

all_significant_list <- list()

for (grp in groups) {
    file_name <- file.path(input_dir, grp, paste0("muscle_", grp, ".cis_qtl.txt.gz"))

    if (!file.exists(file_name)) {
      warning("File not found: ", file_name)
      next
    }

    message("Processing: ", file_name)
    df <- fread(file_name)

    df$is_eGene <- (df$pval_g1_acat < 0.05) & (df$qval_g1 < 0.05)
    df$group <- grp

    all_significant_list[[paste0(grp)]] <- df[is_eGene == TRUE]

    out_file <- file.path(output_dir, paste0(grp, "_significant_beQTLs.tsv"))
    fwrite(df[is_eGene == TRUE], out_file, sep = "\t")
}

all_res <- rbindlist(all_significant_list, fill = TRUE)
fwrite(all_res, file.path(output_dir, "all_significant_beQTLs.tsv"), sep = "\t")
