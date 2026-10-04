# Purpose: export bulk associations.

library(data.table)

groups <- c("bin1", "bin2", "bin3", "bin4")
output_dir <- "./"

if (!dir.exists(output_dir)) dir.create(output_dir)

all_significant_list <- list()

for (grp in groups) {
    file_name <- paste0("muscle_", grp, ".cis_qtl.txt.gz")

    if (!file.exists(file_name)) {
      warning("File not found: ", file_name)
      next
    }

    message("Processing: ", file_name)
    df <- fread(file_name)

    df$is_eGene <- (df$pval_g1_acat < 0.05) & (df$qval_g1 < 0.05)
    df$group <- grp

    all_significant_list[[paste0(grp)]] <- df
}

all_res <- rbindlist(all_significant_list, fill = TRUE)
fwrite(all_res, file.path(output_dir, "all_beQTLs.tsv"), sep = "\t")
