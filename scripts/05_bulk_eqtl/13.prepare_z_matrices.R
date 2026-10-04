# Purpose: prepare z matrices.
analysis_dir <- Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input")

# combine strong set z-score
# R
library(data.table)
mashr_dir <- file.path(analysis_dir, "qtl/mashr")
in_dir <- file.path(mashr_dir, "01_beqtl_mashr_input")
bins <- c("bin1", "bin2", "bin3", "bin4")

z_list <- list()
for (b in bins) {
    f <- file.path(in_dir, paste0(b,"_nominal_pairs.extracted_pairs.txt.gz"))
    dt <- fread(f)
    dt[, zval := beta_g1 / beta_se_g1]

    dt[, pair_id := gsub(",", "_", pair_id)]
    dt <- dt[, .(pair_id, zval)]
    
    setnames(dt, "zval", b)
    z_list[[b]] <- dt
}

zval_strong <- Reduce(function(x,y) merge(x, y, by="pair_id", all=TRUE), z_list)
fwrite(zval_strong, file.path(in_dir, "strong_set_zval.txt"), sep="\t")

## prepare random set
# R
eqtl_dir <- file.path(analysis_dir, "qtl/mapping/bulk_omiga/01_acat")
strong_pairs <- fread(file.path(in_dir,"strong_set_zval.txt"))$pair_id

random_list <- list()
for (b in bins) {
    nom_files <- list.files(path = file.path(eqtl_dir, b),
                            pattern = "cis_qtl_pairs.*.txt.gz",
                            full.names = TRUE)
    dt_all <- rbindlist(lapply(nom_files, fread), use.names=TRUE, fill=TRUE)
    dt_all[, pair_id := paste(variant_id, pheno_id, sep = "_")]
    dt_all <- dt_all[!pair_id %in% strong_pairs]
    dt_all[, zval := beta_g1 / beta_se_g1]
    dt_all <- dt_all[, .(pair_id, zval)]
    setnames(dt_all, "zval", b)
    random_list[[b]] <- dt_all
}

random_set <- Reduce(function(x, y) merge(x, y, by="pair_id", all=TRUE), random_list)

set.seed(123)
random_sample <- random_set[sample(.N, min(50000, .N))]
fwrite(random_sample, file.path(in_dir, "random_set_zval.txt"), sep="\t")
