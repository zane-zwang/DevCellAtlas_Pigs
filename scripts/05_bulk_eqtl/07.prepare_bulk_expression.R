# Purpose: prepare bulk expression.
publication_root <- normalizePath(Sys.getenv("DEVPIGATLAS_ROOT", "."), mustWork = TRUE)
analysis_dir <- Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input")

library(dplyr)
source(file.path(publication_root, "scripts/05_bulk_eqtl/06.write_gct.R"))


all_tpm <- read.csv(file.path(analysis_dir, "qtl/input/pigGTEx/muscle_GTEx_tpm_752.csv"), row.names=1)
all_counts <- read.csv(file.path(analysis_dir, "qtl/input/pigGTEx/muscle_GTEx_counts_752.csv"), row.names=1)

groups <- c('binAll')

for(i in groups){
   sample <- read.table(paste0(i,'/sample_list.txt'))
   counts <- all_counts[, colnames(all_counts)%in%sample$V1]
   tpm <- all_tpm[, colnames(all_tpm)%in%sample$V1]
   write.csv(counts, paste0(i,'/muscle_counts.csv'), quote=F)
   write.csv(tpm, paste0(i,'/muscle_tpm.csv'), quote=F)
   expr_mat <- as.matrix(tpm)
   write_matrix_as_gct(expr_mat, paste0(i,'/muscle.tpm.gct'))
   expr_mat <- as.matrix(counts)
   write_matrix_as_gct(expr_mat, paste0(i,'/muscle.counts.gct'))
  }
