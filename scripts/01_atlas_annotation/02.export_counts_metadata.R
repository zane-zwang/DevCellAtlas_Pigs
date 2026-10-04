# Purpose: export counts metadata.
analysis_dir <- Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input")

library(Seurat)
library(Matrix)

get_meta <- function(obj){
    meta <- obj@meta.data
    return(meta)
}

get_raw_counts <- function(obj){
    counts <- GetAssayData(object = obj, assay = "RNA", layer = "counts")
    return(counts)
}

tissues <- c('adipose', 'cerebrum', 'duodenum', 'heart', 'hypothalamus', 'liver', 'muscle')

for(i in tissues){
    obj <- readRDS(paste0(file.path(analysis_dir, "clock/input/"),i, '_logNorm.rds'))
    counts <- get_raw_counts(obj)
    meta <- get_meta(obj)

    writeMM(counts, paste0(i, '_counts.mtx'))
    write.table(rownames(counts), paste0(i, '_genes.txt'),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(colnames(counts), paste0(i, '_barcodes.txt'),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.csv(meta, paste0(i, '_meta.csv'))
}
