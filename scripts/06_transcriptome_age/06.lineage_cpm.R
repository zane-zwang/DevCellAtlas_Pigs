# Purpose: lineage cpm.
analysis_dir <- Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input")

# average exp from seurat obj

library(Seurat)
library(edgeR)

obj <- readRDS('all_ageGenes_logNrom.rds')

# method 3
dir.create(file.path(analysis_dir, "tai/pseudobulk/03_celllineages"))
setwd(file.path(analysis_dir, "tai/pseudobulk"))
celllineages <- unique(obj@meta.data$celllineage)

# for celllineages
for(i in 1:length(celllineages)){
    ts_sub <- subset(obj, subset=celllineage==celllineages[i])
    current_celllineage <- celllineages[i]
    message("Processing: ", current_celllineage)

    bs <- split(colnames(ts_sub), ts_sub$stage)
    bs <- bs[sapply(bs, length) > 0]

    ct_list <- lapply(names(bs), function(x) {
            kp = colnames(ts_sub) %in% bs[[x]]
            pseudomean = rowSums(as.matrix(ts_sub@assays$RNA@counts[, kp]))
            return(pseudomean)
    })
    ct <- do.call(cbind, ct_list)
    colnames(ct) <- names(bs)

    exprSet <- ct
    exprSet <- exprSet[apply(exprSet, 1, function(x) sum(x > 1) > 1), ]
    exprSet <- cpm(exprSet)

    write.csv(exprSet, paste0('03_celllineages/',celllineages[[i]],'_stage_pseudo_Exp.csv'))
}
