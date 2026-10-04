# Purpose: tissue celltype cpm.
analysis_dir <- Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input")

library(Seurat)
library(edgeR)

obj <- readRDS('all_ageGenes_logNrom.rds')

# method 3
dir.create(file.path(analysis_dir, "tai/pseudobulk"))
dir.create(file.path(analysis_dir, "tai/pseudobulk/01_tissues"))
dir.create(file.path(analysis_dir, "tai/pseudobulk/02_tissues_cells"))
setwd(file.path(analysis_dir, "tai/pseudobulk"))
tissues <- c('adipose','cerebrum','duodenum','heart','hypothalamus','liver','muscle')

# for tissues
for(i in 1:length(tissues)){
    ts_sub <- subset(obj, subset=tissue==tissues[i])
    current_tissue <- tissues[i]
    message("Processing: ", current_tissue)

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

    write.csv(exprSet, paste0('01_tissues/',tissues[[i]],'_stage_pseudo_Exp.csv'))
}

# for cell types
for(i in 1:length(tissues)){
    ts_sub <- subset(obj, subset=tissue==tissues[i])
    celltypes <- unique(ts_sub@meta.data$celltype)
    
    for(j in 1:length(celltypes)){
        cell_sub <- subset(ts_sub, subset = celltype == celltypes[j])

        current_tissue_celltype <- paste0(tissues[i], "_", celltypes[j])
        message("Processing: ", current_tissue_celltype)

        bs <- split(colnames(cell_sub), cell_sub$stage)
        bs <- bs[sapply(bs, length) > 0]

        ct_list <- lapply(names(bs), function(x) {
            kp = colnames(cell_sub) %in% bs[[x]]
            pseudosum = rowSums(as.matrix(cell_sub@assays$RNA@counts[, kp]))
            return(pseudosum)
        })
        ct <- do.call(cbind, ct_list)
        colnames(ct) <- names(bs)

        # filter
        exprSet <- ct
        exprSet <- exprSet[apply(exprSet, 1, function(x) sum(x > 1) > 1), ]
        exprSet <- cpm(exprSet)

        safe_celltype <- gsub("[^a-zA-Z0-9_\\-]", "_", trimws(celltypes[j]))
        safe_tissue <- gsub("[^a-zA-Z0-9_\\-]", "_", trimws(tissues[i]))
        
        output_file <- paste0('02_tissues_cells/', safe_tissue, '_', safe_celltype, '_stage_pseudo_Exp.csv')
        write.csv(exprSet, output_file)
    }
}
