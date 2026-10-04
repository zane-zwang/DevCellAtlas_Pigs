# Purpose: adjacent stage findmarkers.
analysis_dir <- Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input")

library(Seurat)
library(magrittr)
library(future)
plan("multicore", workers = 8)
options(future.globals.maxSize = 100000 * 1024^5)

######data import
all.obj <- readRDS(file.path(analysis_dir, "atlas/all_seurat4.rds"))

######scale data
all.obj <- all.obj %>% NormalizeData()

######diff analysis
celltypes <- unique(all.obj$celltype)
stage <- c('e55d','e90d','0d','30d','90d','180d')

all.obj$celltype.stage <- paste0(all.obj$celltype,"_",all.obj$stage)
Idents(all.obj) <- "celltype.stage"

markers_list=list()

for (j in 1:length(celltypes)) {
    for (i in 1:(length(stage) - 1)) {
        ident1 <- paste0(celltypes[j], "_", stage[i])
        ident2 <- paste0(celltypes[j], "_", stage[i + 1])
        
        num_cells_ident1 <- sum(Idents(all.obj) == ident1)
        num_cells_ident2 <- sum(Idents(all.obj) == ident2)
        
        if (num_cells_ident1 >= 10 & num_cells_ident2 >= 10) {
            subset_name <- paste0(ident1, "_", stage[i + 1])
            markers <- FindMarkers(all.obj, ident.1 = ident1, ident.2 = ident2, verbose = FALSE)
            markers_list[[subset_name]] <- markers

            safe_name <- gsub("[^A-Za-z0-9_\\-]", "_", subset_name)
            safe_name <- gsub(" ", "_", safe_name)

            write.csv(markers_list[[subset_name]], paste0(safe_name, ".csv"))
        } else {
            message(paste("Skipping", ident1, "or", ident2, "cell number < 10"))
        }
    }
}

saveRDS(markers_list, "findmarkers_cell_sta_list.rds")
