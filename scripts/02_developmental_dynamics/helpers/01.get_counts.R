# Export nucleus-level raw counts for pySCENIC.
# Input: Seurat object with an RNA counts layer.
# Output: cell-by-gene RNAcounts.csv.
library(Seurat)

input_file <- Sys.getenv("SCENIC_SEURAT_RDS", "path/to/input/seurat.data.rds")
output_file <- Sys.getenv("SCENIC_COUNTS_CSV", "RNAcounts.csv")

obj <- readRDS(input_file)
counts_matrix <- GetAssayData(obj, layer = "counts")
write.csv(t(as.matrix(counts_matrix)), file = output_file)
