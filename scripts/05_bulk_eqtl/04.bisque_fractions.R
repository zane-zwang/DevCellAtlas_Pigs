# Purpose: bisque fractions.
analysis_dir <- Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input")

### conda activate omnideconv

library("SingleCellExperiment")
library("BisqueRNA")
library("here")
library("sessioninfo")
library("Seurat")

############## BisqueRNA
## bulk input (raw counts)

all_counts <- read.table(file.path(analysis_dir, "qtl/bulk/pigGTEx/PigGTEx_v0.Gene.raw_count.txt"), header=T)
sel_samp <- read.csv(file.path(analysis_dir, "qtl/bulk/pigGTEx/00_muscle_752/Selected_samples_meta_muscle_752.CSV"),header=T)

df <- all_counts[,colnames(all_counts) %in% sel_samp$BioSample]

rownames(df) <- all_counts$X

sc_obj <- readRDS(file.path(analysis_dir, "clock/00_data/muscle_logNorm.rds"))

symbol_ensembl <- read.csv('genes_ensembl.csv', header=T, row.names=1)
id_trans <- function(mx){
    map_sym2ens <- setNames(symbol_ensembl$ENSEMBL, symbol_ensembl$SYMBOL)
    orig_genes <- rownames(mx)
    new_genes <- map_sym2ens[orig_genes]
    new_genes[is.na(new_genes)] <- orig_genes[is.na(new_genes)]
    rownames(mx) <- new_genes
    return(mx)
}

sc_counts <- GetAssayData(object = sc_obj, assay = "RNA", layer = "counts")
sc_counts <- id_trans(sc_counts)
sc_meta <- sc_obj@meta.data

common_genes <- intersect(rownames(sc_counts), rownames(df))
markers <- common_genes

bulk_expr_subset <- df[markers, ]
bulk_meta <- data.frame(SAMPLE_ID = colnames(bulk_expr_subset))
rownames(bulk_meta) <- colnames(bulk_expr_subset)
exp_set_bulk <- ExpressionSet(
  assayData = as.matrix(bulk_expr_subset),
  phenoData = AnnotatedDataFrame(bulk_meta)
)

exp_set_sce <- ExpressionSet(assayData = as.matrix(sc_counts[markers,]), phenoData = AnnotatedDataFrame(sc_meta))


exp_set_sce_temp <- exp_set_sce[markers,]
zero_cell_filter <- colSums(exprs(exp_set_sce_temp)) != 0
message("Exclude ",sum(!zero_cell_filter), " cells")
exp_set_sce_temp <- exp_set_sce_temp[,zero_cell_filter]

est_prop_bisque <- ReferenceBasedDecomposition(bulk.eset = exp_set_bulk[markers,],
                                               sc.eset = exp_set_sce_temp,
                                               cell.types = "celltype",
                                               subject.names = "sample",
                                               use.overlap = FALSE)
saveRDS(est_prop_bisque, 'bisque_full.rds')

bulk_props <- est_prop_bisque$bulk.props
write.csv(bulk_props, 'bulk_props.csv', quote=F)

bulk_prop <- read.csv('bulk_props_bisque.csv', header=T, row.names=1)
df <- t(bulk_prop)
adip_freq <- data.frame(IID = rownames(df), AITERM = df[,1])
write.table(adip_freq, 'adip_freq.txt', sep='\t', row.names=F, quote=F)

celltypes <- c('adip','cap','fap','lec','mac','msc','peri','glial','tc','teno','mf1','mf2ab','mf2x')
for(i in 1:length(celltypes)){
    tmp <- data.frame(IID = rownames(df), AITERM = df[,i])
    write.table(tmp, paste0(celltypes[[i]],'_freq.txt'), sep='\t', row.names=F, quote=F)
}
