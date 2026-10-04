# Purpose: bmind profile.
library(MIND)
library(dplyr)
library(Seurat)
library(tibble)

reference_rds <- Sys.getenv("BMIND_REFERENCE_RDS", "path/to/input/muscle_logNorm.rds")
fractions_csv <- Sys.getenv("BMIND_BISQUE_FRACTIONS_CSV", "path/to/input/bulk_props_bisque.csv")
gene_ids_csv <- Sys.getenv("BMIND_GENE_IDS_CSV", "path/to/input/genes_ensembl.csv")
bulk_counts_csv <- Sys.getenv("BMIND_BULK_COUNTS_CSV", "path/to/input/muscle_GTEx_counts_752.csv")
sc <- readRDS(reference_rds)
prop_raw <- read.csv(fractions_csv, header=T, row.names=1)
prop <- t(prop_raw)

Idents(sc) <- "celltype"
counts <- data.frame(sc@assays$RNA@counts)
counts <- t(counts)
counts <- data.frame(counts)
meta <- data.frame(sc@meta.data)

meta <- droplevels(meta)
counts <- counts %>% mutate(cluster = meta$celltype)
group <- list()
group <- split(counts,counts$cluster)
name<- rownames(table(counts$cluster))

cluMean <- colMeans(group[[1]][, -ncol(group[[1]])])
for (j in 2:length(group)) {
  flmean <- colMeans(group[[j]][, -ncol(group[[j]])])
  cluMean <- rbind(cluMean, flmean)
}

rownames(cluMean) <- name
cluMean <- t(cluMean)
profile <- as.data.frame(cluMean)

profile <- profile[rowSums(profile) > 0, ]

gene_id <- read.csv(gene_ids_csv, row.names=1)
profile <- profile %>%
  rownames_to_column(var = "SYMBOL") %>%
  left_join(gene_id, by = "SYMBOL") %>%
  mutate(RowName = if_else(!is.na(ENSEMBL), ENSEMBL, SYMBOL)) %>%
  select(-SYMBOL, -ENSEMBL) %>%
  column_to_rownames("RowName")


df <- read.csv(bulk_counts_csv, header=T, row.names=1)
bulk <- df

frac <- prop
profile <- profile[,colnames(profile) %in% colnames(frac)]
profile <- data.frame(profile[,order(colnames(profile))])
profile <- profile[! rowSums(profile) == 0,]
profile <- profile[rownames(profile) %in% rownames(bulk),]
bulk <- bulk[rownames(bulk) %in% rownames(profile),]
colnames(frac) <- colnames(profile)
profile <- data.frame(profile[order(rownames(profile)),])
bulk <- data.frame(bulk[order(rownames(bulk)),])
  
colnames(bulk) = rownames(frac) = paste0('s', 1:nrow(frac))
colnames(frac) = colnames(profile) = paste0('c', 1:ncol(frac))

deconv = bMIND(bulk, frac = frac, profile = profile, ncore = 12)

saveRDS(deconv, 'bMIND_deconv_exp_profile.rds')
