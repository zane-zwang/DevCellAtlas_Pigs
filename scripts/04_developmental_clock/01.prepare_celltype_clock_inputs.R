# Purpose: prepare celltype clock inputs.
analysis_dir <- Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input")

library(Seurat)
library(tidyverse)
library(glmnet)
library(ggplot2)
library(factoextra)
library(cluster)
library(DESeq2)

tissues <- c('adipose', 'cerebrum','duodenum', 'heart', 'hypothalamus', 'liver', 'muscle')

all <- readRDS(file.path(analysis_dir, "atlas/all_seurat4.rds"))

for(i in 1:length(tissues)){

obj <- subset(all, subset=tissue==tissues[[i]])
obj@meta.data$celltype <- as.character(obj@meta.data$celltype) %>% tolower() %>%
                str_replace_all(pattern = " ", replacement = "_") %>%
                str_replace_all(pattern = "-", replacement = "_") %>%
                str_replace_all(pattern = ",", replacement = "") %>%
                str_replace_all(pattern = "/", replacement = "_")

obj@meta.data$age <- obj@meta.data$stage

obj@meta.data <- obj@meta.data %>%  
  mutate(age = case_when(  
    stage == "e55d" ~ -59,  
    stage == "e90d" ~ -24,  
    stage == "0d" ~ 0,  
    stage == "30d" ~ 30,  
    stage == "90d" ~ 90,  
    stage == "180d" ~ 180,  
  ))


obj2 <- obj %>% NormalizeData(verbose = F)
saveRDS(obj2, paste0(tissues[[i]],'_logNorm.rds'))

obj1 <- SCTransform(obj)
saveRDS(obj1, paste0(tissues[[i]],'_SCT.rds'))
}
