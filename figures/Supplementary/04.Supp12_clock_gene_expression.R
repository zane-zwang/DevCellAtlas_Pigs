# Purpose: Supplementary clock gene expression.

library(ggplot2)
library(tidyverse)
library(tidyr)

ct_cols = c("Capillary endothelial cell"="#47ad85","Macrophage"="#fc9a00","Fibro-adipogenic progenitor cell"="#d47f4f",
            "Muscle stem cell"="#7c5a33","Pericyte"="#4c2f0d","Type I myonuclei"="#d96666","Type II a/b myonuclei"="#ec7979",
            "Type II x myonuclei"="#ffb3b3")
            
tpm <- read.csv('pseudobulk_data/bulk_tpm.csv', header = TRUE, row.names = 1)
tpm <- t(tpm)
tpm <- as.data.frame(tpm)
ages <- read.table('pseudobulk_data/ages')
ages <- as.numeric(ages$V1)
tpm$age <- ages

long_df <- tpm %>% pivot_longer(cols = -age, names_to = "gene", values_to = "expression")
long_df$age_numeric <- as.numeric(as.character(long_df$age))
long_df$log_expr <- log1p(long_df$expression)

ggplot(long_df %>% filter(gene %in% c("ENSSSCG00000030095")), 
       aes(x = age_numeric, y = log_expr)) +
  geom_point(alpha = 0.6, shape = 15) +
  geom_smooth(method = "loess", se = TRUE, color = 'firebrick', linewidth = 1) +
  labs(x = "Age", y = "log1p(TPM)") +
  theme_bw()+
  theme(
    aspect.ratio = 1,
    axis.text = element_text(colour = 'black', size = 18),
    axis.title = element_text(colour = 'black', size =18),
    panel.grid.minor = element_blank(),
  )

library(Seurat)
obj1 <- readRDS('muscle_annot.rds')
obj2 <- readRDS('muscle_seurat4_hvg.rds')

counts <- GetAssayData(obj1, slot = 'counts')
meta <- obj2@meta.data

obj <- CreateSeuratObject(counts = counts, meta.data = meta)
obj <- NormalizeData(obj)

obj@meta.data <- obj@meta.data %>%  
  mutate(age = case_when(  
    stage == "e55d" ~ -59,  
    stage == "e90d" ~ -24,  
    stage == "0d" ~ 0,  
    stage == "30d" ~ 30,  
    stage == "90d" ~ 90,  
    stage == "180d" ~ 180,  
    #TRUE ~ as.numeric(as.character(age))  
  ))

sel_cts <- c("Muscle stem cell","Capillary endothelial cell",
             "Macrophage","Fibro-adipogenic progenitor cell",
             "Pericyte","Type II x myonuclei","Type I myonuclei",
             "Type II a/b myonuclei")

obj <- subset(obj, subset = celltype %in% sel_cts)

gene_of_interest <- c("ENSSSCG00000029160","ENSSSCG00000056754","LGI2","ZBTB16")
gene_of_interest <- c("ENSSSCG00000029160")
gene_expr <- FetchData(obj, vars = gene_of_interest)
df <- cbind(obj@meta.data, expression = gene_expr[, 1])


df_summary <- df %>%
  group_by(celltype, age) %>%
  summarise(mean_expr = mean(expression, na.rm = TRUE)) %>%
  ungroup()

ggplot(df_summary, aes(x = as.numeric(age), y = mean_expr, color = celltype)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  theme_bw() +
  scale_x_continuous(breaks = unique(df_summary$age))+
  scale_color_manual(values = ct_cols)+
  labs(x = "Stage (days)", y = "log-normalized expression")+
  theme(
    aspect.ratio = 1,
    axis.text = element_text(colour = 'black', size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(colour = 'black', size =18),
    panel.grid.minor = element_blank(),
    legend.position = 'none'
  )

VlnPlot(obj, features = "ZBTB16", group.by = 'celltype', pt.size = 0)
