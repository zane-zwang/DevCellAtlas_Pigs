# Purpose: lineage celltype composition.

library(ggplot2)

meta <- read.csv('all_obs_correct.csv', header = T, row.names = 1)
tissues <- c('adipose','cerebrum','duodenum','heart','hypothalamus','liver','muscle')
stages <- c('e55d','e90d','0d','30d','90d','180d')

celltyeps <- c("Peripheral neuron","Oligodendrocyte progenitor cell","Excitatory neuron","Inhibitory neuron",
               "Oligodendrocyte","Astrocyte","Intermediate progenitor cell","Radial glia",
               "Mesothelial cell","Cholangiocyte","Hepatocyte","Enteroendocrine","BEST4+ epithelial","Enterocyte",
               "Intestinal stem cell","Brunners gland cell","Microfold cell","Tuft cell","Goblet cell",
               "Mesenchymal-like stem cell","Liver endothelial cell","Lymphatic endothelial cell",
               "Endocardial endothelial cell","Arterial endothelial cell","Capillary endothelial cell",
               "Venous endothelial cell","Adipose-derived mesenchymal stem cell",
               "Mast cell","Erythroblast","Microglia","Monocyte","Kupffer cell","Macrophage","Dendritic cell",
               "Pre-dendritic cell","NK cell","NKT cell","T cell","B cell","Plasam cell","Vascular smooth muscle cell",
               "Tenocyte","Fibroblast","Adipose stem and progenitor cell","Fibro-adipogenic progenitor cell",
               "Peripheral glial","Hepatic stellate cell","Pericyte","Interstitial cells of Cajal","Smooth muscle cell",
               "Stromal cell","Adipocyte","Cardiomyocyte","Muscle stem cell","Type II a/b myonuclei","Type I myonuclei",
               "Type II x myonuclei"
               )


cl_cols <- c("Endothelial"="#2d9f5f","Epithelial"="#5f2d9f","Erythroid"="#3f3581","Immune"="#fde725",
             "Muscle"="#9f2d2d","Neural"="#2d5f9f","Stromal"="#9f5f2d")

tis_cols <- c('adipose'='#EEA236FF', 'cerebrum'='#357EBDFF', 'duodenum'='#5CB85CFF', 'heart'='#D43F3A99',
              'hypothalamus'='#46B8DAFF', 'liver'='#20854E99' , 'muscle'='#D43F3AFF')

sta_cols <- c("e55d"="#443a83","e90d"="#31688e", "0d"="#21908c","30d"="#35b779","90d"="#8fd744","180d"="#fde725")

meta$celltype <- factor(meta$celltype, levels = as.factor(celltyeps))

ggplot(meta,aes(x=celltype,fill=tissue))+
  geom_bar(position="fill")+
  scale_fill_manual(values = tis_cols)+
  theme_classic()+
  xlab("")+ylab("Tissue constribution (%)")+
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        legend.title = element_blank())

library(dplyr)

meta <- meta %>%  
  mutate(age = case_when(  
    stage == "e55d" ~ -59,  
    stage == "e90d" ~ -24,  
    stage == "0d" ~ 0,  
    stage == "30d" ~ 30,  
    stage == "90d" ~ 90,  
    stage == "180d" ~ 180,))


df <- meta %>%
  group_by(age, celllineage) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(age) %>%
  mutate(prop = cell_count / sum(cell_count))

ggplot(df, aes(x = age, y = prop, color = celllineage)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, span = 1) +  
  labs(x = "", y = "Proportion", title = "Cell Lineage Changes Over Time") +
  theme_classic()+
  scale_color_manual(values = cl_cols)

ggplot(df, aes(x = age, y = prop, color = celllineage)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE) +
  labs(x = "Stage", y = "Cell Count", title = "Cell Lineage Changes Over Time (Quadratic Model)") +
  theme_minimal()+
  scale_color_manual(values = cl_cols)

ggplot(df, aes(x = age, y = prop, color = celllineage)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 6), se = FALSE) + 
  labs(x = "", y = "Proportion", title = "Cell Lineage Proportion Over Time (GAM)") +
  theme_classic()+
  scale_color_manual(values = cl_cols)

################ tissues
tissue_list <- split(meta, meta$tissue)

plots <- lapply(names(tissue_list), function(t) {
  df <- tissue_list[[t]] %>%
    group_by(age, celllineage) %>%
    summarise(cell_count = n(), .groups = "drop") %>%
    group_by(age) %>%
    mutate(prop = cell_count / sum(cell_count))
  
  ggplot(df, aes(x = age, y = prop, color = celllineage)) +
    geom_point() +
    geom_smooth(method = "loess", se = FALSE, span = 1) +
    labs(x = "", y = "Proportion", title = paste("Cell Lineage Changes in", t)) +
    theme_classic()+
    theme(legend.position = 'none')+
    scale_color_manual(values = cl_cols)
})

library(patchwork)

wrap_plots(plots, ncol = 3)
