# Purpose: clock prediction assessment.

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(stringr)
library(patchwork)

ct_cols = c("Arterial endothelial cell"="#3aa672","Capillary endothelial cell"="#47ad85","Endocardial endothelial cell"="#54b498","Liver endothelial cell"="#61bba1",
            "Lymphatic endothelial cell"="#6ec2b4","Venous endothelial cell"="#7bc9c7",
            
            "BEST4+ epithelial"="#6c3aaa","Cholangiocyte"="#7947b5","Enterochromaffin"="#8654c0","Enterocyte"="#9361cb","Enteroendocrine"="#a06ed6","Goblet cell"="#ad7be1",
            "Hepatocyte"="#ba88ec","Intestinal stem cell"="#c795f7","Mesothelial cell"="#d4a2e2","Microfold cell"="#e1afe3","Tuft cell"="#eebce4",
            
            "Erythroblast"="darkslateblue",
            
            "B cell"="#f53d00","Dendritic cell"="#ff5400","Kupffer cell"="#ff7300","Macrophage"="#fc9a00","Mast cell"="#FFCC00","Microglia"="#FFFF00","Monocyte"="#FFF44F",
            "NK cell"="#FFFF99","NKT cell"="#FFFDD0","Plasam cell"="#fdf7bd","T cell"="#fde3db","Pre-dendritic cell"="khaki",
            
            "Adipocyte"="#7f4a1b", "Adipose stem and progenitor cell"="#bf6a3f", "Fibro-adipogenic progenitor cell"="#d47f4f", "Fibroblast"="#e89a6b", "Fibroblast-like cell"="#fcb281", 
            "Hepatic stellate cell"="#f0b16f", "Interstitial cells of Cajal"="#d59d58", "Stromal cell"="#9e7532", "Smooth muscle cell"="#bb8e57", "Muscle stem cell"="#7c5a33", 
            "Myofibroblast"="#5d4023", "Pericyte"="#4c2f0d","Tenocyte"="#624510", "Vascular leptomeningeal cell"="#785b13","Mesenchymal-like stem cell"="tan", 
            "Vascular smooth muscle cell"="goldenrod","Adipose-derived mesenchymal stem cell"="peru",
            
            "Cardiomyocyte"="#7f1b1b","Type I myonuclei"="#d96666","Type II a/b myonuclei"="#ec7979","Type II x myonuclei"="#ffb3b3",
            
            "Astrocyte"="cornflowerblue","Excitatory neuron"="royalblue","Peripheral glial"="deepskyblue","Inhibitory neuron"="#27408b","Intermediate progenitor cell"="powderblue",
            "Peripheral neuron"="lightskyblue","Oligodendrocyte"="steelblue","Oligodendrocyte progenitor cell"="dodgerblue","Radial glia"="slateblue",
            
            "Brunners gland cell"="#636363",
            
            "Endothelial"="#2d9f5f","Epithelial"="#5f2d9f","Erythroid"="#3f3581","Muscle"="#9f2d2d","Immune"="#fde725","Neural"="#2d5f9f","Stromal"="#9f5f2d",
            'adipose'='#EEA236FF', 'cerebrum'='#357EBDFF', 'duodenum'='#5CB85CFF', 'heart'='#D43F3A99', 'hypothalamus'='#46B8DAFF', 'liver'='#20854E99' , 'muscle'='#D43F3AFF')


names(ct_cols) <- as.character(names(ct_cols)) %>%
  tolower() %>%  
  str_replace_all(pattern = " ", replacement = "_") %>% 
  str_replace_all(pattern = "-", replacement = "_") %>% 
  str_replace_all(pattern = ",", replacement = "") %>%
  str_replace_all(pattern = "/", replacement = "_")

meta <- read.csv('all_obs.csv', row.names = 1)
tissues <- c('adipose','cerebrum','duodenum','heart','hypothalamus','liver','muscle')

# select cells 
cell_counts <- meta %>%
  count(tissue, stage, celltype)

valid_celltypes <- cell_counts %>%
  filter(n > 100) %>%
  group_by(tissue, celltype) %>%
  summarise(n_stages = n_distinct(stage), .groups = "drop") %>%
  filter(n_stages >= 4)

vld_cells <- split(valid_celltypes$celltype, valid_celltypes$tissue)

ts_cells <- lapply(vld_cells, function(cells){
  cells <- as.character(cells) %>% tolower() %>%
    str_replace_all(pattern = " ", replacement = "_") %>%
    str_replace_all(pattern = "-", replacement = "_") %>%
    str_replace_all(pattern = ",", replacement = "") %>%
    str_replace_all(pattern = "/", replacement = "_")
  return(cells)
})


folder_path <- "aging_clock/01_pred_assessment/"
file_names <- list.files(path = folder_path, pattern = "\\.rds$", full.names = TRUE)

pred_list <- lapply(file_names, readRDS)
names(pred_list) <- tissues

  for(i in 1:length(tissues)){
    pred_list[[i]] <- lapply(pred_list[[i]], function(df){
      df_sel <- df[df$Celltype %in% ts_cells[[i]], ]
      return(df_sel)
    })
  }

pred_merge_list <- list()

for(i in 1:length(tissues)){
  for(j in 1:length(pred_list[[i]])){
    pred_list[[i]][[j]]$sample <- paste0("Replicate_", j)
    pred_list[[i]][[j]]$r2 <- as.numeric(pred_list[[i]][[j]]$r2)
  }
  pred_merge_list[[i]] <- bind_rows(pred_list[[i]])
}

names(pred_merge_list) <- tissues

summary_df_list <- lapply(pred_merge_list, function(df){
  summary_df <- df %>%
    group_by(Celltype) %>%
    summarise(
      median_r2 = median(r2),
      sd_r2 = sd(r2),
      se_r2 = sd_r2 / sqrt(n()),
      .groups = "drop"
    ) %>%
    arrange(median_r2)
  summary_df$Celltype <- factor(summary_df$Celltype, levels = summary_df$Celltype)
  return(summary_df)
})


pred_plot_list <- lapply(summary_df_list, function(df){
  p <- ggplot(data = df) +
    geom_point(aes(x = median_r2, y = Celltype, color=Celltype),
               size = 2.5) +
    geom_errorbarh(aes(xmin = median_r2 - se_r2, xmax = median_r2 + se_r2, y = Celltype, color=Celltype),
                   height = 0, linewidth = 1)+
    scale_color_manual(values = ct_cols)+
    theme_bw()+
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13),
          axis.text.y = element_text(size = 13))+
    xlab("")+ylab("")
})

pred_plot_list[[7]]

wrap_plots(pred_plot_list, ncol = 3)


library(purrr)
library(cowplot)

summary_df <- summary_df_list %>%
  imap_dfr(~ .x %>%
             mutate(tissue = .y, 
                    Celltype_tissue = paste0(as.character(Celltype), "_", .y)) %>%
             select(Celltype_tissue, Celltype, tissue, median_r2, sd_r2, se_r2)
  )

# globel
summary_df <- summary_df %>%
  arrange(median_r2) %>%
  mutate(Celltype_tissue = factor(Celltype_tissue, levels = unique(Celltype_tissue)))

# tissue
summary_df <- summary_df %>%
  group_by(tissue) %>%
  arrange(median_r2, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(Celltype_tissue = factor(Celltype_tissue, levels = unique(Celltype_tissue)))

tile_df <- summary_df %>%
  mutate(y_pos = as.numeric(Celltype_tissue))

p_2 <- ggplot(data = tile_df[tile_df$tissue %in% c('muscle'),]) +
  geom_point(aes(x = median_r2, y = Celltype_tissue, color=Celltype),
             size = 1.2) +
  geom_errorbarh(aes(xmin = median_r2 - se_r2, xmax = median_r2 + se_r2, y = y_pos, color=Celltype),
                 height = 0, linewidth = 1)+
  scale_color_manual(values = ct_cols)+
  theme_bw()+
  theme(aspect.ratio = 6,
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13),
        axis.text.y = element_text(size = 13),
        panel.grid = element_blank(),
        )+
  xlab("")+ylab("")

p_2


# selected cell types
tissues <- c('muscle','muscle','hypothalamus','heart','duodenum','adipose')
sel_cells <- c('fibro_adipogenic_progenitor_cell','macrophage','microglia','cardiomyocyte','intestinal_stem_cell','adipose_stem_and_progenitor_cell')

sel_pred_plot_list <- list()

for(i in 1:length(tissues)){
  pred <- readRDS(paste0('aging_clock/00_lobo_pred/',tissues[i],'_lobo_predictions_repeat_2.rds'))
  df <- pred[pred$Celltype==sel_cells[i],]
  
  correlation <- cor(df$age, df$Prediction, method = "pearson")
  r2 <- correlation^2
  n <- length(df$age)
  t_stat <- correlation * sqrt((n - 2) / (1 - correlation^2))
  p_value <- 2 * pt(-abs(t_stat), df = n - 2)
  mae <- mean(abs(df$age - df$Pred))
  
  sel_pred_plot_list[[i]] <- ggplot(df, aes(x = age, y = Prediction)) +
    geom_point(color = "grey") + 
    geom_jitter(color = "grey", width = 2, height = 0)+
    geom_smooth(method = "lm", color = "blue", se = TRUE) +
    theme_bw() +  
    scale_x_continuous(breaks = unique(df$age))+
    labs(x = "Actual Stage", y = "Predicted Stage")+
    annotate("text", x = -8, y = 150, 
             label = as.expression(
               bquote(atop(.(sel_cells[i]), R^2 == .(round(r2, 3))))
             ),
             size = 4.2) +
    theme(
      aspect.ratio = 1,
      axis.text = element_text(colour = 'black', size = 12),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 12)
    )
}

sel_pred_plot_list[[2]]

wrap_plots(sel_pred_plot_list, ncol = 3)


###### cell number & accuracy
obs <- read.csv('all_obs.csv', header = T, row.names = 1)

celltype_counts <- obs %>%
  group_by(tissue, celltype) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  mutate(Celltype_tissue = paste0(celltype, "_", tissue)) %>%
  select(Celltype_tissue, cell_count)

celltype_counts$Celltype_tissue <- as.character(celltype_counts$Celltype_tissue) %>%
  tolower() %>%
  str_replace_all(pattern = " ", replacement = "_") %>%
  str_replace_all(pattern = "-", replacement = "_") %>%
  str_replace_all(pattern = ",", replacement = "") %>%
  str_replace_all(pattern = "/", replacement = "_")

cell_number_r2 <- merge(celltype_counts, summary_df, by = 'Celltype_tissue')

lm_model <- lm(median_r2 ~ cell_count, data = cell_number_r2)

intercept <- coef(lm_model)[1]
slope <- coef(lm_model)[2]

r2 <- summary(lm_model)$r.squared
pval <- summary(lm_model)$coefficients[2, 4]


ggplot(cell_number_r2, aes(cell_count, median_r2))+
  geom_point(color = "steelblue")+
  geom_smooth(method = "lm", color = "firebrick", se = TRUE)+
  theme_bw()+
  theme(
    aspect.ratio = 0.7,
    axis.title = element_text(size = 12),
    axis.text = element_text(colour = "black", size = 12),
    panel.grid = element_blank()
  )+
  labs(
    x = "Number of cells",
    y = expression("Predictive performance (" * R^2 * ")")
  )+
  annotate(
    "text",
    x = Inf, y = -Inf, hjust = 1.1, vjust = -1.2,
    label = paste0(
                   "R² = ", round(r2, 3)),
    size = 4
  )


##### muscle
muscle_rep5 <- pred_merge_list$muscle
muscle_rep5$MedianAbsErr <- as.numeric(muscle_rep5$MedianAbsErr)

p1 <- ggplot(muscle_rep5, aes(x=Celltype, y=MedianAbsErr, colour = Celltype))+
  geom_boxplot()+
  scale_color_manual(values = ct_cols)+
  theme_classic()+
  theme(
    aspect.ratio = 1,
    legend.position = 'none',
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

p2 <- ggplot(muscle_rep5, aes(x=Celltype, y=r2, colour = Celltype))+
  geom_boxplot()+
  scale_color_manual(values = ct_cols)+
  theme_classic()+
  theme(
    aspect.ratio = 1,
    legend.position = 'none',
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

p1 + p2
