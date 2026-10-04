# Purpose: Supplementary tai summary.

library(ggplot2)
library(myTAI)
library(DESeq2)
library(cowplot)
library(dplyr)
library(tidyr)

col_stage = c("e55d"="#443a83","e90d"="#31688e", "0d"="#21908c","30d"="#35b779","90d"="#8fd744","180d"="#fde725")

pseudoExp_cell_tai_list <- readRDS('tai_ensembl/pseudo_exp_raw/cell/pseudoExp_cell_tai_list_log.rds')
pseudoExp_ts_tai_list <- readRDS('tai_ensembl/pseudo_exp_raw/tissue/pseudo_ts_tai_list_log.rds')

#### ct_plot
adipose_cell_tai <- pseudoExp_cell_tai_list$adipose

cell_line_plot <- function(lst){
  plot_data <- bind_rows(
    lapply(names(lst), function(cell_type) {
      data.frame(
        Cell_Type = cell_type,
        Stage = names(lst[[cell_type]]),
        Value = as.numeric(lst[[cell_type]])
      )
    })
  )
  plot_data$Stage <- factor(plot_data$Stage, levels = c("e55d", "e90d", "X0d", "X30d", "X90d", "X180d"))
  ggplot(plot_data, aes(x = Stage, y = Value, color = Cell_Type, group = Cell_Type)) +
    geom_line(size = 1) + 
    geom_point(size = 2) +
    scale_color_brewer(palette = 'Set1')+
    theme_classic() +
    labs(
      x = "",
      y = "TAI",
      color = "cell type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = 'bottom',
          legend.title = element_blank())
}

cell_line_plot(adipose_cell_tai)

####### plot mean & SD
# define function
mean_SD_plot <- function(lst){
  
  # trans to long_df
  long_data <- bind_rows(
    lapply(names(lst), function(cell_type) {
      data.frame(
        Cell_Type = cell_type,
        Stage = names(lst[[cell_type]]),
        Value = as.numeric(lst[[cell_type]])
      )
    })
  )
  
  # calculate
  summary_data <- long_data %>%
    group_by(Stage) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE)
    ) %>%
    ungroup()
  summary_data$Stage <- factor(summary_data$Stage, levels = c("e55d", "e90d", "X0d", "X30d", "X90d", "X180d"))
  
  #plot
  p <- ggplot(summary_data, aes(x = Stage, y = Mean)) +
    geom_line(group = 1, linewidth = 1) +     
    geom_point(size = 3) +              
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), 
                  width = 0.15, linewidth = 0.8) +          
    theme_classic() +
    labs(y = "TAI", x = "") +
    theme(aspect.ratio = 0.55,
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

pseudoExp_cell_tai_plot_list <- lapply(pseudoExp_cell_tai_list, mean_SD_plot)

pseudoExp_cell_tai_plot_list$muscle

tissues <- c('Adipose','Cerebrum','Duodenum','Heart','Hypothalamus','Liver','Skeletal muscle')
plot_grid(plotlist = pseudoExp_cell_tai_plot_list, ncol = 3, 
          labels = tissues, label_size = 12)


## celllineage
library(ggpubr)
library(ggsignif)

cl_cols <- c("Endothelial"="#2d9f5f","Epithelial"="#5f2d9f",
             "Erythroid"="#3f3581","Muscle"="#9f2d2d",
             "Immune"="#fde725","Neural"="#2d5f9f","Stromal"="#9f5f2d")
atlas_metadata_file <- file.path(Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input"), "atlas/all_obs.csv")
meta <- read.csv(atlas_metadata_file, header = T, row.names = 1)
meta$safe_cell <- gsub("[^a-zA-Z0-9_\\-]", "_", trimws(meta$celltype))
cl_ct <- unique(meta[,c('celllineage', 'safe_cell')])


tai_df <- as.data.frame(do.call(rbind, pseudo_ct_tai_list))
rownames(tai_df) <- gsub("_stage_pseudo_Exp\\.csv$", "", rownames(tai_df))
tai_df$safe_cell <- rownames(tai_df)

cl_ct_tai <- merge(tai_df, cl_ct, by='safe_cell')

cl_tai_long <- cl_ct_tai %>%
  pivot_longer(
    cols = e55d:X180d,
    names_to = "Stage",
    values_to = "TAI"
  )

cl_tai_long$Stage <- factor(cl_tai_long$Stage, levels = c("e55d", "e90d", "X0d", "X30d", "X90d", "X180d"))
celllineages <- c("Endothelial","Neural","Stromal","Immune","Muscle","Epithelial","Erythroid")
cl_tai_long$celllineage <- factor(cl_tai_long$celllineage, levels = celllineages)

cl_tai_long <- cl_tai_long %>%
  mutate(Stage_lineage = paste(Stage, celllineage, sep = "_"))

cl_tai_long$Stage_lineage <- factor(
  cl_tai_long$Stage_lineage,
  levels = expand.grid(
    Stage = levels(cl_tai_long$Stage),
    celllineage = levels(cl_tai_long$celllineage)
  ) %>%
    arrange(Stage, celllineage) %>%
    mutate(Stage_lineage = paste(Stage, celllineage, sep = "_")) %>%
    pull(Stage_lineage)
)

cl_tai_long <- cl_tai_long[cl_tai_long$celllineage != "Erythroid", ]

ggplot(cl_tai_long, aes(x = Stage_lineage, y = TAI, color = celllineage)) +
  geom_jitter(width = 0.2, size = 2, shape = 1) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black", size = 0.3) +
  theme_bw() +
  scale_color_manual(values = cl_cols)+
  labs(x = NULL, y = "TAI", title = "TAI variation among lineages and stages") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )


ggplot(cl_tai_long, aes(x = Stage_lineage, y = TAI, colour = celllineage)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.4) +
  geom_point(position = position_jitter(width = 0.15), size = 1.6, shape = 1) +
  scale_color_manual(values = cl_cols)+
  theme_classic() +
  labs(x = NULL, y = "TAI") +
  theme(aspect.ratio = 0.2,
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none")


#### signif

comp_p_cal <- function(df){
  df_A <- df %>% filter(Stage_lineage %in% groupA)
  df_B <- df %>% filter(Stage_lineage %in% groupB)
  
  test_result <- wilcox.test(df_A$TAI, df_B$TAI)
  pval <- test_result$p.value
  return(pval)
}

# e55d
groupA <- c('e55d_Stromal','e55d_Endothelial','e55d_Neural')
groupB <- c('e55d_Epithelial')
comp_p_cal(cl_tai_long)

# e90d
groupA <- c('e90d_Stromal','e90d_Endothelial','e90d_Neural','e90d_Immune','e90d_Muscle')
groupB <- c('e90d_Epithelial')
pval <- comp_p_cal(cl_tai_long)

# 0d
groupA <- c('X0d_Stromal','X0d_Endothelial','X0d_Neural','X0d_Immune')
groupB <- c('X0d_Epithelial','X0d_Muscle')
comp_p_cal(cl_tai_long)

# 30d
groupA <- c('X30d_Stromal','X30d_Endothelial','X30d_Neural','X30d_Immune')
groupB <- c('X30d_Epithelial','X30d_Muscle')
comp_p_cal(cl_tai_long)

# 90d
groupA <- c('X90d_Stromal','X90d_Endothelial','X90d_Neural','X90d_Immune')
groupB <- c('X90d_Epithelial','X90d_Muscle')
comp_p_cal(cl_tai_long)

# 180d
groupA <- c('X180d_Stromal','X180d_Endothelial','X180d_Neural','X180d_Immune')
groupB <- c('X180d_Epithelial','X180d_Muscle')
comp_p_cal(cl_tai_long)


ggplot(cl_tai_long, aes(x = Stage_lineage, y = TAI, colour = celllineage)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.4) +
  geom_point(position = position_jitter(width = 0.15), size = 1.6, shape = 1) +
  scale_color_manual(values = cl_cols)+
  geom_signif(
    annotations = paste0("p = ", signif(pval, 2)),
    y_position = max(cl_tai_long$TAI) + 0.1,
    xmin = 10, xmax = 13
  )+
  theme_classic() +
  labs(x = NULL, y = "TAI") +
  theme(aspect.ratio = 0.2,
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none")


ggplot(cl_tai_long, aes(x = as.numeric(Stage), y = TAI, color = celllineage)) +
  geom_smooth(method = "loess", se = FALSE, size = 1) +
  scale_x_continuous(breaks = 1:6, labels = levels(cl_tai_long$Stage)) +
  scale_color_manual(values = cl_cols) +
  labs(x = "Developmental Stage", y = "TAI", title = "TAI trends across development") +
  theme_bw() +
  theme(
    aspect.ratio = 0.5,
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )

## tissue tai trend
library(reshape2)
ts_tai_df <- as.data.frame(do.call(rbind, pseudoExp_ts_tai_list))
ts_tai_df$tissue <- rownames(ts_tai_df)
ts_tai_long <- melt(ts_tai_df, id.vars = "tissue",
                variable.name = "stage",
                value.name = "value")

ggplot(ts_tai_long, aes(x=tissue, y=value))+
  geom_boxplot()
