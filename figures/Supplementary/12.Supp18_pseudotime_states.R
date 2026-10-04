# Purpose: Supplementary pseudotime states.
publication_root <- normalizePath(Sys.getenv("DEVPIGATLAS_ROOT", "."), mustWork = TRUE)

library(org.Hs.eg.db)
library(clusterProfiler)
library(tibble)
library(dplyr)
library(ggplot2)
source(file.path(publication_root, "scripts/07_cross_species/17.cross_species_go_kegg_mixed.R"))


epr_data <- read.csv('trajectory_classification_exploration_table.csv')
epr_data_clean <- na.omit(epr_data)

epr_data_list <- list(early=subset(epr_data_clean, stage_group=='Early'),
                      mid=subset(epr_data_clean, stage_group=='Mid'),
                      Late=subset(epr_data_clean, stage_group=='Late'))


GO_res_list <- lapply(epr_data_list, perform_enrichment)

bp_list <- lapply(GO_res_list, function(df){
  df1 <- df[df$ONTOLOGY=='BP',]
})

df <- bp_list[[3]][order(bp_list[[3]]$p.adjust),][1 : 10,]
df$RichFactor <- sapply(df$GeneRatio, function(x) eval(parse(text = x)))
df$pval_log <- -log10(df$p.adjust)
df <- df %>%
  arrange(pval_log) %>%
  mutate(Description = factor(Description, levels = unique(Description)))
scale_factor <- max(df$pval_log) / max(df$RichFactor)

cols <- c("#9fd4ed","#028bd8","#7f4ea860","#7f4ea8","#fedeab","#ffa923")

ggplot(df, aes(y = reorder(Description, pval_log))) +
  geom_bar(aes(x = pval_log), stat = "identity", fill = "#fedeab") +
  geom_path(aes(x = RichFactor * scale_factor, y = reorder(Description, pval_log), group = 1), 
            color = "black", linewidth = 0.75) +
  geom_point(aes(x = RichFactor * scale_factor), 
             shape = 21, 
             fill = "#ffa923", 
             color = "black", 
             size = 3) +
  scale_x_continuous(
    name = "-log10 (p.adjust)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Rich Factor")
  ) +
  labs(y = "") +
  theme_minimal() +
  theme(
    aspect.ratio = 2.5,
    panel.background = element_rect(fill = NA, color = NA),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 12),
    axis.ticks = element_line(colour = "black", linewidth = 1),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )


df_raw <- read.csv('3_species_olig_meta.csv')
df <- subset(df_raw, subset=stage %in% stage_order)

stage_order <- c('WPC23','WPC24','E55','E90','P0','P4y','P6y','P14y','P20y','P39y',
                 'E18.5','P4','P14','P32',
                 'P30','P90','P180')

df$stage <- factor(df$stage, levels = stage_order)

ggplot(df, aes(x=stage, y=pseudotime_slingshot, fill=species)) +
  geom_boxplot(outlier.size = 0.1) +
  scale_fill_manual(values = c('#7DAEE0','#edab1c','#725f97')) +
  facet_grid(. ~ species, scales = "free_x", space = "free_x") +
  theme_classic() +
  xlab("stage") + ylab("pseudotime") +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

celltype_order <- c("OPC-1", "OPC-2", "pre-OL", "OL-1", "OL-2")
df$celltype <- factor(df$celltype, levels = celltype_order)

ggplot(df, aes(x=celltype, y=pseudotime_slingshot, fill=species)) +
  geom_boxplot(outlier.size = 0.1) +
  scale_fill_manual(values = c('#7DAEE0','#edab1c','#725f97')) +
  facet_grid(. ~ celltype, scales = "free_x", space = "free_x") +
  theme_classic() +
  xlab("stage") + ylab("pseudotime") +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )


mycolors = c('OL-2'='#1f77b4', 'OPC-2'='#ff7f0e', 'OL-1'='#2ca02c', 'OPC-1'='#d62728', 'pre-OL'='#9467bd')

ggplot(df, aes(x = stage, fill = celltype)) +
  geom_bar(position = "fill") +
  facet_grid(. ~ species, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = alpha(mycolors, .8)) +
  theme_classic() +
  xlab("") + ylab("Abundance") +
  theme(
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
  )
