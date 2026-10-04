# Purpose: Supplementary quality replicates.

library(ggplot2)
library(dplyr)
library(tidyr)

meta <- read.csv('all_obs.csv', header = T, row.names = 1)

tis_cols <- c('adipose'='#EEA236FF', 'cerebrum'='#357EBDFF', 'duodenum'='#5CB85CFF', 'heart'='#D43F3A99',
              'hypothalamus'='#46B8DAFF', 'liver'='#20854E99' , 'muscle'='#D43F3AFF')

sta_cols <- c("e55d"="#443a83","e90d"="#31688e", "0d"="#21908c","30d"="#35b779","90d"="#8fd744","180d"="#fde725")

tissues <- c("adipose", "cerebrum", "duodenum", "heart", "hypothalamus", "liver", "muscle")
stages <- c('e55d','e90d','0d','30d','90d','180d')

samples <- c("adipose_e55d","adipose_e90d","adipose_0d","adipose_30d","adipose_90d","adipose_180d",
             "cerebrum_e55d","cerebrum_e90d","cerebrum_0d","cerebrum_30d","cerebrum_90d","cerebrum_180d",
             "duodenum_e55d","duodenum_e90d","duodenum_0d","duodenum_30d","duodenum_90d","duodenum_180d",
             "heart_e55d","heart_e90d","heart_0d","heart_30d","heart_90d","heart_180d",
             "hypothalamus_e55d","hypothalamus_e90d","hypothalamus_0d","hypothalamus_30d","hypothalamus_90d","hypothalamus_180d",
             "liver_e55d","liver_e90d","liver_0d","liver_30d","liver_90d","liver_90d2","liver_180d",
             "muscle_e55d","muscle_e90d","muscle_0d","muscle_30d","muscle_90d","muscle_180d")

meta$sample <- factor(meta$sample, levels = as.factor(samples))

ggplot(meta, aes(x=sample, y=nCount_RNA, fill=stage))+
  geom_boxplot(outlier.shape = NA, width=0.5, size=0.3)+
  scale_fill_manual(values = sta_cols)+
  xlab("Sample")+ylab("UMI count")+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        aspect.ratio = 0.23)

###### simility of liver
# cell prop
cell_proportions <- meta %>%
  filter(sample %in% c("liver_90d", "liver_90d2")) %>%
  group_by(sample, celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(prop = count / sum(count))


ggplot(cell_proportions, aes(x=celltype, y=prop, fill=sample))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_classic()+
  xlab("Cell type")+ylab("Proportion")+
  theme(#legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3))


df_wide <- cell_proportions %>%
  select(sample, celltype, prop) %>%
  pivot_wider(names_from = sample, values_from = prop, values_fill = 0)

pearson_cor <- cor(df_wide$liver_90d, df_wide$liver_90d2, method = "pearson")
spearman_cor <- cor(df_wide$liver_90d, df_wide$liver_90d2, method = "spearman")


ggplot(df_wide, aes(x = liver_90d, y = liver_90d2)) +
  geom_point(aes(color=celltype), size = 3) +
  scale_color_manual(values = ct_cols)+
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  labs(x = "Sample1 Proportion", y = "Sample2 Proportion",
       title = paste("Cell Type Proportion Correlation (Pearson:", round(pearson_cor, 3), ")")) +
  theme_bw()+
  theme(legend.position = 'none')

# exp
pd_smp <- read.csv('all_pseudobulk_sample.csv', header = T, row.names = 1)
pd_liver90d <- pd_smp[,c('liver_90d','liver_90d2')]
cpm_data <- cpm(pd_liver90d)
log_data <- log2(cpm_data + 1)
log_data <- as.data.frame(log_data)


pearson_cor <- cor(log_data$liver_90d, log_data$liver_90d2, method = "pearson")
spearman_cor <- cor(log_data$liver_90d, log_data$liver_90d2, method = "spearman")

ggplot(log_data, aes(x = liver_90d, y = liver_90d2)) +
  geom_point(size = 1) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  labs(x = "Expression level (Sample1)", y = "Expression level (Sample2)",
       title = paste("Correlation (Pearson:", round(pearson_cor, 3), ")")) +
  theme_bw()+
  theme(legend.position = 'none')
