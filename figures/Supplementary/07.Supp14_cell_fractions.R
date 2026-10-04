# Purpose: Supplementary cell fractions.

library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyverse)

deconv <- read.table('MuscleCIBERSORT-Results.txt', header = T)
sample <- read.table('sample_752.txt')

test <- deconv[deconv$Mixture %in% sample$V1, ]
test$Macrophage <- test$M2.Macrophage + test$M1.Macrophage
test$M1.Macrophage <- NULL
test$M2.Macrophage <- NULL

rownames(test) <- test$Mixture
test$Mixture <- NULL

test$Sample <- rownames(test)

test <- test %>%
  arrange(desc(Type.IIx.myonuclei)) %>%
  mutate(Sample = factor(Sample, levels = Sample))

test_long <- test %>%
  pivot_longer(-Sample, names_to = "CellType", values_to = "Proportion")


ct_cols = c("Macrophage"="#FFFF00","Fibro.adipogenic.progenitor"="#75b947",
            "MuSC"="#c2fea2",
            "Type.I.myonuclei"="#d96666","Type.IIa.b.myonuclei"="#ec7979","Type.IIx.myonuclei"="#ff8c8c",
            "Adipocyte"="#b26a3a","Pericyte"="#ff966e","Tenocyte"="#ffa37a",
            "Endothelial.cell"="#2d9f5f"
)

ggplot(test_long, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = ct_cols)+
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = "Sample", y = "Cell Proportion", fill = "Cell Type")

#### Bisque
sample_bins <- read.table('sample_bins.txt')
colnames(sample_bins) <- c('sample','bin')


bisque <- read.csv('bisque/bulk_props_bisque.csv', header = T, row.names = 1)
df <- as.data.frame(t(bisque))
df$Sample <- rownames(df)

df <- df %>%
  arrange(desc(type_ii_x_myonuclei)) %>%
  mutate(Sample = factor(Sample, levels = Sample))

df_long <- df %>%
  pivot_longer(-Sample, names_to = "CellType", values_to = "Proportion")

base_colors <- brewer.pal(12,"Set3")
extended_colors <- colorRampPalette(base_colors)(13)

ct_cols = c("macrophage"="#FFFF00","fibro_adipogenic_progenitor_cell"="#75b947",
            "muscle_stem_cell"="#c2fea2","t_cell"="#fde3db",
            "type_i_myonuclei"="#d96666","type_ii_a_b_myonuclei"="#ec7979","type_ii_x_myonuclei"="#ff8c8c",
            "adipocyte"="#b26a3a","pericyte"="#ff966e","tenocyte"="#ffa37a",
            "capillary_endothelial_cell"="#47ad85","peripheral_glial"="deepskyblue", "lymphatic_endothelial_cell"="#6ec2b4"
)

ggplot(df_long, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = ct_cols)+
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = "Sample", y = "Cell Proportion", fill = "Cell Type")

