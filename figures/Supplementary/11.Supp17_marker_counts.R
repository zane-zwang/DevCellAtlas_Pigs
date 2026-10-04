# Purpose: Supplementary marker counts.

library(ggplot2)
library(dplyr)
library(purrr)


human_spec_deg <- readRDS('species_specific_genes/human_spec_deg.rds')
mouse_spec_deg <- readRDS('species_specific_genes/mouse_spec_deg.rds')
pig_spec_deg <- readRDS('species_specific_genes/pig_spec_deg.rds')
conserved_deg <- readRDS('species_specific_genes/conserved_deg.rds')
n_gene_conserved_specific <- readRDS('species_specific_genes/n_gene_conserved_specific.rds')


specific_count_by_celltype <- bind_rows(
  human_spec_deg %>% mutate(Category = "Human_specific"),
  mouse_spec_deg %>% mutate(Category = "Mouse_specific"),
  pig_spec_deg   %>% mutate(Category = "Pig_specific")
) %>%
  count(Category, cluster, name = "n_gene")

conserved_count_by_celltype <- conserved_deg %>%
  count(cluster, name = "n_gene") %>%
  mutate(Category = "Conserved")

gene_count_by_celltype <- bind_rows(
  specific_count_by_celltype,
  conserved_count_by_celltype
)


ggplot(gene_count_by_celltype,
       aes(x = cluster, y = n_gene, fill = Category)) +
  geom_col(position = "dodge", width = 0.75) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    aspect.ratio = 0.6
  ) +
  labs(x = NULL, y = "Number of genes")


gene_count_by_celltype_sp <- subset(gene_count_by_celltype, subset = Category %in% c('Human_specific','Mouse_specific','Pig_specific'))
gene_count_by_celltype_con <- subset(gene_count_by_celltype, subset = Category == 'Conserved')

ggplot(gene_count_by_celltype_sp,
       aes(x = cluster, y = n_gene, fill = Category)) +
  geom_col(position = "dodge", width = 0.75) +
  theme_classic() +
  scale_fill_manual(values = c('#7DAEE0','#edab1c','#725f97')) +
  scale_y_continuous(limits = c(0, 450)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    aspect.ratio = 0.6
  ) +
  labs(x = NULL, y = "Number of genes") +
  coord_flip()


ggplot(gene_count_by_celltype_con,
       aes(x = cluster, y = n_gene, fill = Category)) +
  geom_col(position = "dodge", width = 0.75) +
  theme_classic() +
  scale_fill_manual(values = "#9E9E9E") +
  scale_y_continuous(limits = c(0, 450)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    aspect.ratio = 0.6
  ) +
  labs(x = NULL, y = "Number of genes") +
  coord_flip()
