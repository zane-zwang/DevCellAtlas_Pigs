# Purpose: stable clock gene overlap.

library(ggplot2)
library(dplyr)
library(UpSetR)
library(ggupset)


publication_root <- normalizePath(Sys.getenv("DEVPIGATLAS_ROOT", "."), mustWork = TRUE)
clock_input_dir <- Sys.getenv("CLOCK_GENE_INPUT_DIR", "path/to/analysis_data/clock/genes")
clock_output_dir <- Sys.getenv("CLOCK_OVERLAP_OUTPUT_DIR", file.path(publication_root, "results/clock/overlap"))
metadata_file <- Sys.getenv("CLOCK_ATLAS_METADATA_CSV", "path/to/input/all_obs.csv")
clock_public_data_dir <- Sys.getenv("CLOCK_PUBLIC_DATA_DIR", "path/to/input/clock/public_data")
importance_file <- Sys.getenv("CLOCK_IMPORTANCE_CSV", "path/to/input/importance_muscle_stem_cell.csv")
dir.create(clock_output_dir, recursive = TRUE, showWarnings = FALSE)
folder_path <- clock_input_dir
tissues <- c('adipose','cerebrum','duodenum','heart','hypothalamus','liver','muscle')

# Match the manuscript eligibility rule used for the clock figure.
metadata <- read.csv(metadata_file, row.names = 1)
eligible_celltypes <- metadata %>%
  count(tissue, stage, celltype) %>%
  filter(n > 100) %>%
  group_by(tissue, celltype) %>%
  summarise(n_stages = n_distinct(stage), .groups = "drop") %>%
  filter(n_stages >= 4)
ts_cells <- lapply(split(eligible_celltypes$celltype, eligible_celltypes$tissue), function(cells) {
  cells <- tolower(as.character(cells))
  cells <- stringr::str_replace_all(cells, " ", "_")
  cells <- stringr::str_replace_all(cells, "-", "_")
  cells <- stringr::str_replace_all(cells, ",", "")
  stringr::str_replace_all(cells, "/", "_")
})


# Input: repeated clock-gene RDS files and cross-species clock-gene lists.
clock_genes_list <- list()
for(tissue in tissues){
  pattern <- paste0("^", tissue, ".*\\.rds$")
  
  file_names <- list.files(path = folder_path, pattern = pattern, full.names = TRUE)
  
  if(length(file_names) > 0) {
    clock_genes_list[[tissue]] <- lapply(file_names, readRDS)
  } else {
    clock_genes_list[[tissue]] <- NULL
  }
}

# filter cells
clock_genes_filter_list <- list()
for(i in names(clock_genes_list)){
  clock_genes_filter_list[[i]] <- list()
  for(j in 1:length(clock_genes_list[[i]])){
    cell_type_names <- names(clock_genes_list[[i]][[j]])
    selected_cells <- cell_type_names[cell_type_names %in% ts_cells[[i]]]
    clock_genes_filter_list[[i]][[j]] <- clock_genes_list[[i]][[j]][selected_cells]
  }
}

# stable freq >= 3
clock_genes_stable_list <- list()
for(i in names(clock_genes_filter_list)){
  cell_gene_set <- list()
  clock_genes_stable_list[[i]] <- list()
  for(c in ts_cells[[i]]){
    cell_gene_set[[c]] <- list()
    for(j in 1:length(clock_genes_filter_list[[i]])){
      cell_gene_set[[c]][[j]] <- clock_genes_filter_list[[i]][[j]][[c]]
    }
    gene_table <- table(unlist(cell_gene_set[[c]]))
    clock_genes_stable_list[[i]][[c]] <- names(gene_table[gene_table >= 3])
  }
}

clock_module <- do.call(
  rbind,
  lapply(names(clock_genes_stable_list), function(tissue) {
    
    x <- clock_genes_stable_list[[tissue]]
    
    do.call(
      rbind,
      lapply(names(x), function(celltype) {
        
        genes <- setdiff(x[[celltype]], "(Intercept)")
        
        data.frame(
          tissue   = tissue,
          celltype = celltype,
          module   = paste(tissue, celltype, sep = "_"),
          hubs     = paste(genes, collapse = "; "),
          n_genes  = length(genes),
          stringsAsFactors = FALSE
        )
      })
    )
  })
)

rownames(clock_module) <- NULL

write.csv(clock_module, file.path(clock_output_dir, 'clock_module.csv'))

gene_list <- lapply(clock_genes_stable_list, function(org) {
  unlist(org)
})

all_genes <- unlist(gene_list)


gene_frequency <- as.data.frame(table(all_genes))
gene_frequency <- gene_frequency[-1,]
colnames(gene_frequency) <- c("gene", "frequency")

frequency_distribution <- gene_frequency %>%
  group_by(frequency) %>%
  summarise(num_genes = n()) %>%
  arrange(frequency)

ggplot(frequency_distribution, aes(x = frequency, y = num_genes)) +
  geom_line(color = "firebrick", size = 1.3) +
  geom_rect(
    aes(xmin = 1, xmax = 3, ymin = 70, ymax = 550),
    fill = '#9fd4ed', alpha = 0.05
  )+
  geom_rect(
    aes(xmin = 3, xmax = 19, ymin = 0, ymax = 70),
    fill = '#fedeab', alpha = 0.05
  )+
  labs(#title = "Gene Frequency Distribution",
       x = "Number of Cell Types",
       y = "Number of Genes") +
  annotate('text', x = 7, y = 300, size = 5,
           label = '90% cell type specific') +
  annotate('text', x = 10, y = 100, size = 5,
           label = '> 3 cell types') +
  theme_classic()+
  theme(axis.text = element_text(size = 15, colour = 'black'),
        axis.title = element_text(size = 15))

##### upset plot
genes_list_ms <- clock_genes_stable_list$muscle
genes_list_ms <- lapply(genes_list_ms, function(x) x[-1])

binary_matrix <- as.data.frame.matrix(table(unlist(genes_list_ms), 
                                            rep(names(genes_list_ms), 
                                                sapply(genes_list_ms, length))))

upset(binary_matrix, 
      sets = names(genes_list_ms), 
      order.by = "freq",
      main.bar.color = "black", 
      sets.bar.color = "black",
      text.scale = c(2, 2, 2, 2, 2, 2)
      )


######### donut plot
library(caret)
library(glmnet)
library(tidyverse)
library(tidyr)
library(ggthemes)

d <- read.csv(importance_file, header = T, row.names = 1)
ggplot(d, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Sign)) +
  geom_rect(color = "black", size = 0.1) +
  geom_label(x=4.2, aes(y=labelPosition, label=Variable), size=2.5) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_tableau()

########## venn
library(VennDiagram)
library(Seurat)

# -------- AS_2025 human brain -----------
human_clock_data_dir <- Sys.getenv("HUMAN_CLOCK_DATA_DIR", "path/to/input/human_clock_data")
folder_path <- human_clock_data_dir
Oligodendrocytes_gene_h <- read.csv(file.path(folder_path, 'Supplementary_Data_5m.csv'))
Oligodendrocytes_gene_h$cell <- 'Oligodendrocytes'
Astrocytes_gene_h <- read.csv(file.path(folder_path, 'Supplementary_Data_5n.csv'))
Astrocytes_gene_h$cell <- 'Astrocytes'
Microglia_gene_h <- read.csv(file.path(folder_path, 'Supplementary_Data_5o.csv'))
Microglia_gene_h$cell <- 'Microglia'
OPCs_gene_h <- read.csv(file.path(folder_path, 'Supplementary_Data_5p.csv'))
OPCs_gene_h$cell <- 'OPCs'
EN_gene_h <- read.csv(file.path(folder_path, 'Supplementary_Data_5q.csv'))
EN_gene_h$cell <- 'EN'
IN_gene_h <- read.csv(file.path(folder_path, 'Supplementary_Data_5r.csv'))
IN_gene_h$cell <- 'IN'
all_gene_h_list <- list(olig = Oligodendrocytes_gene_h, astro = Astrocytes_gene_h, micro = Microglia_gene_h, 
                        OPCs = OPCs_gene_h, EN = EN_gene_h, IN = IN_gene_h)

flt_any_zero <- function(df){
  df1 <- df[!apply(df == 0, 1, any), ]
  df2 <- df1[df1[1, ] != "intercept", ]
  return(df2)
}

brain_clock_gene_h_list_flt <- lapply(all_gene_h_list, flt_any_zero)
brain_clock_gene_h_df <- do.call(rbind, brain_clock_gene_h_list_flt)

write.csv(brain_clock_gene_h_df, file.path(clock_output_dir, 'brain_clock_gene_h_df.csv'), quote = F)

# -------- nature_aging_2022 mouse brain -----------

# ------- nature_aging 2025 human immune -------
b_gene_h <- read.csv(file.path(clock_public_data_dir, 'b.CSV'))
b_gene_h$cell <- 'b'
cd4t_gene_h <- read.csv(file.path(clock_public_data_dir, 'cd4t.CSV'))
cd4t_gene_h$cell <- 'cd4t'
cd8t_gene_h <- read.csv(file.path(clock_public_data_dir, 'cd8t.CSV'))
cd8t_gene_h$cell <- 'cd8t'
mono_gene_h <- read.csv(file.path(clock_public_data_dir, 'mono.CSV'))
mono_gene_h$cell <- 'mono'
nk_gene_h <- read.csv(file.path(clock_public_data_dir, 'nk.CSV'))
nk_gene_h$cell <- 'nk'

immune_clock_gene_h_df <- rbind(b_gene_h, cd4t_gene_h, cd8t_gene_h, mono_gene_h, nk_gene_h)
write.csv(immune_clock_gene_h_df, file.path(clock_output_dir, 'immune_clock_gene_h_df.csv'), quote = F)

immune_clock_gene_h_df <- read.csv(file.path(clock_output_dir, 'immune_clock_gene_h_df.csv'))
brain_clock_gene_homo <- read.csv(file.path(clock_public_data_dir, 'brain_clock_gene_h_homo.csv'))
brain_clock_gene_homo <- unique(brain_clock_gene_homo$ortholog_name)
brain_clock_gene_homo <- na.omit(brain_clock_gene_homo)
brain_clock_gene_homo <- brain_clock_gene_homo[brain_clock_gene_homo != "N/A"]

brain_clock_gene_mus <- read.csv(file.path(clock_public_data_dir, 'brain_clock_gene_m_homo.csv'))
brain_clock_gene_mus <- unique(brain_clock_gene_mus$ortholog_name)
brain_clock_gene_mus <- na.omit(brain_clock_gene_mus)
brain_clock_gene_mus <- brain_clock_gene_mus[brain_clock_gene_mus != "N/A"]


genes_list_brain <- clock_genes_stable_list$cerebrum
genes_list_brain <- lapply(genes_list_brain, function(x) x[-1])
genes_list_brain_freq <- as.data.frame(sort(table(unlist(genes_list_brain)), decreasing = TRUE))
brain_clock_gene_pig <- as.character(unique(genes_list_brain_freq$Var1))

# -------- fisher test ----------
N <- 12833
K <- 1107
M <- 142
x <- 8

fisher.test(matrix(c(x, K-x, M-x, N-K-M+x), nrow=2),
            alternative = "greater")


library(ggvenn)

brain_clock_gene1 <- list("human"=brain_clock_gene_homo, "pig"=brain_clock_gene_pig)
brain_clock_gene2 <- list("mouse"=brain_clock_gene_mus, "pig"=brain_clock_gene_pig)

ggvenn(brain_clock_gene1, show_percentage = T, show_elements = F, label_sep = ",",
       digits = 3, stroke_color = "white",
       fill_color = c("firebrick","#01a087"))


