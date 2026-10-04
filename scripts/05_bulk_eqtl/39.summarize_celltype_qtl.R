# Purpose: summarize celltype qtl.

library(ggplot2)
library(dplyr)
library(data.table)

######## sample info
samp_info <- read.csv("bisque_prop_bin.CSV", header = T)
bins <- unique(samp_info$group)

sel_cell_list <- list()
for (i in bins){
  df <- samp_info[samp_info$group==i, ]
  column_means <- colMeans(df[, 2:14], na.rm = TRUE)
  sel_cell_list[[i]] <- names(column_means[column_means > 0.05])
}

sel_cell_list1 <- list()
sel_cell_list1[['bin1']] <- c('fap','mf1','mf2x')
sel_cell_list1[['bin2']] <- c('fap','mf1','mf2x','msc','mf2ab')
sel_cell_list1[['bin3']] <- c('fap','mf1','mf2x','msc','mf2ab')
sel_cell_list1[['bin4']] <- c('fap','mf1','mf2x','msc','mf2ab')

bins <- c('bin1', 'bin2', 'bin3', 'bin4')

##################################### ieQTL statis

df_ct <- fread("eQTL_res/02_acat/ieQTL/02_filter/all_significant_ieQTLs.tsv")
df_ct_flt <- data.frame()
for(i in bins) {
  df <- df_ct[df_ct$group==i, ]
  df1 <- df[df$cell %in% sel_cell_list1[[i]], ]
  df_ct_flt <- rbind(df_ct_flt, df1)
}

bar_data <- df_ct_flt[, .N, by = .(group, cell)]

cell_order <- bar_data %>%
  group_by(cell) %>%
  summarise(TotalGeneCount = sum(N, na.rm = TRUE)) %>%
  arrange(desc(TotalGeneCount)) %>%
  pull(cell)

bar_data$cell <- factor(bar_data$cell, levels = cell_order)

bar_data_ie <- bar_data

# bar plot


# eGene number & sample size × cell fraction


egene_count <- df_ct_flt %>%
  distinct(pheno_id, group, cell) %>%
  group_by(group, cell) %>%
  summarise(eGene_count = n(), .groups = "drop")

library(tidyr)

cell_frac_long <- samp_info %>%
  pivot_longer(cols = -c(Sample, group), names_to = "cell", values_to = "cell_fraction")

cell_frac_summary <- cell_frac_long %>%
  group_by(group, cell) %>%
  summarise(n_sample = n(),
            mean_fraction = mean(cell_fraction),
            frac_weighted = n_sample * mean_fraction,
            .groups = "drop")

cell_map <- c(
  "adipocyte" = "adip", "muscle_stem_cell" = "msc",
  "capillary_endothelial_cell" = "cap",
  "fibro_adipogenic_progenitor_cell" = "fap",
  "lymphatic_endothelial_cell" = "lec",
  "macrophage" = "mac", "pericyte" = "peri",
  "peripheral_glial" = "glial",
  "t_cell" = "tc", "tenocyte" = "teno",
  "type_i_myonuclei" = "mf1",
  "type_ii_a_b_myonuclei" = "mf2ab",
  "type_ii_x_myonuclei" = "mf2x"
)

cell_frac_summary <- cell_frac_summary %>%
  mutate(cell = recode(cell, !!!cell_map))

cor_df <- egene_count %>%
  inner_join(cell_frac_summary, by = c("group", "cell"))

cor_test <- cor.test(cor_df$eGene_count, cor_df$frac_weighted, method = "pearson")

ct_cols = c("mac"="#FFFF00","fap"="#75b947", "msc"="#c2fea2","tc"="#fde3db",
            "mf1"="#d96666","mf2ab"="#ec7979","mf2x"="#ff8c8c",
            "adip"="#b26a3a","peri"="#ff966e","teno"="#ffa37a",
            "cap"="#47ad85","glial"="deepskyblue", "lec"="#6ec2b4"
)

ggplot(cor_df, aes(x = frac_weighted, y = eGene_count, colour = cell, size = group)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE)+
  scale_color_manual(values = ct_cols)+
  scale_shape_manual(values = c(15, 19, 17, 18, 7))+
  theme_classic()+
  theme(aspect.ratio = 1,
        )+
  xlab('Sample size × Cell-type fraction')+
  ylab('eGene count')

ggplot(cor_df, aes(x = frac_weighted, y = eGene_count, colour = cell, shape = group)) +
  geom_point(size = 8) +
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE)+
  scale_color_manual(values = ct_cols)+
  scale_shape_manual(values = c(15, 19, 17, 18, 7))+
  theme_classic()+
  theme(aspect.ratio = 1,
  )+
  xlab('Sample size × Cell-type fraction')+
  ylab('eGene count')

##### ieQTL cross bins
cts <- unique(df_ct$cell)

dup_rows_df_ct_list <- list()
for (i in cts){
  dup_rows_df_ct_list[[i]] <- df_ct[df_ct$cell==i, ] %>%
    group_by(pheno_id) %>%
    mutate(count = n()) %>%
    filter(n() > 1) %>%
    ungroup()
}

################################################## beQTL statis

df_bk <- fread("eQTL_res/02_acat/beQTL/02_filter/all_significant_beQTLs.tsv")
bar_data <- df_bk[, .N, by = .(group)]

# bar plot


# intersection plot
library(tibble)
library(UpSetR)
library(VennDiagram)
library(ggvenn)

bin_list <- df_bk %>%
  group_by(group) %>%
  summarise(genes = list(unique(pheno_id))) %>%
  deframe()

gene_bin_matrix <- df_bk %>%
  distinct(pheno_id, group) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = group, values_from = value, values_fill = 0)

gene_bin_matrix <- as.data.frame(gene_bin_matrix)

UpSetR::upset(gene_bin_matrix,
      sets = colnames(gene_bin_matrix)[-1],
      sets.bar.color = "#1f77b4",
      order.by = "freq")

UpSetR::upset(gene_bin_matrix,
      sets = c("bin4", "bin3", "bin2", "bin1"), 
      keep.order = TRUE,                        
      sets.bar.color = "#E69F00",
      matrix.color = "#56B4E9",
      order.by = "freq",
      text.scale = 1.5)

ggvenn(bin_list,  stroke_color = F, show_percentage = F,)+
  scale_fill_viridis_d(option = "D", direction = 1)

############################ cseQTL
df_cse <- fread("eQTL_res/02_acat/cseQTL/02_filter/all_significant_cseQTLs.tsv")
df_cse_flt <- data.frame()
for(i in bins) {
  df <- df_cse[df_cse$group==i, ]
  df1 <- df[df$cell %in% sel_cell_list[[i]], ]
  df_cse_flt <- rbind(df_cse_flt, df1)
}

bar_data <- df_cse_flt[, .N, by = .(group, cell)]
bar_data_cse <- bar_data


egene_count <- df_cse_flt %>%
  distinct(pheno_id, group, cell) %>%
  group_by(group, cell) %>%
  summarise(eGene_count = n(), .groups = "drop")


cell_frac_long <- samp_info %>%
  pivot_longer(cols = -c(Sample, group), names_to = "cell", values_to = "cell_fraction")

cell_frac_summary <- cell_frac_long %>%
  group_by(group, cell) %>%
  summarise(n_sample = n(),
            mean_fraction = mean(cell_fraction),
            frac_weighted = n_sample * mean_fraction,
            .groups = "drop")

cor_df <- egene_count %>%
  inner_join(cell_frac_summary, by = c("group", "cell"))

cor_test <- cor.test(cor_df$eGene_count, cor_df$frac_weighted, method = "pearson")

ct_cols = c("fibro_adipogenic_progenitor_cell"="#75b947", "muscle_stem_cell"="#c2fea2",
            "type_i_myonuclei"="#d96666","type_ii_a_b_myonuclei"="#ec7979","type_ii_x_myonuclei"="#ff8c8c"
)

###### beQTL & ieQTL
library(tidyverse)
library(cowplot)
library(patchwork)

bin_list_bk <- df_bk %>%
  group_by(group) %>%
  summarise(genes = list(unique(pheno_id))) %>%
  deframe()

bin_list_ct <- df_ct_flt %>%
  group_by(group) %>%
  summarise(genes = list(unique(pheno_id))) %>%
  deframe()

bin_list_cse <- df_cse_flt %>%
  group_by(group) %>%
  summarise(genes = list(unique(pheno_id))) %>%
  deframe()

bin_list <- list()
for (i in bins){
  bin_list[[i]] <- list(bk = bin_list_bk[[i]], ct = bin_list_ct[[i]], cse = bin_list_cse[[i]])
}
ggvenn(bin_list[['bin4']])

venn_list <- readRDS('eQTL_res/02_acat/eGenes_venn_list.rds')

venn_list[[1]]

# QTL summary plots.
library(ggbreak)

bar_data_ie$class <- 'ieqtl'
bar_data_cse$class <- 'cseqtl'
bar_data_cse <- bar_data_cse %>%
  mutate(cell = recode(cell, !!!cell_map))

bin_cols <- c("bin3"="#35B779", "bin1"="#440154", "bin2"="#31688E", "bin4"="#FDE725")

plot_df <- rbind(bar_data_ie, bar_data_cse)

plot_data <- copy(plot_df)

bin_order <- paste0("bin", 1:5)
plot_data[, group := factor(group, levels = bin_order)]

plot_data[, N_plot := fifelse(class == "cseqtl", N, -N)]

ggplot(
  plot_data,
  aes(x = group, y = N_plot, fill = group)
) +
  geom_col(width = 0.72) +
  geom_text(
    aes(
      label = N,
      vjust = ifelse(N_plot >= 0, -0.25, 1.25)
    ),
    size = 3
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4,
    colour = "black"
  ) +
  facet_wrap(
    ~ cell,
    scales = "free_y",
    ncol = 3
  ) +
  scale_y_continuous(
    labels = abs,
    expand = expansion(mult = c(0.18, 0.18))
  ) +
  scale_fill_manual(
    values = bin_cols
  ) +
  labs(
    x = "Developmental bin",
    y = "Number of identified QTLs",
    fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_rect(
      fill = "grey95",
      colour = "grey70"
    ),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )+
  coord_flip()


write.csv(plot_data, 'ieqtl_cseqtl_egene_num.csv')


# ----------- upset coloc -------------

coloc_all <- coloc_all %>%
  mutate(cell = recode(cell, !!!cell_map))
coloc_dt <- as.data.table(coloc_all)

sel_cell_list1 <- list(
  bin1 = c("fap", "mf1", "mf2x"),
  bin2 = c("fap", "mf1", "mf2x", "msc", "mf2ab"),
  bin3 = c("fap", "mf1", "mf2x", "msc", "mf2ab"),
  bin4 = c("fap", "mf1", "mf2x", "msc", "mf2ab")
)

sel_bin_cell <- rbindlist(
  lapply(names(sel_cell_list1), function(b) {
    data.table(
      bin = b,
      cell = sel_cell_list1[[b]]
    )
  })
)

sel_bin_cell

coloc_flt <- merge(
  coloc_dt,
  sel_bin_cell,
  by = c("bin", "cell"),
  all = FALSE
)


pp4_cutoff <- 0.75

coloc_sig <- coloc_flt[
  pp4 >= pp4_cutoff
]

coloc_flt[, .(
  total = .N,
  pp4_05 = sum(pp4 >= 0.5, na.rm = TRUE),
  pp4_08 = sum(pp4 >= 0.8, na.rm = TRUE),
  pp4_09 = sum(pp4 >= 0.9, na.rm = TRUE)
), by = comparison]

coloc_sig[, event_id := paste(
  gene_id,
  bin,
  cell,
  sep = "__"
)]

coloc_sig_unique <- unique(
  coloc_sig[, .(
    event_id,
    gene_id,
    bin,
    cell,
    comparison
  )]
)

upset_dt <- dcast(
  coloc_sig_unique,
  event_id + gene_id + bin + cell ~ comparison,
  fun.aggregate = length,
  value.var = "comparison",
  fill = 0
)


set_cols <- c(
  "beQTL_ieQTL",
  "beQTL_cseQTL",
  "ieQTL_cseQTL"
)

missing_cols <- setdiff(set_cols, colnames(upset_dt))

for (x in missing_cols) {
  upset_dt[, (x) := 0L]
}

upset_dt[, (set_cols) := lapply(.SD, function(x) as.integer(x > 0)),
         .SDcols = set_cols]

# ---------- upset eGenes -------
be_hits <- fread('eQTL_res/02_acat/beQTL/02_filter/all_significant_beQTLs.tsv')
ie_hits <- read_hits(cfg$ie_hits, "ieQTL", cfg$ie_cols)
cse_hits <- read_hits(cfg$cse_hits, "cseQTL", cfg$cse_cols)
qtl_df <- rbind(ie_hits, cse_hits)

qtl_df <- qtl_df %>%
  mutate(cell = recode(cell, !!!cell_map))

qtl_flt <- merge(
  qtl_df,
  sel_bin_cell,
  by = c("bin", "cell"),
  all = FALSE
)

be_df <- data.frame(eGene=bulk_hits$pheno_id, variant=bulk_hits$variant_id, bin=bulk_hits$group, type="beQTL")
qtl_flt <- data.frame(eGene=qtl_flt$gene_id, variant=qtl_flt$variant_id, bin=qtl_flt$bin, type=qtl_flt$type)

qtl_all <- rbind(be_df, qtl_flt)


# ----- plot ------
qtl_dt <- as.data.table(qtl_all)

qtl_gene <- unique(
  qtl_dt[!is.na(eGene) & !is.na(type), .(eGene, type)]
)

upset_dt <- dcast(
  qtl_gene,
  eGene ~ type,
  fun.aggregate = length,
  value.var = "type",
  fill = 0
)

set_cols <- c("cseQTL", "ieQTL", "beQTL")

for (x in setdiff(set_cols, names(upset_dt))) {
  upset_dt[, (x) := 0L]
}

upset_dt[, (set_cols) := lapply(
  .SD,
  function(x) as.integer(x > 0)
), .SDcols = set_cols]

upset_plot_df <- as.data.frame(upset_dt[, ..set_cols])

upset(
  upset_plot_df,
  sets = set_cols,
  keep.order = TRUE,
  order.by = "freq",
  decreasing = TRUE,
  nintersects = NA,
  mb.ratio = c(0.65, 0.35),
  mainbar.y.label = "Number of eGenes",
  sets.x.label = "Total number of eGenes",
  text.scale = c(1.4, 1.2, 1.2, 1.1, 1.2, 1.0),
  point.size = 3,
  line.size = 0.8,
  main.bar.color = "grey30",
  sets.bar.color = c("#4C78A8", "#F7C745", "#4FAEA3")
)


# ------- plot by bin ---------
plot_upset_by_bin <- function(dat, bin_use) {
  x <- unique(
    dat[bin == bin_use, .(eGene, type)]
  )
  
  wide <- dcast(
    x,
    eGene ~ type,
    fun.aggregate = length,
    value.var = "type",
    fill = 0
  )
  
  set_cols <- c("beQTL", "ieQTL", "cseQTL")
  
  for (nm in setdiff(set_cols, names(wide))) {
    wide[, (nm) := 0L]
  }
  
  wide[, (set_cols) := lapply(
    .SD,
    function(z) as.integer(z > 0)
  ), .SDcols = set_cols]
  
  upset(
    as.data.frame(wide[, ..set_cols]),
    sets = set_cols,
    keep.order = TRUE,
    order.by = "freq",
    nintersects = NA,
    mainbar.y.label = paste0("Number of eGenes in ", bin_use),
    sets.x.label = "Total number of eGenes"
  )
}

plot_upset_by_bin(qtl_dt, "bin3")


# ---- jaccard ------

qtl_dt <- as.data.table(qtl_all)

type_order <- c("beQTL", "ieQTL", "cseQTL")
bin_order  <- paste0("bin", 1:4)

# Deduplicate eGenes within each bin and QTL type.
qtl_gene <- unique(
  qtl_dt[
    !is.na(eGene) &
      !is.na(bin) &
      !is.na(type) &
      type %in% type_order,
    .(bin, type, eGene)
  ]
)

qtl_gene[, bin := factor(bin, levels = bin_order)]
qtl_gene[, type := factor(type, levels = type_order)]

calc_set_similarity <- function(dat,
                                bin_col = "bin",
                                type_col = "type",
                                gene_col = "eGene",
                                type_order = c("beQTL", "ieQTL", "cseQTL")) {
  
  dat <- as.data.table(dat)
  
  res <- rbindlist(
    lapply(unique(as.character(dat[[bin_col]])), function(b) {
      
      dat_bin <- dat[get(bin_col) == b]
      
      gene_sets <- setNames(
        lapply(type_order, function(tp) {
          unique(dat_bin[get(type_col) == tp, get(gene_col)])
        }),
        type_order
      )
      
      pairs <- combn(type_order, 2, simplify = FALSE)
      
      rbindlist(lapply(pairs, function(pair) {
        
        set1 <- gene_sets[[pair[1]]]
        set2 <- gene_sets[[pair[2]]]
        
        n1 <- length(set1)
        n2 <- length(set2)
        n_intersection <- length(intersect(set1, set2))
        n_union <- length(union(set1, set2))
        
        jaccard <- if (n_union > 0) {
          n_intersection / n_union
        } else {
          NA_real_
        }
        
        overlap_coefficient <- if (min(n1, n2) > 0) {
          n_intersection / min(n1, n2)
        } else {
          NA_real_
        }
        
        data.table(
          bin = b,
          type1 = pair[1],
          type2 = pair[2],
          comparison = paste(pair, collapse = "–"),
          n_type1 = n1,
          n_type2 = n2,
          n_intersection = n_intersection,
          n_union = n_union,
          Jaccard = jaccard,
          overlap_coefficient = overlap_coefficient
        )
      }))
    })
  )
  
  res[, bin := factor(bin, levels = bin_order)]
  res[, comparison := factor(
    comparison,
    levels = c(
      "beQTL–ieQTL",
      "beQTL–cseQTL",
      "ieQTL–cseQTL"
    )
  )]
  
  return(res[])
}

similarity_dt <- calc_set_similarity(qtl_gene)

similarity_dt

library(scales)

