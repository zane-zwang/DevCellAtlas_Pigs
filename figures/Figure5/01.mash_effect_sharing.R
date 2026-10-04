# Purpose: mash effect sharing.

library(data.table)
library(ggplot2)
library(tidyverse)

post_mean <- readRDS("mashr_res/01_beqtl_mashr_output/posterior_means.rds")
lfsr <- readRDS("mashr_res/01_beqtl_mashr_output/lfsr.rds")
pair_annot <- fread("mashr_res/01_beqtl_mashr_output/pair_annot.tsv")

stopifnot(identical(rownames(post_mean), rownames(lfsr)))
Bins <- colnames(post_mean)

lfsr_th <- 0.05
eff_th  <- 0
sig_mat <- (lfsr < lfsr_th) & (abs(post_mean) > eff_th)

sig_cnt_pair <- colSums(sig_mat)

pair_annot2 <- pair_annot %>% distinct(pair_id, gene_id)
stopifnot(all(rownames(sig_mat) %in% pair_annot2$pair_id))
sig_dt <- as.data.table(sig_mat)
sig_dt[, pair_id := rownames(sig_mat)]
sig_long <- melt(sig_dt, id.vars = "pair_id", variable.name = "Bin", value.name = "sig")
sig_long <- sig_long[sig == TRUE][, .(pair_id, Bin)] %>%
  left_join(pair_annot2, by="pair_id")

eGene_cnt <- sig_long %>% distinct(gene_id, Bin) %>% count(Bin, name="n_eGene")

bind_rows(
  data.frame(Bin = names(sig_cnt_pair), n = as.integer(sig_cnt_pair), type="eQTL"),
  data.frame(Bin = eGene_cnt$Bin,        n = eGene_cnt$n_eGene,      type="eGene")
) %>%
  ggplot(aes(Bin, n, fill=type)) +
  geom_col(position = position_dodge(width=0.9)) +
  theme_bw() + 
  theme(aspect.ratio = 0.6,
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x="Bin", y="Count", title="Significant counts per Bin (mash lfsr)")

library(pheatmap)

pm_scaled <- t(scale(t(post_mean)))
cor_mat <- cor(pm_scaled, use="pairwise.complete.obs")
pheatmap(cor_mat, cluster_rows=TRUE, cluster_cols=TRUE,
         main="Bin × Bin correlation (posterior mean)")


library(ComplexUpset)

sig_pair <- (lfsr < lfsr_th) & (abs(post_mean) > eff_th)
sig_dt <- as.data.table(sig_pair)
sig_dt[, pair_id := rownames(sig_pair)]

long <- melt(sig_dt, id.vars="pair_id", variable.name="bin", value.name="sig")
long <- long[sig == TRUE]
long$gene_id <- pair_annot$gene_id[match(long$pair_id, pair_annot$pair_id)]

gene_bin_sig <- long[, .(sig = TRUE), by = .(gene_id, bin)]

gene_bin_mat <- dcast(gene_bin_sig, gene_id ~ bin, value.var = "sig", fill = FALSE)
gene_bin_mat <- gene_bin_mat[, c("gene_id", Bins), with = FALSE]

set_size(8,3)
ComplexUpset::upset(gene_bin_mat, Bins, sort_sets=FALSE)

library(ComplexHeatmap)

sel <- rowSums(sig_mat) > 0
pm_sel <- post_mean[sel,]
pm_sel_scaled <- t(scale(t(pm_sel)))

Heatmap(pm_sel_scaled,
        name = "Z(beta_post)",
        show_row_names = FALSE,
        cluster_rows = TRUE, cluster_columns = TRUE,
        column_title = "Effect patterns across Bins")


library(cowplot)

plot_pair_trend <- function(pid){
  tibble(Bin = Bins,
         beta = as.numeric(post_mean[pid, ]),
         lfsr = as.numeric(lfsr[pid, ])) %>%
    mutate(sig = lfsr < lfsr_th) %>%
    ggplot(aes(x=Bin, y=beta, group=1))+
    geom_line()+
    geom_point(aes(shape=sig))+
    theme_bw()+
    labs(title=paste0("Pair: ", pid),
         y="Posterior mean (beta)", x="Bin")
}

top_pairs <- names(sort(rowSums(sig_mat), decreasing = TRUE))[1:6]
plots <- lapply(top_pairs, plot_pair_trend)
plot_grid(plotlist = plots, ncol = 3)
