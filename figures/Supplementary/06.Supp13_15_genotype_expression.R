# Purpose: Supplementary genotype expression.

library(data.table)
library(ggplot2)
library(cowplot)

bins <- c('bin1','bin2','bin3','bin4')

#### beQTL
geno_exp_plot_list <- list()
for (i in bins) {
  df_plot <- readRDS(paste0('eQTL_res/geno_exp/ENSSSCG00000007108_beQTL/', i, '_17_28087112_T_G_df_plot.rds'))
  gene_id_of_interest <- "ENSSSCG00000007108"
  
  fit <- lm(expression ~ genotype, data = df_plot)
  beta <- round(coef(fit)[2], 3)
  
  geno_exp_plot_list[[i]] <- ggplot(df_plot, aes(x = genotype_label, y = expression)) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 0.5) +
    annotate("text", x = 1.5, y = max(df_plot$expression, na.rm = TRUE), 
             label = paste0("β = ", beta), hjust = 0, size = 4) +
    theme_classic() +
    theme(aspect.ratio = 0.7)+
    labs(x = "17:28087112 G>T", 
         y = "TMM expression",
         subtitle = paste(i, ":", gene_id_of_interest))
}

wrap_plots(geno_exp_plot_list)


#### ieQTL
df_plot <- readRDS('eQTL_res/geno_exp/ENSSSCG00000034313_ieqtl_mf2x/bin3_6_53654971_T_C_df_plot.rds')

mf2x_frac <- data.frame(IID = samp_info$Sample, frac = samp_info$type_ii_x_myonuclei)
inverse_normal_transform <- function(x) {
  r <- rank(x, ties.method = "average")
  p <- (r - 0.5) / length(r)
  qnorm(p)
}
mf2x_frac$frac_int <- inverse_normal_transform(mf2x_frac$frac)

ieqtl_df_plot <- merge(df_plot, mf2x_frac, by='IID')

ENSSSCG00000034313_ieqtl_geno_exp_plot[[2]] <- ggplot(ieqtl_df_plot, aes(x=frac_int, y=expression, color=genotype_label))+
  geom_point(alpha = 0.8, size = 1)+
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE)+
  scale_color_manual(values= c("firebrick", "forestgreen", "steelblue"))+
  theme_classic()+
  theme(aspect.ratio = 1)+
  labs(x='Mf2x enrichment', y='SELENOW TMM expression')

wrap_plots(ENSSSCG00000034313_ieqtl_geno_exp_plot)
