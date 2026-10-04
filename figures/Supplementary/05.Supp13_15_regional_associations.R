# Purpose: Supplementary regional associations.

library(ggplot2)
library(data.table)
library(dplyr)
library(cowplot)

############ manhattan plot
# Supplementary regional association plots.
bins <- c('bin1','bin2','bin3','bin4')


#ENSSSCG00000010258  LY6E  14_73060718_T_A  bin3 bin4
#ENSSSCG00000040535  AIFM2  4_1230263_C_T  bin3 bin4
#ENSSSCG00000007108  KIZ  17_28975311_C_T  bin2 bin3 bin4
#ENSSSCG00000007093  ZNF133  17_26497836_G_A  bin1 bin4


pl_list <- list()
for(i in bins){
  pl_df <- fread(paste0('eQTL_res/04_ld/ENSSSCG00000010258_', i, '_p_ld.csv'))
  pl_df$is_lead <- pl_df$SNP_B == "14_73060718_T_A"
  
  pl_df$R2_group <- cut(pl_df$R2,
                        breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                        labels = c("0.2 > r² ≥ 0.0",
                                   "0.4 > r² ≥ 0.2",
                                   "0.6 > r² ≥ 0.4",
                                   "0.8 > r² ≥ 0.6",
                                   "1.0 > r² ≥ 0.8"),
                        include.lowest = TRUE,
                        right = FALSE)
  
  break_cols <- c("1.0 > r² ≥ 0.8" = "#ff0100", "0.8 > r² ≥ 0.6" = "#ffa304", "0.6 > r² ≥ 0.4" = "#006400",
                  "0.4 > r² ≥ 0.2" = "#8bcfec", "0.2 > r² ≥ 0.0" = "#03018d")
  
  p <- ggplot(pl_df, aes(x = BP_B, y = -log10(pval_g1), shape = is_lead)) +
    geom_point(aes(fill = R2_group), size = 2.5) +
    scale_fill_manual(values = break_cols) +
    scale_shape_manual(values = c(`TRUE` = 23, `FALSE` = 21)) +
    theme_classic() +
    theme(aspect.ratio = 0.5,
          legend.position = "none") +
    labs(y = expression(-log[10](p)), x = "",
         title = i)
  
  pl_list[[i]] <- p
}

#pl_list[[1]] <- p_tmp

wrap_plots(pl_list)
