# Purpose: Supplementary observed expected.

library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)

stages <- c('e55d','e90d','0d','30d','90d','180d')

R_oe.dat <- read.csv('cell_OR/all_R_oe_cl_sta.csv', header = T, row.names = 1, check.names = FALSE)
R_oe.dat <- read.csv('cell_OR/all_R_oe_cl_tis.csv', header = T, row.names = 1, check.names = FALSE)
R_oe.dat <- read.csv('cell_OR/all_R_oe_ct_sta.csv', header = T, row.names = 1 ,check.names = FALSE)
R_oe.dat <- read.csv('cell_OR/all_R_oe_ct_tis.csv', header = T, row.names = 1, check.names = FALSE)

R_oe <- t(R_oe.dat)

max <- ceiling(max(R_oe.dat, na.rm = TRUE))
col_fun = colorRamp2(c(0, max), c("white","firebrick"))
at =  seq(0, max, by = 2)

Heatmap(as.matrix(R_oe),
        show_heatmap_legend = T, 
        cluster_rows = F, 
        cluster_columns = T,
        row_names_side = 'right', 
        show_column_names = T,
        show_row_names = T,
        col = col_fun,
        rect_gp = gpar(col = "dimgray", lwd = 1),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10, rot=45),
        heatmap_legend_param = list(
          title = "R_o/e",
          at =  seq(0, max, by = 2),
          labels =  seq(0, max, by = 2),
          legend_gp = gpar(fill = col_fun(at))
        ),
        row_order = match(stages, rownames(R_oe))
)
