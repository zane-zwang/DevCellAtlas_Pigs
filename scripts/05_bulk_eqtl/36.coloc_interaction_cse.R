# Purpose: coloc interaction cse.

library(data.table)
library(dplyr)
library(coloc)

# beQTL cseQTL coloc
Bins <- c("bin1","bin2","bin3","bin4")
ct <- c("adipocyte", "fibro_adipogenic_progenitor_cell", "lymphatic_endothelial_cell", "muscle_stem_cell", 
        "peripheral_glial", "tenocyte", "type_ii_x_myonuclei", "capillary_endothelial_cell", "macrophage", 
        "pericyte", "t_cell", "type_ii_a_b_myonuclei", "type_i_myonuclei")

cts <- c("adip","fap","lec","msc","glial","teno","mf2x","cap","mac","peri","tc","mf2ab","mf1")

for (i in Bins) {
    for (j in 1:length(ct)) {
        if (file.exists(paste0("cseQTL_for_coloc/muscle_", i, "_", ct[[j]], "_coloc.bed"))) {
            message("Processing ",i, "_", ct[[j]])
            qtl1 <- fread(paste0("ieQTL_for_coloc/muscle_", i, "_", cts[[j]], "_coloc.bed"))
            qtl2 <- fread(paste0("cseQTL_for_coloc/muscle_", i, "_", ct[[j]], "_coloc.bed"))
            names(qtl2) <- c("rs_id","variant_id","gene_id","tss_distance","pval_nominal","slope","slope_se","maf","N")
            gene1 <- unique(qtl1$gene_id)
            gene2 <- unique(qtl2$gene_id)
            gene_list <- intersect(gene1,gene2)
            
            all_summary <- data.frame()
            all_results <- data.frame()
            if (length(gene_list) > 0) {
                for (k in 1:length(gene_list)) {
                    gene_qtl1 <- qtl1[qtl1$gene_id == gene_list[[k]], ]
                    gene_qtl2 <- qtl2[qtl2$gene_id == gene_list[[k]], ]
                    
                    input <- merge(gene_qtl1, gene_qtl2, by="rs_id", all=FALSE, suffixes=c("_celltype_iqtl","_celltype_eqtl"))
                    if (nrow(input) == 0) {
                      next
                      }
                    input <- input[!duplicated(input$rs_id),]
                    input <- input[complete.cases(input[, c("pval_nominal_celltype_iqtl", "pval_nominal_celltype_eqtl", 
                                        "slope_celltype_iqtl", "slope_celltype_eqtl", 
                                        "slope_se_celltype_iqtl", "slope_se_celltype_eqtl", 
                                        "maf_celltype_iqtl", "maf_celltype_eqtl")]), ]
                    
                    if (nrow(input) > 0) {
                        input <- input[!duplicated(input$rs_id),]
                        input$varbeta1 <- input$slope_se_celltype_iqtl ^ 2
                        input$varbeta2 <- input$slope_se_celltype_eqtl ^ 2

                        res <- coloc.abf(
                            dataset1=list(
                                snp = input$rs_id, 
                                pvalues=input$pval_nominal_celltype_iqtl, 
                                beta=input$slope_celltype_iqtl, 
                                varbeta=input$varbeta1, 
                                N=input$N_celltype_iqtl, 
                                type="quant", 
                                MAF=input$maf_celltype_iqtl
                            ),
                            dataset2=list(
                                snp = input$rs_id, 
                                pvalues=input$pval_nominal_celltype_eqtl, 
                                beta=input$slope_celltype_eqtl, 
                                varbeta=input$varbeta2, 
                                N=input$N_celltype_eqtl, 
                                type="quant", 
                                MAF=input$maf_celltype_eqtl
                            )
                        )

                        summary <- data.frame(res$summary)
                        summary$gene_id <- gene_list[[k]]
                        all_summary <- rbind(all_summary, summary)

                        results_df <- res$results
                        results_df$gene_id <- gene_list[[k]]
                        results_df$bin <- i
                        results_df$cell <- ct[[j]]
                        all_results <- rbind(all_results, results_df)
                    }
                }
            }
            write.csv(all_summary, paste0("03_ieQTL_cseQTL_coloc/muscle_", i, "_", ct[[j]],"_ieQTL_cseQTL_coloc.csv"))
            write.csv(all_results, paste0("03_ieQTL_cseQTL_coloc/muscle_", i, "_", ct[[j]],"_ieQTL_cseQTL_coloc_all_snps.csv"))
        }
    }
}
