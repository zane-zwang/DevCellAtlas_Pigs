# Purpose: coloc bulk interaction.

library(data.table)
library(dplyr)
library(coloc)

# beQTL ieQTL coloc
Bins <- c("bin1","bin2","bin3","bin4")
ct <- c("adip","cap","fap","glial","lec","mac","mf1","mf2ab","mf2x","msc","peri","tc","teno")

for (i in Bins) {
    qtl1 <- fread(paste0("beQTL_for_coloc/muscle_", i, "_coloc.bed"))

    for (j in 1:length(ct)) {
        if (file.exists(paste0("ieQTL_for_coloc/muscle_", i, "_", ct[[j]], "_coloc.bed"))) {
            message("Processing ",i, "_", ct[[j]])
            qtl2 <- fread(paste0("ieQTL_for_coloc/muscle_", i, "_", ct[[j]], "_coloc.bed"))
            gene1 <- unique(qtl1$gene_id)
            gene2 <- unique(qtl2$gene_id)
            gene_list <- intersect(gene1,gene2)
            
            all_summary <- data.frame()
            all_results <- data.frame()
            if (length(gene_list) > 0) {
                for (k in 1:length(gene_list)) {
                    gene_qtl1 <- qtl1[qtl1$gene_id == gene_list[[k]], ]
                    gene_qtl2 <- qtl2[qtl2$gene_id == gene_list[[k]], ]
                    
                    input <- merge(gene_qtl1, gene_qtl2, by="rs_id", all=FALSE, suffixes=c("_bulk_eqtl","_celltype_iqtl"))
                    input <- input[!duplicated(input$rs_id),]
                    input <- input[complete.cases(input[, c("pval_nominal_bulk_eqtl", "pval_nominal_celltype_iqtl", 
                                        "slope_bulk_eqtl", "slope_celltype_iqtl", 
                                        "slope_se_bulk_eqtl", "slope_se_celltype_iqtl", 
                                        "maf_bulk_eqtl", "maf_celltype_iqtl")]), ]
                    
                    if (nrow(input) > 0) {
                        input <- input[!duplicated(input$rs_id),]
                        input$varbeta1 <- input$slope_se_bulk_eqtl ^ 2
                        input$varbeta2 <- input$slope_se_celltype_iqtl ^ 2

                        res <- coloc.abf(
                            dataset1=list(
                                snp = input$rs_id, 
                                pvalues=input$pval_nominal_bulk_eqtl, 
                                beta=input$slope_bulk_eqtl, 
                                varbeta=input$varbeta1, 
                                N=input$N_bulk_eqtl, 
                                type="quant", 
                                MAF=input$maf_bulk_eqtl
                            ),
                            dataset2=list(
                                snp = input$rs_id, 
                                pvalues=input$pval_nominal_celltype_iqtl, 
                                beta=input$slope_celltype_iqtl, 
                                varbeta=input$varbeta2, 
                                N=input$N_celltype_iqtl, 
                                type="quant", 
                                MAF=input$maf_celltype_iqtl
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
            write.csv(all_summary, paste0("01_beQTL_ieQTL_coloc/muscle_", i, "_", ct[[j]],"_beQTL_ieQTL_coloc.csv"))
            write.csv(all_results, paste0("01_beQTL_ieQTL_coloc/muscle_", i, "_", ct[[j]],"_beQTL_ieQTL_coloc_all_snps.csv"))
        }
    }
}
