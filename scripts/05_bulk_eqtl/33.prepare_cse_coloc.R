# Purpose: prepare cse coloc.
analysis_dir <- Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input")

library(data.table)
library(dplyr)

########## prepare for cell type eQTL
setwd(file.path(analysis_dir, "qtl/coloc"))
Bins <- c("bin1","bin2","bin3","bin4")
ct <- c("adipocyte", "fibro_adipogenic_progenitor_cell", "lymphatic_endothelial_cell", "muscle_stem_cell", 
        "peripheral_glial", "tenocyte", "type_ii_x_myonuclei", "capillary_endothelial_cell", "macrophage", 
                "pericyte", "t_cell", "type_ii_a_b_myonuclei", "type_i_myonuclei")

for (i in Bins) {
    setwd(paste0(file.path(analysis_dir, "qtl/mapping/cse_omiga/01_acat/"), i))
    snp_freq <- fread(file.path(analysis_dir, "qtl/input/muscle_snp_info.frq"))
    snp_freq <- snp_freq[,c(2,5)]
    names(snp_freq) <- c("variant_id","maf")

    for (j in 1:length(ct)) {
        exp <- fread(paste0(file.path(analysis_dir, "qtl/mapping/cse_omiga/00_data/exp_for_cseqtl/"), ct[[j]], "/", i, "/expr_tmm_inv.bed"))
        sample_num <- ncol(exp) - 4

        eqtl <- fread(paste0(file.path(analysis_dir, "qtl/mapping/cse_omiga/01_acat/"), i, "/", ct[[j]], "/muscle_", i, "_", ct[[j]], ".cis_qtl.txt.gz"))
        eqtl$chr <- sub("_.*", "", eqtl$variant_id)
        chr_num <- unique(eqtl$chr)
        eqtl$is_eGene = eqtl$pval_g1 < eqtl$pval_g1_acat &
                        eqtl$qval_g1 < 0.05
        eGenes <- eqtl[eqtl$is_eGene == TRUE,]
                                                                                    
        all_nom_qtl <- data.frame()
        for (chr in chr_num) {
            chr_qtl <- fread(paste0(file.path(analysis_dir, "qtl/mapping/cse_omiga/01_acat/"), i, "/",ct[[j]],"/muscle_", i, "_",ct[[j]],".cis_qtl_pairs.", chr, ".txt.gz"))
            all_nom_qtl <- rbind(all_nom_qtl, chr_qtl)
        }
        select_qtl <- all_nom_qtl[all_nom_qtl$pheno_id %in% eGenes$pheno_id,]
        select_qtl <- select_qtl[,c(2,2,1,3,7,5,6)]
        names(select_qtl) <- c("rs_id","variant_id","gene_id","tss_distance","pval_nominal","slope","slope_se")
        select_qtl <- select_qtl %>% inner_join(snp_freq, by = "variant_id")
        select_qtl$N <- sample_num
        
        fwrite(select_qtl, paste0(file.path(analysis_dir, "qtl/coloc/cseQTL_for_coloc/muscle_"), i, "_", ct[[j]], "_coloc.bed"))
    }
}
