# Purpose: prepare bulk interaction coloc.
analysis_dir <- Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input")

library(data.table)
library(dplyr)

######## prepare for beQTL
setwd(file.path(analysis_dir, "qtl/coloc"))
Bins <- c("bin1","bin2","bin3","bin4")
path <- NULL

for(i in Bins){
  path[[i]]<-paste0(file.path(analysis_dir, "qtl/mapping/bulk_omiga/01_acat/"), i, "/")
}

for (i in Bins){
    eqtl <- fread(paste0(path[[i]],"muscle_", i, ".cis_qtl.txt.gz"))
    eqtl$is_eGene = eqtl$pval_g1 < eqtl$pval_g1_acat & eqtl$qval_g1 < 0.05
    eGenes <- eqtl[eqtl$is_eGene == TRUE,]

    snp_freq <- fread(file.path(analysis_dir, "qtl/input/muscle_snp_info.frq"))
    snp_freq <- snp_freq[,c(2,5)]
    names(snp_freq) <- c("variant_id","maf")

    exp <- fread(paste0(file.path(analysis_dir, "qtl/input/"), i, "/invTMM/muscle_", i, ".expression.bed.gz"))
    sample_num <- ncol(exp) - 4

    all_nom_qtl <- data.frame()
    for (chr in 1:18) {
        chr_qtl <- fread(paste0(path[[i]],"muscle_", i, ".cis_qtl_pairs.", chr, ".txt.gz"))
        all_nom_qtl <- rbind(all_nom_qtl, chr_qtl)
    }
    select_qtl <- all_nom_qtl[all_nom_qtl$pheno_id %in% eGenes$pheno_id , ]
    select_qtl <- select_qtl[,c(2,2,1,3,7,5,6)]
    names(select_qtl) <- c("rs_id","variant_id","gene_id","tss_distance","pval_nominal","slope","slope_se")
    select_qtl <- select_qtl %>% inner_join(snp_freq, by = "variant_id")
    select_qtl$N <- sample_num

    fwrite(select_qtl, paste0("beQTL_for_coloc/muscle_", i, "_coloc.bed"))
}


########## prepare for ieQTL
setwd(file.path(analysis_dir, "qtl/coloc"))
Bins <- c("bin1","bin2","bin3","bin4")
ct <- c("adip","cap","fap","glial","lec","mac","mf1","mf2ab","mf2x","msc","peri","tc","teno")
path <- NULL

for (i in Bins){
    eqtl <- fread(paste0(file.path(analysis_dir, "qtl/mapping/bulk_omiga/01_acat/"), i, "/muscle_", i, ".cis_qtl.txt.gz"))
    eqtl$is_eGene = eqtl$pval_g1 < eqtl$pval_g1_acat & eqtl$qval_g1 < 0.05
    eGenes <- eqtl[eqtl$is_eGene == TRUE,]

    snp_freq <- fread(file.path(analysis_dir, "qtl/input/muscle_snp_info.frq"))
    snp_freq <- snp_freq[,c(2,5)]
    names(snp_freq) <- c("variant_id","maf")

    exp <- fread(paste0(file.path(analysis_dir, "qtl/input/"), i, "/invTMM/muscle_", i, ".expression.bed.gz"))
    sample_num <- ncol(exp) - 4

    for (j in 1:length(ct)) {
        ieqtl <- fread(paste0(file.path(analysis_dir, "qtl/mapping/bisque_omiga/01_acat/"), i, "/",ct[[j]],"/muscle_", i, "_", ct[[j]], ".cis_qtl.txt.gz"))
        ieqtl$is_ieGene = ieqtl$pval_g2 < ieqtl$pval_g2_acat & ieqtl$qval_g2 < 0.05
        ieGenes <- ieqtl[ieqtl$is_ieGene == TRUE,]
        ieGenes <- ieGenes[ieGenes$pheno_id %in% eGenes$pheno_id,]

        all_nom_qtl <- data.frame()
        for (chr in 1:18) {
            chr_qtl <- fread(paste0(file.path(analysis_dir, "qtl/mapping/bisque_omiga/01_acat/"), i, "/",ct[[j]],"/muscle_", i, "_",ct[[j]],".cis_qtl_pairs.", chr, ".txt.gz"))
            all_nom_qtl <- rbind(all_nom_qtl, chr_qtl)
        }
        select_qtl <- all_nom_qtl[all_nom_qtl$pheno_id %in% ieGenes$pheno_id , ]
        select_qtl <- select_qtl[,c(2,2,1,3,10,8,9)]
        names(select_qtl) <- c("rs_id","variant_id","gene_id","tss_distance","pval_nominal","slope","slope_se")
        select_qtl <- select_qtl %>% inner_join(snp_freq, by = "variant_id")
        select_qtl$N <- sample_num

        fwrite(select_qtl, paste0("ieQTL_for_coloc/muscle_", i, "_", ct[[j]], "_coloc.bed"))
    }
}
