# Purpose: prepare cse phenotypes.
publication_root <- normalizePath(Sys.getenv("DEVPIGATLAS_ROOT", "."), mustWork = TRUE)
analysis_dir <- Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input")

library(rtracklayer)
library(edgeR)
library(preprocessCore)
library(RNOmni)
library(data.table)
library(R.utils)
library(SNPRelate)
library(dplyr)
source(file.path(publication_root, "scripts/05_bulk_eqtl/23.gene_rint.R"))

gtf = rtracklayer::import(file.path(publication_root, "resources/ref/Sus_scrofa.Sscrofa11.1.100.chr.gtf"))
gtf = as.data.frame(gtf);dim(gtf)
exon = gtf[gtf$type=="exon",
           c("start","end","gene_id")]

gle = lapply(split(exon,exon$gene_id),function(x){
  tmp=apply(x,1,function(y){
    y[1]:y[2]
  })
  length(unique(unlist(tmp)))
})

gle=data.frame(gene_id=names(gle),
               length=as.numeric(gle))

ct <- readRDS(file.path(analysis_dir, "qtl/bmind/run_bMIND/cts.rds"))

for(j in 1:length(ct)){
  setwd(file.path(analysis_dir, "qtl/bmind/CT_exp"))
  bulk <- read.csv(paste0(ct[[j]], "_exp.csv"), row.names = 1)
  threshold <- 0.1
  min_frac <- 0.2
  nsamples <- ncol(bulk)
  expr_pass <- rowSums(bulk > threshold) >= (min_frac * nsamples)
  bulk_filtered <- bulk[expr_pass, ]

  expr_for_eqtl <- logCPM_to_eqtl_ready(bulk)

  region_annot <- gtf
  geneid = region_annot$gene_id
  expr_matrix = expr_for_eqtl[rownames(expr_for_eqtl) %in% geneid,]
  bed_annot <- region_annot[region_annot$gene_id %in% rownames(expr_matrix), ]
  expr_matrix_ordered <- expr_matrix[match(bed_annot$gene_id, rownames(expr_matrix)), ]
  bed <- cbind(bed_annot, as.data.frame(expr_matrix_ordered))

  bed <- bed[as.character(bed[[1]]) %in% as.character(1:30), ]
  bed[[1]] <- as.numeric(as.character(bed[[1]]))
  bed <- bed[order(bed[[1]], bed$start), ]
  colnames(bed)[1] <- "#Chr"

  bed <- bed[bed$type == "gene", ]
  bed <- bed[bed$gene_biotype %in% c("protein_coding", "lncRNA"), ]
    
  start = bed$start[bed$strand == "-"]
  end = bed$end[bed$strand == "-"]
  bed$start[bed$strand == "-"] = end
  bed$end[bed$strand == "-"] = start
  bed_unique <- bed[!duplicated(bed$gene_id), ]
  rownames(bed_unique) = bed_unique$gene_id
    
  expr_start_col <- which(colnames(bed_unique) == colnames(bulk)[1])
  bed_unique1 <- as.data.frame(bed_unique[, expr_start_col:ncol(bed_unique)])
  bed_unique <- bed_unique[, 1:(expr_start_col - 1)]
  colnames(bed_unique1) <- colnames(bulk)

  bed <- cbind(bed_unique[,c(1:3,10)], bed_unique1)
  names(bed)[1:4] <- c("#Chr","start","end","gene_id")

  setwd(file.path(analysis_dir, "qtl/bmind/pre_for_cseqtl"))
  if (!dir.exists(ct[[j]])) {
    dir.create(ct[[j]])
  }
  setwd(paste0(file.path(analysis_dir, "qtl/bmind/pre_for_cseqtl/"), ct[[j]]))
  fwrite(bed, file = "expr_tmm_inv.bed", sep = "\t")
}
