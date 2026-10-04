# Purpose: qc doublets example BBZF0D.
analysis_dir <- Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input")

setwd(file.path(analysis_dir, "qc/BBZF0D/output"))
library(magrittr)
library(Seurat)
library(celda)
library(scater)
library(DoubletFinder)

###############data import
BBZF0D <- Read10X(data.dir=file.path(analysis_dir, "qc/BBZF0D/filtered_feature_bc_matrix"))


###############Ambient RNA removal
BBZF0D.sce <- SingleCellExperiment(list(counts = BBZF0D))
BBZF0D.sce <- decontX(BBZF0D.sce)
saveRDS(BBZF0D.sce,'BBZF0D.sce.rds')

BBZF0D_decont.sce <- BBZF0D.sce[,which(BBZF0D.sce$ decontX_contamination < 0.5)]
BBZF0D_decont.so <- CreateSeuratObject(round(decontXcounts(BBZF0D_decont.sce)))
saveRDS(BBZF0D_decont.so, 'BBZF0D_decont.so.rds')
BBZF0D_obj = BBZF0D_decont.so

#################filter
BBZF0D_obj[['percent.mt']] <- PercentageFeatureSet(object = BBZF0D_obj, pattern = 'ENSSSCG00000018060|ENSSSCG00000018061|ENSSSCG00000018062|ENSSSCG00000018063|ENSSSCG00000018064|ENSSSCG00000018065|ENSSSCG00000018066|ENSSSCG00000018067|ENSSSCG00000018068|ENSSSCG00000018069|ENSSSCG00000018070|ENSSSCG00000018071|ENSSSCG00000018072|ENSSSCG00000018073|ENSSSCG00000018074|ENSSSCG00000018075|ENSSSCG00000018076|ENSSSCG00000018077|ENSSSCG00000018078|ENSSSCG00000018079|ENSSSCG00000018080|ENSSSCG00000018081|ENSSSCG00000018082|ENSSSCG00000018083|ENSSSCG00000018084|ENSSSCG00000018085|ENSSSCG00000018086|ENSSSCG00000018087|ENSSSCG00000018088|ENSSSCG00000018089|ENSSSCG00000018090|ENSSSCG00000018091|ENSSSCG00000018092|ENSSSCG00000018093|ENSSSCG00000018094|ENSSSCG00000018095|ENSSSCG00000018096|^ND1$|^ND2$|^COX1$|^COX2$|^ATP8$|^ATP6$|^COX3$|^ND3$|^ND4L$|^ND4$|^ND5$|^ND6$|^CYTB$')
fivenum(BBZF0D_obj[['percent.mt']][,1])
rb.genes <- rownames(BBZF0D_obj)[grep('^RP[S|L]',rownames(BBZF0D_obj))]
C<-GetAssayData(object = BBZF0D_obj, slot = 'counts')
percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
BBZF0D_obj <- AddMetaData(BBZF0D_obj, percent.ribo, col.name = 'percent.ribo')
BBZF0D_obj@meta.data$ log10GenesPerUMI <- log10(BBZF0D_obj@meta.data$ nFeature_RNA)/log10(BBZF0D_obj@meta.data$ nCount_RNA)
saveRDS(BBZF0D_obj,'BBZF0D_decont_qc.so.rds')

## minGene=200 & maxGene=5000 & pctMT=5 & minCount=500 & maxCount=15000
BBZF0D_obj <- subset(x = BBZF0D_obj, subset = nFeature_RNA > 200 & nFeature_RNA<5000 & nCount_RNA>500 & nCount_RNA<15000 & percent.mt < 5 & log10GenesPerUMI > 0.80)

saveRDS(BBZF0D_obj,'BBZF0D_decont_fliter_pctMT5.so.rds')

##################Doublet removal
BBZF0D_obj <- NormalizeData(object = BBZF0D_obj, normalization.method = 'LogNormalize', scale.factor = 10000)
BBZF0D_obj <- FindVariableFeatures(object=BBZF0D_obj,selection.method = 'vst',nfeatures = 3000)
BBZF0D_obj_scal<- ScaleData(object =BBZF0D_obj,features = VariableFeatures(BBZF0D_obj))
BBZF0D_obj_scal<- RunPCA(object = BBZF0D_obj_scal, npcs = 50, pc.genes = VariableFeatures(object = BBZF0D_obj_scal))
##pc Select
pct <- BBZF0D_obj_scal [['pca']]@stdev / sum( BBZF0D_obj_scal [['pca']]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
pcs

##runUMAP
BBZF0D_obj_scal <- RunUMAP(BBZF0D_obj_scal, dims = 1:pcs)
saveRDS(BBZF0D_obj_scal,'BBZF0D_decont_filter_scal_pctMT5.so.rds')

##doubletfinder
sweep.res.list <- paramSweep_v3(BBZF0D_obj_scal, PCs = 1:pcs, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$ pK[which.max(bcmvn$ BCmetric)] %>% as.character() %>% as.numeric()
pK_bcmvn

annotations <- BBZF0D_obj_scal@meta.data$ seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
DoubletRate = 0.05
nExp_poi <- round(DoubletRate*nrow(BBZF0D_obj_scal@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi
nExp_poi.adj

BBZF0D_doublet_obj <- doubletFinder_v3(BBZF0D_obj_scal, PCs = 1:pcs, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi, reuse.pANN = F, sct = F)
BBZF0D_doubletfinder_obj <- doubletFinder_v3(BBZF0D_doublet_obj, PCs = 1:pcs, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = paste0('pANN_0.25_',pK_bcmvn,'_',nExp_poi.adj), sct = F)
saveRDS(BBZF0D_doubletfinder_obj, 'BBZF0D_decont_fliter_doubletfinder_pctMT5.so.rds')

BBZF0D_decont_filter_singlet <- subset(BBZF0D_doubletfinder_obj, cells= rownames(BBZF0D_doubletfinder_obj@meta.data[BBZF0D_doubletfinder_obj@meta.data[[8]]=='Singlet',]))
saveRDS(BBZF0D_decont_filter_singlet, 'BBZF0D_decont_filter_singlet_pctMT5.so.rds')
