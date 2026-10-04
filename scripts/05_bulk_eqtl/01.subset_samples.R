# Purpose: subset samples.
analysis_dir <- Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input")

library(dplyr)

dev <- read.table(file.path(analysis_dir, "qtl/bulk/pigGTEx/00_muscle_752/interaction_files/fap_interaction.txt"), header=F, row.names=1)

colnames(dev) <- c('freq','bin')
bins <- c(1, 2, 3, 4)

for(i in bins){
df <- dev[dev$bin==i,]
df$IID <- rownames(df)
df$AITERM <- df$freq
df2 <- df[,c(3:4)]
write.table(df2,paste0('bin',i,'/fap_freq.txt'), quote=F, row.names=F, sep='\t')
df3 <- df2$IID
write.table(df3,paste0('bin',i,'/sample_list.txt'), quote=F, row.names=F, sep='\t', col.names=F)
df2$AITERM <- "0"
df4 <- data.frame(FID=df2$AITERM, IID= df2$IID)
write.table(df4,paste0('bin',i,'/sample_list_plink.txt'), quote=F, row.names=F, sep='\t', col.names=F)
}
