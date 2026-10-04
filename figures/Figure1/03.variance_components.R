# Purpose: variance components.

library(ggsci)
library(ggplot2)
library(scales)


pvcaObj <- readRDS('pvca_data/tissue_stage_celltype_pvca.rds')
pvcaObj <- data.frame(category=pvcaObj[['label']],value=t(pvcaObj[['dat']]))

pvcaObj$category <- factor(pvcaObj$category, 
                           levels=c('stage','stage:celltype','tissue:stage',
                                    'tissue','tissue:celltype','resid','celltype'))

ggplot(pvcaObj,aes(x=category, y=value))+
  geom_bar(stat = 'identity', fill="#3a6ca8", width = 0.5)+
  geom_text(aes(label = value), vjust = -0.5, color = "black", size=3) +
  theme_classic()+
  xlab("")+ylab("Weighted average proportion variance")+
  theme(axis.text.x = element_text(colour = "black", angle = 45, hjust = 1))
