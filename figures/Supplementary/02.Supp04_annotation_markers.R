# Purpose: Supplementary annotation markers.

library(Seurat)
library(ggplot2)

############ plot feature dotplot
### muscle
annotation_input_dir <- Sys.getenv("ANNOTATION_INPUT_DIR", "path/to/input/annotation")
obj <- readRDS(file.path(annotation_input_dir, "muscle_logNorm.rds"))

marker_list <- list(
    adipocyte = c("ADIPOQ", "PLIN1"),
    capillary_endothelial_cell = c("PECAM1", "VWF"),
    fibro_adipogenic_progenitor_cell = c("PDGFRA", "COL3A1"),
    lymphatic_endothelial_cell = c("LYVE1", "FLT4"),
    macrophage = c("CD163", "CTSS"),
    muscle_stem_cell = c("PAX7", "DPP4"),
    pericyte = c("RGS5", "PDGFRB"),
    peripheral_glial = c("CDH19", "MPZ"),
    t_cell = c("CD247", "CD3E"),
    tenocyte = c("TNMD", "MKX"),
    type_i_myonuclei = c("MYH7", "TNNT1"),
    type_ii_a_b_myonuclei = c("MYH1", "TNNC2"),
    type_ii_x_myonuclei = c("MYH1", "MYH4")
)


features_use <- unique(unlist(marker_list))

DotPlot(obj, features = features_use, group.by = "celltype", cols = c("#f1f1f1", "firebrick")) +
theme_bw()+
theme(panel.grid.major = element_blank()) +
RotatedAxis()+
coord_flip()
  
ggsave('muscle_fea_dot.pdf', width=6, height=6)

### adipose
obj <- readRDS(file.path(annotation_input_dir, "adipose_logNorm.rds"))

marker_list <- list(
    adipocyte = c("ADIPOQ", "PLIN1"),
    adipose_derived_mesenchymal_stem_cell = c("COX1", "LUM"),
    adipose_stem_and_progenitor_cell = c("PDGFRA"),
    arterial_endothelial_cell = c("PECAM1","GJA5"),
    b_cell = c("BANK1"),
    capillary_endothelial_cell = c("PECAM1"),
    dendritic_cell = c("FLT3"),
    erythroblast = c("EPB41", "ALAS2"),
    lymphatic_endothelial_cell = c("LYVE1", "FLT4"),
    macrophage = c("CD163", "C1QB"),
    mast_cell = c("HPGDS", "KIT"),
    mesothelial_cell = c("KRT8"),
    pericyte = c("RGS5", "PDGFRB"),
    peripheral_glial = c("CDH19", "MPZ"),
    stromal_cell = c("KCNN3"),
    t_cell = c("CD247", "CD3E"),
    type_ii_a_b_myonuclei = c("MYH7", "TNNC2"),
    venous_endothelial_cell = c("ACKR1")
)

features_use <- unique(unlist(marker_list))

DotPlot(obj, features = features_use, group.by = "celltype", cols = c("#f1f1f1", "firebrick")) +
theme_bw()+
theme(panel.grid.major = element_blank()) +
RotatedAxis()+
coord_flip()

ggsave('adipose_fea_dot.pdf', width=7, height=7)

### cerebrum
marker_list <- list(
    astrocyte = c("GFAP", "AQP4"),
    capillary_endothelial_cell = c("PECAM1", "VWF"),
    excitatory_neuron = c("SATB2"),
    inhibitory_neuron = c("GAD1", "GAD2"),
    intermediate_progenitor_cell = c("HES6", "EOMES"),
    microglia = c("P2RY12", "C1QB"),
    oligodendrocyte = c("MBP", "MOG"),
    oligodendrocyte_progenitor_cell = c("PDGFRA", "CSPG4"),
    pericyte = c("RGS5", "PDGFRB"),
    radial_glia = c("SOX2", "HES5"),
    vascular_smooth_muscle_cell = c("VIM", "COL1A2")
)

features_use <- unique(unlist(marker_list))

DotPlot(obj, features = features_use, group.by = "celltype", cols = c("#f1f1f1", "firebrick")) +
theme_bw()+
theme(panel.grid.major = element_blank()) +
RotatedAxis()+
coord_flip()

ggsave('cerebrum_fea_dot.pdf', width=6, height=6)

### duodenum
marker_list <- list(
    adipocyte = c("ADIPOQ", "PLIN1"),
    "best4+_epithelial" = c("BEST4", "GUCA2A"),
    brunners_gland_cell = c("A4GNT", "TFF2"),
    ductal_like_epithelial_cell = c("PKHD1", "KRT7"),
    enterocyte = c("APOA4", "FABP1"),
    enteroendocrine = c("CHGA", "INS"),
    erythroblast = c("EPB41", "ALAS2"),
    goblet_cell = c("AGR2", "CLCA1"),
    interstitial_cells_of_cajal = c("KIT", "ETV1"),
    intestinal_stem_cell = c("LGR5", "CCL25"),
    lymphatic_endothelial_cell = c("LYVE1", "FLT4"),
    macrophage = c("CD163", "CTSS"),
    mesothelial_cell = c("KRT8"),
    microfold_cell = c("GP2"),
    peripheral_glial = c("CDH19"),
    peripheral_neuron = c("ELAVL3", "SYT1"),
    plasam_cell = c("MZB1", "JCHAIN"),
    smooth_muscle_cell = c("ACTA2", "MYH11"),
    stromal_cell = c("PDGFRA", "BMP4", "KCNN3"),
    t_cell = c("CD3D", "CD247"),
    tuft_cell = c("POU2F3", "AVIL"),
    venous_endothelial_cell = c("PECAM1", "VWF")
)

features_use <- unique(unlist(marker_list))

DotPlot(obj, features = features_use, group.by = "celltype", cols = c("#f1f1f1", "firebrick")) +
theme_bw()+
theme(panel.grid.major = element_blank()) +
RotatedAxis()+
coord_flip()

ggsave('duodenum_fea_dot.pdf', width=8, height=8)


### heart
marker_list <- list(
    arterial_endothelial_cell = c("PECAM1", "EFNB2", "GJA5"),
    capillary_endothelial_cell = c("KDR", "EMCN"),
    cardiomyocyte = c("MYH7", "MYL2"),
    dendritic_cell = c("HLA-DRA", "IRF8"),
    endocardial_endothelial_cell = c("NPR3", "NRG1"),
    erythroblast = c("EPB41", "ALAS2"),
    fibroblast = c("COL1A1", "DCN"),
    lymphatic_endothelial_cell = c("LYVE1", "FLT4"),
    macrophage = c("CD163", "C1QA"),
    mast_cell = c("KIT", "MS4A2"),
    mesothelial_cell = c("WT1"),
    monocyte = c("CTSS", "LYZ"),
    pericyte = c("RGS5", "PDGFRB"),
    peripheral_glial = c("CDH19"),
    peripheral_neuron = c("ELAVL3", "SYT1"),
    plasam_cell = c("JCHAIN", "PAX5"),
    pre_dendritic_cell = c("P2RY14"),
    t_cell = c("CD3D","CD3E")
)

features_use <- unique(unlist(marker_list))

DotPlot(obj, features = features_use, group.by = "celltype", cols = c("#f1f1f1", "firebrick")) +
theme_bw()+
theme(panel.grid.major = element_blank()) +
RotatedAxis()+
coord_flip()

ggsave('heart_fea_dot.pdf', width=7, height=7)


### hypothalamus
marker_list <- list(
    astrocyte = c("GFAP", "AQP4"),
    capillary_endothelial_cell = c("PECAM1", "EMCN"),
    excitatory_neuron = c("STMN2", "SYT1"),
    inhibitory_neuron = c("GAD1", "GAD2"),
    intermediate_progenitor_cell = c("EGFR"),
    microglia = c("P2RY12", "C1QB"),
    oligodendrocyte = c("MBP", "MOG"),
    oligodendrocyte_progenitor_cell = c("PDGFRA", "CSPG4"),
    pericyte = c("RGS5", "PDGFRB"),
    vascular_smooth_muscle_cell = c("VIM", "COL1A2")
)

features_use <- unique(unlist(marker_list))

DotPlot(obj, features = features_use, group.by = "celltype", cols = c("#f1f1f1", "firebrick")) +
theme_bw()+
theme(panel.grid.major = element_blank()) +
RotatedAxis()+
coord_flip()

ggsave('hypothalamus_fea_dot.pdf', width=6, height=6)


### liver
marker_list <- list(
    b_cell = c("PAX5"),
    cholangiocyte = c("CFTR", "KRT8"),
    dendritic_cell = c("CST3", "SLA-DQB1"),
    erythroblast = c("EPB41", "ALAS2"),
    hepatic_stellate_cell = c("RELN", "COL3A1"),
    hepatocyte = c("ALB", "APOA1"),
    kupffer_cell = c("VSIG4", "C1QB"),
    liver_endothelial_cell = c("WNT2", "PECAM1"),
    lymphatic_endothelial_cell = c("FLT4"),
    nk_cell = c("KLRK1", "KLRB1"),
    nkt_cell = c("CD3E", "IL7R")
)

features_use <- unique(unlist(marker_list))

DotPlot(obj, features = features_use, group.by = "celltype", cols = c("#f1f1f1", "firebrick")) +
theme_bw()+
theme(panel.grid.major = element_blank()) +
RotatedAxis()+
coord_flip()

ggsave('liver_fea_dot.pdf', width=6, height=6)
