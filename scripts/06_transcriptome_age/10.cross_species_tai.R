# Purpose: cross species tai.
publication_root <- normalizePath(Sys.getenv("DEVPIGATLAS_ROOT", "."), mustWork = TRUE)

library(ggplot2)
library(myTAI)
library(cowplot)
library(DESeq2)
library(patchwork)
library(dplyr)
source(file.path(publication_root, "scripts/05_bulk_eqtl/05.counts_to_tmm_logcpm.R"))

# -------- functions ----------
read_exp_csv <- function(name){
  df <- read.csv(name, header = T, row.names = 1, check.names = F)
  df$GeneID <- rownames(df)
  return(df)
}


# -------- cal tai -------
sc.data <- readr::read_tsv('gene18668_ages.tsv')
sc.PhyloMap <- dplyr::select(sc.data, Phylostratum = rank, GeneID = 'SYMBOL')
colnames(sc.PhyloMap) <- c('phylorank', 'GeneID')

tissues <- c('adipose','cerebrum','duodenum','heart','hypothalamus','liver','muscle')

folder_path <- 'pseudobulk_exp_for_tai/01_tissue/'
file_names <- list.files(path = folder_path, pattern = "\\pseudobulk.csv$", full.names = TRUE)
ts_pseudo_exp_list <- lapply(file_names, read_exp_csv)
names(ts_pseudo_exp_list) <- tissues

stage_map <- c(
  "e55d" = "E55",
  "e90d" = "E90",
  "0d"   = "P0",
  "30d"  = "P30",
  "90d"  = "P90",
  "180d" = "P180"
)

target_order <- c("GeneID", "E55", "E90", "P0", "P30", "P90", "P180")

ts_pseudo_exp_list <- lapply(ts_pseudo_exp_list, function(df) {
  
  old_names <- colnames(df)
  
  new_names <- old_names
  idx <- old_names != "GeneID"
  stage <- sub("^.*//", "", old_names[idx])
  new_names[idx] <- stage_map[stage]
  
  colnames(df) <- new_names
  
  df <- df[, target_order, drop = FALSE]
  
  return(df)
})

ts_pseudo_logCPM_list <- convert_pseudobulk_expr(
  x = ts_pseudo_exp_list,
  method = "CPM",
  gene_col = "GeneID",
  log_base = 2,
  pseudocount = 1
)

ts_pseudo_logCPM_list <- lapply(ts_pseudo_exp_list, function(df) {
  convert_counts_to_tmm_logCPM(
    x = df,
    gene_col = "GeneID",
    prior.count = 1
  )
})

# Pig muscle pseudobulk TAI.
sc.PES <- merge(sc.PhyloMap, ts_pseudo_logCPM_list$muscle, by = 'GeneID')
sc.PES.clean <- sc.PES %>%
  dplyr::filter(!is.na(phylorank)) %>%
  dplyr::select(phylorank, GeneID, everything())


phyex_set <- myTAI::BulkPhyloExpressionSet_from_df(
  sc.PES.clean,
  name = "pig developmental pseudobulk",
  species = "Sus scrofa"
)

myTAI::plot_signature(phyex_set, colour = "#813c85") +
  ggplot2::theme(aspect.ratio = 0.55)

myTAI::plot_signature_transformed(phyex_set)

# --------- cross species (cortex) -----------
df <- read.table('animal_age_data/homo_idmapping_2025_05_27.tsv', header = 1)
df_age <- read.csv('animal_age_data/Homo_sapiens_PhyloMap.csv', header = 1)

homo_df <- merge(df, df_age, by='uniport')
homo.PhyloMap <- dplyr::select(homo_df, Phylostratum = Phylostratum, GeneID = 'gene')

df <- read.table('animal_age_data/mouse_idmapping_2025_05_27.tsv', header = 1)
df_age <- read.csv('animal_age_data/Mus_PhyloMap.CSV', header = 1)

mus_df <- merge(df, df_age, by='uniport')
mus.PhyloMap <- dplyr::select(mus_df, Phylostratum = rank, GeneID = 'gene')

hs_pseudobulk_exp <- read_exp_csv('pseudobulk_exp_for_tai/03_cross_species/human_stage_pseudobulk.csv')

hs_pseudo_logCPM <- convert_pseudobulk_expr(
  x = hs_pseudobulk_exp,
  method = "CPM",
  gene_col = "GeneID",
  log_base = 2,
  pseudocount = 1
)

hs_order <- c('GeneID','WPC18','WPC19','WPC23','WPC24','P0','P4y','P6y','P14y','P20y','P39y')
mm_order <- c('GeneID','E18.5','P4','P14','P32','P90')

hs_pseudo_logCPM <- hs_pseudo_logCPM[, hs_order, drop = F]

hs.PES <- merge(homo.PhyloMap, hs_pseudo_logCPM, by = 'GeneID')

hs.PES.clean <- hs.PES %>%
  dplyr::filter(!is.na(Phylostratum)) %>%
  dplyr::select(Phylostratum, GeneID, everything())

hs_phyex_set <- myTAI::BulkPhyloExpressionSet_from_df(
  hs.PES.clean,
  name = "human developmental pseudobulk",
  species = "Homo sapiens"
)

myTAI::plot_signature(hs_phyex_set, colour = "#7DAEE0") + 
  ggplot2::theme(aspect.ratio = 0.55)


mm_pseudobulk_exp <- read_exp_csv('pseudobulk_exp_for_tai/03_cross_species/mouse_stage_pseudobulk.csv')

mm_pseudo_logCPM <- convert_pseudobulk_expr(
  x = mm_pseudobulk_exp,
  method = "CPM",
  gene_col = "GeneID",
  log_base = 2,
  pseudocount = 1
)


mm_order <- c('GeneID','E18.5','P4','P14','P32','P90')

mm_pseudo_logCPM <- mm_pseudo_logCPM[, mm_order, drop = F]

mm.PES <- merge(mus.PhyloMap, mm_pseudo_logCPM, by = 'GeneID')

mm.PES.clean <- mm.PES %>%
  dplyr::filter(!is.na(Phylostratum)) %>%
  dplyr::select(Phylostratum, GeneID, everything())

mm_phyex_set <- myTAI::BulkPhyloExpressionSet_from_df(
  mm.PES.clean,
  name = "human developmental pseudobulk",
  species = "Homo sapiens"
)

myTAI::plot_signature(mm_phyex_set, colour = "#edab1c") + 
  ggplot2::theme_classic() +
  ggplot2::theme(aspect.ratio = 0.55)

myTAI::plot_contribution(hs_phyex_set)

# --------- olig cross species ----------
pb_exp <- read_exp_csv('pseudobulk_exp_for_tai/04_olig/olig_species_stage_match_pseudobulk.csv')

stage_order <- c('S1','S2','S3','S4','S5','S6','S7','S8','S9','S10')

split_pb_exp_by_species <- function(pb_exp) {
  
  gene_col <- "GeneID"
  exp_cols <- setdiff(colnames(pb_exp), gene_col)
  
  species <- sub("//.*$", "", exp_cols)
  stage   <- sub("^.*//", "", exp_cols)
  
  cols_by_species <- split(exp_cols, species)
  
  res <- lapply(names(cols_by_species), function(sp) {
    cols <- cols_by_species[[sp]]
    df <- pb_exp[, c(cols, gene_col), drop = FALSE]
    st <- sub("^.*//", "", cols)
    colnames(df) <- c(paste0("S", st), gene_col)
    df <- df[, c(gene_col, paste0("S", st)), drop = FALSE]
    
    return(df)
  })
  
  names(res) <- names(cols_by_species)
  return(res)
}

pb_exp_list <- split_pb_exp_by_species(pb_exp)
pb_exp_list <- lapply(pb_exp_list, function(df) {
  keep_stage <- intersect(stage_order, colnames(df))
  df <- df[, c("GeneID", keep_stage), drop = FALSE]
  return(df)
})


pb_logCPM_list <- convert_pseudobulk_expr(
  x = pb_exp_list,
  method = "CPM",
  gene_col = "GeneID",
  log_base = 2,
  pseudocount = 1
)

# human
olig_hs.PES <- merge(homo.PhyloMap, pb_logCPM_list$human, by = 'GeneID')

olig_hs.PES.clean <- olig_hs.PES %>%
  dplyr::filter(!is.na(Phylostratum)) %>%
  dplyr::select(Phylostratum, GeneID, everything())

olig_hs_phyex_set <- myTAI::BulkPhyloExpressionSet_from_df(
  olig_hs.PES.clean,
  name = "human developmental pseudobulk",
  species = "Homo sapiens"
)

myTAI::plot_signature(olig_hs_phyex_set, colour = "#7DAEE0") + 
  ggplot2::theme(aspect.ratio = 0.55)

# mouse
olig_mm.PES <- merge(mus.PhyloMap, pb_logCPM_list$mouse, by = 'GeneID')

olig_mm.PES.clean <- olig_mm.PES %>%
  dplyr::filter(!is.na(Phylostratum)) %>%
  dplyr::select(Phylostratum, GeneID, everything())

olig_mm_phyex_set <- myTAI::BulkPhyloExpressionSet_from_df(
  olig_mm.PES.clean,
  name = "human developmental pseudobulk",
  species = "Homo sapiens"
)

myTAI::plot_signature(olig_mm_phyex_set, colour = "#edab1c") + 
  ggplot2::theme(aspect.ratio = 0.55)

# pig
olig_ss.PES <- merge(sc.PhyloMap, pb_logCPM_list$pig, by = 'GeneID')

olig_ss.PES.clean <- olig_ss.PES %>%
  dplyr::filter(!is.na(phylorank)) %>%
  dplyr::select(phylorank, GeneID, everything())

olig_ss_phyex_set <- myTAI::BulkPhyloExpressionSet_from_df(
  olig_ss.PES.clean,
  name = "human developmental pseudobulk",
  species = "Homo sapiens"
)

myTAI::plot_signature(olig_ss_phyex_set, colour = "#725f97") + 
  ggplot2::theme(aspect.ratio = 0.55)


# --------- olig cross species ----------
species <- c('human', 'mouse', 'pig')

folder_path <- 'pseudobulk_exp_for_tai/04_olig/pseudobulk/'
file_names <- list.files(path = folder_path, pattern = "\\pseudobulk.csv$", full.names = TRUE)
pb_exp_list <- lapply(file_names, read_exp_csv)
names(pb_exp_list) <- species

rename_stage_cols <- function(df, gene_col = "GeneID") {
  stage_cols <- setdiff(colnames(df), gene_col)
  new_stage_cols <- paste0("S", sub("^.*//", "", stage_cols))
  colnames(df)[match(stage_cols, colnames(df))] <- new_stage_cols
  df
}

pb_exp_list <- lapply(pb_exp_list, rename_stage_cols)
pb_exp_list <- lapply(pb_exp_list, function(df) {
  keep_stage <- intersect(stage_order, colnames(df))
  df <- df[, c("GeneID", keep_stage), drop = FALSE]
  return(df)
})

pb_logCPM_list <- convert_pseudobulk_expr(
  x = pb_exp_list,
  method = "CPM",
  gene_col = "GeneID",
  log_base = 2,
  pseudocount = 1
)


# ---------------- muscle bulk ----------------
library(tidyverse)
source(file.path(publication_root, "scripts/06_transcriptome_age/09.make_bulk_tai_object.R"))
bulk_input_dir <- Sys.getenv("TAI_BULK_INPUT_DIR", "path/to/input/bulk_pseudobulk")
pheno <- read.csv(file.path(bulk_input_dir, "bulk_pheno.CSV"), stringsAsFactors = FALSE)
gene_age <- read.delim("gene18668_ages.tsv", stringsAsFactors = FALSE)
bk_exp <- read.csv(file.path(bulk_input_dir, "bulk_counts.csv"), row.names = 1)

res_L <- make_bulk_myTAI_object(
  count_mat = bk_exp,
  pheno = pheno,
  gene_age = gene_age,
  breed_name = "Landrace"
)

res_T <- make_bulk_myTAI_object(
  count_mat = bk_exp,
  pheno = pheno,
  gene_age = gene_age,
  breed_name = "Tongcheng"
)


myTAI::plot_signature(res_L$phyex_obj, colour = "#336700") + 
  ggplot2::theme(aspect.ratio = 0.55)

myTAI::plot_signature(res_T$phyex_obj, colour = "#002fa6") + 
  ggplot2::theme(aspect.ratio = 0.55)


breed_cols <- c(
  "Landrace"  = "#0072B2",  # blue
  "Tongcheng" = "#D55E00",  # vermillion
  "Breed3"    = "#009E73"   # bluish green
)
