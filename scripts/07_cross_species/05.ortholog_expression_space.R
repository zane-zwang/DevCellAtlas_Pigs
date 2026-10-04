# Purpose: ortholog expression space.

#### ortholos
library(biomaRt)
library(org.Ss.eg.db)
library(dplyr)
library(stringr)
library(edgeR)

ss <- useMart(biomart = "ensembl", dataset = "sscrofa_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mm <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
hs <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
pig_attributes <- listAttributes(ss)


# human mouse orthologs
human_mouse <- getLDS(
  attributes = c("ensembl_gene_id", "external_gene_name", "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_associated_gene_name", "mmusculus_homolog_orthology_type"),
  filters = "",
  values = "",
  mart = hs,
  attributesL = c("ensembl_gene_id", "external_gene_name"),
  martL = mm
)

colnames(human_mouse) <- c('h.ensembl','h.gene','m.ensembl','m.gene','orthology_type')
human_mouse <- subset(human_mouse, orthology_type == 'ortholog_one2one')
human_mouse <- human_mouse[, 1:4]


human_pig <- getLDS(
  attributes = c("ensembl_gene_id", "external_gene_name", "sscrofa_homolog_ensembl_gene", "sscrofa_homolog_associated_gene_name", "sscrofa_homolog_orthology_type"),
  filters = "",
  values = "",
  mart = hs,
  attributesL = c("ensembl_gene_id", "external_gene_name"),
  martL = ss
)

colnames(human_pig) <- c('h.ensembl','h.gene','s.ensembl','s.gene','orthology_type')
human_pig <- subset(human_pig, orthology_type == 'ortholog_one2one')
human_pig <- human_pig[, 1:4]

hs_mm_ss <- unique(merge(human_mouse, human_pig, by = c("h.ensembl", "h.gene")))

################################################### filter orthologs function

filter_expression_by_ortholog <- function(expr_df, ortholog_df, species = c("pig", "human", "mouse"), gene_col = "gene") {
  
  species <- match.arg(species)
  
  if (!gene_col %in% colnames(expr_df)) {
    stop(paste0("Not find gene_col: '", gene_col, "'"))
  }
  
  expr_df <- expr_df %>%
    as.data.frame() %>%
    rename(gene = all_of(gene_col))
  
  # Define columns corresponding to different species
  species_map <- list(
    pig   = list(symbol = "s.gene",    ensembl = "s.ensembl"),
    human = list(symbol = "h.gene",    ensembl = "h.ensembl"),
    mouse = list(symbol = "m.gene",    ensembl = "m.ensembl")
  )
  sp <- species_map[[species]]
  
  # ---- mapping by gene symbol ----
  matched_by_symbol <- expr_df[expr_df$gene %in% ortholog_df[[sp$symbol]], ]
  matched_by_symbol <- left_join(matched_by_symbol,
                                 ortholog_df[, c("h.gene", "m.gene", "s.gene")],
                                 by = setNames(sp$symbol, "gene"))
  
  # ---- mapping by Ensembl ID ----
  matched_by_ensembl <- expr_df[expr_df$gene %in% ortholog_df[[sp$ensembl]], ]
  matched_by_ensembl <- left_join(matched_by_ensembl,
                                  ortholog_df[, c("h.ensembl", "m.ensembl", "s.ensembl")],
                                  by = setNames(sp$ensembl, "gene"))
  
  # combine & deduplication
  combined <- bind_rows(matched_by_symbol, matched_by_ensembl)
  combined <- combined[!duplicated(combined$gene), ]
  
  colnames(combined)[colnames(combined) == "gene"] <- sp$symbol
  
  combined <- combined %>%
    filter(!is.na(h.gene) & h.gene != "",
           !is.na(m.gene) & m.gene != "",
           !is.na(s.gene) & s.gene != "")
  
  combined <- combined[, -c((ncol(combined) - 1), ncol(combined))]
  
  combined <- combined %>%
    relocate(h.gene, m.gene, s.gene)
  
  return(combined)
}


### hypothalamus
hypo_pig_exp <- read.csv('04_pseudobulk_exp/01_tissues/hypothalamus_stage_pseudo_Exp.csv', header = 1, row.names = 1)
colnames(hypo_pig_exp) <- c('P0','P180','P30','P90','E55','E90')
hypo_pig_exp$s.gene <- rownames(hypo_pig_exp)

hypo_human_exp <- hypo_pseudo_exp_list$human
hypo_human_exp$h.gene <- hypo_human_exp$GeneID
hypo_human_exp$GeneID <- NULL

hypo_mouse_exp <- hypo_pseudo_exp_list$mouse
hypo_mouse_exp$m.gene <- str_to_title(hypo_mouse_exp$GeneID)
hypo_mouse_exp$GeneID <- NULL


pig_exp_orth <- filter_expression_by_ortholog(hypo_pig_exp, ortholog_df = hs_mm_ss, species = "pig", gene_col = "s.gene")
mouse_exp_orth <- filter_expression_by_ortholog(hypo_mouse_exp, ortholog_df = hs_mm_ss, species = "mouse", gene_col = "m.gene")
human_exp_orth <- filter_expression_by_ortholog(hypo_human_exp, ortholog_df = hs_mm_ss, species = "human", gene_col = "h.gene")

common_tmp <- intersect(pig_exp_orth$h.gene, human_exp_orth$h.gene)
common_gene <- intersect(common_tmp, mouse_exp_orth$h.gene)

library(tibble)

exp_for_cross_species_pig <- as_tibble(pig_exp_orth) %>%
  filter(h.gene %in% common_gene) %>%
  dplyr::select(-m.gene, -s.gene) %>%
  tibble::column_to_rownames("h.gene")

exp_for_cross_species_human <- as_tibble(human_exp_orth) %>%
  filter(h.gene %in% common_gene) %>%
  dplyr::select(-m.gene, -s.gene) %>%
  tibble::column_to_rownames("h.gene")

colnames(exp_for_cross_species_human) <- c("PCW5","PCW6","PCW7","PCW9","PCW10","PCW12","PCW13","PCW16","PCW20")

exp_for_cross_species_mouse <- as_tibble(mouse_exp_orth) %>%
  filter(h.gene %in% common_gene) %>%
  dplyr::select(-m.gene, -s.gene) %>%
  tibble::column_to_rownames("h.gene")


### cortex
cortex_pig_exp <- read.csv('04_pseudobulk_exp/01_tissues/cerebrum_stage_pseudo_Exp.csv', header = T, row.names = 1)
colnames(cortex_pig_exp) <- c('P0','P180','P30','P90','E55','E90')
cortex_pig_exp$geneID <- rownames(cortex_pig_exp)

cortex_human_exp <- read.csv('animal_age_data/DevCortex_pseudo_exp/cortex_human_stage_pseudo_Exp_cpm.csv', header = T, row.names = 1)
colnames(cortex_human_exp) <- c('PCW18','PCW19','PCW23','PCW24','P4Y','P6Y','P14Y','P20Y','P39Y','P0')
cortex_human_exp$geneID <- rownames(cortex_human_exp)

cortex_mouse_exp <- read.csv('animal_age_data/DevCortex_pseudo_exp/cortex_mouse_stage_pseudo_Exp_cpm.csv', header = T, row.names = 1)
cortex_mouse_exp$geneID <- rownames(cortex_mouse_exp)


add_ortholog_genes_general <- function(expr_df,
                                       ortholog_df,
                                       gene_col = "geneID",
                                       species  = c("pig", "human", "mouse"),
                                       verbose  = FALSE) {
  species <- match.arg(species)
  if (!gene_col %in% colnames(expr_df)) {
    stop(sprintf("Column '%s' not found in expr_df", gene_col))
  }
  # 1. clean expr_df: rename gene_col to own_id_raw, remove version number, and convert to uppercase
  expr_clean <- expr_df %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    rename(own_id_raw = all_of(gene_col)) %>%
    mutate(
      own_id_nover = str_remove(own_id_raw, "\\.\\d+$"),
      own_id       = toupper(str_trim(own_id_nover))
    ) %>%
    dplyr::select(-own_id_nover)
  
  # 2. Clean the homology table: remove all Ensembl ID version numbers, convert symbols to uppercase and remove spaces
  orth_clean <- ortholog_df %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    mutate(
      h.ensembl_cln = str_remove(h.ensembl, "\\.\\d+$"),
      m.ensembl_cln = str_remove(m.ensembl, "\\.\\d+$"),
      s.ensembl_cln = str_remove(s.ensembl, "\\.\\d+$"),
      h.gene_cln    = toupper(str_trim(h.gene)),
      m.gene_cln    = toupper(str_trim(m.gene)),
      s.gene_cln    = toupper(str_trim(s.gene))
    )
  
  # 3. Determine the species prefix based on species
  if (species == "pig")   { prefix <- "s"; note <- "Pig"   }
  if (species == "human") { prefix <- "h"; note <- "Human" }
  if (species == "mouse") { prefix <- "m"; note <- "Mouse" }
  
  own_ensembl_cln <- paste0(prefix, ".ensembl_cln")
  own_symbol_cln  <- paste0(prefix, ".gene_cln")
  
  # 4. Use "this species Ensembl_cln" and "this species symbol_cln" as keys to build "lookup"
  lookup1 <- orth_clean %>%
    transmute(
      key_id = .data[[own_ensembl_cln]],
      h.gene = h.gene,
      m.gene = m.gene,
      s.gene = s.gene
    )
  lookup2 <- orth_clean %>%
    transmute(
      key_id = .data[[own_symbol_cln]],
      h.gene = h.gene,
      m.gene = m.gene,
      s.gene = s.gene
    )
  lookup <- bind_rows(lookup1, lookup2) %>%
    distinct(key_id, .keep_all = TRUE)
  
  # 5. left_join: align expr_clean$own_id and lookup$key_id
  joined <- expr_clean %>%
    left_join(lookup, by = c("own_id" = "key_id"))
  
  if (verbose) {
    total_rows <- nrow(joined)
    own_sym_col <- paste0(prefix, ".gene")  # pig -> s.gene, human -> h.gene, mouse -> m.gene
    matched_cnt <- sum(!is.na(joined[[own_sym_col]]))
    unm_cnt     <- total_rows - matched_cnt
    message(sprintf("[%s matrix] has %d rows in total, %d rows of homologous species were matched, and %d rows were not matched",
                    note, total_rows, matched_cnt, unm_cnt))
  }
  
  # Retain the input expression columns and append human, mouse and pig gene IDs.
  result <- cbind(
    expr_df,
    data.frame(
      h.gene = joined$h.gene,
      m.gene = joined$m.gene,
      s.gene = joined$s.gene,
      stringsAsFactors = FALSE
    )
  )
  return(result)
}


pig_orth <- add_ortholog_genes_general(
  expr_df     = cortex_pig_exp,
  ortholog_df = hs_mm_ss,
  gene_col    = "geneID",
  species     = "pig",
  verbose     = TRUE
)

human_orth <- add_ortholog_genes_general(
  expr_df     = cortex_human_exp,
  ortholog_df = hs_mm_ss,
  gene_col    = "geneID",
  species     = "human",
  verbose     = TRUE
)

mouse_orth <- add_ortholog_genes_general(
  expr_df     = cortex_mouse_exp,
  ortholog_df = hs_mm_ss,
  gene_col    = "geneID",
  species     = "mouse",
  verbose     = TRUE
)

collapse_to_hgene <- function(annotated_df, gene_col = "h.gene") {
  
  expr_only <- annotated_df %>%
    dplyr::select(-geneID, -h.gene, -m.gene, -s.gene)
  
  hgenes <- annotated_df[[gene_col]]
  
  df_with_key <- cbind(
    data.frame(h.gene = hgenes, stringsAsFactors = FALSE),
    expr_only
  )
  
  df_with_key <- df_with_key %>%
    filter(!is.na(h.gene) & h.gene != "")
  
  collapsed <- df_with_key %>%
    group_by(h.gene) %>%
    summarize(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
    ungroup()
  
  mat <- collapsed %>%
    column_to_rownames(var = "h.gene")
  return(as.data.frame(mat))
}


pig_to_h <- collapse_to_hgene(pig_orth, gene_col = "h.gene")

human_to_h <- collapse_to_hgene(human_orth, gene_col = "h.gene")

mouse_to_h <- collapse_to_hgene(mouse_orth, gene_col = "h.gene")


common_hgenes <- Reduce(intersect, list(
  rownames(pig_to_h),
  rownames(human_to_h),
  rownames(mouse_to_h)
))

exp_for_cross_species_pig <- pig_to_h[common_hgenes, , drop = FALSE]
exp_for_cross_species_human <- human_to_h[common_hgenes, , drop = FALSE]
exp_for_cross_species_mouse <- mouse_to_h[common_hgenes, , drop = FALSE]
