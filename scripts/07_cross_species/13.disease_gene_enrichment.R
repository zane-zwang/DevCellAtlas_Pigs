# Purpose: disease gene enrichment.

library(dplyr)
library(stringr)
library(forcats)
library(purrr)

alz <- read.delim("disease/Alzheimer.tsv", stringsAsFactors = FALSE, na.strings = c("No data", ""))

alz_clean <- alz %>%
  filter(!is.na(symbol), symbol != "") %>%
  mutate(
    symbol = toupper(symbol),
    globalScore = as.numeric(globalScore),
    gwasCredibleSets = as.numeric(gwasCredibleSets)
  ) %>%
  arrange(desc(gwasCredibleSets), desc(globalScore)) %>%
  distinct(symbol, .keep_all = TRUE)

alz_flt <- alz_clean %>%
  filter(gwasCredibleSets >= 0.4)


disease_dir <- "./disease"
disease_files <- list.files(
  disease_dir,
  pattern = "\\.tsv$",
  full.names = TRUE
)

read_disease_gene_set <- function(file, top_n = 150) {
  
  disease_name <- basename(file) %>%
    str_remove("\\.tsv$") %>%
    str_replace_all("_", " ")
  
  df <- read.delim(
    file,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "No data", "no data", "None")
  )
  
  possible_symbol_cols <- c("symbol", "approvedSymbol", "geneSymbol", "targetSymbol")
  symbol_col <- intersect(possible_symbol_cols, colnames(df))[1]
  
  if (is.na(symbol_col)) {
    stop(paste("No gene symbol column found in:", file))
  }
  
  df <- df %>%
    mutate(
      gene = toupper(.data[[symbol_col]])
    ) %>%
    filter(!is.na(gene), gene != "") %>%
    distinct(gene, .keep_all = TRUE)
  
  score_cols <- intersect(c("gwasCredibleSets", "globalScore"), colnames(df))
  
  for (cc in score_cols) {
    df[[cc]] <- suppressWarnings(as.numeric(df[[cc]]))
  }
  
  # Prefer gwasCredibleSets; use globalScore when unavailable.
  if ("gwasCredibleSets" %in% colnames(df) &&
      sum(!is.na(df$gwasCredibleSets)) > 0) {
    
    df_filt <- df %>%
      arrange(desc(gwasCredibleSets), desc(globalScore)) %>%
      slice_head(n = top_n)
    
    filter_method <- paste0("top", top_n, "_gwasCredibleSets")
    
  } else if ("globalScore" %in% colnames(df) &&
             sum(!is.na(df$globalScore)) > 0) {
    
    df_filt <- df %>%
      arrange(desc(globalScore)) %>%
      slice_head(n = top_n)
    
    filter_method <- paste0("top", top_n, "_globalScore")
    
  } else {
    
    warning(paste("No gwasCredibleSets or globalScore found in:", file,
                  "; using all genes."))
    
    df_filt <- df
    filter_method <- "all_genes"
  }
  
  list(
    disease = disease_name,
    genes = unique(df_filt$gene),
    table = df_filt,
    filter_method = filter_method
  )
}

disease_gene_sets <- disease_files %>%
  set_names(basename(.) %>% str_remove("\\.tsv$")) %>%
  map(read_disease_gene_set, top_n = 150)

disease_gene_table <- imap_dfr(
  disease_gene_sets,
  function(x, disease_name) {
    
    x$table %>%
      mutate(
        disease = x$disease,
        rank = seq_len(n())
      ) %>%
      select(
        disease,
        rank,
        gene = symbol,
        gwasCredibleSets,
        globalScore
      ) %>%
      mutate(
        filter_method = x$filter_method
      )
  }
)


target_df <- subset(traj_long_sim, subset = phase =='Early')
target_genes <- unique(target_df$gene)

background_genes <- rownames(pb$counts)


run_fisher_enrichment <- function(target_genes, disease_genes, background_genes) {
  
  target_genes <- intersect(unique(toupper(target_genes)), background_genes)
  disease_genes <- intersect(unique(toupper(disease_genes)), background_genes)
  
  overlap_genes <- intersect(target_genes, disease_genes)
  
  a <- length(overlap_genes)
  b <- length(setdiff(target_genes, disease_genes))
  c <- length(setdiff(disease_genes, target_genes))
  d <- length(setdiff(background_genes, union(target_genes, disease_genes)))
  
  mat <- matrix(
    c(a, b, c, d),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      Target = c("In_target", "Not_in_target"),
      Disease = c("In_disease", "Not_in_disease")
    )
  )
  
  ft <- fisher.test(mat, alternative = "greater")
  
  tibble(
    overlap = a,
    target_size = length(target_genes),
    disease_set_size = length(disease_genes),
    background_size = length(background_genes),
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value,
    overlap_genes = paste(overlap_genes, collapse = ";")
  )
}

enrich_res <- imap_dfr(disease_gene_sets, function(x, nm) {
  
  run_fisher_enrichment(
    target_genes = target_genes,
    disease_genes = x$genes,
    background_genes = background_genes
  ) %>%
    mutate(
      disease = x$disease,
      filter_method = x$filter_method
    )
}) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    minus_log10_padj = -log10(p_adj),
    log2_OR = log2(odds_ratio),
    gene_ratio = overlap / target_size,
    disease_gene_ratio = overlap / disease_set_size
  ) %>%
  arrange(p_adj, desc(odds_ratio))

enrich_res

p_disease_dot <- enrich_res %>%
  mutate(
    disease = str_to_title(disease),
    disease = fct_reorder(disease, log2_OR),
    sig = case_when(
      p_adj < 0.05 ~ "FDR < 0.05",
      p_value < 0.05 ~ "Nominal P < 0.05",
      TRUE ~ "Not significant"
    )
  ) %>%
  ggplot(aes(x = log2_OR, y = disease)) +
  geom_vline(xintercept = 0.585, linetype = "dashed", linewidth = 0.4) +
  geom_point(
    aes(size = overlap, fill = minus_log10_padj),
    shape = 21,
    color = "black",
    alpha = 0.9
  ) +
  scale_fill_gradient(low = "white", high = "firebrick")+
  theme_classic(base_size = 13) +
  labs(
    x = "log2 odds ratio",
    y = NULL,
    size = "Overlap genes",
    fill = "-log10(FDR)",
    title = "Disease gene enrichment in early OPC/OL trajectory genes"
  ) +
  theme(
    aspect.ratio = 3
  )

p_disease_dot
