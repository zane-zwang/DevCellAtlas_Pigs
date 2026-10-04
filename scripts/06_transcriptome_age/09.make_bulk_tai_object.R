# Purpose: make bulk tai object.

make_bulk_myTAI_object <- function(count_mat,
                                   pheno,
                                   gene_age,
                                   breed_name,
                                   gestation_days = 114,
                                   gene_id_col = "ENSEMBL",
                                   ps_col = "rank") {
  stopifnot(all(c("Sample", "age", "breed") %in% colnames(pheno)))
  stopifnot(gene_id_col %in% colnames(gene_age))
  stopifnot(ps_col %in% colnames(gene_age))
  
  common_samples <- intersect(colnames(count_mat), pheno$Sample)
  
  if (length(common_samples) == 0) {
    stop("No shared sample names between count matrix and pheno.")
  }
  
  count_mat <- count_mat[, common_samples, drop = FALSE]
  
  pheno <- pheno %>%
    filter(Sample %in% common_samples) %>%
    arrange(match(Sample, colnames(count_mat)))
  
  stopifnot(identical(pheno$Sample, colnames(count_mat)))
  
  pheno <- pheno %>%
    mutate(
      age = as.numeric(age),
      stage = ifelse(
        age < 0,
        paste0("E", gestation_days + age),
        paste0("D", age)
      )
    )
  
  stage_order <- pheno %>%
    distinct(age, stage) %>%
    arrange(age) %>%
    pull(stage)
  
  pheno$stage <- factor(pheno$stage, levels = stage_order)
  
  dge <- edgeR::DGEList(counts = count_mat)
  
  keep <- edgeR::filterByExpr(
    dge,
    group = interaction(pheno$breed, pheno$age, drop = TRUE)
  )
  
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  
  logCPM <- edgeR::cpm(
    dge,
    normalized.lib.sizes = TRUE,
    log = TRUE,
    prior.count = 1
  )
  
  pheno_b <- pheno %>%
    filter(breed == breed_name) %>%
    arrange(age, Sample)
  
  if (nrow(pheno_b) == 0) {
    stop("No samples found for breed: ", breed_name)
  }
  
  expr_b <- logCPM[, pheno_b$Sample, drop = FALSE]
  
  expr_stage <- sapply(levels(pheno_b$stage), function(st) {
    rowMeans(expr_b[, pheno_b$stage == st, drop = FALSE])
  })
  
  expr_stage <- as.data.frame(expr_stage)
  colnames(expr_stage) <- levels(pheno_b$stage)
  
  gene_age_map <- gene_age %>%
    transmute(
      GeneID = .data[[gene_id_col]],
      Phylostratum = .data[[ps_col]]
    ) %>%
    filter(!is.na(GeneID), !is.na(Phylostratum)) %>%
    distinct(GeneID, .keep_all = TRUE)
  
  phyex_df <- expr_stage %>%
    rownames_to_column("GeneID") %>%
    inner_join(gene_age_map, by = "GeneID") %>%
    select(GeneID, Phylostratum, everything()) %>%
    arrange(Phylostratum, GeneID)

  phyex_df.clean <- phyex_df %>%
    dplyr::filter(!is.na(Phylostratum)) %>%
    dplyr::select(Phylostratum, GeneID, everything())

  phyex_obj <- myTAI::BulkPhyloExpressionSet_from_df(phyex_df.clean)
  
  list(
    phyex_df = phyex_df,
    phyex_obj = phyex_obj,
    pheno = pheno,
    retained_genes = nrow(phyex_df.clean)
  )
}
