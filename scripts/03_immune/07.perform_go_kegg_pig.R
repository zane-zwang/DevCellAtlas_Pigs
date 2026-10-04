# Purpose: perform go kegg pig.

perform_enrichment_ss <- function(df, method = c("GO", "KEGG"), verbose = TRUE) {
  method <- match.arg(method)
  stopifnot(all(c("gene") %in% colnames(df)))
  
  if (nrow(df) == 0) {
    warning("Input data frame is empty.")
    return(NULL)
  }
  
  gen <- unique(as.character(df$gene))
  

  bitr_gene <- clusterProfiler::bitr(gen,
                                     fromType = "SYMBOL",
                                     toType   = "ENTREZID",
                                     OrgDb    = org.Ss.eg.db)
  
  if (verbose) {
    cat("Total input genes:", length(gen), "\n")
    cat("Mapped genes:", nrow(bitr_gene), "\n")
    cat("Mapping success rate:", round(nrow(bitr_gene) / length(gen) * 100, 2), "%\n")
  }
  
  if (nrow(bitr_gene) == 0) {
    warning("No valid gene IDs mapped to ENTREZID.")
    return(NULL)
  }
  
  enrich_result <- tryCatch({
    if (method == "GO") {
      clusterProfiler::enrichGO(
        gene          = unique(bitr_gene$ENTREZID),
        OrgDb         = org.Ss.eg.db,
        keyType       = "ENTREZID",
        ont           = "ALL",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.1,
        qvalueCutoff  = 1,
        readable      = TRUE
      )
    } else {
      clusterProfiler::enrichKEGG(
        gene          = unique(bitr_gene$ENTREZID),
        organism      = "ssc",
        pvalueCutoff  = 0.1
      )
    }
  }, error = function(e) {
    warning("Enrichment failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(enrich_result) || nrow(enrich_result@result) == 0) {
    warning("No enriched terms found (raw result is empty).")
    return(NULL)
  }
  
  if (verbose) cat("Terms before simplify:", nrow(enrich_result@result), "\n")
  
  if (method == "GO") {
    enrich_result <- tryCatch({
      enrich_result <- enrichplot::pairwise_termsim(enrich_result)
      clusterProfiler::simplify(enrich_result, cutoff = 0.7, by = "p.adjust", select_fun = min)
    }, error = function(e) {
      message("GO simplify skipped: ", e$message)
      enrich_result
    })
  }
  
  df_result <- as.data.frame(enrich_result@result)
  
  if (verbose) {
    cat("Terms before final filtering:", nrow(df_result), "\n")
    print(summary(df_result$p.adjust))
  }
  
  has_q <- "qvalue" %in% colnames(df_result) && !all(is.na(df_result$qvalue))
  if (has_q) {
    df_result <- dplyr::filter(df_result, p.adjust <= 0.1, qvalue <= 0.2)
  } else {
    df_result <- dplyr::filter(df_result, p.adjust <= 0.1)
  }
  
  if (nrow(df_result) == 0) {
    warning("Enriched terms existed before filtering, but none passed your cutoff.")
    return(df_result)
  }
  
  parse_ratio <- function(s) {
    sp <- strsplit(s, "/", fixed = TRUE)[[1]]
    as.numeric(sp[1]) / as.numeric(sp[2])
  }
  df_result$Enrichment_fold <- round(mapply(parse_ratio, df_result$GeneRatio) /
                                       mapply(parse_ratio, df_result$BgRatio), 2)
  
  if (verbose) cat("Enriched terms after filtering:", nrow(df_result), "\n\n")
  df_result
}
