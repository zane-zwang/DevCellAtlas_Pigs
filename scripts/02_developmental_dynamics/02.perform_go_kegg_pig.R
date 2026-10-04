# Purpose: perform go kegg pig.

# GO/KEGG enrichment function for pig gene sets.
# Required packages: clusterProfiler, enrichplot, dplyr, org.Ss.eg.db

perform_enrichment <- function(df,
                               method = c("GO", "KEGG"),
                               gene_col = "gene",
                               universe_genes = NULL,
                               OrgDb = org.Ss.eg.db,
                               fromType = "SYMBOL",
                               toType = "ENTREZID",
                               ont = "BP",
                               organism = "ssc",
                               pvalueCutoff = 0.1,
                               qvalueCutoff = 1,
                               pAdjustMethod = "BH",
                               minGSSize = 10,
                               maxGSSize = 500,
                               simplify_go = TRUE,
                               simplify_cutoff = 0.7,
                               final_padj_cutoff = 0.1,
                               final_qvalue_cutoff = NULL,
                               readable = TRUE,
                               verbose = TRUE) {
  method <- match.arg(method)

  if (!gene_col %in% colnames(df)) {
    stop("Column '", gene_col, "' was not found in df.")
  }

  if (nrow(df) == 0) {
    warning("Input data frame is empty.")
    return(data.frame())
  }

  genes <- unique(stats::na.omit(as.character(df[[gene_col]])))
  genes <- genes[genes != ""]

  if (length(genes) == 0) {
    warning("No valid genes after removing NA/empty values.")
    return(data.frame())
  }

  bitr_gene <- tryCatch({
    clusterProfiler::bitr(
      genes,
      fromType = fromType,
      toType   = toType,
      OrgDb    = OrgDb
    )
  }, error = function(e) {
    warning("Gene ID conversion failed: ", e$message)
    return(data.frame())
  })

  bitr_gene <- unique(bitr_gene[, c(fromType, toType), drop = FALSE])
  input_entrez <- unique(as.character(bitr_gene[[toType]]))
  input_entrez <- input_entrez[!is.na(input_entrez) & input_entrez != ""]

  universe_entrez <- NULL
  if (!is.null(universe_genes)) {
    universe_genes <- unique(stats::na.omit(as.character(universe_genes)))
    universe_genes <- universe_genes[universe_genes != ""]

    bitr_bg <- tryCatch({
      clusterProfiler::bitr(
        universe_genes,
        fromType = fromType,
        toType   = toType,
        OrgDb    = OrgDb
      )
    }, error = function(e) {
      warning("Universe gene ID conversion failed; enrichment will run without universe. ", e$message)
      return(data.frame())
    })

    if (nrow(bitr_bg) > 0) {
      bitr_bg <- unique(bitr_bg[, c(fromType, toType), drop = FALSE])
      universe_entrez <- unique(as.character(bitr_bg[[toType]]))
      universe_entrez <- universe_entrez[!is.na(universe_entrez) & universe_entrez != ""]

      # keep only genes that are part of the declared universe
      input_entrez <- intersect(input_entrez, universe_entrez)
    }
  }

  if (verbose) {
    cat("Method:", method, "\n")
    cat("Total input genes:", length(genes), "\n")
    cat("Mapped input genes:", length(input_entrez), "\n")
    cat("Mapping success rate:", round(length(unique(bitr_gene[[fromType]])) / length(genes) * 100, 2), "%\n")
    if (!is.null(universe_entrez)) {
      cat("Mapped universe genes:", length(universe_entrez), "\n")
    }
  }

  if (length(input_entrez) == 0) {
    warning("No valid input genes mapped to ", toType, ".")
    return(data.frame())
  }

  enrich_result <- tryCatch({
    if (method == "GO") {
      clusterProfiler::enrichGO(
        gene          = input_entrez,
        universe      = universe_entrez,
        OrgDb         = OrgDb,
        keyType       = toType,
        ont           = ont,
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff  = pvalueCutoff,
        qvalueCutoff  = qvalueCutoff,
        minGSSize     = minGSSize,
        maxGSSize     = maxGSSize,
        readable      = readable
      )
    } else {
      clusterProfiler::enrichKEGG(
        gene          = input_entrez,
        universe      = universe_entrez,
        organism      = organism,
        keyType       = "kegg",
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff  = pvalueCutoff,
        qvalueCutoff  = qvalueCutoff,
        minGSSize     = minGSSize,
        maxGSSize     = maxGSSize
      )
    }
  }, error = function(e) {
    warning("Enrichment failed: ", e$message)
    return(NULL)
  })

  if (is.null(enrich_result)) {
    return(data.frame())
  }

  raw_df <- as.data.frame(enrich_result)
  if (nrow(raw_df) == 0) {
    warning("No enriched terms found before final filtering.")
    return(raw_df)
  }

  if (verbose) {
    cat("Terms before simplify/filtering:", nrow(raw_df), "\n")
    cat("p.adjust summary before final filtering:\n")
    print(summary(raw_df$p.adjust))
  }

  if (method == "GO" && isTRUE(simplify_go) && nrow(raw_df) > 1) {
    enrich_result <- tryCatch({
      enrich_result2 <- enrichplot::pairwise_termsim(enrich_result)
      clusterProfiler::simplify(
        enrich_result2,
        cutoff = simplify_cutoff,
        by = "p.adjust",
        select_fun = min
      )
    }, error = function(e) {
      message("GO simplify skipped: ", e$message)
      enrich_result
    })
  }

  df_result <- as.data.frame(enrich_result)

  if (!is.null(final_padj_cutoff)) {
    df_result <- dplyr::filter(df_result, .data$p.adjust <= final_padj_cutoff)
  }

  if (!is.null(final_qvalue_cutoff) && "qvalue" %in% colnames(df_result)) {
    df_result <- dplyr::filter(df_result, is.na(.data$qvalue) | .data$qvalue <= final_qvalue_cutoff)
  }

  parse_ratio <- function(x) {
    vapply(strsplit(as.character(x), "/", fixed = TRUE), function(z) {
      if (length(z) != 2) return(NA_real_)
      as.numeric(z[1]) / as.numeric(z[2])
    }, numeric(1))
  }

  if (nrow(df_result) > 0 && all(c("GeneRatio", "BgRatio") %in% colnames(df_result))) {
    df_result$Enrichment_fold <- round(parse_ratio(df_result$GeneRatio) / parse_ratio(df_result$BgRatio), 2)
  }

  if (verbose) {
    cat("Terms after final filtering:", nrow(df_result), "\n\n")
  }

  df_result
}
