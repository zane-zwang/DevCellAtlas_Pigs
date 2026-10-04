# Purpose: counts to tmm logcpm.

convert_counts_to_tmm_logCPM <- function(x,
                                         gene_col = "GeneID",
                                         group = NULL,
                                         min_group_size = NULL,
                                         prior.count = 1,
                                         keep_rownames = TRUE,
                                         verbose = TRUE) {
  convert_one <- function(df, name = NULL) {
    if (!is.data.frame(df)) {
      stop("Each element must be a data.frame.")
    }
    
    if (!gene_col %in% colnames(df)) {
      stop("Gene column '", gene_col, "' not found.")
    }
    
    genes <- as.character(df[[gene_col]])
    
    if (any(is.na(genes) | genes == "")) {
      stop("Missing or empty GeneID detected.")
    }
    
    if (anyDuplicated(genes)) {
      stop("Duplicated GeneID detected.")
    }
    
    expr_cols <- setdiff(colnames(df), gene_col)
    expr <- as.matrix(df[, expr_cols, drop = FALSE])
    suppressWarnings(storage.mode(expr) <- "numeric")
    
    if (anyNA(expr)) {
      stop("NA generated or detected in expression matrix.")
    }
    
    if (any(expr < 0)) {
      stop("Negative values detected. Raw counts should be non-negative.")
    }
    
    if (is.null(group)) {
      group_use <- rep("all", ncol(expr))
    } else {
      if (length(group) != ncol(expr)) {
        stop("Length of group must equal number of expression columns.")
      }
      group_use <- group
    }
    
    dge <- edgeR::DGEList(counts = expr)
    
    keep <- edgeR::filterByExpr(dge, group = group_use)
    
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    
    logCPM <- edgeR::cpm(
      dge,
      normalized.lib.sizes = TRUE,
      log = TRUE,
      prior.count = prior.count
    )
    
    out <- data.frame(
      GeneID = rownames(logCPM),
      logCPM,
      check.names = FALSE
    )
    
    colnames(out)[1] <- gene_col
    
    if (keep_rownames) {
      rownames(out) <- out[[gene_col]]
    }
    
    if (verbose) {
      msg_name <- ifelse(is.null(name), "input", name)
      message("[", msg_name, "] edgeR TMM-logCPM; ",
              nrow(out), " genes × ", length(expr_cols), " samples.")
    }
    
    return(out)
  }
  
  if (is.list(x) && !is.data.frame(x)) {
    if (is.null(names(x))) {
      names(x) <- paste0("element_", seq_along(x))
    }
    
    out_list <- lapply(seq_along(x), function(i) {
      convert_one(x[[i]], name = names(x)[i])
    })
    
    names(out_list) <- names(x)
    return(out_list)
  } else {
    return(convert_one(x))
  }
}
