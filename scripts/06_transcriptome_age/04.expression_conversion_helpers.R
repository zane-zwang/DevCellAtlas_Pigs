# Purpose: expression conversion helpers.

# Convert raw pseudobulk counts to log2(CPM + 1) or log2(TPM + 1)
# Input:
#   - a data.frame with GeneID + raw count columns
#   - or a named list of such data.frames
# Output:
#   - same structure as input
#   - GeneID retained
#   - expression columns transformed

convert_pseudobulk_expr <- function(x,
                                    method = c("CPM", "TPM"),
                                    gene_col = "GeneID",
                                    gene_length = NULL,
                                    length_gene_col = "GeneID",
                                    length_col = "length",
                                    log_base = 2,
                                    pseudocount = 1,
                                    keep_rownames = TRUE,
                                    verbose = TRUE) {
  method <- match.arg(method)

  convert_one <- function(df, name = NULL) {
    if (!is.data.frame(df)) {
      stop("Each element must be a data.frame.")
    }

    if (!gene_col %in% colnames(df)) {
      stop("Gene column '", gene_col, "' not found",
           if (!is.null(name)) paste0(" in element: ", name) else "", ".")
    }

    genes <- as.character(df[[gene_col]])

    if (any(is.na(genes) | genes == "")) {
      stop("Missing or empty GeneID detected",
           if (!is.null(name)) paste0(" in element: ", name) else "", ".")
    }

    if (anyDuplicated(genes)) {
      dup_genes <- unique(genes[duplicated(genes)])
      stop("Duplicated GeneID detected",
           if (!is.null(name)) paste0(" in element: ", name) else "",
           ". Example duplicated genes: ",
           paste(head(dup_genes, 10), collapse = ", "))
    }

    expr_cols <- setdiff(colnames(df), gene_col)

    if (length(expr_cols) == 0) {
      stop("No expression columns found",
           if (!is.null(name)) paste0(" in element: ", name) else "", ".")
    }

    expr <- as.matrix(df[, expr_cols, drop = FALSE])

    suppressWarnings(storage.mode(expr) <- "numeric")

    if (anyNA(expr)) {
      stop("NA generated or detected in expression matrix",
           if (!is.null(name)) paste0(" in element: ", name) else "",
           ". Check whether all expression columns are numeric raw counts.")
    }

    if (any(expr < 0)) {
      stop("Negative values detected",
           if (!is.null(name)) paste0(" in element: ", name) else "",
           ". Raw counts should be non-negative.")
    }

    col_sums <- colSums(expr)

    if (any(col_sums <= 0)) {
      bad_cols <- names(col_sums)[col_sums <= 0]
      stop("Some samples have total counts <= 0",
           if (!is.null(name)) paste0(" in element: ", name) else "",
           ": ", paste(bad_cols, collapse = ", "))
    }

    if (method == "CPM") {
      norm_expr <- sweep(expr, 2, col_sums, "/") * 1e6
    }

    if (method == "TPM") {
      if (is.null(gene_length)) {
        stop("gene_length must be provided when method = 'TPM'.")
      }

      if (!is.data.frame(gene_length)) {
        stop("gene_length must be a data.frame with gene ID and gene length columns.")
      }

      if (!all(c(length_gene_col, length_col) %in% colnames(gene_length))) {
        stop("gene_length must contain columns: ",
             length_gene_col, " and ", length_col, ".")
      }

      gl <- gene_length[, c(length_gene_col, length_col)]
      colnames(gl) <- c("GeneID_tmp", "length_tmp")
      gl$GeneID_tmp <- as.character(gl$GeneID_tmp)

      if (anyDuplicated(gl$GeneID_tmp)) {
        dup_len_genes <- unique(gl$GeneID_tmp[duplicated(gl$GeneID_tmp)])
        stop("Duplicated genes detected in gene_length. Example: ",
             paste(head(dup_len_genes, 10), collapse = ", "))
      }

      idx <- match(genes, gl$GeneID_tmp)
      if (anyNA(idx)) {
        missing_genes <- genes[is.na(idx)]
        stop("Gene length missing for some genes",
             if (!is.null(name)) paste0(" in element: ", name) else "",
             ". Example missing genes: ",
             paste(head(missing_genes, 10), collapse = ", "))
      }

      lengths <- as.numeric(gl$length_tmp[idx])

      if (anyNA(lengths) || any(lengths <= 0)) {
        stop("Invalid gene lengths detected",
             if (!is.null(name)) paste0(" in element: ", name) else "",
             ". Gene lengths must be positive numeric values.")
      }

      # If length is in bp, convert to kb.
      # If your gene length is already in kb, set length_col accordingly before input
      # or modify this line.
      lengths_kb <- lengths / 1000

      rpk <- sweep(expr, 1, lengths_kb, "/")
      scaling_factor <- colSums(rpk)

      if (any(scaling_factor <= 0)) {
        bad_cols <- names(scaling_factor)[scaling_factor <= 0]
        stop("Some samples have TPM scaling factor <= 0: ",
             paste(bad_cols, collapse = ", "))
      }

      norm_expr <- sweep(rpk, 2, scaling_factor, "/") * 1e6
    }

    if (log_base == 2) {
      log_expr <- log2(norm_expr + pseudocount)
    } else if (log_base == exp(1)) {
      log_expr <- log(norm_expr + pseudocount)
    } else if (log_base == 10) {
      log_expr <- log10(norm_expr + pseudocount)
    } else {
      log_expr <- log(norm_expr + pseudocount, base = log_base)
    }

    out <- data.frame(
      GeneID = genes,
      log_expr,
      check.names = FALSE
    )

    colnames(out)[1] <- gene_col

    if (keep_rownames) {
      rownames(out) <- genes
    }

    if (verbose) {
      msg_name <- ifelse(is.null(name), "input", name)
      message("[", msg_name, "] ",
              method, " -> log", log_base, "(",
              method, " + ", pseudocount, "); ",
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
