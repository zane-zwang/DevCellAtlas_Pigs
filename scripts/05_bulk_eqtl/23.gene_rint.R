# Purpose: gene rint.

#' Convert log2(CPM + 1) expression matrix to eQTL-ready form
#'
#' This function applies per-gene rank-based inverse normal transformation
#' to a log2(CPM + 1) gene expression matrix.
#'
#' @param logCPM_mat A gene × sample matrix in log2(CPM + 1) scale
#' @param return_df Logical. Whether to return as a data.frame (default: TRUE)
#'
#' @return An expression matrix suitable for eQTL mapping
#'
logCPM_to_eqtl_ready <- function(logCPM_mat, return_df = TRUE) {
  stopifnot(is.matrix(logCPM_mat) || is.data.frame(logCPM_mat))
  
  # Ensure it's a numeric matrix
  logCPM_mat <- as.matrix(logCPM_mat)
  mode(logCPM_mat) <- "numeric"

  # Inverse normal transformation per gene
  rank_qnorm <- function(x) {
    qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
  }
  
  expr_inv <- t(apply(logCPM_mat, 1, rank_qnorm))
  
  # Return as data.frame if preferred
  if (return_df) {
    expr_inv <- as.data.frame(expr_inv)
    rownames(expr_inv) <- rownames(logCPM_mat)
    colnames(expr_inv) <- colnames(logCPM_mat)
  }
  
  return(expr_inv)
}
