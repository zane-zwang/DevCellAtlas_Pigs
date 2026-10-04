# Purpose: write gct.

write_matrix_as_gct <- function(expr_mat, outpath) {

  num_genes <- nrow(expr_mat)
  num_samples <- ncol(expr_mat)

  con <- file(outpath, open = "wt")
  writeLines("#1.2", con)  # GCT v1.2
  writeLines(paste(num_genes, num_samples, sep = "\t"), con)
  close(con)

  gct_df <- data.frame(
    Name = rownames(expr_mat),
    Description = rep("na", num_genes),
    expr_mat,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  write.table(
    gct_df,
    file = outpath,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    append = TRUE
  )
}
