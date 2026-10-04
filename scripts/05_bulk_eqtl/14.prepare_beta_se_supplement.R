# Purpose: prepare beta se supplement.
analysis_dir <- Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input")

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# Build beta / SE / Z matrices for mashr from extracted pairs
# Input columns:
#   pair_id pval_g1 beta_g1 beta_se_g1
# pair_id format:
#   variant_id,gene_id
# Example:
#   1_65265_T_C,ENSSSCG00000049216

mashr_dir <- file.path(analysis_dir, "qtl/mashr")
in_dir <- file.path(mashr_dir, "01_beqtl_mashr_input")

bins <- c("bin1", "bin2", "bin3", "bin4")

beta_list  <- vector("list", length(bins))
se_list    <- vector("list", length(bins))
z_list     <- vector("list", length(bins))
annot_list <- vector("list", length(bins))

names(beta_list)  <- bins
names(se_list)    <- bins
names(z_list)     <- bins
names(annot_list) <- bins

for (b in bins) {

  f <- file.path(
    in_dir,
    paste0(b, "_nominal_pairs.extracted_pairs.txt.gz")
  )

  if (!file.exists(f)) {
    stop("Missing input file: ", f)
  }

  message("Reading: ", f)

  dt <- fread(f)

  required_cols <- c(
    "pair_id",
    "pval_g1",
    "beta_g1",
    "beta_se_g1"
  )

  missing_cols <- setdiff(required_cols, colnames(dt))

  if (length(missing_cols) > 0) {
    stop(
      "Missing columns in ", f, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # Parse pair_id:
  #   variant_id,gene_id

  pair_split <- tstrsplit(
    dt$pair_id,
    ",",
    fixed = TRUE,
    keep = 1:2
  )

  if (length(pair_split) != 2) {
    stop("Failed to split pair_id into variant_id and gene_id in: ", f)
  }

  dt[, variant_id := pair_split[[1]]]
  dt[, gene_id := pair_split[[2]]]

  if (anyNA(dt$variant_id) || anyNA(dt$gene_id)) {
    stop("Some pair_id values could not be parsed in: ", f)
  }

  if (any(dt$variant_id == "") || any(dt$gene_id == "")) {
    stop("Empty variant_id or gene_id detected in: ", f)
  }

  # Use a separator absent from the input IDs.
  dt[, pair_uid := paste(variant_id, gene_id, sep = "|")]

  # Validate numeric estimates

  dt[, beta_g1 := as.numeric(beta_g1)]
  dt[, beta_se_g1 := as.numeric(beta_se_g1)]
  dt[, pval_g1 := as.numeric(pval_g1)]

  invalid <- !is.finite(dt$beta_g1) |
    !is.finite(dt$beta_se_g1) |
    dt$beta_se_g1 <= 0

  if (any(invalid)) {
    message(
      b, ": replacing ",
      sum(invalid),
      " invalid beta/SE entries with NA"
    )

    dt[
      invalid,
      `:=`(
        beta_g1 = NA_real_,
        beta_se_g1 = NA_real_
      )
    ]
  }

  dt[, zval := beta_g1 / beta_se_g1]

  # Check duplicate SNP-gene pairs within each bin

  duplicate_pairs <- dt[, .N, by = pair_uid][N > 1]

  if (nrow(duplicate_pairs) > 0) {
    stop(
      "Duplicated SNP-gene pairs detected in ", b,
      ". Examples: ",
      paste(head(duplicate_pairs$pair_uid, 5), collapse = ", ")
    )
  }

  # Store one value table per bin

  beta_dt <- dt[, .(pair_uid, value = beta_g1)]
  se_dt   <- dt[, .(pair_uid, value = beta_se_g1)]
  z_dt    <- dt[, .(pair_uid, value = zval)]

  setnames(beta_dt, "value", b)
  setnames(se_dt, "value", b)
  setnames(z_dt, "value", b)

  beta_list[[b]] <- beta_dt
  se_list[[b]]   <- se_dt
  z_list[[b]]    <- z_dt

  annot_list[[b]] <- dt[, .(
    pair_uid,
    pair_id_original = pair_id,
    variant_id,
    gene_id
  )]
}

# Merge bins

merge_by_pair <- function(x, y) {
  merge(x, y, by = "pair_uid", all = TRUE, sort = FALSE)
}

beta_df <- Reduce(merge_by_pair, beta_list)
se_df   <- Reduce(merge_by_pair, se_list)
z_df    <- Reduce(merge_by_pair, z_list)

pair_annot <- unique(
  rbindlist(annot_list, use.names = TRUE, fill = TRUE)
)

# Confirm one-to-one mapping
annot_check <- pair_annot[, .(
  n_original = uniqueN(pair_id_original),
  n_variant = uniqueN(variant_id),
  n_gene = uniqueN(gene_id)
), by = pair_uid]

bad_annot <- annot_check[
  n_original != 1 |
    n_variant != 1 |
    n_gene != 1
]

if (nrow(bad_annot) > 0) {
  stop("Some pair_uid values map to inconsistent annotations.")
}

pair_annot <- pair_annot[!duplicated(pair_uid)]

# Force identical row order in all matrices

all_pair_ids <- Reduce(
  union,
  list(
    beta_df$pair_uid,
    se_df$pair_uid,
    z_df$pair_uid
  )
)

reorder_table <- function(df, ids, bins) {

  template <- data.table(
    pair_uid = ids,
    row_order = seq_along(ids)
  )

  out <- merge(
    template,
    df,
    by = "pair_uid",
    all.x = TRUE,
    sort = FALSE
  )

  setorder(out, row_order)

  out[, row_order := NULL]

  out[, c("pair_uid", bins), with = FALSE]
}

beta_df <- reorder_table(beta_df, all_pair_ids, bins)
se_df   <- reorder_table(se_df, all_pair_ids, bins)
z_df    <- reorder_table(z_df, all_pair_ids, bins)

# Convert to matrices

beta_mat <- as.matrix(beta_df[, ..bins])
se_mat   <- as.matrix(se_df[, ..bins])
z_mat    <- as.matrix(z_df[, ..bins])

rownames(beta_mat) <- beta_df$pair_uid
rownames(se_mat)   <- se_df$pair_uid
rownames(z_mat)    <- z_df$pair_uid

storage.mode(beta_mat) <- "double"
storage.mode(se_mat)   <- "double"
storage.mode(z_mat)    <- "double"

# Check alignment of effect, standard-error and z-score matrices.

stopifnot(
  is.matrix(beta_mat),
  is.matrix(se_mat),
  is.matrix(z_mat),

  identical(dim(beta_mat), dim(se_mat)),
  identical(dim(beta_mat), dim(z_mat)),

  identical(rownames(beta_mat), rownames(se_mat)),
  identical(rownames(beta_mat), rownames(z_mat)),

  identical(colnames(beta_mat), bins),
  identical(colnames(se_mat), bins),
  identical(colnames(z_mat), bins)
)

if (any(se_mat[!is.na(se_mat)] <= 0)) {
  stop("Non-positive standard errors remain in se_mat.")
}

na_mismatch <- xor(is.na(beta_mat), is.na(se_mat))

if (any(na_mismatch)) {
  warning(
    sum(na_mismatch),
    " entries have beta/SE missingness mismatch."
  )
}

# Reorder annotation to matrix rows
pair_annot <- pair_annot[
  match(rownames(beta_mat), pair_uid)
]

stopifnot(
  identical(pair_annot$pair_uid, rownames(beta_mat))
)

# Save outputs

saveRDS(
  beta_mat,
  file.path(in_dir, "strong_set_beta.rds")
)

saveRDS(
  se_mat,
  file.path(in_dir, "strong_set_se.rds")
)

saveRDS(
  z_mat,
  file.path(in_dir, "strong_set_zval.rds")
)

fwrite(
  data.table(
    pair_uid = rownames(beta_mat),
    beta_mat
  ),
  file.path(in_dir, "strong_set_beta.tsv.gz"),
  sep = "\t"
)

fwrite(
  data.table(
    pair_uid = rownames(se_mat),
    se_mat
  ),
  file.path(in_dir, "strong_set_se.tsv.gz"),
  sep = "\t"
)

fwrite(
  data.table(
    pair_uid = rownames(z_mat),
    z_mat
  ),
  file.path(in_dir, "strong_set_zval.tsv.gz"),
  sep = "\t"
)

fwrite(
  pair_annot,
  file.path(in_dir, "pair_annot.tsv"),
  sep = "\t"
)

message("Done.")
message("Number of SNP-gene pairs: ", nrow(beta_mat))
message("Number of bins: ", ncol(beta_mat))
message("Output directory: ", in_dir)
