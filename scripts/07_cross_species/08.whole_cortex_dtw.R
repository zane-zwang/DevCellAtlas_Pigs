# Purpose: whole cortex dtw.
publication_root <- normalizePath(Sys.getenv("DEVPIGATLAS_ROOT", "."), mustWork = TRUE)

library(dtw)
library(ggplot2)
library(dplyr)
library(tidyr)
library(edgeR)
library(limma)
library(MatrixGenerics)

source(file.path(publication_root, "scripts/07_cross_species/07.dtw_helpers.R"))

# expression processing

edger_exp_df_log <- function(df, group = NULL) {
  df <- as.matrix(df)
  dge <- edgeR::DGEList(counts = df)
  
  if (is.null(group)) {
    group <- colnames(df)
  }
  
  keep <- edgeR::filterByExpr(dge, group = group)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  out <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  return(out)
}

make_pb_meta <- function(colnames_vec) {
  tibble(gr = colnames_vec) |>
    tidyr::separate(gr, c("species", "stage"), sep = "//", remove = FALSE)
}

run_split <- function(df, pb_meta) {
  expr_ss <- df[, pb_meta$species == "pig", drop = FALSE]
  colnames(expr_ss) <- sub("^pig//", "", colnames(expr_ss))
  
  expr_hs <- df[, pb_meta$species == "human", drop = FALSE]
  colnames(expr_hs) <- sub("^human//", "", colnames(expr_hs))
  
  expr_mm <- df[, pb_meta$species == "mouse", drop = FALSE]
  colnames(expr_mm) <- sub("^mouse//", "", colnames(expr_mm))
  
  exp_list <- list(
    human = expr_hs,
    mouse = expr_mm,
    pig = expr_ss
  )
  return(exp_list)
}


# load data

df_pb <- read.csv("pseudobulk/pseudobulk_species_sample_stage.csv", row.names = 1, check.names = FALSE)

hs_order <- c('WPC18','WPC19','WPC23','WPC24','P0','P4y','P6y','P14y','P20y','P39y')

mm_order <- c('E18.5','P4','P14','P32','P90')

ss_order <- c('E55','E90','P0','P30','P90','P180')

stage_order_list <- list(
  human = hs_order,
  mouse = mm_order,
  pig = ss_order
)

# aggregate smaples

df_pb_stage <- aggregate_stage_pseudobulk(
  df = df_pb,
  sep = "//",
  agg_fun = "sum"
)

# normalize

pb_df.logCPM <- edger_exp_df_log(df_pb_stage)
pb_meta <- make_pb_meta(colnames(pb_df.logCPM))

# Main analysis retains species effects.
expr_list <- run_split(pb_df.logCPM, pb_meta = pb_meta)


# optional anchor genes

# No anchor genes were supplied.
anchor_genes <- NULL

# run DTW

res <- run_dtw_all_pairs_to_list(
  pseudobulk_list = expr_list,
  stage_order_list = stage_order_list,
  top_n_hvgs = 2000,
  anchor_genes = anchor_genes,
  similarity_method = "cosine_distance",
  dynamic_method = "spearman",
  gene_combine = "freq",
  min_species = 2
)

dir.create("dtw_results", showWarnings = FALSE, recursive = TRUE)

saveRDS(res, file = "dtw_results/dtw_all_pairs_results.rds")
write.csv(res$scores, file = "dtw_results/dtw_scores.csv", row.names = FALSE)
write.csv(
  dplyr::bind_rows(res$alignments, .id = "pair_name"),
  file = "dtw_results/dtw_alignments.csv",
  row.names = FALSE
)
write.table(
  res$feature_genes,
  file = "dtw_results/feature_genes.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# bootstrap summary
if (length(res$bootstrap) > 0) {
  boot_sum <- bind_rows(
    lapply(names(res$bootstrap), function(nm) {
      res$bootstrap[[nm]]$summary |> mutate(pair = nm)
    })
  )
  write.csv(
    boot_sum,
    file = "dtw_results/dtw_bootstrap_summary.csv",
    row.names = FALSE
  )
}

# save plots

pdf("dtw_results/dtw_heatmaps.pdf", width = 6, height = 6)
for (nm in names(res$plots)) {
  print(res$plots[[nm]])
}
dev.off()
