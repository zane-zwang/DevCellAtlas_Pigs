# Purpose: Supplementary celltype expression heatmap.

# Cross-species cell-type correlation heatmap from pseudobulk expression
# Input:
#   - pseudobulk matrix generated from single-cell/snRNA-seq data
#   - rows are ortholog-mapped genes
#   - columns are pseudobulk samples named as:
#       species//sample_id//stage//celltype
# Main purpose:
#   Quantify cross-species cell-type similarity using pseudobulk expression and
#   visualize the cell-type correlation structure with a publication-quality
#   ComplexHeatmap figure.
# Key methodological choices:
#   1. Use mean-based pseudobulk expression as input, not sum-based counts.
#   2. Filter weakly detected genes before feature selection.
#   3. Select informative genes using pseudobulk-level variance.
#   4. Apply log1p transformation.
#   5. Apply sample-level centering and gene-wise z-score normalization.
#   6. Optionally apply species-level centering to remove global species shift.
#   7. Aggregate sample/stage-level pseudobulks to species-by-celltype profiles.
#   8. Compute Spearman or Pearson correlation between species-celltype profiles.
# Author note:
#   This script assumes the input genes have already been mapped to a shared
#   ortholog space. If multiple species-specific gene IDs remain, harmonize them
#   before running this script.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(matrixStats)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(RColorBrewer)
})

# 0. User-defined parameters

input_file <- "pseudobulk/pseudobulk_species_sample_stage_celltype_mean.csv"
outdir <- "celltype_correlation_heatmap"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Column name separator in pseudobulk matrix.
sep <- "//"

# Minimum expression threshold after mean-pseudobulk construction.
# For mean expression, do NOT use min_expr = 1 unless values are known to be high.
# A small threshold is safer. Set to 0 to keep any gene expressed above zero.
min_expr <- 0

# Gene must be expressed above min_expr in at least this fraction of pseudobulk samples.
# This removes genes dominated by zeros and species-specific dropout.
min_sample_fraction <- 0.10

# Number of pseudobulk-level highly variable genes to use.
# For cross-species cell-type correlation, 3000-5000 is often more stable than 2000.
n_hvg <- 3000

# Correlation method. Spearman is usually more robust across species.
cor_method <- "spearman"   # options: "spearman", "pearson"

# Whether to remove species-level global expression shifts before aggregation.
# This is often the key step when samples cluster by species instead of cell type.
do_species_centering <- TRUE

# Whether to remove stage-level shifts. Use with caution.
# Recommended only when matched developmental stages are not the main signal of interest.
do_stage_centering <- FALSE

# Whether to aggregate sample/stage pseudobulks to species x celltype profiles.
# For a clean cell-type conservation heatmap, TRUE is recommended.
do_celltype_aggregation <- TRUE

# Optional fixed cell type order. Set NULL to infer from data.
celltype_order <- NULL

# Optional fixed species order. Edit according to your study design.
species_order <- c("human", "mouse", "pig")

# Output file names.
output_prefix <- file.path(outdir, "cross_species_celltype_correlation")

# 1. Utility functions

clean_matrix <- function(mat) {
  # Convert data frame to numeric matrix and remove invalid values.
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- NA_real_
  return(mat)
}

safe_zscore_rows <- function(mat) {
  # Gene-wise z-score. Genes with zero variance are removed.
  # Input: gene x sample matrix.
  sds <- matrixStats::rowSds(mat, na.rm = TRUE)
  keep <- is.finite(sds) & sds > 0
  mat <- mat[keep, , drop = FALSE]

  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- NA_real_
  return(z)
}

center_by_group <- function(mat, group) {
  # Subtract group-specific gene mean from each sample.
  # This removes broad group-specific baseline shifts while preserving relative
  # cell-type or stage variation within each group.
  stopifnot(ncol(mat) == length(group))

  out <- mat
  group <- as.character(group)

  for (g in unique(group)) {
    idx <- which(group == g)
    if (length(idx) <= 1) next
    out[, idx] <- sweep(out[, idx, drop = FALSE], 1,
                        rowMeans(out[, idx, drop = FALSE], na.rm = TRUE),
                        FUN = "-")
  }

  return(out)
}

make_named_palette <- function(values, palette = "Set2") {
  # Generate a named discrete color vector for annotation.
  values <- unique(as.character(values))
  n <- length(values)

  if (n <= 8) {
    cols <- RColorBrewer::brewer.pal(max(3, n), palette)[seq_len(n)]
  } else {
    cols <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, palette))(n)
  }

  stats::setNames(cols, values)
}

parse_pseudobulk_metadata <- function(sample_names, sep = "//") {
  # Parse pseudobulk sample names into metadata.
  meta <- data.frame(sample = sample_names, stringsAsFactors = FALSE)
  meta <- tidyr::separate(
    meta,
    col = "sample",
    into = c("species", "sample_id", "stage", "celltype"),
    sep = sep,
    remove = FALSE,
    extra = "merge",
    fill = "right"
  )

  if (any(is.na(meta$species)) || any(is.na(meta$sample_id)) ||
      any(is.na(meta$stage)) || any(is.na(meta$celltype))) {
    stop(
      "Some column names could not be parsed into species/sample_id/stage/celltype.\n",
      "Expected format: species//sample_id//stage//celltype"
    )
  }

  rownames(meta) <- meta$sample
  return(meta)
}

# 2. Read input matrix and metadata

message("[1/8] Reading pseudobulk matrix: ", input_file)
expr_raw <- read.csv(input_file, row.names = 1, check.names = FALSE)
expr_raw <- clean_matrix(expr_raw)

if (any(duplicated(rownames(expr_raw)))) {
  message("[Info] Duplicated gene names detected. Collapsing duplicated genes by mean expression.")
  rs <- rowsum(expr_raw, group = rownames(expr_raw), reorder = FALSE)
  gene_n <- as.numeric(table(factor(rownames(expr_raw), levels = rownames(rs))))
  expr_raw <- sweep(rs, 1, gene_n, FUN = "/")
}

meta <- parse_pseudobulk_metadata(colnames(expr_raw), sep = sep)

# Keep only samples with complete metadata and expression columns.
meta <- meta[colnames(expr_raw), , drop = FALSE]
stopifnot(all(rownames(meta) == colnames(expr_raw)))

# Apply user-defined factor order where possible.
if (!is.null(species_order)) {
  species_order <- intersect(species_order, unique(meta$species))
  meta$species <- factor(meta$species,
                         levels = c(species_order,
                                    setdiff(unique(meta$species), species_order)))
} else {
  meta$species <- factor(meta$species)
}

if (!is.null(celltype_order)) {
  celltype_order <- intersect(celltype_order, unique(meta$celltype))
  meta$celltype <- factor(meta$celltype,
                          levels = c(celltype_order,
                                     setdiff(unique(meta$celltype), celltype_order)))
} else {
  meta$celltype <- factor(meta$celltype)
}

message("[Info] Input matrix: ", nrow(expr_raw), " genes x ", ncol(expr_raw), " pseudobulk samples")
message("[Info] Species: ", paste(levels(meta$species), collapse = ", "))
message("[Info] Cell types: ", paste(levels(meta$celltype), collapse = ", "))

# 3. Gene filtering and log transformation

message("[2/8] Filtering weakly detected genes")

# For mean pseudobulk values, expression > 0 is usually a reasonable detection call.
keep_detected <- rowMeans(expr_raw > min_expr, na.rm = TRUE) >= min_sample_fraction
expr_filt <- expr_raw[keep_detected, , drop = FALSE]

message("[Info] Genes retained after detection filter: ",
        nrow(expr_filt), " / ", nrow(expr_raw))

if (nrow(expr_filt) < 500) {
  warning("Fewer than 500 genes retained. Consider lowering min_sample_fraction or min_expr.")
}

message("[3/8] Applying log1p transformation")
expr_log <- log1p(expr_filt)

# 4. Pseudobulk-level HVG selection

message("[4/8] Selecting pseudobulk-level highly variable genes")

# Select HVGs after log transformation, before z-score.
rv <- matrixStats::rowVars(expr_log, na.rm = TRUE)
n_hvg_use <- min(n_hvg, sum(is.finite(rv) & rv > 0))

if (n_hvg_use < n_hvg) {
  message("[Info] Requested ", n_hvg, " HVGs, but only ", n_hvg_use,
          " genes have non-zero finite variance.")
}

hvg <- names(sort(rv, decreasing = TRUE))[seq_len(n_hvg_use)]
expr_hvg <- expr_log[hvg, , drop = FALSE]

write.table(
  data.frame(gene = hvg, variance = rv[hvg]),
  file = paste0(output_prefix, "_selected_HVGs.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

message("[Info] HVGs used for correlation: ", length(hvg))

# 5. Normalization and bias reduction

message("[5/8] Normalizing expression matrix")

# 5.1 Sample-level centering:
# Remove sample-wide mean expression differences. This helps reduce residual
# library size or global expression shifts after pseudobulk generation.
expr_norm <- sweep(expr_hvg, 2, colMeans(expr_hvg, na.rm = TRUE), FUN = "-")

# 5.2 Gene-wise z-score:
# Put genes on a comparable scale so correlation is driven by relative patterns.
expr_norm <- safe_zscore_rows(expr_norm)

# 5.3 Optional species-level centering:
# This is especially useful if the correlation heatmap clusters by species.
if (isTRUE(do_species_centering)) {
  message("[Info] Applying species-level gene-wise centering")
  expr_norm <- center_by_group(expr_norm, meta[colnames(expr_norm), "species"])
}

# 5.4 Optional stage-level centering:
# Use this only if stage-specific signal is a nuisance for the heatmap.
if (isTRUE(do_stage_centering)) {
  message("[Info] Applying stage-level gene-wise centering")
  expr_norm <- center_by_group(expr_norm, meta[colnames(expr_norm), "stage"])
}

# Remove genes that became invalid after normalization or centering.
valid_gene <- rowSums(is.finite(expr_norm)) == ncol(expr_norm)
expr_norm <- expr_norm[valid_gene, , drop = FALSE]

message("[Info] Matrix after normalization: ",
        nrow(expr_norm), " genes x ", ncol(expr_norm), " samples")

# 6. Aggregate to species x celltype profiles

message("[6/8] Building correlation input matrix")

if (isTRUE(do_celltype_aggregation)) {
  # Convert to long format for transparent aggregation.
  expr_df <- as.data.frame(expr_norm)
  expr_df$gene <- rownames(expr_df)

  expr_long <- tidyr::pivot_longer(
    expr_df,
    cols = -gene,
    names_to = "sample",
    values_to = "expr"
  )

  expr_long <- expr_long %>%
    left_join(meta, by = "sample")

  # Average across samples and stages within each species-celltype combination.
  # If developmental stage is the target of analysis, do not aggregate across stage;
  # instead use species x stage x celltype profiles.
  expr_celltype <- expr_long %>%
    group_by(gene, species, celltype) %>%
    summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") %>%
    mutate(profile_id = paste(species, celltype, sep = "__"))

  expr_wide <- expr_celltype %>%
    select(gene, profile_id, expr) %>%
    tidyr::pivot_wider(names_from = profile_id, values_from = expr)

  expr_mat <- as.matrix(expr_wide[, -1, drop = FALSE])
  rownames(expr_mat) <- expr_wide$gene

  profile_meta <- data.frame(profile_id = colnames(expr_mat), stringsAsFactors = FALSE) %>%
    tidyr::separate(profile_id,
                    into = c("species", "celltype"),
                    sep = "__",
                    remove = FALSE)
  rownames(profile_meta) <- profile_meta$profile_id

} else {
  # Use each pseudobulk sample as one profile.
  expr_mat <- expr_norm
  profile_meta <- meta[colnames(expr_mat), , drop = FALSE]
  profile_meta$profile_id <- rownames(profile_meta)
}

# Remove profiles with any NA after aggregation.
valid_profile <- colSums(is.finite(expr_mat)) == nrow(expr_mat)
expr_mat <- expr_mat[, valid_profile, drop = FALSE]
profile_meta <- profile_meta[colnames(expr_mat), , drop = FALSE]

# Optional ordering before heatmap. Clustering can still override this visually,
# but the metadata order is useful when cluster_rows/columns are disabled.
profile_meta$species <- factor(profile_meta$species, levels = levels(meta$species))
profile_meta$celltype <- factor(profile_meta$celltype, levels = levels(meta$celltype))

ord <- order(profile_meta$celltype, profile_meta$species)
expr_mat <- expr_mat[, ord, drop = FALSE]
profile_meta <- profile_meta[colnames(expr_mat), , drop = FALSE]

write.csv(expr_mat, paste0(output_prefix, "_correlation_input_matrix.csv"))
write.csv(profile_meta, paste0(output_prefix, "_profile_metadata.csv"), row.names = FALSE)

message("[Info] Correlation input: ",
        nrow(expr_mat), " genes x ", ncol(expr_mat), " profiles")

# 7. Correlation analysis

message("[7/8] Computing ", cor_method, " correlation")

cor_mat <- cor(expr_mat, method = cor_method, use = "pairwise.complete.obs")

write.csv(cor_mat, paste0(output_prefix, "_", cor_method, "_correlation_matrix.csv"))

# 8. Publication-quality heatmap

message("[8/8] Drawing heatmap")

# Annotation colors.

# Keep color vectors named and aligned with actual annotation values.

species_cols = c("human"="#7DAEE0", "mouse"="#edab1c", "pig"="#725f97")
celltype_cols = c("astrocyte"="#FED43999", "endothelial_cell"="#709AE199", "excitatory_neuron"="#8A919799", 
                  "inhibitory_interneuron"="#D2AF8199", "microglia"="#FD744699", "mural"="#D5E4A299",
                  "oligodendrocyte"="#197EC099", "oligodendrocyte_progenitor_cell"="#F05C3B99", 
                  "choroid_plexus"="#46732E99", "fibroblast"="#71D0F599", 
                  "intermediate_progenitor_cell"="#37033599","mixed_neuron"="#07514999", "radial_glia"="#C8081399")

ha_top <- HeatmapAnnotation(
  species = profile_meta$species,
  celltype = profile_meta$celltype,
  col = list(
    species = species_cols,
    celltype = celltype_cols
  ),
  annotation_name_side = "left",
  annotation_legend_param = list(
    species = list(title = "Species"),
    celltype = list(title = "Cell type")
  )
)

ha_left <- rowAnnotation(
  species = profile_meta$species,
  celltype = profile_meta$celltype,
  col = list(
    species = species_cols,
    celltype = celltype_cols
  ),
  annotation_name_side = "top",
  annotation_legend_param = list(
    species = list(title = "Species"),
    celltype = list(title = "Cell type")
  )
)

# Correlation color scale. The lower bound is adaptive but capped for readability.
cor_min <- max(-1, floor(min(cor_mat, na.rm = TRUE) * 10) / 10)
cor_col_fun <- circlize::colorRamp2(
  c(cor_min, 0, 1),
  c("#2166AC", "#F7F7F7", "#B2182B")
)

ht <- Heatmap(
  cor_mat,
  name = paste0(tools::toTitleCase(cor_method), "\ncorrelation"),
  col = cor_col_fun,
  top_annotation = ha_top,
  left_annotation = ha_left,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_method_rows = "average",
  clustering_method_columns = "average",
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 8),
  column_names_gp = grid::gpar(fontsize = 8),
  column_names_rot = 45,
  heatmap_legend_param = list(
    title = paste0(tools::toTitleCase(cor_method), "\ncorrelation"),
    at = c(round(cor_min, 1), 0, 0.5, 1),
    labels = c(round(cor_min, 1), 0, 0.5, 1)
  )
)

pdf_file <- paste0(output_prefix, "_", cor_method, "_heatmap.pdf")

pdf(pdf_file, width = 12, height = 7, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

message("[Done] Output files:")
message("  - ", pdf_file)
message("  - ", png_file)
message("  - ", paste0(output_prefix, "_", cor_method, "_correlation_matrix.csv"))
message("  - ", paste0(output_prefix, "_correlation_input_matrix.csv"))
message("  - ", paste0(output_prefix, "_profile_metadata.csv"))
message("  - ", paste0(output_prefix, "_selected_HVGs.tsv"))
