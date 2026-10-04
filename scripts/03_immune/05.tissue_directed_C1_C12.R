# Purpose: tissue directed C1 C12.
analysis_dir <- Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input")

library(Seurat)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(dynamicTreeCut)
library(cluster)

set.seed(123)

# Input / output
obj <- readRDS(file.path(analysis_dir, "immune/immune_anno_scVI.rds"))

outdir <- "result4_tissue_directed_programs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Step 0. Define immune lineage
lineage_map <- c(
  "CNS-associated macrophage"          = "Macrophage",
  "Cycling macrophage"                = "Macrophage",
  "IFN-response macrophage"           = "Macrophage",
  "Kupffer-like macrophage"           = "Macrophage",
  "Macrophage-like mixed"             = "Macrophage",
  "Microglia"                         = "Macrophage",
  "Monocyte-derived macrophage"       = "Macrophage",
  "Resident macrophage"               = "Macrophage",

  "CD4/activated T cell"              = "T_NK",
  "CD8 T cell"                        = "T_NK",
  "Cycling T cell"                    = "T_NK",
  "Duodenum-enriched activated T cell" = "T_NK",
  "Innate-like T cell"                = "T_NK",
  "NK cell"                           = "T_NK",

  "B cell"                            = "B_Plasma",
  "Plasma cell"                       = "B_Plasma",

  "pDC-like/pre-DC"                   = "DC",
  "cDC1-like DC"                      = "DC",
  "cDC2/activated DC"                 = "DC",

  "Mast cell"                         = "Mast"
)

obj$immune_lineage <- unname(lineage_map[as.character(obj$immune_subtype)])

unmapped_subtypes <- setdiff(unique(as.character(obj$immune_subtype)), names(lineage_map))
if (length(unmapped_subtypes) > 0) {
  stop("Unmapped immune_subtype detected: ", paste(unmapped_subtypes, collapse = ", "))
}
if (any(is.na(obj$immune_lineage))) {
  stop("NA values detected in obj$immune_lineage.")
}

# Main lineages for Result 4. DC is retained only if it passes cell number filters.
lineages.use <- c("Macrophage", "T_NK", "B_Plasma", "DC")

DefaultAssay(obj) <- "RNA"

# Step 1. Define tissue-restricted / tissue-shared subtypes
min_cells_subtype_tissue <- 100

subtype_tissue_stat <- obj@meta.data %>%
  count(immune_subtype, tissue, name = "n") %>%
  group_by(immune_subtype) %>%
  mutate(
    total = sum(n),
    frac = n / total,
    max_frac = max(frac),
    n_tissues_ge_min = sum(n >= min_cells_subtype_tissue)
  ) %>%
  ungroup()

subtype_class <- subtype_tissue_stat %>%
  group_by(immune_subtype) %>%
  summarise(
    max_frac = max(max_frac),
    n_tissues_ge_min = max(n_tissues_ge_min),
    .groups = "drop"
  ) %>%
  mutate(
    tissue_class = case_when(
      max_frac >= 0.70 | n_tissues_ge_min <= 2 ~ "tissue_restricted",
      n_tissues_ge_min >= 3 ~ "tissue_shared",
      TRUE ~ "intermediate"
    )
  )

write.table(
  subtype_tissue_stat,
  file.path(outdir, "subtype_tissue_distribution.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  subtype_class,
  file.path(outdir, "subtype_tissue_class.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

print(table(subtype_class$tissue_class))

# Step 2. Lineage-level tissue-associated markers
# Tissue-restricted cells are intentionally retained here.
# This step identifies lineage-level tissue programs, not strict replicate-supported DEGs.
run_lineage_tissue_markers <- function(obj,
                                       lineage,
                                       min_cells = 200,
                                       max_cells_per_group = 5000,
                                       logfc.threshold = 0.25,
                                       min.pct = 0.10,
                                       seed = 123) {
  message("Lineage: ", lineage)

  obj.lin <- subset(obj, subset = immune_lineage == lineage)
  obj.lin$tissue <- droplevels(factor(obj.lin$tissue))

  tissue_counts <- table(obj.lin$tissue)
  tissues.use <- names(tissue_counts)[tissue_counts >= min_cells]

  if (length(tissues.use) < 2) {
    message("  skipped: fewer than two tissues with sufficient cells")
    return(NULL)
  }

  res <- map_dfr(seq_along(tissues.use), function(i) {
    tt <- tissues.use[i]
    set.seed(seed + i)

    cells.1.all <- colnames(obj.lin)[obj.lin$tissue == tt]
    cells.2.all <- colnames(obj.lin)[obj.lin$tissue %in% setdiff(tissues.use, tt)]

    if (length(cells.1.all) < min_cells || length(cells.2.all) < min_cells) {
      return(NULL)
    }

    n1 <- min(length(cells.1.all), max_cells_per_group)
    n2 <- min(length(cells.2.all), max_cells_per_group)

    cells.1 <- sample(cells.1.all, n1)
    cells.2 <- sample(cells.2.all, n2)

    cells.use <- c(cells.1, cells.2)

    obj.test <- subset(obj.lin, cells = cells.use)

    obj.test$de_group <- ifelse(
      colnames(obj.test) %in% cells.1,
      "target",
      "other"
    )

    Idents(obj.test) <- "de_group"

    deg <- FindMarkers(
      obj.test,
      ident.1 = "target",
      ident.2 = "other",
      test.use = "wilcox",
      only.pos = TRUE,
      logfc.threshold = logfc.threshold,
      min.pct = min.pct
    )

    if (nrow(deg) == 0) return(NULL)

    deg %>%
      rownames_to_column("gene") %>%
      mutate(
        immune_lineage = lineage,
        tissue = tt,
        n_tissue_total = length(cells.1.all),
        n_other_total = length(cells.2.all),
        n_tissue_used = length(cells.1),
        n_other_used = length(cells.2)
      )
  })

  if (nrow(res) == 0) return(NULL)
  res
}

lineage_marker_res <- map_dfr(
  lineages.use,
  ~ run_lineage_tissue_markers(
    obj = obj,
    lineage = .x,
    min_cells = 200,
    max_cells_per_group = 5000,
    logfc.threshold = 0.25,
    min.pct = 0.10,
    seed = 123
  )
)

if (is.null(lineage_marker_res) || nrow(lineage_marker_res) == 0) {
  stop("No lineage-level tissue markers were detected. Check cell number thresholds.")
}

lineage_marker_filt <- lineage_marker_res %>%
  filter(
    avg_log2FC >= 0.35,
    pct.1 >= 0.20,
    pct.1 - pct.2 >= 0.10,
    p_val_adj < 1e-3
  )

if (nrow(lineage_marker_filt) == 0) {
  stop("No markers passed filtering. Consider relaxing avg_log2FC/min.pct thresholds.")
}

write.table(
  lineage_marker_res,
  file.path(outdir, "lineage_level_tissue_associated_markers_all.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  lineage_marker_filt,
  file.path(outdir, "lineage_level_tissue_associated_markers_filtered.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Step 3. Cluster lineage-level tissue-associated genes
top_n_per_lineage_tissue <- 100
min_cells_lineage_tissue <- 200

obj.use <- subset(obj, subset = immune_lineage %in% lineages.use)
obj.use$lineage_tissue <- paste(obj.use$immune_lineage, obj.use$tissue, sep = "|||")

lineage_tissue_counts <- obj.use@meta.data %>%
  count(immune_lineage, tissue, lineage_tissue, name = "n_cells") %>%
  arrange(immune_lineage, tissue)

write.table(
  lineage_tissue_counts,
  file.path(outdir, "lineage_tissue_cell_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

valid_lineage_tissue <- lineage_tissue_counts %>%
  filter(n_cells >= min_cells_lineage_tissue) %>%
  pull(lineage_tissue)

if (length(valid_lineage_tissue) < 2) {
  stop("Too few lineage_tissue groups with sufficient cells for gene program clustering.")
}

obj.use <- subset(obj.use, subset = lineage_tissue %in% valid_lineage_tissue)

top_lineage_markers <- lineage_marker_filt %>%
  mutate(lineage_tissue = paste(immune_lineage, tissue, sep = "|||")) %>%
  filter(lineage_tissue %in% valid_lineage_tissue) %>%
  group_by(immune_lineage, tissue) %>%
  arrange(desc(avg_log2FC), desc(pct.1 - pct.2)) %>%
  slice_head(n = top_n_per_lineage_tissue) %>%
  ungroup()

genes.use <- intersect(unique(top_lineage_markers$gene), rownames(obj))

if (length(genes.use) < 50) {
  stop("Too few candidate genes for clustering. Check marker filtering.")
}

avg_expr <- AverageExpression(
  obj.use,
  assays = "RNA",
  features = genes.use,
  group.by = "lineage_tissue",
  slot = "data"
)$RNA

avg_mat <- as.matrix(avg_expr)

# Remove uninformative genes across lineage-tissue groups.
gene_sd <- apply(avg_mat, 1, sd, na.rm = TRUE)
avg_mat <- avg_mat[gene_sd > quantile(gene_sd, 0.25, na.rm = TRUE), , drop = FALSE]

avg_mat_z <- t(scale(t(avg_mat)))
avg_mat_z[is.na(avg_mat_z)] <- 0

gene_dist <- dist(avg_mat_z, method = "euclidean")
gene_hc <- hclust(gene_dist, method = "ward.D2")
# )

k_candidates <- 4:15

k_eval <- lapply(k_candidates, function(k) {
  cl <- cutree(gene_hc, k = k)
  tab <- table(cl)
  sil <- silhouette(cl, gene_dist)
  
  data.frame(
    k = k,
    n_clusters = length(tab),
    min_cluster_size = min(tab),
    median_cluster_size = median(tab),
    max_cluster_size = max(tab),
    n_clusters_lt_10_genes = sum(tab < 10),
    mean_silhouette = mean(sil[, "sil_width"]),
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows()

write.table(
  k_eval,
  file.path(outdir, "gene_cluster_k_evaluation.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

print(k_eval)

# Final k should be selected after checking k_eval.
k_gene <- 12
gene_cluster <- cutree(gene_hc, k = k_gene)

gene_cluster_df <- data.frame(
  gene = names(gene_cluster),
  gene_cluster = paste0("C", gene_cluster),
  stringsAsFactors = FALSE
)

# Filter very small clusters before module scoring.
min_genes_per_cluster <- 5
cluster_size <- table(gene_cluster_df$gene_cluster)
valid_clusters <- names(cluster_size)[cluster_size >= min_genes_per_cluster]

gene_cluster_df <- gene_cluster_df %>%
  filter(gene_cluster %in% valid_clusters)

write.table(
  gene_cluster_df,
  file.path(outdir, "lineage_level_tissue_gene_clusters.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cluster_sig_df <- gene_cluster_df %>%
  group_by(gene_cluster) %>%
  summarise(genes = list(unique(gene)), .groups = "drop") %>%
  arrange(as.integer(sub("^C", "", gene_cluster)))

cluster_signatures <- setNames(cluster_sig_df$genes, cluster_sig_df$gene_cluster)

saveRDS(
  cluster_signatures,
  file.path(outdir, "lineage_level_tissue_gene_cluster_signatures.rds")
)

# Step 4. Map gene clusters back to immune_subtype
# Remove existing score columns before recalculating them.
old_cols <- grep("^tissue_gene_cluster_", colnames(obj@meta.data), value = TRUE)
if (length(old_cols) > 0) {
  obj@meta.data[, old_cols] <- NULL
}

obj <- AddModuleScore(
  obj,
  features = cluster_signatures,
  name = "tissue_gene_cluster_"
)

cluster_score_cols <- grep("^tissue_gene_cluster_", colnames(obj@meta.data), value = TRUE)
cluster_score_cols <- cluster_score_cols[order(match(cluster_score_cols, colnames(obj@meta.data)))]

if (length(cluster_score_cols) != length(cluster_signatures)) {
  stop("Number of module score columns does not match number of gene clusters.")
}

cluster_score_map <- data.frame(
  score_col = cluster_score_cols,
  gene_cluster = names(cluster_signatures),
  stringsAsFactors = FALSE
)

write.table(
  cluster_score_map,
  file.path(outdir, "tissue_gene_cluster_score_map.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cluster_score_df <- obj@meta.data %>%
  select(immune_lineage, immune_subtype, tissue, stage, all_of(cluster_score_cols)) %>%
  pivot_longer(
    cols = all_of(cluster_score_cols),
    names_to = "score_col",
    values_to = "score"
  ) %>%
  left_join(cluster_score_map, by = "score_col") %>%
  group_by(gene_cluster, immune_lineage, immune_subtype, tissue) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  left_join(subtype_class, by = "immune_subtype")

write.table(
  cluster_score_df,
  file.path(outdir, "tissue_gene_cluster_scores_by_subtype_tissue.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Step 5. Identify subtype drivers and classify program source
program_driver <- cluster_score_df %>%
  filter(n_cells >= 50) %>%
  group_by(gene_cluster) %>%
  arrange(desc(mean_score)) %>%
  slice_head(n = 10) %>%
  ungroup()

write.table(
  program_driver,
  file.path(outdir, "top_subtype_drivers_of_tissue_gene_cluster_programs.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

program_source <- program_driver %>%
  group_by(gene_cluster) %>%
  summarise(
    top_tissue_class = tissue_class[which.max(mean_score)],
    frac_top10_restricted = mean(tissue_class == "tissue_restricted", na.rm = TRUE),
    frac_top10_shared = mean(tissue_class == "tissue_shared", na.rm = TRUE),
    top_driver = paste0(immune_subtype[which.max(mean_score)], " | ", tissue[which.max(mean_score)]),
    .groups = "drop"
  ) %>%
  mutate(
    program_source = case_when(
      top_tissue_class == "tissue_restricted" & frac_top10_restricted >= 0.5 ~ "restricted_subtype_driven",
      top_tissue_class == "tissue_shared" & frac_top10_shared >= 0.5 ~ "shared_subtype_associated",
      TRUE ~ "mixed_or_intermediate"
    )
  )

write.table(
  program_source,
  file.path(outdir, "tissue_gene_cluster_program_source_classification.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Step 6. Within-subtype validation in tissue-shared subtypes
shared_subtypes <- subtype_class %>%
  filter(tissue_class == "tissue_shared") %>%
  pull(immune_subtype)

validate_within_subtype <- cluster_score_df %>%
  filter(
    immune_subtype %in% shared_subtypes,
    n_cells >= 50
  ) %>%
  group_by(gene_cluster, immune_subtype) %>%
  summarise(
    max_score = max(mean_score, na.rm = TRUE),
    min_score = min(mean_score, na.rm = TRUE),
    delta_score = max_score - min_score,
    tissue_with_max_score = tissue[which.max(mean_score)],
    n_tissues_tested = n_distinct(tissue),
    .groups = "drop"
  ) %>%
  filter(n_tissues_tested >= 2) %>%
  arrange(desc(delta_score))

write.table(
  validate_within_subtype,
  file.path(outdir, "within_shared_subtype_cluster_validation.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Optional marker-level validation in tissue-shared subtypes.
run_within_subtype_tissue_markers <- function(obj,
                                              subtype,
                                              min_cells = 100,
                                              max_cells_per_group = 3000,
                                              logfc.threshold = 0.25,
                                              min.pct = 0.10,
                                              seed = 123) {
  message("Subtype: ", subtype)

  obj.sub <- subset(obj, subset = immune_subtype == subtype)
  obj.sub$tissue <- droplevels(factor(obj.sub$tissue))

  tissue_counts <- table(obj.sub$tissue)
  tissues.use <- names(tissue_counts)[tissue_counts >= min_cells]

  if (length(tissues.use) < 2) return(NULL)

  Idents(obj.sub) <- "tissue"

  res <- map_dfr(seq_along(tissues.use), function(i) {
    tt <- tissues.use[i]
    set.seed(seed + i)

    cells.1.all <- WhichCells(obj.sub, idents = tt)
    cells.2.all <- WhichCells(obj.sub, idents = setdiff(tissues.use, tt))

    if (length(cells.1.all) < min_cells || length(cells.2.all) < min_cells) return(NULL)

    n1 <- min(length(cells.1.all), max_cells_per_group)
    n2 <- min(length(cells.2.all), max_cells_per_group)

    cells.1 <- sample(cells.1.all, n1)
    cells.2 <- sample(cells.2.all, n2)

    cells.use <- c(cells.1, cells.2)
    obj.test <- subset(obj.sub, cells = cells.use)

    obj.test$de_group <- ifelse(
      colnames(obj.test) %in% cells.1,
      "target",
      "other"
    )

    Idents(obj.test) <- "de_group"

    deg <- FindMarkers(
      obj.test,
      ident.1 = "target",
      ident.2 = "other",
      test.use = "wilcox",
      only.pos = TRUE,
      logfc.threshold = logfc.threshold,
      min.pct = min.pct
    )

    if (nrow(deg) == 0) return(NULL)

    deg %>%
      rownames_to_column("gene") %>%
      mutate(
        immune_subtype = subtype,
        tissue = tt,
        n_tissue_total = length(cells.1.all),
        n_other_total = length(cells.2.all),
        n_tissue_used = length(cells.1),
        n_other_used = length(cells.2)
      )
  })

  if (nrow(res) == 0) return(NULL)

  res
}

within_subtype_markers <- map_dfr(
  shared_subtypes,
  ~ run_within_subtype_tissue_markers(
    obj = obj,
    subtype = .x,
    min_cells = 100,
    max_cells_per_group = 3000,
    logfc.threshold = 0.25,
    min.pct = 0.10,
    seed = 123
  )
)

if (!is.null(within_subtype_markers) && nrow(within_subtype_markers) > 0) {
  within_subtype_markers_filt <- within_subtype_markers %>%
    filter(
      avg_log2FC >= 0.30,
      pct.1 >= 0.20,
      pct.1 - pct.2 >= 0.10,
      p_val_adj < 1e-3
    )

  write.table(
    within_subtype_markers,
    file.path(outdir, "within_shared_subtype_tissue_markers_all.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  write.table(
    within_subtype_markers_filt,
    file.path(outdir, "within_shared_subtype_tissue_markers_filtered.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
} else {
  warning("No within-subtype tissue markers were detected.")
}

# Step 7. Check stage confounding / stage support
prog_stage_df <- obj@meta.data %>%
  select(immune_lineage, immune_subtype, tissue, stage, all_of(cluster_score_cols)) %>%
  pivot_longer(
    cols = all_of(cluster_score_cols),
    names_to = "score_col",
    values_to = "score"
  ) %>%
  left_join(cluster_score_map, by = "score_col") %>%
  group_by(gene_cluster, immune_lineage, immune_subtype, tissue, stage) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  left_join(subtype_class, by = "immune_subtype")

write.table(
  prog_stage_df,
  file.path(outdir, "tissue_gene_cluster_scores_by_subtype_tissue_stage.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

stage_support <- prog_stage_df %>%
  filter(n_cells >= 30) %>%
  group_by(gene_cluster, immune_subtype, tissue) %>%
  summarise(
    n_stages_observed = n_distinct(stage),
    mean_score_across_stages = mean(mean_score, na.rm = TRUE),
    .groups = "drop"
  )

write.table(
  stage_support,
  file.path(outdir, "tissue_gene_cluster_stage_support.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Save object with module scores for downstream Result 5.
saveRDS(
  obj,
  file.path(outdir, "immune_with_tissue_gene_cluster_scores.rds")
)

message("Result 4 tissue-directed program pipeline finished successfully.")
