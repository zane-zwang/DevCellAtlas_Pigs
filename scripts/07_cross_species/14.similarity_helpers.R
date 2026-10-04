# Purpose: similarity helpers.

library(Seurat)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(Matrix)
library(edgeR)
library(ComplexHeatmap)
library(circlize)

# define functions
make_sample_celltype_pseudobulk <- function(obj,
                                            species_col = "species",
                                            sample_col = "orig.ident",
                                            celltype_col = "celltype",
                                            assay = "RNA",
                                            layer = "counts",
                                            min_cells = 30) {
  counts <- GetAssayData(obj, assay = assay, slot = layer)
  meta <- obj@meta.data
  meta$cell_id <- rownames(meta)
  
  group_count <- meta %>%
    count(
      species = .data[[species_col]],
      sample = .data[[sample_col]],
      celltype = .data[[celltype_col]],
      name = "n_cells"
    ) %>%
    filter(n_cells >= min_cells)
  
  meta_use <- meta %>%
    mutate(
      species = .data[[species_col]],
      sample = .data[[sample_col]],
      celltype = .data[[celltype_col]]
    ) %>%
    inner_join(
      group_count,
      by = c("species", "sample", "celltype")
    ) %>%
    mutate(group = paste(species, sample, celltype, sep = "|"))
  
  groups <- unique(meta_use$group)
  
  pb_counts <- sapply(groups, function(g) {
    cells <- meta_use$cell_id[meta_use$group == g]
    Matrix::rowSums(counts[, cells, drop = FALSE])
  })
  
  pb_counts <- as.matrix(pb_counts)
  pb_logcpm <- edgeR::cpm(pb_counts, log = TRUE, prior.count = 1)
  
  pb_info <- tibble(group = colnames(pb_logcpm)) %>%
    separate(group, into = c("species", "sample", "celltype"), sep = "\\|", remove = FALSE)
  
  list(
    counts = pb_counts,
    logcpm = pb_logcpm,
    info = pb_info
  )
}


calc_pairwise_human_similarity <- function(
    pb_mat,
    pb_info,
    celltype_use,
    species_col = "species",
    sample_col = "sample",
    celltype_col = "celltype",
    gene_use = NULL,
    method = "spearman"
) {
  if (is.null(gene_use)) {
    gene_use <- rownames(pb_mat)
  }
  
  gene_use <- intersect(gene_use, rownames(pb_mat))
  mat <- pb_mat[gene_use, , drop = FALSE]
  
  info_ct <- pb_info %>%
    filter(.data[[celltype_col]] == celltype_use)
  
  human_pb <- info_ct %>% filter(.data[[species_col]] == "human")
  pig_pb   <- info_ct %>% filter(.data[[species_col]] == "pig")
  mouse_pb <- info_ct %>% filter(.data[[species_col]] == "mouse")
  
  hp <- expand.grid(
    human_pb = human_pb$pb_id,
    animal_pb = pig_pb$pb_id,
    stringsAsFactors = FALSE
  ) %>%
    mutate(pair = "human_pig")
  
  hm <- expand.grid(
    human_pb = human_pb$pb_id,
    animal_pb = mouse_pb$pb_id,
    stringsAsFactors = FALSE
  ) %>%
    mutate(pair = "human_mouse")
  
  pair_df <- bind_rows(hp, hm)
  
  res <- pair_df %>%
    rowwise() %>%
    mutate(
      correlation = suppressWarnings(
        cor(
          mat[, human_pb],
          mat[, animal_pb],
          method = method,
          use = "pairwise.complete.obs"
        )
      )
    ) %>%
    ungroup() %>%
    left_join(
      info_ct %>%
        select(
          human_pb = pb_id,
          human_sample = all_of(sample_col)
        ),
      by = "human_pb"
    ) %>%
    left_join(
      info_ct %>%
        select(
          animal_pb = pb_id,
          animal_species = all_of(species_col),
          animal_sample = all_of(sample_col)
        ),
      by = "animal_pb"
    ) %>%
    mutate(celltype = celltype_use)
  
  res
}

calc_human_marker_program_similarity <- function(
    pb_mat,
    pb_info,
    human_deg,
    celltype_use,
    top_n = 200,
    species_col = "species",
    sample_col = "sample",
    celltype_col = "celltype"
) {
  fc_col <- if ("avg_log2FC" %in% colnames(human_deg)) "avg_log2FC" else "avg_logFC"
  
  genes <- human_deg %>%
    filter(cluster == celltype_use, p_val_adj < 0.05) %>%
    arrange(desc(.data[[fc_col]])) %>%
    pull(gene) %>%
    unique() %>%
    head(top_n)
  
  calc_pairwise_human_similarity(
    pb_mat = pb_mat,
    pb_info = pb_info,
    celltype_use = celltype_use,
    species_col = species_col,
    sample_col = sample_col,
    celltype_col = celltype_col,
    gene_use = genes,
    method = "spearman"
  ) %>%
    mutate(
      program = paste0(celltype_use, "_human_marker_top", top_n),
      n_marker = length(intersect(genes, rownames(pb_mat)))
    )
}

calc_pairwise_celltype_cor <- function(pb_mat,
                                       pb_info,
                                       celltypes,
                                       gene_use = NULL,
                                       method = "spearman") {
  if (is.null(gene_use)) {
    gene_use <- rownames(pb_mat)
  }
  
  gene_use <- intersect(gene_use, rownames(pb_mat))
  mat <- pb_mat[gene_use, , drop = FALSE]
  
  res <- lapply(celltypes, function(ct) {
    info_ct <- pb_info %>% filter(celltype == ct)
    
    h_groups <- info_ct %>% filter(species == "human")
    m_groups <- info_ct %>% filter(species == "mouse")
    p_groups <- info_ct %>% filter(species == "pig")
    
    if (nrow(h_groups) == 0 || nrow(m_groups) == 0 || nrow(p_groups) == 0) {
      return(NULL)
    }
    
    hp <- expand.grid(
      human_group = h_groups$group,
      animal_group = p_groups$group,
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        celltype = ct,
        comparison = "human_pig",
        human_sample = sub("^human\\|([^|]+)\\|.*$", "\\1", human_group),
        animal_sample = sub("^pig\\|([^|]+)\\|.*$", "\\1", animal_group),
        cor = map2_dbl(human_group, animal_group, ~ cor(mat[, .x], mat[, .y], method = method))
      )
    
    hm <- expand.grid(
      human_group = h_groups$group,
      animal_group = m_groups$group,
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        celltype = ct,
        comparison = "human_mouse",
        human_sample = sub("^human\\|([^|]+)\\|.*$", "\\1", human_group),
        animal_sample = sub("^mouse\\|([^|]+)\\|.*$", "\\1", animal_group),
        cor = map2_dbl(human_group, animal_group, ~ cor(mat[, .x], mat[, .y], method = method))
      )
    
    bind_rows(hp, hm)
  }) %>%
    bind_rows()
  
  res
}

calc_human_specific_program_cor <- function(pb_mat,
                                            pb_info,
                                            human_spec_list,
                                            method = "spearman",
                                            min_genes = 10) {
  res <- lapply(seq_len(nrow(human_spec_list)), function(i) {
    ct <- human_spec_list$cluster[i]
    genes <- human_spec_list$genes[[i]]
    genes <- intersect(genes, rownames(pb_mat))
    
    if (length(genes) < min_genes) {
      return(NULL)
    }
    
    calc_pairwise_celltype_cor(
      pb_mat = pb_mat,
      pb_info = pb_info,
      celltypes = ct,
      gene_use = genes,
      method = method
    ) %>%
      mutate(
        gene_set = paste0(ct, "_human_specific_program"),
        n_genes = length(genes)
      )
  }) %>%
    bind_rows()
  
  res
}
