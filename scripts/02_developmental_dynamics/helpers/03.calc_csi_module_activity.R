# Calculate mean CSI-module regulon AUC by cell group.
# Input: module membership, SCENIC regulon AUC, and cell metadata.
# Output: module-by-cell-group activity matrix.
calc_csi_module_activity <- function(clusters_df, regulonAUC, metadata, cell_type_column) {
  metadata <- metadata[!is.na(metadata[[cell_type_column]]), ]

  metadata$cell_type <- metadata[, cell_type_column]
  cell_types <- unique(metadata$cell_type)
  regulons <- unique(clusters_df$regulon)

  regulonAUC_sub <- regulonAUC@assays@data@listData$AUC
  regulonAUC_sub <- regulonAUC_sub[regulons, , drop = FALSE]

  csi_activity_matrix_list <- list()
  csi_cluster_activity <- data.frame("csi_cluster" = c(),
                                     "mean_activity" = c(),
                                     "cell_type" = c())

  for (ct in cell_types) {
    cell_indices <- rownames(subset(metadata, cell_type == ct))
    if (length(cell_indices) > 0 && ncol(regulonAUC_sub[, cell_indices, drop = FALSE]) > 1) {
      cell_type_aucs <- rowMeans(regulonAUC_sub[, cell_indices, drop = FALSE])
    } else {
      cell_type_aucs <- numeric(nrow(regulonAUC_sub))
      names(cell_type_aucs) <- rownames(regulonAUC_sub)
    }

    cell_type_aucs_df <- data.frame("regulon" = names(cell_type_aucs),
                                   "activity" = cell_type_aucs,
                                   "cell_type" = ct)
    csi_activity_matrix_list[[ct]] <- cell_type_aucs_df
  }

  for (ct in names(csi_activity_matrix_list)) {
    for (cluster in unique(clusters_df$csi_cluster)) {
      csi_regulon <- subset(clusters_df, csi_cluster == cluster)
      csi_regulon_activity <- subset(csi_activity_matrix_list[[ct]], regulon %in% csi_regulon$regulon)
      csi_activity_mean <- mean(csi_regulon_activity$activity, na.rm = TRUE)
      this_cluster_ct_activity <- data.frame("csi_cluster" = cluster,
                                             "mean_activity" = csi_activity_mean,
                                             "cell_type" = ct)
      csi_cluster_activity <- rbind(csi_cluster_activity, this_cluster_ct_activity)
    }
  }

  csi_cluster_activity[is.na(csi_cluster_activity)] <- 0

  csi_cluster_activity_wide <- tidyr::spread(csi_cluster_activity, cell_type, mean_activity)
  rownames(csi_cluster_activity_wide) <- csi_cluster_activity_wide$csi_cluster
  csi_cluster_activity_wide <- as.matrix(csi_cluster_activity_wide[2:ncol(csi_cluster_activity_wide)])

  return(csi_cluster_activity_wide)
}
