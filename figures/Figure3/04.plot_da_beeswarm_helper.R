# Purpose: plot da beeswarm helper.

# Plot differential-abundance results as a beeswarm.
plotDAbeeswarm <- function(da.res, group.by=NULL, alpha=0.1, subset.nhoods=NULL){
  if (!is.null(group.by)) {
    if (!group.by %in% colnames(da.res)) {
      stop(group.by, " is not a column in da.res. Have you forgot to run annotateNhoods(x, da.res, ", group.by,")?")
    }
    if (is.numeric(da.res[,group.by])) {
    }
    da.res <- mutate(da.res, group_by = da.res[,group.by])
  } else {
    da.res <- mutate(da.res, group_by = "g1")
  }

  if (!is.factor(da.res[,"group_by"])) {
    message("Converting group_by to factor...")
    da.res <- mutate(da.res, group_by = factor(group_by, levels=unique(group_by)))
  }

  if (!is.null(subset.nhoods)) {
    da.res <- da.res[subset.nhoods,]
  }

  # Get position with ggbeeswarm
  beeswarm_pos <- ggplot_build(
    da.res %>%
      mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
      arrange(group_by) %>%
      ggplot(aes(group_by, logFC)) +
      geom_quasirandom()
  )

  pos_x <- beeswarm_pos$data[[1]]$x
  pos_y <- beeswarm_pos$data[[1]]$y

  n_groups <- unique(da.res$group_by) %>% length()

  da.res %>%
    mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
    mutate(logFC_color = ifelse(is_signif==1, logFC, NA)) %>%
    arrange(group_by) %>%
    mutate(Nhood=factor(Nhood, levels=unique(Nhood))) %>%
    mutate(pos_x = pos_x, pos_y=pos_y) %>%
    ggplot(aes(pos_x, pos_y, color=logFC_color)) +
    scale_color_gradient2() +
    guides(color="none") +
    xlab(group.by) + ylab("Log Fold Change") +
    scale_x_continuous(
      breaks = seq(1,n_groups),
      labels = setNames(levels(da.res$group_by), seq(1,n_groups))
      ) +
    geom_point() +
    coord_flip() +
    theme_bw(base_size=12) +
    theme(strip.text.y =  element_text(angle=0),
          aspect.ratio = 2)
}
