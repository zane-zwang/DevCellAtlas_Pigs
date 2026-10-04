# Purpose: tai model decomposition.

# TAI model effects and tissue-level decomposition.
# Match filename-derived cell-type identifiers to metadata cell types before
# constructing tissue-stage-cell-type analysis units.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(car)
  library(effectsize)
  library(emmeans)
  library(scales)
})

# 0. Parameters
tai_rds  <- "tai_ensembl/pseudo_exp_raw/cell/pseudoExp_cell_tai_list_log.rds"
meta_csv <- file.path(Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input"), "atlas/all_obs.csv")

stage_order <- c("E55", "E90", "P0", "P30", "P90", "P180")

min_cells_unit <- 30L
min_common_mass <- 0.80
weight_mode <- "sqrt_n"   # "none", "sqrt_n", or "n"

sanitize_celltype <- function(x) {
  x <- trimws(as.character(x))
  gsub("[^a-zA-Z0-9_\\-]", "_", x)
}

standardize_stage <- function(x) {
  dplyr::recode(
    as.character(x),
    "e55d" = "E55",
    "E55" = "E55",
    "e90d" = "E90",
    "E90" = "E90",
    "X0d" = "P0",
    "0d" = "P0",
    "P0" = "P0",
    "X30d" = "P30",
    "30d" = "P30",
    "P30" = "P30",
    "X90d" = "P90",
    "90d" = "P90",
    "P90" = "P90",
    "X180d" = "P180",
    "180d" = "P180",
    "P180" = "P180",
    .default = as.character(x)
  )
}

standardize_tissue <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[[:space:]_-]+", " ", x)

  dplyr::recode(
    x,
    "skeletal muscle" = "muscle",
    "skeletal muscles" = "muscle",
    "muscle" = "muscle",
    "cerebral cortex" = "cerebrum",
    "cortex" = "cerebrum",
    .default = x
  )
}

# Critical repair:
# Strip BOTH the tissue prefix and filename suffix.
extract_celltype_from_filename <- function(filename, tissue) {
  x <- basename(as.character(filename))

  # Remove suffix; accept both Exp and aveExp variants.
  x <- sub(
    "_stage_pseudo_(ave)?Exp\\.csv$",
    "",
    x,
    ignore.case = TRUE
  )

  # Remove exact tissue prefix.
  tissue_pattern <- paste0("^", stringr::str_replace_all(tissue, "([.\\+*?\\[\\](){}^$|\\\\])", "\\\\\\1"), "_")
  x <- sub(tissue_pattern, "", x, ignore.case = TRUE)

  sanitize_celltype(x)
}

check_required_columns <- function(df, cols, object_name) {
  missing <- setdiff(cols, colnames(df))
  if (length(missing) > 0) {
    stop(
      object_name, " is missing columns: ",
      paste(missing, collapse = ", ")
    )
  }
}

# 1. Read and inspect nested TAI list
pseudoExp_cell_tai_list <- readRDS(tai_rds)

if (!is.list(pseudoExp_cell_tai_list) ||
    is.null(names(pseudoExp_cell_tai_list))) {
  stop("The TAI RDS must be a named nested list: tissue -> celltype file -> TAI vector.")
}

message("TAI tissues: ", paste(names(pseudoExp_cell_tai_list), collapse = ", "))

# 2. Convert nested list to long format
tai_long <- imap_dfr(
  pseudoExp_cell_tai_list,
  function(cell_list, tissue_name) {
    imap_dfr(
      cell_list,
      function(tai_vec, file_name) {
        stage_names <- names(tai_vec)

        if (is.null(stage_names)) {
          stop("TAI vector has no stage names: ", file_name)
        }

        tibble(
          tissue = standardize_tissue(tissue_name),
          source_name = as.character(file_name),
          safe_cell = extract_celltype_from_filename(
            filename = file_name,
            tissue = tissue_name
          ),
          stage = standardize_stage(stage_names),
          TAI = as.numeric(tai_vec)
        )
      }
    )
  }
) %>%
  mutate(
    stage = factor(stage, levels = stage_order, ordered = TRUE)
  ) %>%
  filter(!is.na(TAI), !is.na(stage), safe_cell != "")

message("TAI rows after list expansion: ", nrow(tai_long))
message("TAI cell types: ", n_distinct(tai_long$safe_cell))

# 3. Read metadata and count cells
meta <- read.csv(
  meta_csv,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

check_required_columns(
  meta,
  c("tissue", "stage", "celltype", "celllineage"),
  "meta"
)

count_df <- meta %>%
  as_tibble() %>%
  transmute(
    tissue = standardize_tissue(tissue),
    stage = standardize_stage(stage),
    safe_cell = sanitize_celltype(celltype),
    celllineage = as.character(celllineage)
  ) %>%
  filter(
    !is.na(tissue),
    !is.na(stage),
    !is.na(safe_cell),
    !is.na(celllineage),
    safe_cell != ""
  ) %>%
  count(
    tissue,
    stage,
    safe_cell,
    celllineage,
    name = "n_cells"
  ) %>%
  mutate(
    stage = factor(stage, levels = stage_order, ordered = TRUE)
  )

# 4. Diagnose the join BEFORE filtering
tai_keys <- tai_long %>%
  distinct(tissue, stage, safe_cell)

meta_keys <- count_df %>%
  distinct(tissue, stage, safe_cell)

unmatched_tai <- anti_join(
  tai_keys,
  meta_keys,
  by = c("tissue", "stage", "safe_cell")
)

unmatched_meta <- anti_join(
  meta_keys,
  tai_keys,
  by = c("tissue", "stage", "safe_cell")
)

write_tsv(unmatched_tai, "TAI_unmatched_keys.tsv")
write_tsv(unmatched_meta, "metadata_unmatched_keys.tsv")

joined <- tai_long %>%
  left_join(
    count_df,
    by = c("tissue", "stage", "safe_cell")
  )

join_summary <- joined %>%
  summarise(
    n_tai_rows = n(),
    n_matched_rows = sum(!is.na(n_cells) & !is.na(celllineage)),
    match_rate = n_matched_rows / n_tai_rows
  )

print(join_summary)

if (join_summary$n_matched_rows == 0) {
  stop(
    "Zero TAI rows matched metadata. Inspect TAI_unmatched_keys.tsv ",
    "and metadata_unmatched_keys.tsv."
  )
}

if (join_summary$match_rate < 0.80) {
  warning(
    "TAI-metadata match rate is below 80%: ",
    round(100 * join_summary$match_rate, 1), "%."
  )
}

# 5. Construct the modeling table
tai_unit <- joined %>%
  filter(
    !is.na(celllineage),
    !is.na(n_cells),
    n_cells >= min_cells_unit
  ) %>%
  group_by(tissue, stage, safe_cell, celllineage) %>%
  summarise(
    TAI = weighted.mean(TAI, w = n_cells, na.rm = TRUE),
    n_cells = sum(n_cells),
    .groups = "drop"
  ) %>%
  mutate(
    tissue = factor(tissue),
    stage = factor(stage, levels = stage_order, ordered = FALSE),
    safe_cell = factor(safe_cell),
    celllineage = factor(celllineage)
  ) %>%
  droplevels()

write_tsv(tai_unit, "TAI_celltype_tissue_stage_long.tsv")

if (nrow(tai_unit) == 0) {
  stop(
    "tai_unit is empty after applying n_cells >= ",
    min_cells_unit,
    ". Lower min_cells_unit or inspect the join diagnostics."
  )
}

factor_check <- tibble(
  variable = c("tissue", "celllineage", "stage"),
  n_levels = c(
    nlevels(tai_unit$tissue),
    nlevels(tai_unit$celllineage),
    nlevels(tai_unit$stage)
  ),
  values = c(
    paste(levels(tai_unit$tissue), collapse = ", "),
    paste(levels(tai_unit$celllineage), collapse = ", "),
    paste(levels(tai_unit$stage), collapse = ", ")
  )
)

print(factor_check)

bad_factors <- factor_check %>%
  filter(n_levels < 2) %>%
  pull(variable)

if (length(bad_factors) > 0) {
  stop(
    "Cannot fit the model; fewer than two levels remain for: ",
    paste(bad_factors, collapse = ", ")
  )
}

# Require adequate stage coverage for estimating lineage x stage.
eligible_lineages <- tai_unit %>%
  distinct(celllineage, stage) %>%
  count(celllineage, name = "n_stages") %>%
  filter(n_stages >= 4) %>%
  pull(celllineage)

tai_model_df <- tai_unit %>%
  filter(celllineage %in% eligible_lineages) %>%
  droplevels()

if (nlevels(tai_model_df$celllineage) < 2) {
  stop("Fewer than two eligible lineages remain after stage-coverage filtering.")
}

# D. Model effects
tai_model_df <- tai_model_df %>%
  mutate(
    model_weight = case_when(
      weight_mode == "none" ~ 1,
      weight_mode == "sqrt_n" ~ sqrt(n_cells),
      weight_mode == "n" ~ as.numeric(n_cells),
      TRUE ~ NA_real_
    )
  )

if (anyNA(tai_model_df$model_weight)) {
  stop("Unsupported weight_mode: ", weight_mode)
}

old_contrasts <- options("contrasts")
on.exit(options(contrasts = old_contrasts$contrasts), add = TRUE)
options(contrasts = c("contr.sum", "contr.poly"))

fit_D <- lm(
  TAI ~ tissue + celllineage * stage,
  data = tai_model_df,
  weights = model_weight
)

anova_D <- car::Anova(fit_D, type = 3)

eta_D <- effectsize::eta_squared(
  anova_D,
  partial = TRUE
) %>%
  as.data.frame() %>%
  as_tibble()

effect_D <- eta_D %>%
  filter(Parameter != "(Intercept)", Parameter != "Residuals") %>%
  transmute(
    term_raw = Parameter,
    term = dplyr::recode(
      Parameter,
      "tissue" = "Tissue",
      "celllineage" = "Lineage",
      "stage" = "Stage",
      "celllineage:stage" = "Lineage × Stage",
      .default = Parameter
    ),
    partial_eta2 = Eta2_partial,
    ci_low = CI_low,
    ci_high = CI_high
  ) %>%
  mutate(
    term = factor(
      term,
      levels = c("Tissue", "Lineage", "Stage", "Lineage × Stage")
    )
  )

write_tsv(effect_D, "TAI_model_effects.tsv")

# E. Exact symmetric two-component decomposition
prep_tissue_stage <- tai_unit %>%
  group_by(tissue, stage) %>%
  mutate(p_all = n_cells / sum(n_cells)) %>%
  ungroup()

decompose_pair <- function(df_tissue, stage0, stage1) {
  d0 <- df_tissue %>%
    filter(stage == stage0) %>%
    select(safe_cell, TAI0 = TAI, p0_all = p_all)

  d1 <- df_tissue %>%
    filter(stage == stage1) %>%
    select(safe_cell, TAI1 = TAI, p1_all = p_all)

  common <- inner_join(d0, d1, by = "safe_cell")

  if (nrow(common) == 0) {
    return(tibble())
  }

  coverage0 <- sum(common$p0_all)
  coverage1 <- sum(common$p1_all)

  common <- common %>%
    mutate(
      p0 = p0_all / coverage0,
      p1 = p1_all / coverage1,
      within = ((p0 + p1) / 2) * (TAI1 - TAI0),
      composition = ((TAI0 + TAI1) / 2) * (p1 - p0)
    )

  tibble(
    stage_from = stage0,
    stage_to = stage1,
    n_common_celltypes = nrow(common),
    coverage_from = coverage0,
    coverage_to = coverage1,
    TAI_from_common = sum(common$p0 * common$TAI0),
    TAI_to_common = sum(common$p1 * common$TAI1),
    delta_TAI = TAI_to_common - TAI_from_common,
    within = sum(common$within),
    composition = sum(common$composition),
    reconstruction_error = delta_TAI - within - composition
  )
}

stage_pairs <- tibble(
  stage_from = stage_order[-length(stage_order)],
  stage_to = stage_order[-1]
)

decomp_pairwise <- prep_tissue_stage %>%
  split(.$tissue) %>%
  imap_dfr(
    function(df_tissue, tissue_name) {
      map2_dfr(
        stage_pairs$stage_from,
        stage_pairs$stage_to,
        ~ decompose_pair(df_tissue, .x, .y)
      ) %>%
        mutate(tissue = tissue_name, .before = 1)
    }
  ) %>%
  mutate(
    transition = paste(stage_from, stage_to, sep = "→"),
    transition = factor(
      transition,
      levels = paste(
        stage_order[-length(stage_order)],
        stage_order[-1],
        sep = "→"
      ),
      ordered = TRUE
    ),
    pass_coverage = coverage_from >= min_common_mass &
      coverage_to >= min_common_mass
  )

write_tsv(decomp_pairwise, "TAI_decomposition.tsv")

summary_df <- decomp_pairwise %>%
  filter(pass_coverage) %>%
  group_by(tissue) %>%
  summarise(
    within_total = sum(abs(within)),
    composition_total = sum(abs(composition))
  ) %>%
  mutate(
    total = within_total + composition_total,
    within_fraction =
      within_total / total,
    composition_fraction =
      composition_total / total
  )


# plot E
plot_E <- summary_df %>%
  select(
    tissue,
    within_fraction,
    composition_fraction
  ) %>%
  pivot_longer(
    cols = c(
      within_fraction,
      composition_fraction
    ),
    names_to = "component",
    values_to = "fraction"
  ) %>%
  mutate(
    component = dplyr::recode(
      component,
      within_fraction =
        "Within-cell-type remodeling",
      composition_fraction =
        "Cell composition change"
    ),
    component = factor(
      component,
      levels=c(
        "Within-cell-type remodeling",
        "Cell composition change"
      )
    )
  )

plot_E

message("Finished successfully.")


decomp_use <- decomp_pairwise %>%
  filter(
    coverage_from >= 0.8,
    coverage_to >=0.8
  )

summary_contribution <- decomp_use %>%
  group_by(tissue) %>%
  summarise(
    within_total = sum(abs(within),na.rm=T),
    composition_total = sum(abs(composition),na.rm=T)
  ) %>%
  mutate(
    total = within_total + composition_total,
    within_fraction = within_total/total,
    composition_fraction = composition_total/total
    )


plot_summary <- summary_contribution %>%
  pivot_longer(
    cols=c(
      within_fraction,
      composition_fraction
    ),
    names_to="component",
    values_to="fraction"
  ) %>%
  mutate(
    
    component =
      dplyr::recode(
        component,
        
        within_fraction=
          "Within-cell-type",
        
        composition_fraction=
          "Cell composition"
      )
    
  )

plot_summary$tissue <-
  factor(
    plot_summary$tissue,
    levels=rev(unique(plot_summary$tissue))
  )

heat_df <- decomp_use %>%
  select(tissue, transition, within, composition)%>%
  pivot_longer(
    cols=c(within,composition),
    names_to="component",
    values_to="value"
  )

plot_df <- decomp_pairwise %>%
  filter(pass_coverage) %>%
  mutate(
    total_change = abs(within) + abs(composition)
  )

p <- ggplot(plot_df, aes(x = composition, y = within)) +
  geom_point(aes(size = total_change, color = transition), alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~tissue) +
  theme_bw()+
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    aspect.ratio = 0.75
    )

p

# ----- new ----
emm_ls <- emmeans(
  fit_D,
  specs = ~ stage | celllineage,
  weights = "proportional"
)

emm_ls_df <- as.data.frame(emm_ls) %>%
  as_tibble() %>%
  mutate(
    stage = factor(
      stage,
      levels = c("E55", "E90", "P0", "P30", "P90", "P180"),
      ordered = TRUE
    )
  )

print(emm_ls_df)


stage_order <- c("E55", "E90", "P0", "P30", "P90", "P180")

delta_df <- emm_ls_df %>%
  select(celllineage, stage, emmean) %>%
  mutate(stage = factor(stage, levels = stage_order, ordered = TRUE)) %>%
  arrange(celllineage, stage) %>%
  group_by(celllineage) %>%
  mutate(
    previous_stage = lag(stage),
    delta_TAI = emmean - lag(emmean)
  ) %>%
  ungroup() %>%
  filter(!is.na(delta_TAI)) %>%
  mutate(
    transition = paste0(previous_stage, " \u2192 ", stage),
    transition = factor(
      transition,
      levels = paste0(
        stage_order[-length(stage_order)],
        " \u2192 ",
        stage_order[-1]
      ),
      ordered = TRUE
    )
  )


max_abs_delta <- max(abs(delta_df$delta_TAI), na.rm = TRUE)

pD_new <- ggplot(
  delta_df,
  aes(x = transition, y = celllineage, fill = delta_TAI)
) +
  geom_tile(
    colour = "white",
    linewidth = 0.8,
    width = 0.96,
    height = 0.94
  ) +
  scale_fill_gradient2(
    low = "#3B6FB6",
    mid = "#F7F7F7",
    high = "#B24745",
    midpoint = 0,
    limits = c(-max_abs_delta, max_abs_delta),
    oob = squish,
    name = expression(Delta * "TAI")
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      angle = 40,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      face = "italic"
    ),
    legend.position = "right",
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

pD_new

stage_contrasts <- contrast(
  emm_ls,
  method = "consec",
  reverse = FALSE,
  adjust = "BH"
)

stage_contrast_df <- as.data.frame(stage_contrasts) %>%
  dplyr::as_tibble()


stage_contrast_df2 <- stage_contrast_df %>%
  mutate(
    transition = gsub(" - ", " \u2192 ", contrast),
    
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )

effect_D_plot <- effect_D %>%
  mutate(
    term = factor(
      term,
      levels = rev(
        c("Tissue", "Lineage", "Stage", "Lineage × Stage")
      )
    )
  )


