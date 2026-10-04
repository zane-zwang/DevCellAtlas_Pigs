# Purpose: qtl validation.

# Developmental bulk-eQTL / ieQTL / cseQTL validation workflow
# Purpose:
#   1) Validate atlas-guided deconvolution
#   2) Quantify stage/bin-specific eQTL effect heterogeneity
#   3) Evaluate ieQTL/cseQTL calibration and power
#   4) Compare ieQTL and cseQTL effect direction and overlap
#   5) Summarize coloc evidence using matched loci
#   6) Evaluate enrichment of xeGenes among bin-shared bulk eGenes
# Input: association, mashR, coloc and cell-fraction tables.
# Output: validation plots and summary tables.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(broom)
})

# 0. CONFIG

cfg <- list(
  outdir = "eQTL_validation_figures",
  bins = c("bin1", "bin2", "bin3", "bin4"),
  bin_labels = c(
    bin1 = "Prenatal",
    bin2 = "Early postnatal / juvenile",
    bin3 = "Adolescent",
    bin4 = "Adult"
  ),

  # Sample-level metadata and cell fractions
  sample_info = "bisque_prop_bin.CSV",
  sample_col = "Sample",
  bin_col = "group",

  # Optional second deconvolution result for method-concordance analysis.
  # Set to NA_character_ to skip.
  bisque_prop = "bisque/bulk_props_bisque.csv",
  cibersort_prop = "MuscleCIBERSORT-Results.txt",
  cibersort_sample_col = "Mixture",

  # Bulk eQTL mashr output
  mash_post_mean = "mashr_res/01_beqtl_mashr_output/posterior_means.rds",
  mash_lfsr = "mashr_res/01_beqtl_mashr_output/lfsr.rds",
  mash_pair_annot = "mashr_res/01_beqtl_mashr_output/pair_annot.tsv",

  # Optional frequentist effect and SE matrices for formal heterogeneity tests.
  # Rows must be pair_id, columns must be bins.
  # If absent, the heterogeneity section is skipped.
  bulk_beta_mat = "mashr_res/01_beqtl_mashr_input/strong_set_beta.rds",
  bulk_se_mat = "mashr_res/01_beqtl_mashr_input/strong_set_se.rds",
  
  # Significant hit tables
  bulk_hits = "eQTL_res/02_acat/beQTL/02_filter/all_significant_beQTLs.tsv",
  ie_hits = "eQTL_res/02_acat/ieQTL/02_filter/all_significant_ieQTLs.tsv",
  cse_hits = "eQTL_res/02_acat/cseQTL/02_filter/all_significant_cseQTLs.tsv",

  # Full association tables for QQ plots and lambda.
  # Each file should contain at least p-value, bin/group, cell and gene columns.
  # Multiple files can be listed. Set character(0) to skip.
  ie_full_files = "eQTL_res/02_acat/ieQTL/01_raw/all_ieQTLs.tsv",
  cse_full_files = "eQTL_res/02_acat/cseQTL/01_raw/all_cseQTLs.tsv",

  # Optional permutation association tables for calibration.
  ie_perm_files = character(0),
  cse_perm_files = character(0),

  # Parsed, merged coloc tables.
  coloc_bulk_ie = "eQTL_res/03_coloc/01_beQTL_ieQTL_coloc/merged_be_ie_coloc.tsv",
  coloc_bulk_cse = "eQTL_res/03_coloc/02_beQTL_cseQTL_coloc/merged_be_cse_coloc.tsv",
  coloc_ie_cse = "eQTL_res/03_coloc/03_ieQTL_cseQTL_coloc/merged_ie_cse_coloc.tsv",

  # Column mapping
  hit_gene_col = "pheno_id",
  hit_variant_col = "variant_id",
  hit_bin_col = "group",
  hit_cell_col = "cell",

  # OmiGA term-specific columns:
  # bulk eQTL and cseQTL use the genotype main effect (g1);
  # ieQTL uses the genotype × cell-fraction interaction effect (g2).
  # Raw SNP-level P values are used for QQ/calibration and effect reporting.
  # q-values are retained as gene-level discovery statistics, but are not
  # substituted for raw P values in SNP-level QQ plots.
  bulk_cols = list(
    beta = "beta_g1",
    se = "beta_se_g1",
    p = "pval_g1",
    q = "qval_g1"
  ),
  cse_cols = list(
    beta = "beta_g1",
    se = "beta_se_g1",
    p = "pval_g1",
    q = "qval_g1"
  ),
  ie_cols = list(
    beta = "beta_g2",
    se = "beta_se_g2",
    p = "pval_g2",
    q = "qval_g2"
  ),

  # Cell-fraction filters
  min_mean_fraction = 0.05,
  min_nonzero_fraction = 0.01,
  min_nonzero_n = 30,
  min_fraction_sd = 0.02,

  # Significance / coloc thresholds
  mash_lfsr = 0.05,
  heterogeneity_fdr = 0.05,
  coloc_pp4_thresholds = c(0.75, 0.90),

  # Shared bulk-eGene definition
  min_bins_shared = 2,

  seed = 123
)

dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(cfg$seed)

# 1. Utility functions

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

file_exists_all <- function(x) {
  length(x) > 0 && all(file.exists(x))
}

read_any <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") return(readRDS(path))
  as.data.frame(data.table::fread(path), stringsAsFactors = FALSE)
}

save_plot <- function(p, stem, width = 7, height = 5) {
  pdf_file <- file.path(cfg$outdir, paste0(stem, ".pdf"))
  png_file <- file.path(cfg$outdir, paste0(stem, ".png"))
  ggsave(pdf_file, p, width = width, height = height, useDingbats = FALSE)
  ggsave(png_file, p, width = width, height = height, dpi = 300)
}

safe_cor_test <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(list(estimate = NA_real_, p.value = NA_real_, n = sum(ok)))
  out <- suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE))
  list(estimate = unname(out$estimate), p.value = out$p.value, n = sum(ok))
}

std_names <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

inverse_normal_transform <- function(x) {
  r <- rank(x, ties.method = "average", na.last = "keep")
  p <- (r - 0.5) / sum(!is.na(r))
  qnorm(p)
}

bin_cols <- c("bin3"="#35B779", "bin1"="#440154", "bin2"="#31688E", "bin4"="#FDE725")

# 2. Sample metadata and deconvolution QC

sample_info <- as.data.frame(data.table::fread(cfg$sample_info), stringsAsFactors = FALSE)
stopifnot(cfg$sample_col %in% colnames(sample_info))
stopifnot(cfg$bin_col %in% colnames(sample_info))

names(sample_info)[names(sample_info) == cfg$sample_col] <- "sample"
names(sample_info)[names(sample_info) == cfg$bin_col] <- "bin"

bin_raw <- as.character(sample_info$bin)
unexpected_bins <- setdiff(unique(bin_raw), cfg$bins)
if (length(unexpected_bins) > 0) {
  stop(
    "Unexpected values in the bin column: ",
    paste(unexpected_bins, collapse = ", "),
    ". Expected: ", paste(cfg$bins, collapse = ", ")
  )
}
sample_info$bin <- factor(bin_raw, levels = cfg$bins)

cell_cols <- setdiff(colnames(sample_info), c("sample", "bin"))
cell_is_numeric <- vapply(
  sample_info[cell_cols],
  is.numeric,
  logical(1)
)
cell_cols <- cell_cols[cell_is_numeric]

if (length(cell_cols) == 0) {
  stop("No numeric cell-fraction columns found in sample_info.")
}

frac_long <- sample_info %>%
  pivot_longer(cols = all_of(cell_cols), names_to = "cell", values_to = "fraction")

frac_summary <- frac_long %>%
  group_by(bin, cell) %>%
  summarise(
    n_sample = n(),
    mean_fraction = mean(fraction, na.rm = TRUE),
    median_fraction = median(fraction, na.rm = TRUE),
    sd_fraction = sd(fraction, na.rm = TRUE),
    nonzero_n = sum(fraction > cfg$min_nonzero_fraction, na.rm = TRUE),
    n_eff_sum = sum(fraction, na.rm = TRUE),
    n_x_mean = n_sample * mean_fraction,
    n_x_var = n_sample * var(fraction, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    eligible = mean_fraction >= cfg$min_mean_fraction &
      nonzero_n >= cfg$min_nonzero_n &
      sd_fraction >= cfg$min_fraction_sd
  )

data.table::fwrite(frac_summary, file.path(cfg$outdir, "cell_fraction_summary.tsv"), sep = "\t")

# 2A. Fraction distributions by developmental bin
p_frac_box <- ggplot(frac_long, aes(x = bin, y = fraction, color = bin)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  facet_wrap(~cell, scales = "free_y", ncol = 4) +
  scale_x_discrete(labels = cfg$bin_labels) +
  scale_color_manual(values = bin_cols) +
  labs(x = NULL, y = "Estimated cell fraction") +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.background = element_blank()
  )
p_frac_box
save_plot(p_frac_box, "S1_cell_fraction_by_bin", width = 6, height = 4.7)

# 2B. Heatmap-style eligibility plot
p_elig <- ggplot(frac_summary, aes(x = bin, y = cell, fill = mean_fraction)) +
  geom_tile() +
  geom_point(
    aes(shape = eligible),
    size = 2.2,
    fill = "white"
  ) +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 4)) +
  scale_x_discrete(labels = cfg$bin_labels) +
  scale_fill_viridis_c(option = "E") +
  labs(x = NULL, y = NULL, fill = "Mean fraction", shape = "Eligible") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
p_elig
save_plot(p_elig, "S2_cell_fraction_eligibility", width = 4.4, height = 5.7)

# 2C. Bisque vs CIBERSORT concordance
plot_method_concordance <- function(bisque_path, cibersort_path) {
  if (!file_exists_all(c(bisque_path, cibersort_path))) return(NULL)

  bisque <- as.data.frame(data.table::fread(bisque_path), stringsAsFactors = FALSE)
  if (!("sample" %in% std_names(colnames(bisque)))) {
    first <- colnames(bisque)[1]
    names(bisque)[names(bisque) == first] <- "sample"
  } else {
    names(bisque) <- std_names(names(bisque))
  }

  cib <- as.data.frame(data.table::fread(cibersort_path), stringsAsFactors = FALSE)
  if (!(cfg$cibersort_sample_col %in% colnames(cib))) {
    stop("CIBERSORT sample column not found: ", cfg$cibersort_sample_col)
  }
  names(cib)[names(cib) == cfg$cibersort_sample_col] <- "sample"
  names(cib) <- std_names(names(cib))

  # Standardize sample column name after all column-name cleaning
  if (!("sample" %in% colnames(bisque))) {
    names(bisque)[1] <- "sample"
  }
  if (!("sample" %in% colnames(cib))) stop("Could not standardize CIBERSORT sample column.")

  common_cells <- intersect(setdiff(colnames(bisque), "sample"), setdiff(colnames(cib), "sample"))
  if (length(common_cells) == 0) {
    message("No identically named cell types between Bisque and CIBERSORT; skipping concordance.")
    return(NULL)
  }

  b_long <- bisque %>%
    select(sample, all_of(common_cells)) %>%
    pivot_longer(-sample, names_to = "cell", values_to = "bisque")

  c_long <- cib %>%
    select(sample, all_of(common_cells)) %>%
    pivot_longer(-sample, names_to = "cell", values_to = "cibersort")

  cc <- inner_join(b_long, c_long, by = c("sample", "cell"))

  stat <- cc %>%
    group_by(cell) %>%
    summarise(
      rho = safe_cor_test(bisque, cibersort, "spearman")$estimate,
      p = safe_cor_test(bisque, cibersort, "spearman")$p.value,
      n = safe_cor_test(bisque, cibersort, "spearman")$n,
      .groups = "drop"
    )

  data.table::fwrite(stat, file.path(cfg$outdir, "deconvolution_method_concordance.tsv"), sep = "\t")

  p <- ggplot(cc, aes(x = bisque, y = cibersort)) +
    geom_point(alpha = 0.35, size = 0.7) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
    facet_wrap(~cell, scales = "free", ncol = 4) +
    labs(x = "Bisque fraction", y = "CIBERSORT fraction") +
    theme_classic(base_size = 10) +
    theme(strip.background = element_blank())

  save_plot(p, "S3_deconvolution_method_concordance", width = 10, height = 7)
  stat
}

try(plot_method_concordance(cfg$bisque_prop, cfg$cibersort_prop), silent = TRUE)

# 3. mashr sharing and formal bulk-eQTL heterogeneity

post_mean <- readRDS(cfg$mash_post_mean)
lfsr <- readRDS(cfg$mash_lfsr)
pair_annot <- as.data.frame(data.table::fread(cfg$mash_pair_annot), stringsAsFactors = FALSE)

stopifnot(identical(dim(post_mean), dim(lfsr)))
stopifnot(identical(rownames(post_mean), rownames(lfsr)))

bins_use <- intersect(cfg$bins, colnames(post_mean))
if (length(bins_use) < 2) stop("Fewer than two configured bins found in mashr matrices.")

sig_mat <- lfsr[, bins_use, drop = FALSE] < cfg$mash_lfsr

# 3A. Posterior-effect correlation
cor_mat <- cor(post_mean[, bins_use, drop = FALSE], use = "pairwise.complete.obs")

cor_df <- as.data.frame(as.table(cor_mat))
colnames(cor_df) <- c("bin1", "bin2", "correlation")
cor_df$bin1 <- factor(cor_df$bin1, levels = bins_use)
cor_df$bin2 <- factor(cor_df$bin2, levels = bins_use)

# 3B. Number of significant bins per gene
# Do not derive bin-specific versus bin-shared genes from the mashr strong set.
# The strong set is preselected for robust signals, and mashr borrowing across
# conditions can make nearly every included gene significant in all bins.
# Classify genes from the per-bin significant bulk-eQTL calls.

bulk_sig_for_sharing <- as.data.frame(
  data.table::fread(cfg$bulk_hits),
  stringsAsFactors = FALSE
)

sharing_required <- c(cfg$hit_gene_col, cfg$hit_bin_col)
sharing_missing <- setdiff(sharing_required, colnames(bulk_sig_for_sharing))

if (length(sharing_missing) > 0) {
  stop(
    "bulk_hits is missing columns required for sharing classification: ",
    paste(sharing_missing, collapse = ", ")
  )
}

gene_bin_sig <- bulk_sig_for_sharing %>%
  transmute(
    gene_id = .data[[cfg$hit_gene_col]],
    bin = as.character(.data[[cfg$hit_bin_col]])
  ) %>%
  filter(
    !is.na(gene_id),
    bin %in% bins_use
  ) %>%
  distinct(gene_id, bin)

gene_sharing <- gene_bin_sig %>%
  count(gene_id, name = "n_sig_bins") %>%
  mutate(
    class_gene = case_when(
      n_sig_bins == 1 ~ "bin-specific",
      n_sig_bins >= cfg$min_bins_shared ~ "bin-shared",
      TRUE ~ "none"
    )
  )

all_genes <- unique(gene_bin_sig$gene_id)

gene_sharing <- tibble(gene_id = all_genes) %>%
  left_join(gene_sharing, by = "gene_id") %>%
  mutate(
    n_sig_bins = replace_na(n_sig_bins, 0L),
    class_gene = replace_na(class_gene, "none")
  )

sharing_counts <- gene_sharing %>%
  count(class_gene, name = "n")

print(sharing_counts)

if (!all(c("bin-specific", "bin-shared") %in% sharing_counts$class_gene)) {
  warning(
    "The bulk-eQTL calls still do not contain both bin-specific and bin-shared genes. ",
    "Check that cfg$bulk_hits includes significant results from all four bins and ",
    "that the file was not prefiltered to cross-bin shared signals."
  )
}

data.table::fwrite(
  gene_sharing,
  file.path(cfg$outdir, "bulk_eGene_sharing_class.tsv"),
  sep = "\t"
)

# 3C. Cochran Q heterogeneity test from beta and SE
heterogeneity_test <- function(beta_mat, se_mat, pair_annot, bins) {
  stopifnot(identical(dim(beta_mat), dim(se_mat)))
  stopifnot(identical(rownames(beta_mat), rownames(se_mat)))

  b <- beta_mat[, bins, drop = FALSE]
  s <- se_mat[, bins, drop = FALSE]

  calc_one <- function(i) {
    bi <- as.numeric(b[i, ])
    sei <- as.numeric(s[i, ])
    ok <- is.finite(bi) & is.finite(sei) & sei > 0
    if (sum(ok) < 2) return(c(Q = NA, df = NA, p_het = NA, beta_range = NA))
    wi <- 1 / sei[ok]^2
    bbar <- sum(wi * bi[ok]) / sum(wi)
    Q <- sum(wi * (bi[ok] - bbar)^2)
    df <- sum(ok) - 1
    p <- pchisq(Q, df = df, lower.tail = FALSE)
    c(Q = Q, df = df, p_het = p, beta_range = diff(range(bi[ok])))
  }

  out <- t(vapply(seq_len(nrow(b)), calc_one,
                  FUN.VALUE = c(Q = 0, df = 0, p_het = 0, beta_range = 0)))
  out <- as.data.frame(out)
  out$pair_id <- rownames(b)
  out$fdr_het <- p.adjust(out$p_het, method = "BH")
  out <- left_join(out, pair_annot %>% distinct(pair_id, gene_id), by = "pair_id")
  out
}

if (file_exists_all(c(cfg$bulk_beta_mat, cfg$bulk_se_mat))) {
  beta_mat <- readRDS(cfg$bulk_beta_mat)
  se_mat <- readRDS(cfg$bulk_se_mat)

  het <- heterogeneity_test(beta_mat, se_mat, pair_annot, bins_use)
  data.table::fwrite(het, file.path(cfg$outdir, "bulk_eQTL_heterogeneity.tsv"), sep = "\t")

  p_het <- ggplot(het, aes(x = beta_range, y = -log10(p_het))) +
    geom_point(alpha = 0.35, size = 0.7) +
    geom_hline(yintercept = -log10(0.05 / max(1, nrow(het))),
               linetype = 2) +
    labs(x = "Range of estimated effects across bins",
         y = expression(-log[10](P[heterogeneity]))) +
    theme_classic(base_size = 11)
  save_plot(p_het, "M3_bulk_eQTL_effect_heterogeneity", width = 5.8, height = 4.8)

  # Select robust examples: significant heterogeneity and observed in >=3 bins
  candidate <- het %>%
    filter(fdr_het < cfg$heterogeneity_fdr) %>%
    arrange(fdr_het, desc(beta_range)) %>%
    slice_head(n = 20)

  data.table::fwrite(candidate, file.path(cfg$outdir, "candidate_stage_heterogeneous_bulk_eQTL.tsv"),
         sep = "\t")
} else {
  message("bulk_beta_mat or bulk_se_mat missing; formal heterogeneity test skipped.")
}

p_het

# 4. ieQTL/cseQTL power diagnostics

read_hits <- function(path, type, col_map) {
  x <- as.data.frame(data.table::fread(path), stringsAsFactors = FALSE)

  required <- c(
    cfg$hit_gene_col,
    cfg$hit_bin_col,
    cfg$hit_cell_col,
    col_map$beta,
    col_map$se,
    col_map$p
  )
  miss <- setdiff(required, colnames(x))
  if (length(miss)) {
    stop(type, " hit table missing columns: ", paste(miss, collapse = ", "))
  }

  x %>%
    transmute(
      gene_id = .data[[cfg$hit_gene_col]],
      variant_id = if (cfg$hit_variant_col %in% colnames(x)) {
        .data[[cfg$hit_variant_col]]
      } else {
        NA_character_
      },
      bin = .data[[cfg$hit_bin_col]],
      cell = .data[[cfg$hit_cell_col]],
      beta = .data[[col_map$beta]],
      se = .data[[col_map$se]],
      p = .data[[col_map$p]],
      q = if (!is.null(col_map$q) && col_map$q %in% colnames(x)) {
        .data[[col_map$q]]
      } else {
        NA_real_
      },
      type = type
    )
}

ie_hits <- read_hits(cfg$ie_hits, "ieQTL", cfg$ie_cols)
cse_hits <- read_hits(cfg$cse_hits, "cseQTL", cfg$cse_cols)


cell_name_map <- c(
  "adip"  = "adipocyte",
  "cap"   = "capillary_endothelial_cell",
  "fap"   = "fibro_adipogenic_progenitor_cell",
  "glial" = "peripheral_glial",
  "lec"   = "lymphatic_endothelial_cell",
  "mac"   = "macrophage",
  "mf1"   = "type_i_myonuclei",
  "mf2ab" = "type_ii_a_b_myonuclei",
  "mf2x"  = "type_ii_x_myonuclei",
  "msc"   = "muscle_stem_cell",
  "peri"  = "pericyte",
  "tc"    = "t_cell",
  "teno"  = "tenocyte"
)

standardize_cell_names <- function(x) {
  x <- as.character(x)
  
  mapped <- unname(cell_name_map[x])
  mapped[is.na(mapped)] <- x[is.na(mapped)]
  
  mapped
}

ie_hits <- ie_hits %>%
  mutate(cell = standardize_cell_names(cell))

cse_hits <- cse_hits %>%
  mutate(cell = standardize_cell_names(cell))


power_table <- function(hits, frac_summary) {
  eg <- hits %>%
    distinct(gene_id, bin, cell) %>%
    count(bin, cell, name = "eGene_count")

  eg %>%
    left_join(frac_summary, by = c("bin", "cell")) %>%
    mutate(log1p_eGene = log1p(eGene_count))
}

ie_power <- power_table(ie_hits, frac_summary)
cse_power <- power_table(cse_hits, frac_summary)

data.table::fwrite(ie_power, file.path(cfg$outdir, "ieQTL_power_metrics.tsv"), sep = "\t")
data.table::fwrite(cse_power, file.path(cfg$outdir, "cseQTL_power_metrics.tsv"), sep = "\t")

plot_power <- function(df, prefix) {
  metrics <- c(
    n_x_mean = "N × mean(cell fraction)",
    n_x_var = "N × variance(cell fraction)",
    sd_fraction = "SD of cell fraction",
    nonzero_n = "Samples with fraction > threshold"
  )

  long <- df %>%
    pivot_longer(cols = all_of(names(metrics)), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = names(metrics), labels = metrics))

  stat <- long %>%
    group_by(metric) %>%
    summarise(
      rho = safe_cor_test(value, eGene_count, "spearman")$estimate,
      p = safe_cor_test(value, eGene_count, "spearman")$p.value,
      n = safe_cor_test(value, eGene_count, "spearman")$n,
      .groups = "drop"
    )

  data.table::fwrite(stat, file.path(cfg$outdir, paste0(prefix, "_power_metric_correlations.tsv")), sep = "\t")

  lab <- stat %>%
    mutate(label = paste0("rho=", sprintf("%.2f", rho), "\nP=", format.pval(p, digits = 2)))

  p <- ggplot(long, aes(x = value, y = eGene_count)) +
    geom_point(aes(fill = bin), size = 2, alpha = 0.85, shape = 21, colour = "white") +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.6) +
    facet_wrap(~metric, scales = "free_x", ncol = 2) +
    scale_fill_manual(values = bin_cols) +
    geom_text(
      data = lab,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.05, vjust = 1.1, size = 3.4
    ) +
    labs(x = NULL, y = "Number of eGenes") +
    theme_classic(base_size = 10) +
    theme(strip.background = element_blank(), aspect.ratio = 1)

  save_plot(p, paste0(prefix, "_power_diagnostics"), width = 5.5, height = 4)
}

plot_power(ie_power, "M4_ieQTL")
plot_power(cse_power, "M5_cseQTL")

# 5. QQ plots, lambda and optional permutation calibration

read_assoc_files <- function(files, type) {
  if (length(files) == 0) return(NULL)
  bind_rows(lapply(files, function(f) {
    x <- as.data.frame(data.table::fread(f), stringsAsFactors = FALSE)
    x$source_file <- basename(f)
    x$type <- type
    x
  }))
}

calc_lambda <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (length(p) < 10) return(NA_real_)
  chisq <- qchisq(1 - p, df = 1)
  median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
}

qq_table <- function(df, p_col, group_cols) {
  df %>%
    filter(is.finite(.data[[p_col]]), .data[[p_col]] > 0, .data[[p_col]] <= 1) %>%
    group_by(across(all_of(group_cols))) %>%
    arrange(.data[[p_col]], .by_group = TRUE) %>%
    mutate(
      rank = row_number(),
      n = n(),
      expected = -log10((rank - 0.5) / n),
      observed = -log10(.data[[p_col]])
    ) %>%
    ungroup()
}

plot_qq_and_lambda <- function(files, type, p_col) {
  x <- read_assoc_files(files, type)
  if (is.null(x)) return(NULL)

  required <- c(p_col, cfg$hit_bin_col, cfg$hit_cell_col)
  miss <- setdiff(required, colnames(x))
  if (length(miss)) {
    message(type, " full association files missing columns: ", paste(miss, collapse = ", "))
    return(NULL)
  }

  x <- x %>%
    rename(
      bin = all_of(cfg$hit_bin_col),
      cell = all_of(cfg$hit_cell_col)
    )

  lam <- x %>%
    group_by(bin, cell) %>%
    summarise(
      lambda_gc = calc_lambda(.data[[p_col]]),
      n_test = sum(is.finite(.data[[p_col]])),
      .groups = "drop"
    )

  data.table::fwrite(lam, file.path(cfg$outdir, paste0(type, "_lambda_gc.tsv")), sep = "\t")

  qq <- qq_table(x, p_col, c("bin", "cell"))

  p <- ggplot(qq, aes(expected, observed)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point(alpha = 0.25, size = 0.4) +
    facet_grid(bin ~ cell, scales = "free") +
    labs(x = expression(Expected~~-log[10](P)),
         y = expression(Observed~~-log[10](P))) +
    theme_classic(base_size = 8) +
    theme(strip.background = element_blank())

  save_plot(p, paste0(type, "_QQ_plots"), width = 12, height = 8)
  lam
}

# Use raw SNP-level P values for QQ plots and lambdaGC.
# qval_g1/qval_g2 are gene-level multiple-testing statistics and should not
# replace pval_g1/pval_g2 in these SNP-level calibration plots.
plot_qq_and_lambda(cfg$ie_full_files, "ieQTL", cfg$ie_cols$p)
plot_qq_and_lambda(cfg$cse_full_files, "cseQTL", cfg$cse_cols$p)

plot_observed_vs_permuted <- function(obs_files, perm_files, type, p_col) {
  if (length(obs_files) == 0 || length(perm_files) == 0) return(NULL)
  obs <- read_assoc_files(obs_files, "Observed")
  prm <- read_assoc_files(perm_files, "Permuted")

  required <- c(p_col, cfg$hit_bin_col, cfg$hit_cell_col)
  if (!all(required %in% colnames(obs)) || !all(required %in% colnames(prm))) return(NULL)

  x <- bind_rows(obs, prm) %>%
    rename(bin = all_of(cfg$hit_bin_col), cell = all_of(cfg$hit_cell_col))

  qq <- qq_table(x, p_col, c("type", "bin", "cell"))

  p <- ggplot(qq, aes(expected, observed, colour = type)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_line(linewidth = 0.5) +
    facet_grid(bin ~ cell, scales = "free") +
    labs(x = expression(Expected~~-log[10](P)),
         y = expression(Observed~~-log[10](P)),
         colour = NULL) +
    theme_classic(base_size = 8) +
    theme(strip.background = element_blank())

  save_plot(p, paste0(type, "_observed_vs_permuted_QQ"), width = 12, height = 8)
}

plot_observed_vs_permuted(
  cfg$ie_full_files, cfg$ie_perm_files, "ieQTL", cfg$ie_cols$p
)
plot_observed_vs_permuted(
  cfg$cse_full_files, cfg$cse_perm_files, "cseQTL", cfg$cse_cols$p
)

# 6. ieQTL vs cseQTL overlap and effect concordance

ie_key <- ie_hits %>%
  mutate(key = paste(gene_id, variant_id, bin, cell, sep = "|")) %>%
  distinct(key, .keep_all = TRUE)

cse_key <- cse_hits %>%
  mutate(key = paste(gene_id, variant_id, bin, cell, sep = "|")) %>%
  distinct(key, .keep_all = TRUE)

joint <- inner_join(
  ie_key %>% select(key, gene_id, variant_id, bin, cell, beta_ie = beta, p_ie = p),
  cse_key %>% select(key, beta_cse = beta, p_cse = p),
  by = "key"
)

if (nrow(joint) > 0) {
  joint <- joint %>%
    mutate(direction_concordant = sign(beta_ie) == sign(beta_cse))

  conc_summary <- joint %>%
    group_by(bin, cell) %>%
    summarise(
      n_overlap = n(),
      concordance = mean(direction_concordant, na.rm = TRUE),
      rho = safe_cor_test(beta_ie, beta_cse, "spearman")$estimate,
      p = safe_cor_test(beta_ie, beta_cse, "spearman")$p.value,
      .groups = "drop"
    )

  data.table::fwrite(joint, file.path(cfg$outdir, "ieQTL_cseQTL_matched_effects.tsv"), sep = "\t")
  data.table::fwrite(conc_summary, file.path(cfg$outdir, "ieQTL_cseQTL_effect_concordance.tsv"), sep = "\t")

  p_joint <- ggplot(joint, aes(beta_ie, beta_cse)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_point(alpha = 0.35, size = 0.7) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
    facet_wrap(~cell, scales = "free", ncol = 4) +
    labs(x = "ieQTL effect", y = "cseQTL effect") +
    theme_classic(base_size = 9) +
    theme(strip.background = element_blank())

  save_plot(p_joint, "M6_ieQTL_cseQTL_effect_concordance", width = 10, height = 7)
} else {
  message("No exact gene-variant-bin-cell matches between ieQTL and cseQTL tables.")
}

# Gene-level overlap by bin and cell
overlap_gene <- full_join(
  ie_hits %>% distinct(gene_id, bin, cell) %>% mutate(ie = TRUE),
  cse_hits %>% distinct(gene_id, bin, cell) %>% mutate(cse = TRUE),
  by = c("gene_id", "bin", "cell")
) %>%
  mutate(
    ie = replace_na(ie, FALSE),
    cse = replace_na(cse, FALSE),
    class = case_when(
      ie & cse ~ "Shared",
      ie ~ "ieQTL only",
      cse ~ "cseQTL only"
    )
  )

overlap_summary <- overlap_gene %>%
  count(bin, cell, class, name = "n_gene")

# 7. Colocalization summary using matched tested loci

standardize_coloc <- function(path, comparison) {
  if (!file.exists(path)) return(NULL)
  x <- as.data.frame(data.table::fread(path), stringsAsFactors = FALSE)

  pp4_candidates <- c("PP.H4.abf", "PP_H4_abf", "pp4", "PP4")
  pp3_candidates <- c("PP.H3.abf", "PP_H3_abf", "pp3", "PP3")
  pp4_col <- intersect(pp4_candidates, colnames(x))[1]
  pp3_col <- intersect(pp3_candidates, colnames(x))[1]

  if (is.na(pp4_col)) stop("No PP4 column found in ", path)

  x %>%
    transmute(
      gene_id = if ("gene_id" %in% colnames(x)) gene_id else NA_character_,
      bin = if ("bin" %in% colnames(x)) bin else NA_character_,
      cell = if ("cell" %in% colnames(x)) cell else NA_character_,
      pp4 = .data[[pp4_col]],
      pp3 = if (!is.na(pp3_col)) .data[[pp3_col]] else NA_real_,
      comparison = comparison
    )
}

coloc_all <- bind_rows(
  standardize_coloc(cfg$coloc_bulk_ie, "Bulk eQTL vs ieQTL"),
  standardize_coloc(cfg$coloc_bulk_cse, "Bulk eQTL vs cseQTL"),
  standardize_coloc(cfg$coloc_ie_cse, "ieQTL vs cseQTL")
)

if (nrow(coloc_all) > 0) {
  coloc_summary <- coloc_all %>%
    group_by(comparison) %>%
    summarise(
      n_tested = n(),
      median_pp4 = median(pp4, na.rm = TRUE),
      mean_pp4 = mean(pp4, na.rm = TRUE),
      n_pp4_075 = sum(pp4 >= 0.75, na.rm = TRUE),
      prop_pp4_075 = mean(pp4 >= 0.75, na.rm = TRUE),
      n_pp4_090 = sum(pp4 >= 0.90, na.rm = TRUE),
      prop_pp4_090 = mean(pp4 >= 0.90, na.rm = TRUE),
      prop_pp3_gt_pp4 = mean(pp3 > pp4, na.rm = TRUE),
      .groups = "drop"
    )

  data.table::fwrite(coloc_summary, file.path(cfg$outdir, "coloc_summary.tsv"), sep = "\t")

  p_coloc_density <- ggplot(coloc_all, aes(x = pp4, fill = comparison)) +
    geom_density(alpha = 0.35) +
    geom_vline(xintercept = cfg$coloc_pp4_thresholds, linetype = 2) +
    facet_wrap(~comparison, ncol = 1) +
    labs(x = "Posterior probability of a shared causal variant (PP4)",
         y = "Density",
         fill = NULL) +
    theme_classic(base_size = 10) +
    theme(legend.position = "none")
  save_plot(p_coloc_density, "M8_coloc_PP4_density", width = 6.5, height = 7)

  coloc_thr <- coloc_all %>%
    group_by(comparison) %>%
    summarise(
      `PP4 ≥ 0.75` = mean(pp4 >= 0.75, na.rm = TRUE),
      `PP4 ≥ 0.90` = mean(pp4 >= 0.90, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(-comparison, names_to = "threshold", values_to = "proportion")

  p_coloc_bar <- ggplot(coloc_thr, aes(x = comparison, y = proportion, fill = threshold)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    scale_y_continuous(labels = percent_format()) +
    labs(x = NULL, y = "Proportion of tested loci", fill = NULL) +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  save_plot(p_coloc_bar, "M9_coloc_threshold_summary", width = 7, height = 4.8)
} else {
  message("No merged coloc tables found; coloc summary skipped.")
}


# 8. Figure index

index <- tibble::tribble(
  ~file_stem, ~purpose,
  "S1_cell_fraction_by_bin", "Distribution of inferred cell fractions across developmental bins",
  "S2_cell_fraction_eligibility", "Eligibility of each bin × cell-type combination",
  "S3_deconvolution_method_concordance", "Bisque versus CIBERSORT concordance",
  "M1_mash_posterior_effect_correlation", "Global sharing of bulk-eQTL effects across bins",
  "M2_bulk_eGene_number_of_shared_bins", "Number of developmental bins in which each bulk eGene is significant",
  "M3_bulk_eQTL_effect_heterogeneity", "Formal Cochran-Q test of effect heterogeneity",
  "M4_ieQTL_power_diagnostics", "Relationship between ieGene counts and cell-information metrics",
  "M5_cseQTL_power_diagnostics", "Relationship between cseGene counts and cell-information metrics",
  "M6_ieQTL_cseQTL_effect_concordance", "Matched ieQTL–cseQTL effect concordance",
  "M7_ieQTL_cseQTL_gene_overlap", "Gene-level overlap of ieQTL and cseQTL discoveries",
  "M8_coloc_PP4_density", "Full PP4 distributions",
  "M9_coloc_threshold_summary", "Proportion of tested loci with PP4 above prespecified thresholds",
)

data.table::fwrite(index, file.path(cfg$outdir, "figure_index.tsv"), sep = "\t")

message("Finished. Outputs written to: ", normalizePath(cfg$outdir, mustWork = FALSE))
