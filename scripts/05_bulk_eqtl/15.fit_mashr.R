# Purpose: fit mashr.
analysis_dir <- Sys.getenv("DEVPIGATLAS_ANALYSIS_DIR", "path/to/input")

## mashR
library(mashr)
library(data.table)

in_dir  <- file.path(analysis_dir, "qtl/mashr/01_beqtl_mashr_input")
out_dir <- file.path(analysis_dir, "qtl/mashr/01_beqtl_mashr_output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

bins <- c("bin1","bin2","bin3","bin4")
vhat_sub_n <- 20000


strong <- fread(file.path(in_dir, "strong_set_zval.txt"))
random <- fread(file.path(in_dir, "random_set_zval.txt"))

stopifnot(all(bins %in% names(strong)),
          all(bins %in% names(random)))

S <- as.matrix(strong[, bins, with=FALSE])
R <- as.matrix(random[, bins, with=FALSE])
rownames(S) <- strong$pair_id
rownames(R) <- random$pair_id

cat(sprintf("Strong: %d pairs x %d bins\n", nrow(S), ncol(S)))
cat(sprintf("Random: %d pairs x %d bins\n", nrow(R), ncol(R)))

drop_cols <- which(colSums(!is.na(S)) == 0 | colSums(!is.na(R)) == 0)
if (length(drop_cols) > 0) {
  stop(sprintf("These cols are all NA in strong or random: %s", paste(bins[drop_cols], collapse=",")))
}

set.seed(9823)
sub_n <- min(vhat_sub_n, nrow(R))
sub_idx <- sample.int(nrow(R), sub_n)
R_sub <- R[sub_idx, , drop=FALSE]

data.temp <- mash_set_data(Bhat = R_sub, alpha = 1, zero_Bhat_Shat_reset = 1e6)
Vhat <- estimate_null_correlation_simple(data.temp)
cat("Estimated Vhat from random subset.\n")

data.random <- mash_set_data(Bhat = R, alpha = 1, V = Vhat, zero_Bhat_Shat_reset = 1e6)
data.strong <- mash_set_data(Bhat = S, alpha = 1, V = Vhat, zero_Bhat_Shat_reset = 1e6)

U.c <- cov_canonical(data.random)
U.pca <- cov_pca(data.strong, 4)
U.ed <- cov_ed(data.strong, U.pca)
saveRDS(U.ed, file.path(out_dir, "U.ed.rds"))

cat("Fitting mixture on random set...\n")
m.r <- mash(data.random, Ulist = c(U.ed, U.c), outputlevel = 1)

cat("Computing posterior on strong set with fixed g...\n")
m.s <- mash(data.strong, g = get_fitted_g(m.r), fixg = TRUE, outputlevel = 3)

saveRDS(m.r, file.path(out_dir, "mash_random_fit.rds"))
saveRDS(m.s, file.path(out_dir, "mash_strong_posterior.rds"))

saveRDS(get_pm(m.s), file.path(out_dir, "posterior_means.rds"))
saveRDS(get_lfsr(m.s), file.path(out_dir, "lfsr.rds"))
saveRDS(get_estimated_pi(m.s), file.path(out_dir, "mixture_weights_pi.rds"))
saveRDS(get_pairwise_sharing(m.s), file.path(out_dir, "pairwise_sharing_lfsr05.rds"))
saveRDS(get_pairwise_sharing(m.s, factor=0), file.path(out_dir, "pairwise_sharing_any_effect.rds"))
