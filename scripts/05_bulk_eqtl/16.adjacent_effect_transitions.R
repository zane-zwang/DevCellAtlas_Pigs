# Purpose: adjacent effect transitions.

library(dplyr)
library(tidyr)
library(ggplot2)

bulk_hits <- bulk_sig_for_sharing

bulk_sig <- bulk_hits %>%
  transmute(
    pair_id = paste(variant_id, pheno_id, sep = "|"),
    gene_id = pheno_id,
    bin = group
  ) %>%
  distinct(pair_id, gene_id, bin)

bulk_share <- bulk_sig %>%
  count(pair_id, gene_id, name = "n_sig_bins")

heterogeneity_test_simple <- function(beta_mat, se_mat) {
  
  calc_one <- function(i) {
    
    beta <- as.numeric(beta_mat[i, ])
    se   <- as.numeric(se_mat[i, ])
    
    ok <- is.finite(beta) & is.finite(se) & se > 0
    
    if (sum(ok) < 2) {
      return(c(Q = NA, p_het = NA, beta_range = NA))
    }
    
    w <- 1 / se[ok]^2
    beta_mean <- sum(w * beta[ok]) / sum(w)
    
    Q <- sum(w * (beta[ok] - beta_mean)^2)
    df <- sum(ok) - 1
    
    c(
      Q = Q,
      p_het = pchisq(Q, df, lower.tail = FALSE),
      beta_range = diff(range(beta[ok]))
    )
  }
  
  res <- t(vapply(
    seq_len(nrow(beta_mat)),
    calc_one,
    FUN.VALUE = c(Q = 0, p_het = 0, beta_range = 0)
  ))
  
  res <- as.data.frame(res)
  res$pair_id <- rownames(beta_mat)
  res$fdr_het <- p.adjust(res$p_het, method = "BH")
  
  res
}


bulk_class <- bulk_share %>%
  left_join(het, by = "pair_id") %>%
  mutate(
    class = case_when(
      n_sig_bins >= 2 & fdr_het >= 0.05 ~ "Stable",
      n_sig_bins >= 2 & fdr_het < 0.05  ~ "Developmentally modulated",
      n_sig_bins == 1 & fdr_het < 0.05  ~ "Stage-restricted candidate",
      TRUE ~ "Uncertain"
    )
  )


# --------- res2 ----------
adjacent_pairs <- list(
  c("bin1", "bin2"),
  c("bin2", "bin3"),
  c("bin3", "bin4")
)

compare_two_bins <- function(beta_mat, se_mat, bin_a, bin_b) {
  
  beta_a <- beta_mat[, bin_a]
  beta_b <- beta_mat[, bin_b]
  
  se_a <- se_mat[, bin_a]
  se_b <- se_mat[, bin_b]
  
  z_diff <- (beta_b - beta_a) / sqrt(se_a^2 + se_b^2)
  p_diff <- 2 * pnorm(-abs(z_diff))
  
  data.frame(
    pair_id = rownames(beta_mat),
    bin_a = bin_a,
    bin_b = bin_b,
    beta_a = beta_a,
    beta_b = beta_b,
    delta_beta = beta_b - beta_a,
    p_diff = p_diff,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      fdr_diff = p.adjust(p_diff, method = "BH")
    )
}

classify_transition <- function(df) {
  
  df %>%
    mutate(
      transition = case_when(
        fdr_diff < 0.05 &
          sign(beta_a) == sign(beta_b) &
          abs(beta_b) > abs(beta_a) ~ "Strengthened",
        
        fdr_diff < 0.05 &
          sign(beta_a) == sign(beta_b) &
          abs(beta_b) < abs(beta_a) ~ "Weakened",
        
        fdr_diff < 0.05 &
          sign(beta_a) != sign(beta_b) ~ "Direction switched",
        
        TRUE ~ "Stable"
      )
    )
}

library(purrr)

transition_res <- map_dfr(
  adjacent_pairs,
  function(x) {
    compare_two_bins(
      beta_mat,
      se_mat,
      bin_a = x[1],
      bin_b = x[2]
    ) %>%
      classify_transition()
  }
)

transition_summary <- transition_res %>%
  count(bin_a, bin_b, transition) %>%
  group_by(bin_a, bin_b) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

ggplot(
  transition_summary,
  aes(
    x = paste(bin_a, bin_b, sep = " → "),
    y = prop,
    fill = transition
  )
) +
  geom_col(width = 0.7, alpha=0.8) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#F5A623","#2EA67E","#E44B4B","#347EC7")) +
  labs(
    x = NULL,
    y = "Proportion of tested eQTLs",
    fill = NULL
  ) +
  theme_classic()

# ----------- res3 -----------
bulk_gene <- be_hits %>%
  transmute(
    gene_id = pheno_id,
    bin = group,
    bulk = TRUE
  ) %>%
  distinct()

ie_gene <- qtl_flt[qtl_flt$type=='ieQTL', ] %>%
  transmute(
    gene_id = eGene,
    bin = bin,
    ie = TRUE
  ) %>%
  distinct()

cse_gene <- qtl_flt[qtl_flt$type=='cseQTL', ] %>%
  transmute(
    gene_id = eGene,
    bin = bin,
    cse = TRUE
  ) %>%
  distinct()

regulatory_class <- full_join(
  bulk_gene,
  ie_gene,
  by = c("gene_id", "bin")
) %>%
  full_join(
    cse_gene,
    by = c("gene_id", "bin")
  ) %>%
  mutate(
    bulk = replace_na(bulk, FALSE),
    ie   = replace_na(ie, FALSE),
    cse  = replace_na(cse, FALSE)
  )

regulatory_class <- regulatory_class %>%
  mutate(
    class = case_when(
      bulk & !ie & !cse ~ "beQTL only",
      bulk & ie & !cse  ~ "beQTL + ieQTL",
      bulk & !ie & cse  ~ "beQTL + cseQTL",
      bulk & ie & cse   ~ "beQTL + ieQTL + cseQTL",
      !bulk & ie & !cse ~ "ieQTL only",
      !bulk & !ie & cse ~ "cseQTL only",
      !bulk & ie & cse  ~ "ieQTL + cseQTL"
    )
  )

class_summary <- regulatory_class %>%
  count(bin, class) %>%
  group_by(bin) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# ---------- res5 -----------
