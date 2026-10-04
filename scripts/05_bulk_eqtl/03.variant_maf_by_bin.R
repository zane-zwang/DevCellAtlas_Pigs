# Purpose: variant maf by bin.

library(data.table)
library(dplyr)

bins <- c("bin1", "bin2", "bin3", "bin4")

variant_maf_by_bin <- lapply(bins, function(b) {
  x <- fread(
    paste0("./", b, ".frq")
  )

  x %>%
    transmute(
      variant_id = SNP,
      bin = b,
      allele1 = A1,
      allele2 = A2,
      maf = MAF,
      n_genotyped = NCHROBS / 2
    )
}) %>%
  bind_rows()

fwrite(
  variant_maf_by_bin,
  "variant_maf_by_bin.tsv",
  sep = "\t"
)
