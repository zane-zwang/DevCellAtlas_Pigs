#!/usr/bin/env bash
# Purpose: calculate ld.
ANALYSIS_DIR="${DEVPIGATLAS_ANALYSIS_DIR:-path/to/input}"

for i in bin1 bin2 bin3 bin4
  do
    plink --bfile ${ANALYSIS_DIR}/qtl/input/${i}/genotype/Muscle_subset \
    --r2 --ld-snp-list snp_list.txt \
    --ld-window-kb 1000 \
    --ld-window 99999 \
    --ld-window-r2 0 \
    --out ${i}_ld_r2_result
  done
