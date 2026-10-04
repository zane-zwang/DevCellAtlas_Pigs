#!/usr/bin/env bash
# Purpose: map bulk cis.
ANALYSIS_DIR="${DEVPIGATLAS_ANALYSIS_DIR:-path/to/input}"

cd ${ANALYSIS_DIR}/qtl/mapping/bulk_omiga/01_acat

for i in bin1 bin2 bin3 bin4
  do
    mkdir ${i}
    omiga --mode cis \
          --genotype ${ANALYSIS_DIR}/qtl/input/${i}/genotype/Muscle_subset \
          --phenotype ${ANALYSIS_DIR}/qtl/input/${i}/invTMM/muscle_${i}.expression.bed.gz \
          --covariates ${ANALYSIS_DIR}/qtl/input/${i}/covFile/muscle_${i}_cov.tsv \
          --prefix muscle_${i} \
          --verbose  --dprop-pc-covar 0.001 \
          --rm-collinear-covar 0.95 \
          --calcu-variant-threshold \
          --output-dir ${ANALYSIS_DIR}/qtl/mapping/bulk_omiga/01_acat/${i}
  done
