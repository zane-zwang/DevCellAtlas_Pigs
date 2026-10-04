#!/usr/bin/env bash
# Purpose: map interaction bins1 4.
ANALYSIS_DIR="${DEVPIGATLAS_ANALYSIS_DIR:-path/to/input}"

cd ${ANALYSIS_DIR}/qtl/mapping/bisque_omiga/01_acat

for i in bin1 bin4
  do
    for c in adip cap fap lec mac msc peri glial tc teno mf1 mf2ab mf2x
      do
        mkdir ${i}
        mkdir ${i}/${c}
        omiga --mode cis_interaction \
        --genotype ${ANALYSIS_DIR}/qtl/input/${i}/genotype/Muscle_subset \
        --phenotype ${ANALYSIS_DIR}/qtl/input/${i}/invTMM/muscle_${i}.expression.bed.gz \
        --covariates ${ANALYSIS_DIR}/qtl/input/${i}/covFile/muscle_${i}_cov.tsv \
        --interaction ${ANALYSIS_DIR}/qtl/input/bisque/${c}_freq.txt \
        --verbose  --dprop-pc-covar 0.001 \
        --rm-collinear-covar 0.95 \
        --calcu-variant-threshold \
        --prefix muscle_${i}_${c} \
        --output-dir ${ANALYSIS_DIR}/qtl/mapping/bisque_omiga/01_acat/${i}/${c}
      done
  done
