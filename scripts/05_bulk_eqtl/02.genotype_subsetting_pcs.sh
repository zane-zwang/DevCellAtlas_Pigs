#!/usr/bin/env bash
# Purpose: genotype subsetting pcs.
ANALYSIS_DIR="${DEVPIGATLAS_ANALYSIS_DIR:-path/to/input}"

cd ${ANALYSIS_DIR}/qtl/input
for i in binAll
  do
    mkdir ${i}/genotype
    mkdir ${i}/geno_pca

    plink --bfile ${ANALYSIS_DIR}/qtl/input/pigGTEx/PigGTEx_v0.ALL_Tissues_Genotype/Muscle --keep ${i}/sample_list_plink.txt --make-bed --out ${i}/genotype/Muscle_subset
    plink --bfile ${i}/genotype/Muscle_subset --pca 30 --out ${i}/geno_pca/genotypePC_30
  done
