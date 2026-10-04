#!/usr/bin/env bash
# Purpose: map cse bin2.
ANALYSIS_DIR="${DEVPIGATLAS_ANALYSIS_DIR:-path/to/input}"

cd ${ANALYSIS_DIR}/qtl/mapping/cse_omiga/01_acat

for i in bin2
  do
    for c in adipocyte fibro_adipogenic_progenitor_cell lymphatic_endothelial_cell muscle_stem_cell peripheral_glial tenocyte type_ii_x_myonuclei capillary_endothelial_cell macrophage pericyte t_cell type_ii_a_b_myonuclei type_i_myonuclei
      do
        mkdir ${i}
        mkdir ${i}/${c}
        omiga --mode cis \
        --genotype ${ANALYSIS_DIR}/qtl/input/${i}/genotype/Muscle_subset \
        --phenotype ${ANALYSIS_DIR}/qtl/mapping/cse_omiga/00_data/exp_for_cseqtl/${c}/${i}/expr_tmm_inv.bed \
        --covariates ${ANALYSIS_DIR}/qtl/mapping/cse_omiga/00_data/muscle_752_cov.tsv \
        --verbose  --dprop-pc-covar 0.001 \
        --rm-collinear-covar 0.95 \
        --calcu-variant-threshold \
        --prefix muscle_${i}_${c} \
        --output-dir ${ANALYSIS_DIR}/qtl/mapping/cse_omiga/01_acat/${i}/${c}
      done
  done
