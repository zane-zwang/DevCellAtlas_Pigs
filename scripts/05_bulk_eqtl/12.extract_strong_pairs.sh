#!/usr/bin/env bash
# Purpose: extract strong pairs.
ANALYSIS_DIR="${DEVPIGATLAS_ANALYSIS_DIR:-path/to/input}"


eqtl_dir="${ANALYSIS_DIR}/qtl"
mashr_dir="${ANALYSIS_DIR}/qtl/mashr"
scripts_dir="${CATTLECELL_PIPELINE_DIR:-path/to/CattleCell-GTEx-Pipeline-v0-main}"
bulk_filter_dir="${QTL_BULK_FILTER_OUTPUT_DIR:-${eqtl_dir}/mapping/bulk_omiga/01_acat/beQTL/02_filter}"
mkdir -p ${mashr_dir}/01_beqtl_mashr_input

# step 1 combine_signif_pairs_tjy.py
find "${bulk_filter_dir}" -name 'bin*_significant_beQTLs.tsv' | sort > "${mashr_dir}/01_beqtl_mashr_input/signifpair_list.txt"
python "${scripts_dir}/6.Tissue and Celltype sharing/combine_signif_pairs_tjy.py" "${mashr_dir}/01_beqtl_mashr_input/signifpair_list.txt" nominal_pairs -o "${mashr_dir}/01_beqtl_mashr_input"
# nominal_pairs.combined_signifpairs.txt.gz

for bin in bin1 bin2 bin3 bin4
  do
    ls ${eqtl_dir}/mapping/bulk_omiga/01_acat/${bin}/muscle_${bin}.cis_qtl_pairs.*.txt.gz > "${mashr_dir}/01_beqtl_mashr_input/${bin}.nominal_files.txt"
    python "${scripts_dir}/6.Tissue and Celltype sharing/extract_pairs_tjy.py" \
           ${mashr_dir}/01_beqtl_mashr_input/${bin}.nominal_files.txt \
           ${mashr_dir}/01_beqtl_mashr_input/nominal_pairs.combined_signifpairs.txt.gz \
           ${bin}_nominal_pairs \
           -o ${mashr_dir}/01_beqtl_mashr_input
  done
