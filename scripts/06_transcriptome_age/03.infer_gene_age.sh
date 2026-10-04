#!/usr/bin/env bash
# Purpose: infer gene age.
ANALYSIS_DIR="${DEVPIGATLAS_ANALYSIS_DIR:-path/to/input}"


genEra -q Sus_scrofa.Sscrofa11.1.pep.all.fa -t 9823 -b ${ANALYSIS_DIR}/tai/reference/nr_db/nr -d ${ANALYSIS_DIR}/tai/reference/taxdump -r ${ANALYSIS_DIR}/tai/reference/ncbi_lineages_2025-05-20.csv -n 80 -i true
