#!/usr/bin/env bash
# Purpose: build taxonomy lineages.
ANALYSIS_DIR="${DEVPIGATLAS_ANALYSIS_DIR:-path/to/input}"


echo "Start......"
ncbitax2lin --nodes-file ${ANALYSIS_DIR}/tai/reference/taxdump/nodes.dmp \
            --names-file ${ANALYSIS_DIR}/tai/reference/taxdump/names.dmp \
            --output ${ANALYSIS_DIR}/tai/reference/ncbi_lineages_2025-05-20.csv.gz

date
echo "Done......"
