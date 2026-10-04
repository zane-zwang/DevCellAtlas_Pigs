#!/usr/bin/env bash
# Purpose: run pyscenic.
PUBLICATION_ROOT="${DEVPIGATLAS_ROOT:-$PWD}"
REFERENCE_DIR="${SCENIC_REFERENCE_DIR:-path/to/reference/cisTarget_databases_pig}"
SCENIC_OUTPUT_DIR="${SCENIC_OUTPUT_DIR:-${PUBLICATION_ROOT}/results/scenic}"
export SCENIC_OUTPUT_DIR
mkdir -p "${SCENIC_OUTPUT_DIR}"
export SCENIC_COUNTS_CSV="${SCENIC_OUTPUT_DIR}/RNAcounts.csv"
export SCENIC_COUNTS_LOOM="${SCENIC_OUTPUT_DIR}/RNAcounts.loom"
export SCENIC_RANKING_FEATHER="${SCENIC_RANKING_FEATHER:-${REFERENCE_DIR}/10KbUP_10KbDOWN.regions_vs_motifs.rankings.feather}"


Rscript "${PUBLICATION_ROOT}/scripts/02_developmental_dynamics/helpers/01.get_counts.R"
python "${PUBLICATION_ROOT}/scripts/02_developmental_dynamics/helpers/02.counts_loom.py"


tfs="${REFERENCE_DIR}/pig_tfs.lst"
feather="${SCENIC_RANKING_FEATHER}"
tbl="${REFERENCE_DIR}/motifs-v10nr_clust-nr.pig-m0.001-o0.0.tbl"

##grn
pyscenic grn "${SCENIC_COUNTS_LOOM}" \
--method grnboost2 \
--num_workers 20 \
-o "${SCENIC_OUTPUT_DIR}/adj.tsv" \
"${tfs}"

##cistarget
pyscenic ctx "${SCENIC_OUTPUT_DIR}/adj.tsv" \
"${feather}" \
--annotations_fname "${tbl}" \
--expression_mtx_fname "${SCENIC_COUNTS_LOOM}" \
--mode "dask_multiprocessing" \
--output "${SCENIC_OUTPUT_DIR}/reg.csv" \
--mask_dropouts \
--num_workers 20

##AUCell
pyscenic aucell "${SCENIC_COUNTS_LOOM}" \
"${SCENIC_OUTPUT_DIR}/reg.csv" \
--output "${SCENIC_OUTPUT_DIR}/out_SCENIC.loom" \
--num_workers 20
