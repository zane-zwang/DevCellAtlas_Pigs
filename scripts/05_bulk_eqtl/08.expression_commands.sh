#!/usr/bin/env bash
# Purpose: expression commands.
GTEX_PIPELINE_DIR="${GTEX_PIPELINE_DIR:-path/to/gtex-pipeline}"

for i in binAll
  do
    mkdir ${i}/invTMM
    python "${GTEX_PIPELINE_DIR}/qtl/src/eqtl_prepare_expression_tss_gene.py" \
    ${i}/exp_process/muscle.tpm.gct \
    ${i}/exp_process/muscle.counts.gct \
    Sus_scrofa.Sscrofa11.1.100.gtf \
    sample_to_participant.tsv \
    ${i}/invTMM/muscle_${i} \
    --chrs chrAuto.list \
    --tpm_threshold 0.1 --count_threshold 6 --sample_frac_threshold 0.2
  done
