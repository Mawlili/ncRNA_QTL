#!/bin/bash
#BSUB -W 0:10                           # Wall time
#BSUB -J gbat
#BSUB -cwd /rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/PRAD_gbat   # CHANGE ME
#BSUB -q short                          # Queue
#BSUB -o log/gbat_%J.out
#BSUB -e log/gbat_%J.error

#mkdir -p ./log ./results ./results/best_models

module load R/4.3.1

# 1. GReX
echo "Submitting predictive models."
grex_id=$(bsub < sub_grex.sh)
# Check if job was submitted successfully
if [[ -z "$grex_id" ]]; then
    echo "Failed to submit GReX job."
    exit 1
fi
# Extract job ID (remove angle brackets)
grex_id=$(echo $grex_id | awk '{print $2}' | tr -d '<>')
echo "GReX job ID: ${grex_id}"

# 2. Matrix eQTL
mqtl_id=$(bsub -q short \
               -J mqtl\
               -o log/mqtl_%J.out \
               -e log/mqtl_%J.error \
               -M 10G \
               -W 1:00 \
               -w "ended(${grex_id})" \
               "Rscript ./matrixeqtl.R")

# Check if job was submitted successfully
if [[ -z "$mqtl_id" ]]; then
    echo "Failed to submit Matrix eQTL job."
    exit 1
fi
# Extract job ID (remove angle brackets)
mqtl_id=$(echo $mqtl_id | awk '{print $2}' | tr -d '<>')
echo "Matrix eQTL job ID: ${mqtl_id}"
