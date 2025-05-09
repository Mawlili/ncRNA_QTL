#!/bin/bash
#BSUB -W 36:00                           # Wall time
#BSUB -q long                            # Queue
#BSUB -R "rusage[mem=25G]"              # Memory per core
#BSUB -M 165G                            # Total memory
#BSUB -N                                 # Email on job completion
#BSUB -u whwu1@mdanderson.org           # Email address
#BSUB -J grex[1-5]                       # Job array with 5 workers (index from 1 to 5)
#BSUB -o log/grex_%J_%I.out                 # Output file, %I is array index
#BSUB -e log/grex_%J_%I.error               # Error file

NUM_WORKERS=5

# Total ncRNA (excluding header row)
NUM_NCRNA=$(tail -n +2 /rsrch5/home/epi/bhattacharya_lab/projects/ncRNA_QTL/BRCA_gbat/1_dataprep/data/nc.bed | wc -l)
echo "Total ncRNA: $NUM_NCRNA"

# nCRNA per worker
LINES_PER_WORKER=$(( (NUM_NCRNA + NUM_WORKERS - 1) / NUM_WORKERS ))  # ceiling division
echo "Each worker processes ~${LINES_PER_WORKER} lines"

#mkdir -p ./log ./results ./results/best_models

module load R/4.3.1

Rscript ./grex.R $NUM_WORKERS $LSB_JOBINDEX
