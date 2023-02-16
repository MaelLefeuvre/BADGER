#!/usr/bin/env bash

# Usage: './loop-pipeline.sh 10 80' will run the snakemake pipeline 
#  - 10 times 
#  - using 80 cores 
# ...and archive the results between each iteration

set -euo pipefail
N=$1
CORES="$2"

GREEN='\u001b[32m'
RED='\u001b[31m'
RESET='\u001b[0m'

RESULTS_DIR="./results"
[[ -d $RESULTS_DIR ]] && rm -ri $RESULTS_DIR

mkdir -p $RESULTS_DIR
for ((i=1; i<=$N; i++)); do
    printf "${GREEN}[INFO]: ($i) Removing previous results${RESET}\n" &&
    rm -r $RESULTS_DIR &&
    printf "${GREEN}[INFO]: ($i) Running Snakemake pipeline${RESET}\n" &&
    snakemake all --cores $CORES --use-conda --conda-frontend mamba --printshellcmds --rerun-incomplete &&
    printf "${GREEN}[INFO]: ($i) Archiving results${RESET}" &&
    snakemake archive --cores $CORES --use-conda --conda-frontend mamba --printshellcmds --rerun-triggers 'mtime' ||
    { echo -e "${RED}[ERROR]: ($i) Failed Run.${RESET}"; exit 1; }
done
