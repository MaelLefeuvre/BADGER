#!/usr/bin/env bash

# Usage: './loop-pipeline.sh 10 80' will run the snakemake pipeline 
#  - 10 times 
#  - using 80 cores 
# ...and archive the results between each iteration
#
# To rerun the pipeline without cleaning up results on the first iteration, add '--no-remove' as a third argument, i.e
# ./loop-pipeline.sh 10 80 --no-remove

set -euo pipefail
N=$1
CORES="$2"

NO_REMOVE_ARG="--no-remove"
if [ $# -lt 3 ]; then
    DO_REMOVE=true
else
    [[ ! "$3" = "$NO_REMOVE_ARG" ]] && DO_REMOVE=true || DO_REMOVE=false
fi

echo "REMOVING REQUESTED?: $DO_REMOVE"

GREEN='\u001b[32m'
RED='\u001b[31m'
RESET='\u001b[0m'

RESULTS_DIR="./results"

if $DO_REMOVE; then
    [[ -d $RESULTS_DIR ]] && rm -ri $RESULTS_DIR
fi

mkdir -p $RESULTS_DIR
for ((i=1; i<=$N; i++)); do
    (if [ $i -eq 1 ] && ! $DO_REMOVE; then
    printf "${GREEN}[INFO]: ($i) Skipping removal for the first run, as $NO_REMOVE_ARG requested${RESET}\n"
    else
        printf "${GREEN}[INFO]: ($i) Removing previous results${RESET}\n" &&
        rm -r $RESULTS_DIR
    fi) &&
    printf "${GREEN}[INFO]: ($i) Running Snakemake pipeline${RESET}\n" &&
    snakemake READ KIN TKGWV2 GRUPS --cores $CORES --resources mem_mb=100000 --use-conda --conda-frontend mamba --printshellcmds --rerun-incomplete --keep-going --rerun-triggers mtime &&
    printf "${GREEN}[INFO]: ($i) Archiving results${RESET}" &&
    snakemake archive --cores $CORES --use-conda --conda-frontend mamba --printshellcmds --rerun-triggers 'mtime' ||
    { echo -e "${RED}[ERROR]: ($i) Failed Run.${RESET}"; exit 1; }
done
