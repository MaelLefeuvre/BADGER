#!/usr/bin/bash
set -euxo pipefail

# ---- Usage: 'rerun-command.sh path/to/archive/dir/ '
# - run this from the root of BADGER
# Example: ./rerun-commands.sh /data/mlefeuvre/Projects/20230206-aDNA-kinship-simulations/Chan-Meso-CEU-004X/

BASE_DIR=$1
CORES=`nproc`

# Find out where this pipeline is expected to store its archive
ARCHIVE_DIR=$(grep "archive-dir:" config/config.yml | cut -d":" -f2 | sed 's/ //g' | sed 's/"//g')


touch_output() {
    SLEEP="1"
    TIME="202304240000.00"
    touch -amt $TIME results/meta/pipeline-metadata.yml && sleep $SLEEP
    touch -amt $TIME results/00-ped-sim/pedigree.def && sleep $SLEEP
    touch -amt $TIME results/00-ped-sim/CEU-pedigrees.* && sleep $SLEEP
    touch -amt $TIME results/00-ped-sim/CEU-pedigrees-everyone.fam && sleep $SLEEP
    touch -amt $TIME results/meta/pipeline-metadata.yml && sleep $SLEEP
    touch -amt $TIME results/00-ped-sim/CEU-pedigrees.* && sleep $SLEEP
    touch -amt $TIME results/00-ped-sim/CEU-pedigrees-everyone.fam  && sleep $SLEEP
    touch -amt $TIME results/00-ped-sim/CEU-pedigrees-twins-merged.vcf.gz* && sleep $SLEEP
    touch -amt $TIME results/00-ped-sim/sample_names.tsv && sleep $SLEEP
    touch -amt $TIME results/01-gargammel/contaminants/contaminants.tsv && sleep $SLEEP

    find results/02-preprocess/05-dedup -name *.bam -exec touch -am {} \; 
}



for run in $BASE_DIR/run-* ; do
    if [ -d "$ARCHIVE_DIR/$(basename $run)" ]; then
        echo "$(basename $run) already found in $ARCHIVE_DIR"
	echo "Skipping..."
	sleep 1
        continue	
    fi
    resources/utilities/unpack.sh metadata contaminants pedigree uncram $run/
    find results -type f -exec chmod 644 {} \;
    touch_output
    rm -rf .snakemake/incomplete/ .snakemake/locks/ .snakemake/metadata/
    snakemake READ --cores 80 --use-conda --conda-frontend mamba --touch --keep-going --quiet rules 
    snakemake READ GRUPS KIN TKGWV2 --cores ${CORES} --resources mem_mb=80000 --use-conda --conda-frontend mamba --printshellcmds --rerun-triggers mtime
    snakemake archive --cores ${CORES} --use-conda --conda-frontend mamba --printshellcmds --rerun-triggers mtime
    rm -rf results
done

