#!/usr/bin/env bash

# Usage: ./unpack uncram /data/mlefeuvre/Projects/20230206-aDNA-kinship-simulations/archive/run-008/
# Dependencies: samtools

RED='\e[31m'
GREEN='\e[32m'
YELLOW='\e[33m'
RESET='\e[0m'
BOLD='\e[1;'

function log() {
    local log_level="$1"
    case $log_level in
        "INFO")
	    log_color="${GREEN}"
	    ;;
	"WARN")
	    log_color="${YELLOW}"
	    ;;
	"ERROR")
	    log_color="${RED}"
	    ;;
    esac

    printf "${BOLD}${log_color}[${log_level}]:${RESET} ${@:2}\n"
}

function uncram() {
    local cramfile="${1}"
    local target_dir="${2}"

    # Keep track of the pedigree number
    rep=$(basename ${cramfile} | grep -Po "ped[0-9]+")

    # Go into the directory where those bam came from and unpack the cram
    cd ${target_dir}
    samtools split -@ ${THREADS} -f "%!".bam --output-fmt BAM  ${cramfile}

    # Put each unpacked bam in its own directory 
    for bam in ./${rep}*.bam; do
        local bam_dir="${bam%.bam}"
	log INFO "Moving ${bam} in ${target_dir}/${bam_dir}"
	mkdir $bam_dir
	chmod 444 $bam
	mv $bam $bam_dir
    done
    cd - 
}

function uncram_all() {
    local OUT_DIR="${1}"
    for cramfile in $(find ${OUT_DIR} -type f | grep ".cram$"); do
        TARGET_OUTPUT=$(dirname ${cramfile#${OUT_DIR}})
        log INFO "Unpacking $cramfile in ${TARGET_OUTPUT}"
        mkdir -p ${TARGET_OUTPUT} && uncram "${cramfile}" "${TARGET_OUTPUT}"
    done
}


function copy_file() {
    local out_dir="${1}"
    local file_pattern="${2}"
    declare -a metadata=(`find ${OUT_DIR} -type f | grep "${file_pattern}$"`)
    
    [[ "${#metadata[@]}" -gt 1 ]] && { log ERROR "Found multiple candidate files: ${metadata[*]}"; exit 1; }
    local target_out_dir="$(dirname ${metadata#${out_dir}})"

    log INFO "Copying ${metadata[0]} in ${target_out_dir}"
    mkdir -p ${target_out_dir} && cp ${metadata[0]} ${target_out_dir}
}

function unpack_pedigree() {
    local out_dir="${1}"
    local expected_tar="results/00-ped-sim.tar.xz"
    log INFO "Attempting to untar ${OUT_DIR}/${expected_tar} in current directory..."
    XZ_OPT="-T ${THREADS}" tar -xJvf "${out_dir}/${expected_tar}" -C"./"
}


function unpack_xz_file() {
   local out_dir="${1}"
   local expected_file="${2}"

   declare -a tarfile=(`find ${out_dir} -type f | grep "${expected_file}$"`)
   [[ "${#tarfile[@]}" -gt 1 ]] && { log ERROR "Found multiple candidate files: ${tarfile[*]}"; exit 1; }

   local target_out_dir="$(dirname ${tarfile[0]#${out_dir}})"
   log INFO "Attempting to unzip ${out_dir}/${expected_file} in ${target_out_dir}"
   out_file="${target_out_dir}/${expected_file%.xz}"
   mkdir -p ${target_out_dir} && xz --threads ${THREADS} -d -c ${tarfile} > "${out_file}"
   chmod 444 ${out_file}
}

# ---- Main
function main() {
    FILETYPE="${1}"
    OUT_DIR="${2}"
    BASE_DIR=`pwd`
    
    THREADS=16

    case "${FILETYPE}" in
        "uncram")
	    uncram_all "${OUT_DIR}"
            ;;
	"metadata")
	    copy_file "${OUT_DIR}" "meta/pipeline-metadata.yml"
	    ;;
	"contaminants")
	    copy_file "${OUT_DIR}" "contaminants.tsv"
	    ;;
	"pedigree")
	    unpack_pedigree "${OUT_DIR}"
	    ;;
	"variants")
	    unpack_xz_file "${OUT_DIR}" "variants-intersect-[A-Z]{3}-maf[.0-9]+.ucscbed.xz"
    esac
}


if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    set -euo pipefail
    main "$@"
fi

