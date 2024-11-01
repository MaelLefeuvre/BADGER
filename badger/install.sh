#!/usr/bin/env bash
set -euo pipefail

VERSION="0.4.0"
CONDA="${CONDA_EXE:-$(which conda)}"
TERM="${TERM:-dumb}"

THIS_SCRIPT_DIR=$(dirname $(readlink -f "${BASH_SOURCE[0]}"))

BADGER_YAML="$THIS_SCRIPT_DIR/envs/badger-${VERSION}.yml"
BADGER_PLOT_YAML="$THIS_SCRIPT_DIR/src/badger-plots/envs/badger-plots.yml"


hr(){
    printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | sed 's/ /─/g'
}

log(){
    GREEN='\e[1;32m'
    RESET='\e[0m'
    echo -e "${GREEN}${@}${RESET}"
}


abort(){
    RED='\e[1;31m'
    RESET='\e[0m'
    echo -e "${RED}[ERROR] ${@}${RESET}"
    exit 1
}


_conda_optargs(){    
    [[ $# -gt 0 ]] && optargs="--name ${1}" || optargs=""
    echo "${optargs}"
}

_install_badger(){
    log "Installing ${1}"
    hr
    $CONDA env create $(_conda_optargs $@) --force -f "${BADGER_YAML}"
}

_install_badger_plots(){
    log "Installing badger-plots"
    hr
    $CONDA env update $(_conda_optargs $@) --file "${BADGER_PLOT_YAML}"
}

_post_deploy(){
    bash ${1%.yml}.post-deploy.sh
}

_test(){
    what=${1}

    SMK_VERSION="7.20.0"

    _search_path_ok(){
        which $1 | grep -q "${what}"
    }
    
    eval "$($CONDA shell.bash hook)"

    log "Running quick test-suite for ${what}"
    hr

    conda activate ${what} \
    && log "  - ${what} conda environment activated" \
    || abort "cannot activate ${what} conda environment"

    [[ "$(snakemake --version)" =  "${SMK_VERSION}" ]] \
    && log "  - snakemake version ${SMK_VERSION} found" \
    && _search_path_ok snakemake \
    || abort "  - snakemake ${SMK_VERSION} not found"

    _search_path_ok badger \
    && log "  - badger command line program found" \
    || abort "failed to find badger command line program not found"

    R --slave -e 'library(badger.plots)' \
    && log "  - badger.plots R library found" \
    && _search_path_ok R \
    || abort "badger.plots R library not found"
 
    _search_path_ok badger-plots \
    && badger-plots > /dev/null 2>&1 \
    && log "  - badger-plots command line program found" \
    || abort "badger-plots command line program not found."

    # ---- Run pytest tests
    hr
    log "Running pytest ..."
    pytest ./badger --verbose || abort "Some integration tests failed (See above)."

    hr
    log "Done. All tests successful!"

}

main(){
    local what="${1}"
    case "${what}" in
        badger)
            _install_badger "${what}-${VERSION}"
            _install_badger_plots "${what}-${VERSION}"
            _post_deploy ${BADGER_YAML}
            ;;

        badger-minimal)
            _install_badger "${what}-${VERSION}"
            ;;

        badger-plots)
            _install_badger_plots
            _post_deploy ${BADGER_PLOT_YAML}
            ;;
        test)
            _test "badger-${VERSION}"
            ;;
    esac
    log "Done"

}

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    [[ $# -eq 0 ]] && what="badger" || what="${1}"
    main "${what}"
fi

