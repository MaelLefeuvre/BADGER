#!/usr/bin/env bash
set -euo pipefail 

THIS_DIR="$(dirname $(readlink -f "${BASH_SOURCE[0]}"))"

BADGER_PLOTS_POSTDEPLOY_SCRIPT="$THIS_DIR/../src/badger-plots/envs/badger-plots.post-deploy.sh"
source $BADGER_PLOTS_POSTDEPLOY_SCRIPT

# ---- rpy2 vs r-archive shenanigans...
# See: https://github.com/rpy2/rpy2/issues/1121
CONDA_ETC_DIR="$CONDA_PREFIX/etc/conda"
cd $CONDA_PREFIX/etc/conda
mkdir -p "$CONDA_ETC_DIR/activate.d" "$CONDA_ETC_DIR/deactivate.d"


tee -a $CONDA_ETC_DIR/activate.d/env_vars.sh  > /dev/null << 'EOF'
#!/bin/sh

export BADGER_LIBR_PATH="$CONDA_PREFIX/lib/R/lib/libR.so"
export LD_PRELOAD="${BADGER_LIBR_PATH}${LD_PRELOAD+:${LD_PRELOAD}}"
EOF

tee -a $CONDA_ETC_DIR/deactivate.d/env_vars.sh  > /dev/null << 'EOF'
#!/bin/sh

export LD_PRELOAD="$(echo "${LD_PRELOAD#$BADGER_LIBR_PATH}" | sed 's/^://')"
[[ ${#LD_PRELOAD} -eq 0 ]] && unset LD_PRELOAD
EOF



# ---- End rpy2 vs r-archive shenanigans
pip install -U pip setuptools~=74.1.1 packaging~=24.1
pip install "$THIS_DIR/.."
