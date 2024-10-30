#!/usr/bin/env bash
set -euo pipefail

# ---- SOURCE_DIR is the location of this script, regardless of who is sourcing/executing it.
SOURCE_DIR=$(dirname $(readlink -f "${BASH_SOURCE[0]}"))

# ---- THIS_SCRIPT is the path of the caller script.
THIS_SCRIPT=$(readlink -f "${0}")
THIS_DIR=$(dirname $THIS_SCRIPT)

echo This script is : $THIS_SCRIPT
echo Source dir is : $SOURCE_DIR


eval "$(conda shell.bash hook)"

NAME=$(grep -P 'name:' "${THIS_DIR}/$(basename ${THIS_SCRIPT%.post-deploy.sh}).yml"  | cut -d: -f2 | sed 's/ //g')

echo "Attempting to activate conda environment: ${NAME}"
conda activate ${NAME}
echo Done. Current environment is: $(conda env list | grep "*" | cut -f1 -d' ')

echo "Installing badger.plots R library..."
BADGER_PATH="$SOURCE_DIR/../badger.plots"
R --slave -e "devtools::install(\"$BADGER_PATH\", upgrade='never')"

echo "Installing badger-plots command line interface..."
cp $SOURCE_DIR/../badger-plots.R $CONDA_PREFIX/bin/badger-plots

echo "Configuring r-reticulate..."
R --slave -e "reticulate::use_condaenv(\"${NAME}\", conda=Sys.getenv(\"CONDA_EXE\"))"
R_ENVIRON=$(R --slave -e 'cat(R.home(component="home"))')/etc/Renviron
echo "BADGER_PLOTS_CONDA_ENV=\"${NAME}\"" >> $R_ENVIRON
echo "BADGER_PLOTS_CONDA_EXE=\"${CONDA_EXE}\"" >> $R_ENVIRON 
echo "Done."
