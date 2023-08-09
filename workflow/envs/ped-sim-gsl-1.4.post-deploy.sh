#!/usr/bin/env bash

#eval "$(conda shell.bash hook)"
#conda activate "$(basename ${0%.post-deploy.sh})"


#set -euo pipefail

OPT_LEVEL=2

GIT_RELEASE="https://github.com/williamslab/ped-sim/archive/refs/tags/v1.4.tar.gz"
BUILD_DIR="${CONDA_PREFIX}/ped-sim-1.4"


cd $CONDA_PREFIX &&
wget -O- ${GIT_RELEASE} | tar -xvzf- &&
cd ${BUILD_DIR} &&
sed -i 's/-O2/-O$(OPT_LEVEL)/' Makefile-gsl &&

make -f Makefile-gsl OPT_LEVEL=$OPT_LEVEL DEFINES="-march=native -DUSEGSL -I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib" && 
cp ped-sim fam2def.py plot-fam.R $CONDA_PREFIX/bin && rm -r ${BUILD_DIR}
