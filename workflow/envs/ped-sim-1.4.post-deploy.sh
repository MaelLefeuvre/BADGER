#!/usr/bin/env bash

set -euxo pipefail

GIT_RELEASE="https://github.com/williamslab/ped-sim/archive/refs/tags/v1.4.tar.gz"
BUILD_DIR="${CONDA_PREFIX}/ped-sim-1.4"


cd $CONDA_PREFIX &&
wget -O- ${GIT_RELEASE} | tar -xvzf- &&
cd ${BUILD_DIR} &&
make -f Makefile-gsl CFLAGS="-Wall -DUSEGSL -I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib" && 
cp ped-sim fam2def.py plot-fam.R $CONDA_PREFIX/bin && rm -r ${BUILD_DIR}
