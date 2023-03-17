#!/usr/bin/env bash

set -euxo pipefail

OPT_LEVEL=3

GIT_RELEASE="https://github.com/williamslab/ped-sim/archive/refs/tags/v1.4.tar.gz"
BUILD_DIR="${CONDA_PREFIX}/ped-sim-1.4"


cd $CONDA_PREFIX &&
wget -O- ${GIT_RELEASE} | tar -xvzf- &&
cd ${BUILD_DIR} &&
sed -i 's/-O2/-O$(OPT_LEVEL)/' Makefile &&
make -f Makefile OPT_LEVEL=$OPT_LEVEL CFLAGS="\$(CFLAGS) -I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib" && 
cp ped-sim fam2def.py plot-fam.R $CONDA_PREFIX/bin && rm -r ${BUILD_DIR}
