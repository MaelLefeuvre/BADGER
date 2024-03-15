#!/usr/bin/env bash

set -euxo pipefail

#eval "$(conda shell.bash hook)"
#conda activate "$(basename ${0%.post-deploy.sh})"

cd $CONDA_PREFIX

REPO="https://github.com/zmaroti/correctKin.git"
BUILD_DIR="$(basename ${REPO%.git})"

git clone $REPO
pushd ${BUILD_DIR} && make && cp ./bin/* $CONDA_PREFIX/bin/ && popd && rm -rf ${BUILD_DIR}