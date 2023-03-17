#!/usr/bin/env bash

BUILD_DIR="workflow/scripts/grups/target"
[[ -d ${BUILD_DIR} ]] && rm -rf ${BUILD_DIR}

export CARGO_INCREMENTAL=0
export CC="x86_64-conda-linux-gnu-gcc"
CMAKE_C_FLAGS="-I./ -L./ -I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib" \
cargo install --path $(pwd)/workflow/scripts/grups --root $CONDA_PREFIX

