#!/usr/bin/env bash
set -eo pipefail

REPO_URL="https://github.com/grenaud/gargammel.git"
REPO_DIR="$(basename ${REPO_URL%.git})"

#if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
#    eval "$(conda shell.bash hook)"
#    conda activate "$(basename ${0%.post-deploy.sh})"
#fi


cd $CONDA_PREFIX
git clone --recursive $REPO_URL
cd $REPO_DIR
#git checkout 1.1.4
make
sed -i -e 's/$pathdir\.//' -e 's/\/src\///' -e 's/\/art_src_MountRainier\///' -e 's/fileExists(.*);//' gargammel.pl

cp gargammel.pl $CONDA_PREFIX/bin/gargammel

cp -r art_src_MountRainier/art_illumina src/* $CONDA_PREFIX/bin/

cd $CONDA_PREFIX
rm gargammel -rf
