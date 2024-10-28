#!/usr/bin/env bash

cd $CONDA_PREFIX/bin

git clone https://bitbucket.org/tguenther/read.git

cd read
git checkout v1.0
cp READ.py READscript.R ../
cd - 

chmod +x READ.py
chmod +x READscript.R

rm -rf ./read
