#!/usr/bin/env bash

TK_URL="https://github.com/danimfernandes/tkgwv2/archive/refs/heads/master.zip"

cd $CONDA_PREFIX/bin

wget "${TK_URL}"
unzip $(basename ${TK_URL}) && rm $(basename "${TK_URL}")

chmod +x tkgwv2-master/TKGWV2.py
chmod +x tkgwv2-master/scripts/*
chmod +x tkgwv2-master/helpers/*.R

ln -s tkgwv2-master/helpers/* .

# Permanently add to path
#conda develop `pwd`/tkgwv2-master/scripts
#conda develop `pwd`/tkgwv2-master/helpers
