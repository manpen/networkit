#!/bin/bash
python3 -mvenv pyenv
. pyenv/bin/activate

git submodule update --init --recursive
echo "Build and install local networkit source"
pip3 install --upgrade pip
pip3 install cython
pip3 install -e .

