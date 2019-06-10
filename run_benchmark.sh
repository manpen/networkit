#!/bin/bash
PYENV="pyenv"
PYCMD="python3"

if [ ! -d $PYENV  ]; then
	$PYCMD -mvenv $PYENV
fi

. $PYENV/bin/activate

test=`$PYCMD -c "import networkit; print('found')"` 2> /dev/null || true
if [[ $string == *"found"* ]]; then
   git submodule update --init --recursive
	echo "Build and install local networkit source"
	pip3 install --upgrade pip
	pip3 install cython ipython
	pip3 install -e .
fi

./benchmark.py

