#!/bin/bash
PYENV="pyenv"
PYCMD="python3"


if [ ! -d $PYENV  ]; then
	$PYCMD -mvenv $PYENV
fi

. $PYENV/bin/activate

$PYCMD -c "import networkit"
if [ $? -ne 0 ]; then
	echo "Build and install local networkit source"
	pip3 install --upgrade pip
	pip3 install cython
	pip3 install -e .
fi

./benchmark.py

