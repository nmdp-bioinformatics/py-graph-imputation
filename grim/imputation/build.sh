#!/usr/bin/env bash

if [[ -z "$1" ]]
then
  echo "Usage: `basename $0` <config.json>"
  exit 1
fi

# Install package
pip uninstall -y imputegl
pip uninstall -y imputegl

python setup.py install
pip install .

# Build graph
cd graph_generation

make nemo CONFIG_FILE=${1}
echo "Graph data:"
ls -lah output/*

cd ..
