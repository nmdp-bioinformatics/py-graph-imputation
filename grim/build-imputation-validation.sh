#!/usr/bin/env bash

CONFIG_FILE=${1:-conf/minimal-configuration.json}
CONFIG_FULL_PATH=${PWD}/${CONFIG_FILE}

echo "Using config file:" ${CONFIG_FULL_PATH}

# Build imputation files
cd imputation
./build.sh ${CONFIG_FULL_PATH}
cd ..

# validation
cd validation
./runimputation.sh ${CONFIG_FULL_PATH}
