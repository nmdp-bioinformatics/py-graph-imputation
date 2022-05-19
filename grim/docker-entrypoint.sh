#!/usr/bin/env bash

CMD=${1:-minimal}

cd /grimm/

if [[ "${CMD}" = "minimal" ]]; then
    echo "Running minimal build"
    ./build-imputation-validation.sh
elif [[ "${CMD}" = "clean" ]]; then
    echo "Clean output directories"
    cd validation/haplogic
    make clean_all
elif [[ "${CMD}" = "wmda" ]]; then
    echo "Running wmda validation"
    cd validation/wmda
    ./wmda-pipeline.sh
elif [[ "${CMD}" = "haplogic" ]]; then
    echo "**Running Haplogic Validation**"
    shift # remove haplogic
    SUB_CMD=$1
    shift # remove 5loci/6loci
    if [[ "${SUB_CMD}" = "6loci" ]]; then
        echo "6 Loci validation"
        cd validation/haplogic
        make PAT_CONFIG_FILE=$PWD/../../conf/haplogic/6loci-pat-configuration.json\
             DON_CONFIG_FILE=$PWD/../../conf/haplogic/6loci-don-configuration.json \
             "$@"
    elif [[ "${SUB_CMD}" = "5loci" ]]; then
        echo "5 Loci validation"
        cd validation/haplogic
        make PAT_CONFIG_FILE=$PWD/../../conf/haplogic/5loci-pat-configuration.json\
             DON_CONFIG_FILE=$PWD/../../conf/haplogic/5loci-don-configuration.json \
             "$@"
    else
        echo "Unknown option:" "$CMD" "$SUB_CMD" "$@"
    fi
elif [[ "${CMD}" = "make" ]]; then
    echo "Running Haplogic validation"
    cd validation/haplogic
    exec "$@"
else
    exec "$@"
fi

