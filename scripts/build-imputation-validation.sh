#!/usr/bin/env bash
CONFIG_FILE=${1:-conf/minimal-configuration.json}
echo "Using config file:" ${CONFIG_FILE}

# prep and load the data
python graph_generation/generate_hpf.py -c ${CONFIG_FILE}

# generate the graph
python graph_generation/generate_neo4j_multi_hpf.py -c ${CONFIG_FILE}

# run imputation
python scripts/runfile.py -c ${CONFIG_FILE}
