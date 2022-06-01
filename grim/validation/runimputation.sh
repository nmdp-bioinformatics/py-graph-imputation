#!/bin/sh

echo "Current Directory: "  `pwd`

# Clean the imputation output directory everytime
rm -rf ./output

python runfile.py -c $1
