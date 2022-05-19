# imputegl

## install networkx
```
pip install networkx
```

## run a test file
```
python runfile.py
```

This will look for a file data/testwmda.100.txt and generate output files: 
```
output/testwmda.100.txt_out     - imputation output in HPF
output/testwmda.100.txt_miss    - no calls 
output/testwmda.100.txt_val     - more detail on no calls
```

# Directory Structure
```
.
├── README.md
├── cypher_query.py
├── data
│   └── testwmda.100.txt
├── impute.py
├── neo4j.json
├── output
└── runfile.py
```

