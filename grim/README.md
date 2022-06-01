![Build Status](https://codebuild.us-east-1.amazonaws.com/badges?uuid=eyJlbmNyeXB0ZWREYXRhIjoiVjVVUmpITGROeStDNzBrLzkreDZWVEhWblBuYjVIRzc1a3QxY1gyTE9aMjdMeEc1UFdGVFVHWGdrWWVpNURIcG44Z1FOdHBHTWZYdytYZkxyZ3haQWI4PSIsIml2UGFyYW1ldGVyU3BlYyI6Ik1aVFQ1dFJtQVVsaGl0ZW4iLCJtYXRlcmlhbFNldFNlcmlhbCI6MX0%3D&branch=master)

# graph-imputation-match

![Build Status](https://codebuild.us-east-1.amazonaws.com/badges?uuid=eyJlbmNyeXB0ZWREYXRhIjoiVjVVUmpITGROeStDNzBrLzkreDZWVEhWblBuYjVIRzc1a3QxY1gyTE9aMjdMeEc1UFdGVFVHWGdrWWVpNURIcG44Z1FOdHBHTWZYdytYZkxyZ3haQWI4PSIsIml2UGFyYW1ldGVyU3BlYyI6Ik1aVFQ1dFJtQVVsaGl0ZW4iLCJtYXRlcmlhbFNldFNlcmlhbCI6MX0%3D&branch=master)

Graph-based imputation and matching.

| Subdirectory    | Description |
| :-------------- | :---------- |
| imputation      | graph-based imputation |
| validation      | validation of components of this repository |
| matching        | graph-based matching |

# Walkthrough (With Docker)

Make sure you have [Docker](https://www.docker.com/products/docker-desktop) installed. If you're running Docker Desktop,
you will need to bump up the memory given to Docker. See Docker -> Preferences -> Resources.

A Docker image is available with all the pre-requisite software installed to run the whole pipeline.
It includes JDK, Neo4j, Perl + Libraries and Python.

Make sure you're in the `graph-imputation-match` directory when running `docker` commands.

You may need to run `sudo docker` depending on how docker was installed.

## Run minimal build

This will build the imputation library and run validation on 20 AAFA simulated sample data.

```
docker run --name grimm --rm -v $PWD:/grimm nmdpbioinformatics/grimm-networkx
```

## Run WMDA Validation

Run imputation and matching on the sample WMDA validation data. 
This requires at least 16Gb of memory. Suggest EC2 size: t2.xlarge

```
docker run --name grimm --rm -v $PWD:/grimm nmdpbioinformatics/grimm-networkx wmda
```

This will produce a `matching/output/cmp.txt` file which lists the difference in the grimm version versus Consensus rsults.

```
wc -l output/cmp.txt
       0 output/cmp.txt
```

## Run Haplogic Validation

### 6 Loci

This will run the 6-Loci version of validation. (includes DRB4/5/6)

```
docker run --name grimm --rm -v $PWD:/grimm nmdpbioinformatics/grimm-networkx haplogic 6loci
```

### 5 Loci

This will run the 5-Loci version of validation. 

```
docker run --name grimm --rm -v $PWD:/grimm nmdpbioinformatics/grimm-networkx haplogic 5loci
```

### Number of parallel processes

By default Haplogic processes will use `n - 1` cpu available. You can specify number of parallel processes by passing in `NUM_OF_PROCESSES` argument.

```
docker run --name grimm --rm -v $PWD:/grimm nmdpbioinformatics/grimm-networkx haplogic 5loci NUM_OF_PROCESSES=6
```

## Clean output directory

Remove the output directories.

```
docker run --name grimm --rm -v $PWD:/grimm nmdpbioinformatics/grimm-networkx clean
```

## Manually run inside the container

The following gives you a `bash` shell inside the container where all the pre-req software is installed.
This mounts the current directory in `/grimm/`

```
docker run --name grimm --rm -it -v $PWD:/grimm nmdpbioinformatics/grimm-networkx bash
```

# Walkthrough (Without Docker)

## 0. Pre-requisites

- Clone the repository 
```
git clone https://github.com/nmdp-bioinformatics/graph-imputation-match
```

- Python 3
The project requires Python version 3. Install from <https://www.python.org/> If  on MacOS python3 can be installed with homebrew as well.
```
brew install python3
```

- Setup Python3 virtual environment
```
cd graph-imputation-match
python3 -m venv venv
source venv/bin/activate 
```
# Validate 
`./build-imputation-validation.sh` script runs through graph generation and imputation with minimal configuration. 

The configuration file is a json file that's used by graph generation and imputation process to select parameters.
Supply a new configuration file to test new population set/freq threshold/priors etc.

`conf/minimal-configuration.json` is used for testing that any code changes does not break the build. Make a copy of this file to create your own config.

# Individual Steps
## 1. Generate Imputation Graph
 Generate haplotype frequency CSV files.

### Build nemo nodes/csv files

```
cd imputation
./build.sh $PWD/../conf/minimal-configuration.json
cd ..
 ```

Output directory is created and *edges*, *nodes*, *top links* and *info node* CSV files are generated.
```
imputation/graph_generation/output
├── [ 192]  csv
│   ├── [969M]  edges.csv
│   ├── [ 221]  info_node.csv
│   ├── [ 92M]  nodes.csv
│   └── [ 71M]  top_links.csv
└── [ 25M]  hpf.csv
```

## 2. Perform Imputation
Perform imputation on simulated data.

```
cd validation
```

The input data lives in `data` folder. 
```
data
└── [110K]  AAFA_NAMER_DONOR_20.in
```

Start Simulated Imputation
 ```
 ./runimputation.sh ../conf/minimal-configuration.json 
 ```

The imputed results are in the `output` directory:

```
output
├── [1.9M]  AAFA_NAMER_DONOR_20.in_outdeltaFunc_hap
├── [ 848]  AAFA_NAMER_DONOR_20.in_outdeltaFunc_pop
├── [   0]  AAFA_NAMER_DONOR_20.miss
└── [   0]  AAFA_NAMER_DONOR_20.problem
```

## 3. Perform Matching
  Follow [Match Validation](matching/README.md) to perform matching and compare against consensus results.
