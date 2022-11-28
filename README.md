# py-graph-imputation
[![PyPi Version](https://img.shields.io/pypi/v/py-graph-imputation.svg)](https://pypi.python.org/pypi/py-graph-imputation)

## Graph Imputation

`py-graph-imputation` is the successor of [GRIMM](https://github.com/nmdp-bioinformatics/grimm) written in Python and based on [NetworkX](https://networkx.org/)

![GRIM Dependencies](images/py-graph-imputation.png)

How to develop on the project locally.

1. Make sure the following pre-requites are installed.
   1. `git`
   2. `python >= 3.8`
   3. build tools eg `make`
2. Clone the repository locally
    ```shell
    git clone git@github.com:nmdp-bioinformatics/py-graph-imputation.git
    cd py-graph-imputation
    ```
3. Make a virtual environment and activate it, run `make venv`
   ```shell
    > make venv
      python3 -m venv venv --prompt py-graph-imputation-venv
      =====================================================================
    To activate the new virtual environment, execute the following from your shell
    source venv/bin/activate
   ```
4. Source the virtual environment
   ```shell
   source venv/bin/activate
   ```
5. Development workflow is driven through `Makefile`. Use `make` to list show all targets.
   ```
    > make
    clean                remove all build, test, coverage and Python artifacts
    clean-build          remove build artifacts
    clean-pyc            remove Python file artifacts
    clean-test           remove test and coverage artifacts
    lint                 check style with flake8
    behave               run the behave tests, generate and serve report
    pytest               run tests quickly with the default Python
    test                 run all(BDD and unit) tests
    coverage             check code coverage quickly with the default Python
    dist                 builds source and wheel package
    docker-build         build a docker image for the service
    docker               build a docker image for the service
    install              install the package to the active Python's site-packages
    venv                 creates a Python3 virtualenv environment in venv
    activate             activate a virtual environment. Run `make venv` before activating.
   ```
6. Install all the development dependencies. Will install packages from all `requirements-*.txt` files.
   ```shell
    make install
   ```
7. The Gherkin Feature files, step files and pytest files go in `tests` directory:
    ```
    tests
    |-- features
    |   |-- algorithm
    |   |   `-- SLUG\ Match.feature
    |   `-- definition
    |       `-- Class\ I\ HLA\ Alleles.feature
    |-- steps
    |   |-- HLA_alleles.py
    |   `-- SLUG_match.py
    `-- unit
        `-- test_grim.py
    ```
8. Package Module files go in the `grim` directory.
    ```
    grim
    |-- __init__.py
    |-- algorithm
    |   `-- match.py
    |-- model
    |   |-- allele.py
    |   `-- slug.py
    `-- grim.py
    ```
9. Run all tests with `make test` or different tests with `make behave` or `make pytest`. `make behave` will generate report files and open the browser to the report.
10. Use `python app.py` to run the Flask service app in debug mode. Service will be available at http://localhost:8080/
11. Use `make docker-build` to build a docker image using the current `Dockerfile`.
12. `make docker` will build and run the docker image with the service.  Service will be available at http://localhost:8080/


# Runing a minimal configuration example

From the main directory of the repo run:
```
scripts//build-imputation-validation.sh
```

This will pepare and load frequency data into the graph and run imputation on a sample set of subjects.

The execution is driven by the configuration file:
`conf/minimal-configuration.json`

It takes input from this file:
```
data/subjects/donor.csv
```


And genrates an `output` directory with these contents:

```
output
├── don.miss
├── don.pmug
├── don.pmug.pops
├── don.problem
├── don.umug
└── don.umug.pops
```

The .miss and .problem files are cases contain cases that failed due to errors.

The .pmug file contains the Phased Multi-locus Unambiguous Genotypes.

The .umug file contains the Unphased Multi-locus Unambiguous Genotypes.


The format of both files is:

* id
* genotype - in glstring format
* frequency
* rank


The .pmug.pops and .umug.pops contain the corresponding population assignments.

The format of the .pops files is:

* id
* pop1
* pop2
* frequency
* rank
