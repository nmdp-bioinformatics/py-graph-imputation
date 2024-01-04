# py-graph-imputation
[![PyPi Version](https://img.shields.io/pypi/v/py-graph-imputation.svg)](https://pypi.python.org/pypi/py-graph-imputation)

* [Graph Imputation](#graph-imputation)
* [Development](#develop)
* [Running A Minimal Example Imputation](#running-a-minimal-configuration-example)

### Graph Imputation

`py-graph-imputation` is the successor of [GRIMM](https://github.com/nmdp-bioinformatics/grimm) written in Python and based on [NetworkX](https://networkx.org/)

![GRIM Dependencies](images/py-graph-imputation.png)

### Development
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
5. Run the following commands:
    cd grim/imputation
    python setup.py build_ext --inplace
6. Development workflow is driven through `Makefile`. Use `make` to list show all targets.
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
7. Install all the development dependencies. Will install packages from all `requirements-*.txt` files.
   ```shell
    make install
   ```
8. Package Module files go in the `grim` directory.
    ```
    grim
    |-- __init__.py
    |-- grim.py
    `-- imputation
        |-- __init__.py
        |-- cutils.pyx
        |-- cypher_plan_b.py
        |-- cypher_query.py
        |-- impute.py
        `-- networkx_graph.py
    ```
9. Run all tests with `make test` or different tests with `make behave` or `make pytest`.
10. Run `make lint` to run the linter and black formatter.
11. Use `python app.py` to run the Flask service app in debug mode. Service will be available at http://localhost:8080/
12. Use `make docker-build` to build a docker image using the current `Dockerfile`.
13. `make docker` will build and run the docker image with the service.  Service will be available at http://localhost:8080/


### Running a minimal configuration example

From the main directory of the repo run:
```
scripts/build-imputation-validation.sh
```

This will prepare and load frequency data into the graph and run imputation on a sample set of subjects.

The execution is driven by the configuration file: `conf/minimal-configuration.json`

It takes input from this file:
```
data/subjects/donor.csv
```

And generates an `output` directory with these contents:

```
output
├── don.miss
├── don.pmug
├── don.pmug.pops
├── don.problem
├── don.umug
└── don.umug.pops
```

The `.problem` file contains cases that failed due to serious errors (e.g., invalid HLA).

The `.miss` file contains cases where there was no output possible given the input, frequencies and configuration options.

The `.pmug` file contains the Phased Multi-locus Unambiguous Genotypes.

The `.umug` file contains the Un-phased Multi-locus Unambiguous Genotypes.


The format of both files is (csv):

* id
* genotype - in glstring format
* frequency
* rank


The `.pmug.pops` and `.umug.pops` contain the corresponding population assignments.

The format of the `.pops` files is (csv):

* id
* pop1
* pop2
* frequency
* rank
