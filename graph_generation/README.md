# Graph Data Generation

## Pre-requisites

### Mac OS

- JDK 8
	- Install JDK 1.8 from Oracle
	- Set JAVA_HOME
	     On MacOS, add the follwoing to .bash_profile:
		```
		export JAVA_HOME="$(/usr/libexec/java_home -v 1.8)"
		```

- Python 3
	- On MacOS install with
		```
		brew install python3
		```

- Install Neo4J
	- On MacOS install with
		```
		brew install neo4j
		```

	- Setup NEO4J_HOME

        Point NEO4J_HOME to the root of the NEO4J directory.
		```
		export NEO4J_HOME=/usr/local/Cellar/neo4j/3.2.2/libexec
		```

### Linux
- JDK 8
	- Install JDK 1.8 from Oracle
	- add JAVA_HOME to ~/.bash_profile
		```
		export JAVA_HOME=$(readlink -f /usr/bin/java | sed "s:/bin/java::")
		```


- Install Neo4J
	- Download latest compressed version of Neo4j, for example
		```
		wget https://neo4j.com/artifact.php?name=neo4j-community-3.5.7-unix.tar.gz
		```
    - uncompress tar file
        ```
        tar -xf artifact.php?name=neo4j-community-3.5.7-unix.tar.gz
        ```

	- Point NEO4J_HOME to the root of the uncompressed NEO4J directory and add the following line to ~/.bash_profile

		```
		export NEO4J_HOME=path/to/neo4j-community-3.5.7
		```



# Using Makefile

For WMDA
```
make wmda
```

For Nemo
```
make nemo
```

# Manually Step by Step

# WMDA Dataset
- Download and prepare wmda data. Python script downloads reference wmda data and untars it in wmda directory
```
	cd graph_generation
	python wmda_download.py
```

- Generate nodes/edges/toplinks from the reference wmda data. The freqs file is cnverted to HPF format first.
```
	python wmda_to_hpf_csv.py
	python generate_neo4j_single_hpf.py
```

- Generated nodes and edges (creates edges; and nodes with “Frequency” float attribute, not array) files are in the output/csv directory.

```
	output
	└── csv
	    ├── edges.csv
	    ├── nodes.csv
	    └── top_links.csv
```

# NEMO Dataset

To use a different set of frequencies use the following procedure:

- Starting in the graph generator directory, convert the data from frequency format to hpf (haplotype, population, frequency).
```
   python nemo_to_hpf_csv.py
```

- This program looks for a data/NEMO2011 directory and reads the individual frequency files and generates this csv:
```
    output/hpf.csv
```
- To generate the graph files run this script
```
    python generate_neo4j_multi_hpf.py
```
    Note: the script reads from data/NEMO2011/hpf.csv and writes the three graph csv files to the directory
```
	output/csv
	├── edges.csv
	├── nodes.csv
	└── top_links.csv
```
    Note: there is an option to trim the frequency set below a frequency threshold.  If the trimming threshold is 1e-6 it will take 9m35s to generate the graph csv files on a mid-2015 MacBook Pro (2.5 GHz Intel Core i7) and will result in 1,088,817 nodes (159MB), 14,868,976 edges (2.0GB)and 5,947,591 top links (108MB).
