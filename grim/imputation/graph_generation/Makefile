CONFIG_FILE = "../../minimal-configuration.json"
wmda:
	rm -rf output
	python wmda_download.py
	python wmda_to_hpf_csv.py
	python generate_neo4j_multi_hpf.py -c $(CONFIG_FILE)


nemo:
	rm -rf output
	mkdir -p output
	python nemo_to_hpf_csv.py -c $(CONFIG_FILE)
	python generate_neo4j_multi_hpf.py -c $(CONFIG_FILE)

em2014:
	rm -rf output
	mkdir -p output
	python nemo_to_hpf_csv.py -c $(CONFIG_FILE)
	python generate_neo4j_multi_hpf.py -c $(CONFIG_FILE)
