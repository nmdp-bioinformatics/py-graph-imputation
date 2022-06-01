import argparse
import cProfile
import json
import math
import os
import pathlib
import string
from multiprocessing import Process

from imputegl import Imputation
from imputegl.networkx_graph import Graph

# Profiler start
pr = cProfile.Profile()
pr.enable()

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--config",
                    required=False,
                    default="../minimal-configuration.json",
                    help="Configuration JSON file",
                    type=str)

parser.add_argument("-d", "--processes",
                    required=False,
                    default=4,
                    help="number of processes",
                    type=int)

args = parser.parse_args()
configuration_file = args.config
num_processes = args.processes

project_dir = "../"
output_dir = "output/"

# Read configuration file and load properties
with open(configuration_file) as f:
    json_conf = json.load(f)

config = {
    "planb": json_conf.get('planb', True),
    "pops": json_conf.get('populations'),
    "priority": json_conf.get('priority'),
    "epsilon": json_conf.get('epsilon', 1e-3),
    "number_of_results": json_conf.get('number_of_results', 1000),
    "number_of_pop_results": json_conf.get('number_of_pop_results', 100),
    "output_MUUG": json_conf.get("output_MUUG", True),
    "output_haplotypes": json_conf.get("output_haplotypes", False),
    "node_file": project_dir + json_conf.get("node_csv_file"),
    "top_links_file": project_dir + json_conf.get("top_links_csv_file"),
    "edges_file": project_dir + json_conf.get("edges_csv_file"),
    "imputation_input_file": project_dir + json_conf.get("imputation_in_file"),
    "imputation_out_umug_freq_file": output_dir + json_conf.get("imputation_out_umug_freq_filename"),
    "imputation_out_umug_pops_file": output_dir + json_conf.get("imputation_out_umug_pops_filename"),
    "imputation_out_hap_freq_file": output_dir + json_conf.get("imputation_out_hap_freq_filename"),
    "imputation_out_hap_pops_file": output_dir + json_conf.get("imputation_out_hap_pops_filename"),
    "imputation_out_miss_file": output_dir + json_conf.get("imputation_out_miss_filename"),
    "imputation_out_problem_file": output_dir + json_conf.get("imputation_out_problem_filename")
}

# Display the configurations we are using
print('****************************************************************************************************')
print("Performing imputation based on:")
print("\tPopulation: {}".format(config["pops"]))
print("\tPriority: {}".format(config["priority"]))
print("\tEpsilon: {}".format(config["epsilon"]))
print("\tPlan B: {}".format(config["planb"]))
print("\tNumber of Results: {}".format(config["number_of_results"]))
print("\tNumber of Population Results: {}".format(config["number_of_pop_results"]))
print("\tNodes File: {}".format(config["node_file"]))
print("\tTop Links File: {}".format(config["edges_file"]))
print("\tInput File: {}".format(config["imputation_input_file"]))
print("\tOutput UMUG Format: {}".format(config["output_MUUG"]))
print("\tOutput UMUG Freq Filename: {}".format(config["imputation_out_umug_freq_file"]))
print("\tOutput UMUG Pops Filename: {}".format(config["imputation_out_umug_pops_file"]))
print("\tOutput Haplotype Format: {}".format(config["output_haplotypes"]))
print("\tOutput HAP Freq Filename: {}".format(config["imputation_out_hap_freq_file"]))
print("\tOutput HAP Pops Filename: {}".format(config["imputation_out_hap_pops_file"]))
print("\tOutput Miss Filename: {}".format(config["imputation_out_miss_file"]))
print("\tOutput Problem Filename: {}".format(config["imputation_out_problem_file"]))
print('****************************************************************************************************')

# Perform imputation
graph = Graph()
graph.build_graph(config["node_file"], config["top_links_file"], config["edges_file"])
imputation = Imputation(graph, config["pops"])

# Create output directory if it doesn't exist
pathlib.Path(output_dir).mkdir(parents=False, exist_ok=True)

# Write out the results from imputation

input_file = config["imputation_input_file"]
in_dir = os.path.dirname(input_file)
in_file_basename = os.path.basename(input_file)

num_subjects = os.popen('wc -l ' + config["imputation_input_file"]).read()
num_subjects = int(num_subjects.strip().split(" ")[0])
print(num_subjects)

split_cmd = 'split  -l ' + str(
    int(math.ceil(num_subjects / num_processes))) + " " + input_file + " " + in_dir + "/" + in_file_basename[0:4]

os.system(split_cmd)

# Write out the results from imputation
alpha = string.ascii_lowercase

processes = []

for i in range(num_processes):
    in_file = in_dir + "/" + in_file_basename[0:4] + "a" + alpha[i]
    print(in_file)
    config["imputation_input_file"] = in_file
    out_file = output_dir + os.path.basename(in_file)

    config["imputation_out_umug_freq_file"] = out_file + '.umug.freqs'
    config["imputation_out_umug_pops_file"] = out_file + '.umug.pops'
    config["imputation_out_hap_freq_file"] = out_file + '.hap.freqs'
    config["imputation_out_hap_pops_file"] = out_file + '.hap.pops'
    config["imputation_out_miss_file"] = out_file + '.miss'
    config["imputation_out_problem_file"] = out_file + '.problem'
    t = Process(target=imputation.impute_file, args=(config,))
    t.start()
    processes.append(t)

for p in processes:
    p.join()

# Profiler end
pr.disable()
pr.print_stats(sort="time")
