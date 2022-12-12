import argparse
import cProfile
import copy
import csv
import gzip
import json
import os
import pathlib
from multiprocessing.pool import Pool

from grim.imputation.impute import Imputation
from imputegl.impute import write_best_prob, write_best_prob_genotype

from grim.imputation.networkx_graph import Graph

# Profiler start
pr = cProfile.Profile()
pr.enable()

parser = argparse.ArgumentParser()
parser.add_argument(
    "-c",
    "--config",
    required=False,
    default="../minimal-configuration-script.json",
    help="Configuration JSON file",
    type=str,
)

parser.add_argument(
    "-p",
    "--processes",
    required=False,
    default=0,
    help="Number of processes to use. Defaults to total cores - 1",
    type=int,
)

args = parser.parse_args()
configuration_file = args.config
num_of_cpus = args.processes

project_dir = "../"
output_dir = "output/"

# Read configuration file and load properties
with open(configuration_file) as f:
    json_conf = json.load(f)

config = {
    "planb": json_conf.get("planb", True),
    "pops": json_conf.get("populations"),
    "priority": json_conf.get("priority"),
    "epsilon": json_conf.get("epsilon", 1e-3),
    "number_of_results": json_conf.get("number_of_results", 1000),
    "number_of_pop_results": json_conf.get("number_of_pop_results", 100),
    "output_MUUG": json_conf.get("output_MUUG", True),
    "output_haplotypes": json_conf.get("output_haplotypes", False),
    "node_file": project_dir + json_conf.get("graph_files_path") + json_conf.get("node_csv_file"),
    "top_links_file": project_dir + json_conf.get("graph_files_path") + json_conf.get("top_links_csv_file"),
    "edges_file": project_dir + json_conf.get("graph_files_path") + json_conf.get("edges_csv_file"),
    "imputation_input_file": project_dir + json_conf.get("imputation_in_file"),
    "imputation_out_umug_freq_file": output_dir
    + json_conf.get("imputation_out_umug_freq_filename"),
    "imputation_out_umug_pops_file": output_dir
    + json_conf.get("imputation_out_umug_pops_filename"),
    "imputation_out_hap_freq_file": output_dir
    + json_conf.get("imputation_out_hap_freq_filename"),
    "imputation_out_hap_pops_file": output_dir
    + json_conf.get("imputation_out_hap_pops_filename"),
    "imputation_out_miss_file": output_dir
    + json_conf.get("imputation_out_miss_filename"),
    "imputation_out_problem_file": output_dir
    + json_conf.get("imputation_out_problem_filename"),
    "factor_missing_data": json_conf.get("factor_missing_data", 0.01),
    "loci_map": json_conf.get(
        "loci_map", {"A": 1, "B": 3, "C": 2, "DQB1": 4, "DRB1": 5}
    ),
    "matrix_planb": json_conf.get(
        "Plan_B_Matrix",
        [
            [[1, 2, 3, 4, 5]],
            [[1, 2, 3], [4, 5]],
            [[1], [2, 3], [4, 5]],
            [[1, 2, 3], [4], [5]],
            [[1], [2, 3], [4], [5]],
            [[1], [2], [3], [4], [5]],
        ],
    ),
    "pops_count_file": project_dir + json_conf.get("pops_count_file", ""),
    "use_pops_count_file": json_conf.get("pops_count_file", False),
    "number_of_options_threshold": json_conf.get("number_of_options_threshold", 100000),
    "max_haplotypes_number_in_phase": json_conf.get(
        "max_haplotypes_number_in_phase", 100
    ),
    "bin_imputation_input_file": project_dir
    + json_conf.get("bin_imputation_in_file", "None"),
    "nodes_for_plan_A": json_conf.get("Plan_A_Matrix", []),
    "save_mode": json_conf.get("save_space_mode", False),
}

all_loci_set = set()
for _, val in config["loci_map"].items():
    all_loci_set.add(str(val))

config["full_loci"] = "".join(sorted(all_loci_set))

# Display the configurations we are using
print(
    "****************************************************************************************************"
)
print("\tLoci: {}".format(config["full_loci"]))
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
print("\tFactor Missing Data: {}".format(config["factor_missing_data"]))
print("\tLoci Map: {}".format(config["loci_map"]))
print("\tPlan B Matrix: {}".format(config["matrix_planb"]))
print("\tPops Count File: {}".format(config["pops_count_file"]))
print("\tUse Pops Count File: {}".format(config["use_pops_count_file"]))
print("\tNumber of Options Threshold: {}".format(config["number_of_options_threshold"]))
print(
    "\tMax Number of haplotypes in phase: {}".format(
        config["max_haplotypes_number_in_phase"]
    )
)
print("\tbin_imputation_input_file: {}".format(config["bin_imputation_input_file"]))
if config["nodes_for_plan_A"]:
    print("\tNodes in plan A: {}".format(config["nodes_for_plan_A"]))
print("\tSave space mode: {}".format(config["save_mode"]))
print(
    "****************************************************************************************************"
)

if num_of_cpus == 0:
    # Use all the CPUs but 1
    num_of_cpus = os.cpu_count() - 1
print(
    "Using {} processes with parent PID {}".format(num_of_cpus, os.getpid()), flush=True
)

priority = config["priority"]
epsilon = config["epsilon"]
planb = config["planb"]
number_of_results = config["number_of_results"]
number_of_pop_results = config["number_of_pop_results"]
output_umug = config["output_MUUG"]
output_haplotypes = config["output_haplotypes"]
MUUG_output = config["output_MUUG"]
haps_output = config["output_haplotypes"]
imputation_infile = config["imputation_input_file"]

# Perform imputation
graph = Graph(config)
graph.build_graph(config["node_file"], config["top_links_file"], config["edges_file"])


def perform_impute_one(*args):
    subject_id, gl, race1, race2 = args[0]
    subject_bin = [1] * (len(config["full_loci"]) - 1)
    try:
        imputation = Imputation(graph, copy.deepcopy(config))
        impute_results = imputation.impute_one(
            subject_id,
            gl,
            subject_bin,
            race1,
            race2,
            priority,
            epsilon,
            1000,
            output_umug,
            output_haplotypes,
            planb,
            em=False,
        )
        imputation.print_options_count(subject_id)
    except MemoryError as e:
        print("Memory Error for subject_id: ", subject_id)
        print("Error: ", e)
        return subject_id, None, None

    return impute_results


def imputation_record_from(imputation_file):
    if not imputation_file.endswith(".gz"):
        with open(imputation_file, "r") as imputation_file_unzipped:
            for line in imputation_file_unzipped:
                subject_id, gl, race1, race2 = line.rstrip().split(",")
                yield subject_id, gl, race1, race2
    else:
        with gzip.open(imputation_file, "rt") as imputation_gzip:
            for line in imputation_gzip:
                subject_id, gl, race1, race2 = line.rstrip().split(",")
                yield subject_id, gl, race1, race2


# Create output directory if it doesn't exist
pathlib.Path(output_dir).mkdir(parents=False, exist_ok=True)

problem_ids = []
miss_ids = []
if MUUG_output:
    fout_hap_muug = open(config["imputation_out_umug_freq_file"], "w")
    fout_pop_muug = open(config["imputation_out_umug_pops_file"], "w")
if haps_output:
    fout_hap_haplo = open(config["imputation_out_hap_freq_file"], "w")
    fout_pop_haplo = open(config["imputation_out_hap_pops_file"], "w")

with Pool(processes=num_of_cpus) as pool:
    subject_records = imputation_record_from(imputation_infile)
    for imputation_results in pool.imap_unordered(perform_impute_one, subject_records):
        subject_id, umug_results, haps_results = imputation_results
        if umug_results is None:
            problem_ids.append(subject_id)
            continue

        if (len(haps_results["Haps"]) == 0 or haps_results["Haps"] == "NaN") and len(
            umug_results["Haps"]
        ) == 0:
            miss_ids.append(subject_id)

        if haps_output:
            haps = haps_results["Haps"]
            probs = haps_results["Probs"]
            pops = haps_results["Pops"]
            print(
                "Subject: {id} {hap_length} haplotypes".format(
                    id=subject_id, hap_length=len(haps)
                )
            )
            write_best_prob(subject_id, haps, probs, number_of_results, fout_hap_haplo)
            write_best_prob(
                subject_id, pops, probs, number_of_pop_results, fout_pop_haplo
            )

        if MUUG_output:
            haps = umug_results["Haps"]
            pops = umug_results["Pops"]
            print(
                "Subject: {id} {hap_length} haplotypes".format(
                    id=subject_id, hap_length=len(haps)
                )
            )
            write_best_prob_genotype(subject_id, haps, number_of_results, fout_hap_muug)
            write_best_prob_genotype(
                subject_id, pops, number_of_pop_results, fout_pop_muug
            )

num_of_misses = len(miss_ids)
print("Missing: {}".format(num_of_misses))
if num_of_misses > 0:
    with open(config["imputation_out_miss_file"], "w") as miss_csv_file:
        miss_writer = csv.writer(miss_csv_file)
        for miss_id in miss_ids:
            miss_writer.writerow([miss_id])

num_of_problems = len(problem_ids)
print("Problems: {}".format(num_of_problems))
if num_of_problems > 0:
    with open(config["imputation_out_problem_file"], "w") as problem_csv_file:
        problem_writer = csv.writer(problem_csv_file)
        for problem_id in problem_ids:
            problem_writer.writerows([problem_id])
