import csv
import gzip
import json
import pathlib
import argparse

project_dir = "../../"
output_dir = "output/"

parser = argparse.ArgumentParser()
parser.add_argument(
    "-c",
    "--config",
    required=False,
    default="../../conf/minimal-configuration.json",
    help="Configuration JSON file",
    type=str,
)

args = parser.parse_args()
configuration_file = args.config

# Read configuration file and load properties
with open(configuration_file) as f:
    conf = json.load(f)

pops = conf.get("populations")
freq_data_dir = project_dir + conf.get("freq_data_dir")
pop_ratio_dir = project_dir + conf.get(
    "pops_count_file", "imputation/graph_generation/output/pop_ratio.txt"
)

# Output in HaplotypePopulationFrequency (hpf) csv file
hpf_file = output_dir + "hpf.csv"

# Create output directory if it doesn't exist
pathlib.Path(output_dir).mkdir(parents=False, exist_ok=True)

# Display the configurations we are using
print(
    "****************************************************************************************************"
)
print("Conversion to HPF file based on following configuration:")
print("\tPopulation: {}".format(pops))
print("\tFrequency File Directory: {}".format(freq_data_dir))
print("\tOutput File: {}".format(hpf_file))
print(
    "****************************************************************************************************"
)

haplist_overall = {}  # list of haplotypes across all populations
pop_hap_combos = {}

list_pop_count = []
#### Load initial frequency files
for pop in pops:
    freq_file = freq_data_dir + pop + ".freqs.gz"
    print("Reading Frequency File:\t {}".format(freq_file))
    with gzip.open(freq_file, "rb") as zf:
        count_pop = 0
        lines = [x.decode("utf8").strip() for x in zf.readlines()]
        for hap_line in lines:
            haplotype, count, freq = hap_line.split(",")
            if haplotype == "Haplo":
                continue
            freq = float(freq)
            # Ignore lines with 0 freq
            if freq == 0.0:
                continue

            pop_haplotype = pop + "-" + haplotype
            haplist_overall[haplotype] = 1
            pop_hap_combos[pop_haplotype] = freq

            count_pop += float(count)
        list_pop_count.append(count_pop)

sum_pops = sum(list_pop_count)
pop_ratio_file = open(pop_ratio_dir, "w")
for pop, ratio in zip(pops, list_pop_count):
    pop_ratio_file.write("{},{},{}\n".format(pop, ratio, (ratio / sum_pops)))


header = ["hap", "pop", "freq"]


print("Writing hpf File:\t {}".format(hpf_file))
with open(hpf_file, mode="w") as csv_file:
    csv_writer = csv.writer(csv_file, delimiter=",", quoting=csv.QUOTE_NONE)
    csv_writer.writerow(header)
    for pop_haplotype in pop_hap_combos:
        (pop, haplotype) = pop_haplotype.split("-")
        freq = pop_hap_combos[pop_haplotype]
        csv_writer.writerow([haplotype, pop, freq])
