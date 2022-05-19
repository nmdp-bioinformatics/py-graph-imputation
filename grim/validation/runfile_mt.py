import argparse
import cProfile
import json
import os
import pathlib
import threading
import string,math

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

parser.add_argument("-d", "--threads",
                    required=False,
                    default=4,
                    help="number of threads", 
                    type=int)


args = parser.parse_args()
configuration_file = args.config
num_threads = args.threads

project_dir = "../"
output_dir = "output/"

# Read configuration file and load properties
with open(configuration_file) as f:
    conf = json.load(f)

    
number_of_pop_results = conf.get('number_of_pop_results', 100)
number_of_results = conf.get('number_of_results', 1000)
planb = conf.get('planb')  
epsilon = conf.get('epsilon', 1e-3)

pops = conf.get("populations")
priority = conf.get('priority')

output_haplotypes = conf.get("output_haplotypes", False)
output_MUUG = conf.get("output_MUUG", True)

node_file = project_dir + conf.get("node_csv_file")
top_linkks_file = project_dir + conf.get("top_links_csv_file")
edges_file =project_dir + conf.get("edges_csv_file")

input_file = project_dir + conf.get("imputation_in_file")
# Output filename is based on input filename
#output_file = output_dir + os.path.basename(input_file) + '_out'

# Display the configurations we are using
print('****************************************************************************************************')
print("Performing imputation based on:")
print("\tPopulation: {}".format(pops))
print("\tPriority: {}".format(priority))
print("\tEpsilon: {}".format(epsilon))
print("\tPlan B: {}".format(planb))
print("\tNumber of Results: {}".format(number_of_results))
print("\tNumber of Population Results: {}".format(number_of_pop_results))
print("\tNodes File: {}".format(node_file))
print("\tTop Links File: {}".format(edges_file))
print("\tInput File: {}".format(input_file))
print('****************************************************************************************************')



# Perform imputation
graph = Graph()
graph.build_graph(node_file,top_linkks_file, edges_file)

imputation = Imputation(graph, pops)

# Create output directory if it doesn't exist
pathlib.Path(output_dir).mkdir(parents=False, exist_ok=True)

# Write out the results from imputation

in_dir = os.path.dirname(input_file)
in_file_basename = os.path.basename(input_file)

num_subjects = os.popen('wc -l ' + input_file ).read()
num_subjects = int(num_subjects.strip().split(" ")[0])
print(num_subjects)

split_cmd ='split  -l ' + str(int(math.ceil(num_subjects/num_threads)))+ " "+ input_file +" "+in_dir+"/"  + in_file_basename[0:4]

os.system(split_cmd)

alpha = string.ascii_lowercase 
for i in range(num_threads):
    in_file = in_dir +"/"+in_file_basename[0:4] +"a" +  alpha[i] 
    print(in_file)
    output_file = output_dir + os.path.basename(in_file) + '_out'
    t = threading.Thread(target=imputation.impute_file, args = (in_file, output_file, priority, output_MUUG, 
                                                       output_haplotypes,1000,epsilon,number_of_results,
                                                      number_of_pop_results,planb, ))
    t.start()

# Profiler end
pr.disable()
pr.print_stats(sort="time")
