#
# read frequency file in hpf format and create 3 csv files for the graph
#
import argparse
import csv
import json
import pathlib
from collections import defaultdict, namedtuple
from itertools import combinations, count
from operator import indexOf, add

import sys
import os

# sys.path.insert(0, os.path.join(".."))


# set haplotype frequency trimming threshold
#
# The NEMO frequency generating sample is 6.93M and the largest single
# cohort is CAU 3,912,440 so a haplotype with a frequency of < 1/7,824,880
# would have an expected value < 1 in all  populations
# freq_trim = 1e-5  # 0.00001          1 in 100K
# freq_trim = 1.1278e-7  # 0.0000001278... 1 in 7.8M
# freq_trim = 1.66e-7  # 0.000000166...   1 in 6M
# freq_trim = 1e-6  # 0.000001         1 in 1M
# freq_trim = 1.66e-07


# If the trimming threshold is 1e-6 it will take 9m35s to generate the graph
# csv files  on a mid-2015 MacBook Pro (2.5 GHz Intel Core i7) and will result
# in 1,088,817 nodes (159MB),
# 14,868,976 edges (2.0GB)and 5,947,591 top links (108MB).

# If the trimming threshold is 1e-7 it will take 19m37 to generate the graph
# csv files on a mid-2015 MacBook Pro (2.5 GHz Intel Core i7) and will result
# in 1.8M nodes (262MB),
# 26M edges (3.6B)and 107M top links (200MB).
# FULL_LOCI = 'ABCQR'


##############################################################################
# functions
##############################################################################
# nCr: n Choose r
def nCr(nchars, r):
    return ["".join(x) for x in combinations(nchars, r)]


# dividebyzero: safe division
def dividebyzero(a, b):
    if b == 0:
        return 0
    else:
        return a / b


# make_allele_list: convert haplotype into sorted list of alleles
def make_allele_list(haplotype, full_name_index_dict, num_of_alleles):
    hl = haplotype.split("~")
    # Remove g from the allele names
    hl = list(map(lambda h: h[:-1] if (h[-1] == "g") else h, hl))
    sorted_h1 = ["0"] * num_of_alleles
    for allele in hl:
        locus = allele.split("*")[0]
        # print(allele, hl)
        sorted_h1[full_name_index_dict[locus] - 1] = allele
    return sorted_h1


# Find the index of a loci name in `FULL_LOCI`. ie. `BC` in `ABCRQ` is `[1,2]`
def find_indices(FULL_LOCI, loci):
    return list(map(lambda x: indexOf(FULL_LOCI, x), loci))


# Given a list of indices find the corresponding loci combo
def find_loci_name(FULL_LOCI, loci_indices):
    return "".join([FULL_LOCI[i] for i in sorted(loci_indices)])


# Find all the parents of the loci with loci_indices against loci_allele_list.
def find_parents(
    FULL_LOCI, loci_indices, loci_allele_list, p, ParentLoci, top_nodes_plan_b
):
    fi = find_indices(FULL_LOCI, FULL_LOCI)
    new_loci = set(fi).difference(loci_indices)
    parents = []
    for i in list(new_loci):
        new_loci = loci_indices.copy()
        new_loci.append(i)
        # print(new_loci)
        loci_name = find_loci_name(FULL_LOCI, new_loci)
        if loci_name in top_nodes_plan_b:
            # print(loci_name)
            haplotype = "~".join([loci_allele_list[i] for i in sorted(new_loci)])
            parents.append(ParentLoci(loci=loci_name, haplotype=haplotype, p=p))
    return parents


# Find the wanted nodes for plan A and plan B
def labels_for_grap(conf, full_loci, csvdir):

    nodes_graph = conf.get("Plan_A_Matrix", [])

    # Created all the node combinations if there is no limit on the requested nodes
    if not nodes_graph:
        list_all_combo = [full_loci]
        for i in range(len(full_loci) - 1, 0, -1):
            list_all_combo.extend(nCr(full_loci, i))
        return list_all_combo, list_all_combo, list_all_combo, list_all_combo

    nodes_plan_a = []
    nodes_plan_b = []  # nodes from which edges to the parents came out in plan b
    top_nodes_plan_b = (
        []
    )  # nodes to which edges will reach in plan b or single locus without requested parents.

    matrix_plan_b = conf.get(
        "Plan_B_Matrix",
        [
            [[1, 2, 3, 4, 5]],
            [[1, 2, 3], [4, 5]],
            [[1], [2, 3], [4, 5]],
            [[1, 2, 3], [4], [5]],
            [[1], [2, 3], [4], [5]],
            [[1], [2], [3], [4], [5]],
        ],
    )[1]
    # convert the first recombination event from loci indexes to loci names,
    # and add the event subcomponents to top_nodes_plan_b
    for subcomponent in matrix_plan_b:
        for idx, locus in enumerate(subcomponent):
            subcomponent[idx] = str(locus)
        top_nodes_plan_b.append(("").join(subcomponent))

    list_complement = []
    full_loci_list = []
    full_loci_list[:0] = full_loci
    for node_label in nodes_graph:
        # add nodes of plan a
        label = ""
        list_label = []
        for locus_idx in node_label:
            label += str(locus_idx)
            list_label.append(str(locus_idx))

        nodes_plan_a.append(label)
        complement_node = ("").join(
            sorted([item for item in full_loci_list if item not in list_label])
        )
        if complement_node != "":
            list_complement.append(complement_node)
        # find the sub-nodes of label for plan b
        for subcomponent in matrix_plan_b:
            node = ""
            for locus in label:
                if locus in subcomponent:
                    node += locus
            if node != "" and not node in top_nodes_plan_b:
                nodes_plan_b.append(node)

    # add the full-locus to plan a nodes
    if not full_loci in nodes_plan_a:
        nodes_plan_a.append(full_loci)

    # add the single locus to top_nodes_plan_b
    all_loci = list(full_loci)
    for locus in all_loci:
        if not locus in nodes_plan_b:
            top_nodes_plan_b.append(locus)

    # add the complementary nodes if they are not in one of the graphs
    for node in list_complement:
        if (
            (not node in nodes_plan_a)
            and (not node in nodes_plan_b)
            and (not node in top_nodes_plan_b)
        ):
            top_nodes_plan_b.append(node)

    nodes_plan_b = list(set(nodes_plan_b))
    all_combo_list = list(dict.fromkeys(nodes_plan_a + nodes_plan_b + top_nodes_plan_b))
    with open(csvdir + "nodes_for_plan_a.txt", "w") as f:
        for item in nodes_plan_a:
            f.write("%s\n" % item)
    with open(csvdir + "nodes_for_plan_b.txt", "w") as f:
        for item in nodes_plan_b + top_nodes_plan_b:
            f.write("%s\n" % item)
    # pickle.dump(nodes_plan_a, open(csvdir + '/nodes_for_plan_a.pkl', "wb"))
    # pickle.dump(nodes_plan_b + top_nodes_plan_b, open(csvdir + '/nodes_for_plan_b.pkl', "wb"))

    return all_combo_list, nodes_plan_a, nodes_plan_b, top_nodes_plan_b


# create full loci order,
# map beteen the loci neme and his index
def loci_order(loc_values):
    full_name_index_dict = {}
    all_loci_set = set()
    for locus, val in loc_values.items():
        full_name_index_dict[locus] = val
        all_loci_set.add(str(val))

    FULL_LOCI = "".join(sorted(all_loci_set))

    return FULL_LOCI, full_name_index_dict


def generate_graph(
    config_file="../conf/minimal-configuration-script.json",
    em_pop=None,
    em=False,
    use_default_path=False,
):
    ##############################################################################
    # Configure
    ##############################################################################
    # set output directory and create it if it doesn't exist
    # csvdir = "output/csv"

    # Input file
    # freq_file = path  + freq_file
    path = ""
    if use_default_path:
        path = os.path.dirname(os.path.realpath(__file__)) + "/"
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--config",
        required=False,
        default=config_file,
        help="Configuration JSON file",
        type=str,
    )

    args = parser.parse_args()
    configuration_file = args.config

    # Read configuration file and load properties
    with open(configuration_file) as f:
        conf = json.load(f)

    csvdir = conf.get("graph_files_path")
    pathlib.Path(csvdir).mkdir(parents=True, exist_ok=True)
    if csvdir[-1] != "/":
        csvdir += "/"

    pops = conf.get("populations")
    if em_pop:
        pops = em_pop
    freq_trim = conf.get("freq_trim_threshold")

    freq_file = path + conf.get("freq_file")
    dict_count_of_pop = {}

    pop_ratio_dir = path + conf.get("pops_count_file", "")
    path = pathlib.Path(pop_ratio_dir)

    if em or not path.is_file():
        for pop in pops:
            dict_count_of_pop[pop] = freq_trim
    else:
        with open(pop_ratio_dir) as f_count:
            for line in f_count:
                pop, count_pop, ratio = line.strip().split(",")
                dict_count_of_pop[pop] = freq_trim / float(count_pop)

    # Display the configurations we are using
    print(
        "****************************************************************************************************"
    )
    print("Performing graph generation based on following configuration:")
    print("\tPopulation: {}".format(pops))
    print("\tFreq File: {}".format(freq_file))
    print("\tFreq Trim Threshold: {}".format(freq_trim))
    print(
        "****************************************************************************************************"
    )

    ##############################################################################
    # setup graph elements
    ##############################################################################
    HaplotypeNode = namedtuple(
        "HaplotypeNode", "node_id, freq_list, allele_list, parents, top_links"
    )
    TopLink = namedtuple("TopLink", "node_id, haplotype")
    ParentLoci = namedtuple("ParentLoci", "loci, haplotype, p")

    header = [":START_ID(HAPLOTYPE)", ":END_ID(HAPLOTYPE)", "CP:DOUBLE[]", ":TYPE"]

    FULL_LOCI, full_name_index_dict = loci_order(conf.get("loci_map"))

    all_combo_list, nodes_plan_a, nodes_plan_b, top_nodes_plan_b = labels_for_grap(
        conf, FULL_LOCI, csvdir
    )
    sequence = count(0)

    loci_combo_map = {
        combo_list: defaultdict(
            lambda: HaplotypeNode(
                node_id=next(sequence),
                freq_list=[],
                allele_list=[],
                parents=[],
                top_links=set(),
            )
        )
        for combo_list in all_combo_list
    }

    # for key in loci_combo_map:
    #    print (key)

    # Load the frequency file into the FULL_LOCI dictionary with the cleaned up haplotype as the key.

    haplist_overall = {}  # list of haplotypes across all populations
    pop_hap_combos = {}

    num_of_alleles = len(FULL_LOCI)
    with open(freq_file) as f:
        for hap_line in f:
            if not hap_line:
                continue
            haplotype, pop, freq = hap_line.split(",")
            if haplotype == "hap":
                continue
            freq = float(freq)
            # Ignore lines with 0 freq
            if freq == 0.0:
                continue
            # trimming frequencies
            if freq < dict_count_of_pop[pop]:
                continue

            hap_list = make_allele_list(haplotype, full_name_index_dict, num_of_alleles)
            haplotype = "~".join(hap_list)
            pop_haplotype = pop + "-" + haplotype
            haplist_overall[haplotype] = 1
            pop_hap_combos[pop_haplotype] = freq

    for haplotype in haplist_overall:
        # make frequency array
        freqs = []
        for pop in pops:
            # check to see if population + haplotype combo exists
            pop_haplotype = pop + "-" + haplotype
            if pop_haplotype in pop_hap_combos:
                freqs.append(pop_hap_combos[pop_haplotype])
            else:
                freqs.append(0)
        hap_list = make_allele_list(haplotype, full_name_index_dict, num_of_alleles)
        loci_combo_map[FULL_LOCI][haplotype] = HaplotypeNode(
            node_id=next(sequence),
            freq_list=freqs,
            allele_list=hap_list,
            parents=None,
            top_links=set(),
        )

    # list(loci_combo_map[FULL_LOCI].items())[99]

    # For each loci block, get all the allele combinations. Create parents link, and accumulate the freqs to create new Haplotype Nodes.

    for loci in all_combo_list:
        if loci != FULL_LOCI:  # FULL_LOCI is already in the dictionary
            loci_map = loci_combo_map[loci]
            loci_indices = find_indices(FULL_LOCI, loci)
            for full_hap_name, full_hap_node in loci_combo_map[FULL_LOCI].items():
                allele_list = full_hap_node.allele_list
                haplotype = "~".join([allele_list[i] for i in loci_indices])
                haplotype_node = loci_map[haplotype]

                if (
                    loci in nodes_plan_b
                ):  # find parents just for internal nodes of plan B
                    parents = find_parents(
                        FULL_LOCI,
                        loci_indices,
                        allele_list,
                        full_hap_node.freq_list,
                        ParentLoci,
                        top_nodes_plan_b,
                    )
                    parents_list = haplotype_node.parents
                    for parent in parents:
                        parents_list.append(parent)
                else:
                    parents_list = []

                if loci in nodes_plan_a:  # top link just for nodes of plan A
                    top_links = haplotype_node.top_links
                    top_links.add(
                        TopLink(node_id=full_hap_node.node_id, haplotype=full_hap_name)
                    )
                else:
                    top_links = None

                # have to make freq_list with all zeros otherwise map add function returns empty list
                new_freq_list = []
                if not haplotype_node.freq_list:
                    for pop in pops:
                        new_freq_list.append(0)
                else:
                    new_freq_list = haplotype_node.freq_list
                freq_sum = list(map(add, new_freq_list, full_hap_node.freq_list))

                new_node = HaplotypeNode(
                    node_id=haplotype_node.node_id,
                    freq_list=freq_sum,
                    allele_list=None,
                    parents=parents_list,
                    top_links=top_links,
                )

                loci_map[haplotype] = new_node

    # #### Build Nodes file

    header = ["haplotypeId:ID(HAPLOTYPE)", "name", "loci:LABEL", "frequency:DOUBLE[]"]
    node_file = csvdir + conf.get("node_csv_file")
    with open(node_file, mode="w") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(header)
        for loci in all_combo_list:
            loci_map = loci_combo_map[loci]
            for haplotype, haplotype_node in loci_map.items():
                freq_array = ";".join(map(str, haplotype_node.freq_list))
                csv_writer.writerow(
                    [haplotype_node.node_id, haplotype, loci, freq_array]
                )

    # #### Build Edges File

    edgeheader = [":START_ID(HAPLOTYPE)", ":END_ID(HAPLOTYPE)", "CP:DOUBLE[]", ":TYPE"]
    edge_file = csvdir + conf.get("edges_csv_file")
    with open(edge_file, mode="w") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(edgeheader)

        for loci in nodes_plan_b:
            if loci != FULL_LOCI:  # FULL_LOCI is already in the dictionary
                loci_map = loci_combo_map[loci]
                for haplotype, haplotype_node in list(loci_map.items()):
                    for parent in haplotype_node.parents:
                        # assert(haplotype_node.freq > 0.0)
                        # avoid division by zero

                        cp = list(map(dividebyzero, parent.p, haplotype_node.freq_list))
                        prob_array = ";".join(map(str, cp))
                        loci_combo = parent.loci
                        hap = parent.haplotype
                        parent_id = loci_combo_map[loci_combo][hap].node_id
                        csv_writer.writerow(
                            [haplotype_node.node_id, parent_id, prob_array, "CP"]
                        )

    # #### Generate Top Links file

    topheader = [":START_ID(HAPLOTYPE)", ":END_ID(HAPLOTYPE)", ":TYPE"]
    top_links_file = csvdir + conf.get("top_links_csv_file")
    with open(top_links_file, mode="w") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(topheader)
        for loci in nodes_plan_a:
            if loci != FULL_LOCI:  # FULL_LOCI is already in the dictionary
                #
                loci_map = loci_combo_map[loci]
                for haplotype, haplotype_node in list(loci_map.items()):
                    top_links = haplotype_node.top_links
                    for top_link in top_links:
                        csv_writer.writerow(
                            [haplotype_node.node_id, top_link.node_id, "TOP"]
                        )

    # #### Generate Info Node file

    infonode_header = [
        "INFO_NODE_ID:ID(INFO_NODE)",
        "populations:STRING[]",
        "INFO_NODE:LABEL",
    ]
    top_links_file = csvdir + conf.get("info_node_csv_file")
    with open(top_links_file, mode="w") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(infonode_header)
        csv_writer.writerow([1, ";".join(pops), "INFO_NODE"])


if __name__ == "__main__":
    generate_graph()
