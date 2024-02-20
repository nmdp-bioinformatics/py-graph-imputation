import copy
import logging
import math
import operator
import timeit
from collections import defaultdict
import os.path
import json


import numpy as np
from .cutils import open_ambiguities, create_hap_list, deepcopy_list
from .cypher_plan_b import CypherQueryPlanB
from .cypher_query import CypherQuery

comp_cand_epsilon = 1e-15


def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i : i + n]


def write_best_prob(name_gl, res, probs, numOfResult, fout, sign=","):
    sumProbsDict = defaultdict(list)
    # loop over the result and sum the prob by populations/haplotype
    for k in range(len(res)):
        key = res[k][0] + sign + res[k][1]
        if key in sumProbsDict:
            sumProb = probs[k] + sumProbsDict[key]
            sumProbsDict[key] = sumProb
        else:
            key2 = res[k][1] + sign + res[k][0]
            if key2 in sumProbsDict:
                sumProb = probs[k] + sumProbsDict[key2]
                sumProbsDict[key2] = sumProb
            else:
                sumProbsDict[key] = probs[k]

    multProbs = []
    for k in sumProbsDict:
        multProbs.append([sumProbsDict[k], [k, sumProbsDict[k]]])

    multProbs.sort(key=lambda x: x[0], reverse=True)

    # write the output to file
    minBestResult = min(numOfResult, len(multProbs))
    for k in range(minBestResult):
        fout.write(
            name_gl
            + ","
            + str(multProbs[k][1][0])
            + ","
            + str(multProbs[k][0])
            + ","
            + str(k)
            + "\n"
        )


def write_best_prob_genotype(name_gl, res, numOfResult, fout):
    sorted_by_value = sorted(res.items(), key=lambda kv: kv[1], reverse=True)

    # write the output to file
    minBestResult = min(numOfResult, len(sorted_by_value))
    for k in range(minBestResult):
        fout.write(
            name_gl
            + ","
            + str(sorted_by_value[k][0])
            + ","
            + str(sorted_by_value[k][1])
            + ","
            + str(k)
            + "\n"
        )


def write_best_hap_race_pairs(name_gl, haps, pops, probs, numOfResult, fout):
    all_res = []

    for i in range(len(probs)):
        pair = haps[i][0] + ";" + pops[i][0] + "," + haps[i][1] + ";" + pops[i][1]
        all_res.append([probs[i], pair])
    all_res.sort(key=lambda x: x[0], reverse=True)

    # write the output to file
    minBestResult = min(numOfResult, len(all_res))
    for k in range(minBestResult):
        fout.write(
            name_gl
            + ","
            + str(all_res[k][1])
            + ","
            + str(all_res[k][0])
            + ","
            + str(k)
            + "\n"
        )


#
# given a gl string, remove g,L and UUUUs from it
#
def clean_up_gl(gl):
    gl = gl.replace("g", "")
    # keep nulls
    # gl = gl.replace('N', '')
    gl = gl.replace("L", "")
    unknowns = []
    locus_gl = gl.split("^")
    for j in range(len(locus_gl)):
        gen = locus_gl[j]
        if not gen.strip("UUUU") == gen:
            unknowns.append(gen)
    for miss in unknowns:
        locus_gl.remove(miss)
    return "^".join(locus_gl)


class Imputation(object):
    __slots__ = (
        "logger",
        "verbose",
        "populations",
        "netGraph",
        "priorMatrix",
        "full_hapl",
        "index_dict",
        "full_loci",
        "factor",
        "_factor_missing_data",
        "cypher",
        "cypher_plan_b",
        "matrix_planb",
        "count_by_prob",
        "number_of_options_threshold",
        "plan",
        "option_1",
        "option_2",
        "haplotypes_number_in_phase",
        "save_space_mode",
        "nodes_for_plan_A",
        "unk_priors",
    )

    def __init__(self, net=None, config=None, count_by_prob=None, verbose=False):
        """Constructor
        Intialize an instance of `Imputation` with a py2neo graph
        and a `CypherQuery` object.
        """
        # graph is a link to the graph and populations is the list of populations
        self.logger = logging.getLogger("Logger." + __name__)
        self.verbose = verbose
        #    self.graph = graph
        if not config is None:

            # can not let populations be set on construction
            # because the code assumes it is always the same length as the
            # arrays in the graph
            self.populations = config["pops"]
            self.netGraph = net
            self.priorMatrix = np.ones((len(self.populations), len(self.populations)))
            self.unk_priors = config["UNK_priors"]

            # For plan b
            # self.full_loci = config["full_loci"]
            # loc_values = list(set(config["loci_map"].values()))
            # loc_values.sort()
            # indxes_loc_dict = {}
            # indxes_loc_dict_planb = {}
            # for i, loc in enumerate(loc_values):
            # indxes_loc_dict[loc] = i+1
            # indxes_loc_dict_planb[i+1] = loc
            self.full_hapl = []
            # self.index_dict = {}

            # for loci in config["loci_map"]:
            # self.full_hapl.append(loci)
            # self.index_dict[loci] = indxes_loc_dict[config["loci_map"][loci]]

            loc_values = config["loci_map"]
            #####################
            indxes_loc_dict_planb = {}
            self.index_dict = {}
            # short_name_index_dict = {}
            for locus, val in loc_values.items():
                self.full_hapl.append(locus)
                self.index_dict[locus] = val  # full_name_index_dict
                # indxes_loc_dict_planb[val[1]] = val[0] #indxes_loc_short_name_dict

            self.full_loci = config["full_loci"]

            ##############################

            self.factor = 0.0001
            self._factor_missing_data = config["factor_missing_data"]
            self.cypher = CypherQuery(config["loci_map"])
            self.cypher_plan_b = CypherQueryPlanB(
                config["loci_map"]
            )  # , indxes_loc_dict_planb)

            self.matrix_planb = config["matrix_planb"]

            if count_by_prob is None:
                self.count_by_prob = np.ones(len(self.populations))
                if config["use_pops_count_file"]:
                    with open(config["pops_count_file"]) as f_count:
                        for i, line in enumerate(f_count):
                            self.count_by_prob[i] = float(line.strip().split(",")[2])
            else:
                self.count_by_prob = count_by_prob
            self.number_of_options_threshold = config["number_of_options_threshold"]

            # The following is not reliable in parallel mode.
            self.plan = "a"
            self.option_1 = 0
            self.option_2 = 0
            self.haplotypes_number_in_phase = config["max_haplotypes_number_in_phase"]
            self.save_space_mode = config["save_mode"]
            self.nodes_for_plan_A = config["nodes_for_plan_A"]

    def print_options_count(self, subject_id):
        message = (
            "Subject: {id} plan: {plan}, open_phases - count of open regular option: {option1}, "
            "count of alternative opening: {option2} "
        )
        print(
            message.format(
                id=subject_id,
                plan=self.plan,
                option1=self.option_1,
                option2=self.option_2,
            )
        )

    def power_find(self, n):
        """produces all powers of 2"""
        result = []
        binary = bin(n)[:1:-1]
        for x in range(len(binary)):
            if int(binary[x]):
                result.append(x)
        return result

    def gl2haps(self, GL_String):
        # Receives a GL string adn produces a genotype in list structure
        if GL_String == "" or GL_String == " ":
            return []
        split_hap = GL_String.split("^")
        N_Loci = len(split_hap)
        t1 = []
        t2 = []
        count = 0
        for i in range(N_Loci):
            # if not split_hap[i]=='+':
            if split_hap[i][0] == "+":
                split_hap[i] = split_hap[i][1:]
            curr_locus = split_hap[i].split("+")
            if len(curr_locus) == 1:
                if curr_locus == [""]:
                    count = count + 1
                    continue
                else:
                    return []
            t1.append(curr_locus[0])
            t2.append(curr_locus[1])
        for i in range(count):
            split_hap.remove("")
            N_Loci = N_Loci - 1
        Gen = [sorted(t1), sorted(t2)]
        return {"Genotype": Gen, "N_Loc": N_Loci}

    def gen_phases(self, gen, n_loci, b_phases):
        # haplotypes in gen is expected to be sorted
        # Generates all phases, but does not handle locus ambiguities
        if b_phases is not None:
            phase_indices = [i for i, e in enumerate(b_phases) if e == 1]
        Phases = []
        N_Phases = 2 ** (n_loci - 1)  # Total Number of phases
        exists = {}
        for i in range(0, N_Phases):
            H1 = []  # Hap lists
            H2 = []
            M1 = self.power_find(i)  # find all the powers of 2 in i
            L = [0] * n_loci  # Initialize to 0 for all loci
            for m in M1:
                L[m] = 1
                if (b_phases is not None) and (m not in phase_indices):
                    L[m] = 0
                # take a phase and set it to 1,
                # all others go to the other phase.
            for k in range(n_loci):
                H1.append(gen[L[k]][k])
                H2.append(gen[1 - L[k]][k])
            geno1 = "^".join(["~".join(H1), "~".join(H2)])
            geno2 = "^".join(["~".join(H2), "~".join(H1)])
            if (geno1 not in exists) or (geno2 not in exists):
                exists[geno1] = 1
                exists[geno2] = 1
                Phases.append([H1, H2])

        return Phases

    def open_gl_string(self, gl_string, cutoff):
        # receives a list of phases and computes haps and
        # probabilties and accumulate cartesian productEpsilon=0.0001
        chr = self.gl2haps(gl_string)
        if chr == []:
            return

        # generate the 2^(N-1) phases
        chr1 = self.gen_phases(chr["Genotype"], chr["N_Loc"], None)

        # return if the result is empty (why would that be?)
        if chr1 == []:
            return

        phases = self.open_phases_for_em(chr1, chr["N_Loc"], cutoff)
        return phases

    def open_phases_for_em(self, haps, N_Loc, cutoff):
        phases = []
        for j in range(len(haps)):
            H1 = []
            H2 = []

            for k in range(2):
                hap_list = [haps[j][k]]
                hap_list_splits = [tuple(allele.split("/")) for allele in hap_list[0]]
                # compute the number of options:
                options = 1
                for i in range(N_Loc):
                    options *= len(hap_list_splits[i])

                # if the number of options is smaller than the total number of nodes:
                if options < cutoff:  # open ambiguities regularly:
                    for i in range(N_Loc):
                        hap_list = self.open_ambiguities(
                            hap_list, i, hap_list_splits[i]
                        )
                    if k == 0:
                        H1.append(hap_list)
                    else:
                        H2.append(hap_list)

                # if there are more options than actual haplotypes possible:

            if H1 and H2:
                phases.append([H1, H2])
        return phases

    def comp_hap_prob(self, Hap, N_Loc, epsilon, n):
        haplo_probs = self.get_haplo_freqs(Hap, epsilon, n)
        probs = list(haplo_probs.values())
        # probs=haplo_probs.values()
        haplos = list(haplo_probs.keys())
        if not haplo_probs:
            return {"Haps": "", "Probs": ""}
        return {"Haps": haplos, "Probs": probs}

    # get_haplo_freqs
    # given a set of possible haplotypes, look them up
    # def get_haplo_freqs(self, haplos, epsilon, n=25000):
    #     haplo_probs = {}
    #     haplos_joined = ["~".join(sorted(item)) for sublist in haplos for item in sublist]
    #
    #     # Break up the haplotype list into
    #     # chunks of 25,000 if the total
    #     # number of haplotypes is greater
    #     # than 25,000.
    #     #
    #     # - Is it faster to break up into many
    #     # smaller chunks, or few large
    #
    #     if len(haplos_joined) > n:
    #         i=0
    #         for chunk in chunks(haplos_joined, n):
    #             haplo_query = self.cypher.buildQuery(chunk)
    #             fq = pa.DataFrame(self.graph.data(haplo_query))
    #             i = i+1
    #             if not fq.empty:
    #                 freq1_dic = fq.set_index('abcqr.name')['abcqr.frequency'].to_dict()
    #                 haplo_probs.update(freq1_dic)
    #     else:
    #         haplo_query1 = self.cypher.buildQuery(haplos_joined)
    #         fq = pa.DataFrame(self.graph.data(haplo_query1))
    #         if not fq.empty:
    #             freq1_dic = fq.set_index('abcqr.name')['abcqr.frequency'].to_dict()
    #             haplo_probs.update(freq1_dic)
    #     return haplo_probs

    def get_haplo_freqs(self, haplos, epsilon, n=25000):
        haplos_joined = ["~".join(item) for item in haplos[0]]  ###
        # haplos_joined = [item for sublist in haplos for item in sublist]  ###
        # haplos_joined = ["~".join(sorted(item)) for sublist in haplos for item in sublist]
        return self.netGraph.adjs_query(haplos_joined)

    # def get_haplo_freqs_miss(self, haplos, epsilon):
    #     haplo_probs = {}
    #     haplos_joined = ["~".join(sorted(item)) for sublist in haplos for item in sublist]
    #     haplo_query1 = self.cypher.buildQuery(haplos_joined)
    # #    fq = pa.DataFrame(self.graph.data(haplo_query1))
    #     if not fq.empty:
    #         freq1_dic = fq.set_index('abc.name')['abc.frequency'].to_dict()
    #         haplo_probs.update(freq1_dic)
    #     return haplo_probs

    def cal_prob(self, probs1, probs2, epsilon):
        # This is the part where we loop over all race combinations.
        # N*N loop, where N is the number of races/
        places = []
        # print ("cal_prob\n\n")
        # print ("probs1 = ", probs1)
        # print ("probs2 = ", probs2)
        for i in range(len(probs1)):
            for j in range(len(probs2)):
                if probs1[i] * probs2[j] >= epsilon and probs1[i] * probs2[j] != 0:
                    places.append([i, j])
        # places are the indices of the positions in the population array.
        return places

    # move to one dim sorted list, and remove prob=0
    def convert_list_to_one_dim(self, prob):
        ProbWithIndexes = []

        prob_with_indexes_by_prior = []
        for k in range(len(prob)):
            for j in range(len(prob[k])):
                if prob[k][j] > 0:
                    prob_with_indexes_by_prior.append(
                        [(prob[k][j] * self.priorMatrix[j][j]), [prob[k][j], [k, j]]]
                    )
        prob_with_indexes_by_prior.sort(key=lambda x: x[0], reverse=True)

        minBestResult = min(
            self.haplotypes_number_in_phase, len(prob_with_indexes_by_prior)
        )
        for k in range(minBestResult):
            ProbWithIndexes.append(prob_with_indexes_by_prior[k][1])

        return ProbWithIndexes

    def calc_haps_pairs(
        self,
        Haps1,
        Haps2,
        Prob1WithIndexes,
        Prob2WithIndexes,
        epsilon,
        hap_total,
        pop_res,
        maxProb,
        geno_seen,
    ):

        for h in range(len(Prob1WithIndexes)):
            x = epsilon / Prob1WithIndexes[h][0]
            x_homozygote = x * 2
            prob2Len = len(Prob2WithIndexes)
            #  k=0
            # while  k < prob2Len and Prob2WithIndexe s[k][0] >= x:
            for k in range(prob2Len):
                if Prob2WithIndexes[k][0] >= x:
                    if (
                        self.priorMatrix[Prob1WithIndexes[h][1][1]][
                            Prob2WithIndexes[k][1][1]
                        ]
                        > 0
                    ):
                        hap1 = Haps1[Prob1WithIndexes[h][1][0]]
                        hap2 = Haps2[Prob2WithIndexes[k][1][0]]
                        if (
                            hap1 != hap2
                            and (
                                self.priorMatrix[Prob1WithIndexes[h][1][1]][
                                    Prob2WithIndexes[k][1][1]
                                ]
                                * Prob2WithIndexes[k][0]
                            )
                            >= x
                        ) or (
                            hap1 == hap2
                            and (
                                self.priorMatrix[Prob1WithIndexes[h][1][1]][
                                    Prob2WithIndexes[k][1][1]
                                ]
                                * Prob2WithIndexes[k][0]
                            )
                            >= x_homozygote
                        ):

                            # hap1_list = sorted(hap1.split('~'))
                            # hap2_list = sorted(hap2.split('~'))
                            # genotype = ("~".join(sorted(hap1.split("~") + hap2.split("~"))))
                            # if hap1_list[-1].split("*")[0] == hap2_list[-1].split("*")[0]:
                            genotype = "^".join(
                                map(
                                    lambda x: "+".join(sorted(x)),
                                    zip(
                                        sorted(hap1.split("~")), sorted(hap2.split("~"))
                                    ),
                                )
                            )
                            race1 = self.populations[Prob1WithIndexes[h][1][1]]
                            race2 = self.populations[Prob2WithIndexes[k][1][1]]

                            h1_id = hap1 + "," + race1
                            h2_id = hap2 + "," + race2
                            geno_id = "-".join(sorted([h1_id, h2_id]))
                            # geno_id = tuple(sorted([h1_id, h2_id]))
                            if geno_id not in geno_seen:
                                geno_seen.add(geno_id)

                                prob = (
                                    Prob1WithIndexes[h][0]
                                    * Prob2WithIndexes[k][0]
                                    * self.priorMatrix[Prob1WithIndexes[h][1][1]][
                                        Prob2WithIndexes[k][1][1]
                                    ]
                                )

                                if hap1 != hap2:
                                    prob = prob * 2

                                if prob > maxProb:
                                    maxProb = prob

                                if genotype in hap_total:
                                    sumProb = hap_total[genotype] + prob
                                    hap_total[genotype] = sumProb
                                else:
                                    hap_total[genotype] = prob

                                races = sorted([race1, race2])

                                races = races[0] + "," + races[1]

                                if races in pop_res:
                                    sumProb = pop_res[races] + prob
                                    pop_res[races] = sumProb
                                else:
                                    pop_res[races] = prob

                else:
                    break

        return maxProb

    def calc_haps_pairs_haplotype(
        self,
        Haps1,
        Haps2,
        Prob1WithIndexes,
        Prob2WithIndexes,
        epsilon,
        hap_total,
        pop_res,
        maxProb,
        haps_pairs,
        geno_seen,
        pop_res_haplo,
        p_total,
    ):

        for h in range(len(Prob1WithIndexes)):
            x = epsilon / Prob1WithIndexes[h][0]
            x_homozygote = x * 2
            prob2Len = len(Prob2WithIndexes)
            #  k=0
            # while  k < prob2Len and Prob2WithIndexe s[k][0] >= x:
            for k in range(prob2Len):
                if Prob2WithIndexes[k][0] >= x:
                    if (
                        self.priorMatrix[Prob1WithIndexes[h][1][1]][
                            Prob2WithIndexes[k][1][1]
                        ]
                        > 0
                    ):
                        hap1 = Haps1[Prob1WithIndexes[h][1][0]]
                        hap2 = Haps2[Prob2WithIndexes[k][1][0]]
                        if (
                            hap1 != hap2
                            and (
                                self.priorMatrix[Prob1WithIndexes[h][1][1]][
                                    Prob2WithIndexes[k][1][1]
                                ]
                                * Prob2WithIndexes[k][0]
                            )
                            >= x
                        ) or (
                            hap1 == hap2
                            and (
                                self.priorMatrix[Prob1WithIndexes[h][1][1]][
                                    Prob2WithIndexes[k][1][1]
                                ]
                                * Prob2WithIndexes[k][0]
                            )
                            >= x_homozygote
                        ):

                            # hap1_list = sorted(hap1.split('~'))
                            # hap2_list = sorted(hap2.split('~'))
                            # genotype = ("~".join(sorted(hap1.split("~") + hap2.split("~"))))
                            # if hap1_list[-1].split("*")[0] == hap2_list[-1].split("*")[0]:
                            genotype = "~".join(
                                sorted(hap1.split("~") + hap2.split("~"))
                            )
                            race1 = self.populations[Prob1WithIndexes[h][1][1]]
                            race2 = self.populations[Prob2WithIndexes[k][1][1]]

                            h1_id = hap1 + "," + race1
                            h2_id = hap2 + "," + race2
                            geno_id = "-".join(sorted([h1_id, h2_id]))
                            # geno_id = tuple(sorted([h1_id, h2_id]))
                            # if len(geno_seen) > 25000000:
                            # k = 1 / 0
                            if geno_id not in geno_seen:
                                geno_seen.add(geno_id)

                                prob = (
                                    Prob1WithIndexes[h][0]
                                    * Prob2WithIndexes[k][0]
                                    * self.priorMatrix[Prob1WithIndexes[h][1][1]][
                                        Prob2WithIndexes[k][1][1]
                                    ]
                                )

                                if hap1 != hap2:
                                    prob = prob * 2

                                if prob > maxProb:
                                    maxProb = prob

                                if genotype in hap_total:
                                    sumProb = hap_total[genotype] + prob
                                    hap_total[genotype] = sumProb
                                else:
                                    hap_total[genotype] = prob

                                races = sorted([race1, race2])

                                races = races[0] + "," + races[1]

                                if races in pop_res:
                                    sumProb = pop_res[races] + prob
                                    pop_res[races] = sumProb
                                else:
                                    pop_res[races] = prob

                                haps_pairs.append([hap1, hap2])
                                pop_res_haplo.append([race1, race2])
                                p_total.append(prob)

                else:
                    break

        return maxProb

    def comp_phase_prob_haplotype(self, phases, N_Loc, epsilon, n):
        # receives a list of phases and computes haps and
        # probabilties and accumulate cartesian product
        pop_res_haplo = []
        p_total = []
        geno_seen = set([])
        genotypes_total = {}
        hap_total = []
        pop_res = {}
        Prob2 = []
        Haps2 = []
        maxProb = 0
        for i in range(len(phases)):
            P1 = self.comp_hap_prob(phases[i][0], N_Loc, epsilon, n)
            # This will open locus ambiguities and comp probabilities for Hap1
            Haps1 = P1["Haps"]
            Prob1 = P1["Probs"]
            # if the first haplotype returned something then look at 2nd
            if len(Prob1) > 0:
                P2 = self.comp_hap_prob(phases[i][1], N_Loc, epsilon, n)
                # This will do the same for Hap 2;
                Haps2 = P2["Haps"]
                Prob2 = P2["Probs"]

            Prob1WithIndexes = self.convert_list_to_one_dim(Prob1)

            Prob2WithIndexes = self.convert_list_to_one_dim(Prob2)

            maxProb = self.calc_haps_pairs_haplotype(
                Haps1,
                Haps2,
                Prob1WithIndexes,
                Prob2WithIndexes,
                epsilon,
                genotypes_total,
                pop_res,
                maxProb,
                hap_total,
                geno_seen,
                pop_res_haplo,
                p_total,
            )

        # p_total returns an array of N*N (N is number of populations), hap_total - pairs of haplotypes.
        # pop_res are the names of the populations
        # return {'MaxProb': maxProb, 'Haps': hap_total, 'Pops': pop_res}

        return {
            "MaxProb": maxProb,
            "Haps": hap_total,
            "Probs": p_total,
            "Pops": pop_res_haplo,
        }

    def comp_phase_prob_genotype(self, phases, N_Loc, epsilon, n):
        # receives a list of phases and computes haps and
        # probabilties and accumulate cartesian product
        geno_seen = set([])
        hap_total = {}
        p_total = {}
        pop_res = {}
        Prob2 = []
        Haps2 = []
        maxProb = 0
        for i in range(len(phases)):
            P1 = self.comp_hap_prob(phases[i][0], N_Loc, epsilon, n)
            # This will open locus ambiguities and comp probabilities for Hap1
            Haps1 = P1["Haps"]
            Prob1 = P1["Probs"]
            # if the first haplotype returned something then look at 2nd
            if len(Prob1) > 0:
                P2 = self.comp_hap_prob(phases[i][1], N_Loc, epsilon, n)
                # This will do the same for Hap 2;
                Haps2 = P2["Haps"]
                Prob2 = P2["Probs"]

            Prob1WithIndexes = self.convert_list_to_one_dim(Prob1)

            Prob2WithIndexes = self.convert_list_to_one_dim(Prob2)

            maxProb = self.calc_haps_pairs(
                Haps1,
                Haps2,
                Prob1WithIndexes,
                Prob2WithIndexes,
                epsilon,
                hap_total,
                pop_res,
                maxProb,
                geno_seen,
            )

        # p_total returns an array of N*N (N is number of populations), hap_total - pairs of haplotypes.
        # pop_res are the names of the populations
        return {"MaxProb": maxProb, "Haps": hap_total, "Pops": pop_res}

    def comp_phase_prob(self, phases, N_Loc, epsilon, n):
        # receives a list of phases and computes haps and
        # probabilties and accumulate cartesian product
        geno_seen = set([])
        hap_total = []
        p_total = []
        pop_res = []
        Prob2 = []
        for i in range(len(phases)):
            P1 = self.comp_hap_prob(phases[i][0], N_Loc, epsilon, n)
            # This will open locus ambiguities and comp probabilities for Hap1
            Haps1 = P1["Haps"]
            Prob1 = P1["Probs"]
            # if the first haplotype returned something then look at 2nd
            if len(Prob1) > 0:
                P2 = self.comp_hap_prob(phases[i][1], N_Loc, epsilon, n)
                # This will do the same for Hap 2;
                Haps2 = P2["Haps"]
                Prob2 = P2["Probs"]

            Prob1WithIndexes = self.convert_list_to_one_dim(Prob1)

            Prob2WithIndexes = self.convert_list_to_one_dim(Prob2)

            #   for h in range(len(Prob1)):
            #       for k in range(len(Prob2)):
            for h in range(len(Prob1WithIndexes)):
                x = epsilon / Prob1WithIndexes[h][0]
                prob2Len = len(Prob2WithIndexes)
                #  k=0
                # while  k < prob2Len and Prob2WithIndexes[k][0] >= x:
                for k in range(prob2Len):
                    if Prob2WithIndexes[k][0] >= x:
                        if (
                            self.priorMatrix[Prob1WithIndexes[h][1][1]][
                                Prob2WithIndexes[k][1][1]
                            ]
                            > 0
                        ):
                            if (
                                self.priorMatrix[Prob1WithIndexes[h][1][1]][
                                    Prob2WithIndexes[k][1][1]
                                ]
                                > 0
                            ):
                                if (
                                    self.priorMatrix[Prob1WithIndexes[h][1][1]][
                                        Prob2WithIndexes[k][1][1]
                                    ]
                                    * Prob2WithIndexes[k][0]
                                ) >= x:
                                    # NOTE: there is a bug in places
                                    #  places = self.cal_prob(Prob1[h], Prob2[k], epsilon)
                                    # print ("places= ", places)
                                    #  for i in range(len(places)):
                                    # Prob1WithIndexes[h][0]
                                    # p=(places[i])
                                    # avoid reporting the same haplotype pair more than once
                                    #
                                    # crashing at the next line
                                    # p[0] = 3 but there only one population
                                    # "IndexError: list index out of range"
                                    #
                                    # this is due to an issue where the code expects
                                    # the population list to be the same as the graph
                                    #
                                    #
                                    # print ("p[0] = ", p[0])
                                    # print ("p[1] = ", p[1])
                                    h1_id = (
                                        Haps1[Prob1WithIndexes[h][1][0]],
                                        self.populations[Prob1WithIndexes[h][1][1]],
                                    )
                                    h2_id = (
                                        Haps2[Prob2WithIndexes[k][1][0]],
                                        self.populations[Prob2WithIndexes[k][1][1]],
                                    )
                                    geno_id = tuple(sorted([h1_id, h2_id]))
                                    if geno_id not in geno_seen:
                                        geno_seen.add(geno_id)
                                        hap_total.append([geno_id[0][0], geno_id[1][0]])
                                        pop_res.append([geno_id[0][1], geno_id[1][1]])
                                        # record prob in same order as associated hap & pop
                                        # for example, in the WMDA imputation output, this:
                                        # D000001,A*02:01~B*44:03~C*02:02~DQB1*02:01~DRB1*07:01,0.00009000,CAU,A*11:01~B*13:02~C*06:02~DQB1*05:01~DRB1*01:01,0.00006000,CAU
                                        # should be this:
                                        # D000001,A*02:01~B*44:03~C*02:02~DQB1*02:01~DRB1*07:01,0.00006000,CAU,A*11:01~B*13:02~C*06:02~DQB1*05:01~DRB1*01:01,0.00009000,CAU
                                        if geno_id[0] == h1_id:
                                            p_total.append(
                                                [
                                                    Prob1WithIndexes[h][0],
                                                    Prob2WithIndexes[k][0],
                                                ]
                                            )
                                        else:
                                            p_total.append(
                                                [
                                                    Prob2WithIndexes[k][0],
                                                    Prob1WithIndexes[h][0],
                                                ]
                                            )
                        #   k+=1
                    else:
                        break
        # p_total returns an array of N*N (N is number of populations), hap_total - pairs of haplotypes.
        # pop_res are the names of the populations
        return {"Haps": hap_total, "Probs": p_total, "Pops": pop_res}

    def reduce_phase_to_valid_allels(self, haps, N_Loc, planc=False):
        for j in range(len(haps)):
            for k in range(2):
                hap_list = []
                hap_list.append(haps[j][k])

                options = 1
                for i in range(N_Loc):
                    options = options * (len(hap_list[0][i].split("/")))
                if options >= self.number_of_options_threshold or planc:
                    for hap_k in hap_list:
                        for i, g in enumerate(hap_k):
                            gen = g.split("/")
                            probs = self.check_if_alleles_exist(gen)
                            if not probs == {}:
                                haps[j][k][i] = "/".join(list(probs.keys()))

    def reduce_phase_to_commons_alleles(
        self, haps, N_Loc, commons_number=1, planc=False
    ):
        for j in range(len(haps)):
            for k in range(2):
                hap_list = []
                hap_list.append(haps[j][k])

                options = 1
                for i in range(N_Loc):
                    options = options * (len(hap_list[0][i].split("/")))
                if options >= self.number_of_options_threshold or planc:
                    for hap_k in hap_list:
                        for i, g in enumerate(hap_k):
                            gen = g.split("/")
                            probs = self.check_if_alleles_exist(gen)
                            if not probs == {}:
                                dict_allels = {}
                                for allele in probs:
                                    list_allel = probs[allele]
                                    sum = 0
                                    for p, prob in enumerate(list_allel):
                                        sum += prob * self.priorMatrix[p, p]
                                    dict_allels[allele] = sum
                                commons_alleles = dict(
                                    sorted(
                                        dict_allels.items(),
                                        key=lambda kv: kv[1],
                                        reverse=True,
                                    )[:commons_number]
                                )
                                haps[j][k][i] = ("/").join(list(commons_alleles.keys()))

    def open_phases(self, haps, N_Loc, gl_string):
        phases = []
        for j in range(len(haps)):
            H1 = []
            H2 = []
            ##  fq = pa.DataFrame()
            fq = []

            for k in range(2):
                hap_list = [haps[j][k]]
                hap_list_splits = [tuple(allele.split("/")) for allele in hap_list[0]]

                # compute the number of options:
                options = 1
                for i in range(N_Loc):
                    options *= len(hap_list_splits[i])

                # if the number of options is smaller than the total number of nodes:
                if options < self.number_of_options_threshold:
                    # open ambiguities regularly:
                    for i in range(N_Loc):
                        hap_list = self.open_ambiguities(
                            hap_list, i, hap_list_splits[i]
                        )

                    if k == 0:
                        H1.append(hap_list)
                    else:
                        H2.append(hap_list)

                    self.option_1 += 1

                # if there are more options than actual haplotypes possible:
                else:
                    self.option_2 += 1
                    optionDict = {}  # set()
                    if len(fq) == 0:
                        _list = []
                        for gen, name in self.cypher.loc_map.items():
                            count = 0
                            for i in range(len(hap_list[0])):
                                if hap_list[0][i].split("*", 1)[0] == gen:
                                    count = count + 1
                            if count > 0:
                                _list.append(name)
                        # we'll get all the options possible
                        # (query,lc)=self.cypher.buildQuery(["~".join(_list)])

                        # fq = pa.DataFrame(self.graph.data(query))
                        label = "".join(sorted(_list))
                        fq = self.netGraph.haps_by_label(label)
                        # fq = pa.DataFrame(self.netGraph.abcdq_allele(), )
                        # fq = self.netGraph.abcdq_allele()
                        # we'll find which of all the options are compatable with the donor
                        """gl_list=gl_string.split("^")
                        for gen in gl_list:
                            gens = gen.split('+')
                            for g in gens:
                                gen=g.split('/')
                                for option in gen:
                                    optionDict[option]=True"""

                    for gen in hap_list_splits:
                        for option in gen:
                            optionDict[option] = True
                            # optionDict.add(option)
                    ##all_haps = fq.values.tolist()
                    hap_list = create_hap_list(fq, optionDict, N_Loc)

                    if k == 0:
                        H1.append(hap_list)
                    else:
                        H2.append(hap_list)
            if H1[0] and H2[0]:
                phases.append([H1, H2])
        return phases

    def find_type_loci(self, loci):
        return (loci.split("*"))[0]

    def find_missing_indexes(self, haplo):

        locus_in_haplo = []
        index = []
        for locus in haplo:
            locus_in_haplo.append(self.index_dict[self.find_type_loci(locus)])
        for locus in self.full_hapl:
            locus = self.index_dict[locus]
            if locus not in locus_in_haplo:
                # Find what index this is
                if not locus in index:
                    index.append(locus)
        return index

    def open_dict_data(self, haplo_probs):
        probs = list(haplo_probs.values())
        haplos = list(haplo_probs.keys())
        if not haplo_probs:
            return {"Haps": "", "Probs": ""}
        return {"Haps": haplos, "Probs": probs}

    def create_haplos_string(self, haplos, division, missing):
        string_haplos = []
        wanted_haplos = haplos[0]
        # Go over all the hap
        for hap in wanted_haplos:
            # The first item
            hap_string = ""
            # Create the wanted string
            for i in range(0, len(division)):
                place = division[i]
                # If this is not a missing option
                if division[i] not in missing:
                    # Normalize the places
                    for miss in missing:
                        if division[i] > miss:
                            place = place - 1

                    # [division[i]]- to get the index for the item in haplox we want
                    if hap_string == "":
                        hap_string = str(hap[place - 1])
                    else:
                        hap_string = hap_string + "~" + str(hap[place - 1])
            if not hap_string == "":
                string_haplos.append(hap_string)
        return string_haplos

    def open_option_(self, dict2, dict1, planc=False, num_of_options=10):
        # Will contain the two option combine
        dict_all = {}
        # The size of the vector of freq will always be equal to size of the population
        size = len(self.populations)
        if planc:
            size = 1
        if self.save_space_mode:
            for dict_ in [dict1, dict2]:
                if len(dict_) > num_of_options:
                    dict_tmp = {}
                    for hap in dict_:
                        dict_tmp[hap] = sum(dict_[hap])

                    dict_tmp = sorted(dict_tmp.items(), key=lambda kv: kv[1])
                    while len(dict_) > num_of_options:
                        del dict_[dict_tmp[0][0]]
                        del dict_tmp[0]
                    del dict_tmp
        # Go over all the option of the combination
        for key1 in dict1:
            for key2 in dict2:
                freq1 = dict1[key1]
                freq2 = dict2[key2]
                list_prob = [freq1[i] * freq2[i] * self.factor for i in range(size)]
                if max(list_prob) > 0:
                    key = ("~").join(sorted(key1.split("~") + key2.split("~")))
                    dict_all[key] = list_prob
        return dict_all

    # Find the freq of option from the matrix
    def find_option_freq(self, option, haplos, missing):
        # one = True
        # Here we calculate the  freq of first the first division of the option
        # The first division
        division = option[0]

        # The string of the haplos
        string_option = self.create_haplos_string(haplos, division, missing)

        # Find the dict of the haplos and their freqs
        dict_all = self.get_haplo_freqs_pan_b(string_option, division)

        # If we found the freq of the first optin- if no-no need to continue
        if dict_all != {}:
            # Go over all the parts of the haplo according to option
            # Skip the first options
            for i in range(1, len(option)):
                # Create the string of the option
                # Example: option [1,2,3] will get :A*01:01~B*07:02~C*07:02
                # Example: option [5] will get DRB1*13:02
                division = option[i]
                string_option = self.create_haplos_string(haplos, division, missing)
                if string_option == []:
                    t = 0
                div_dict = self.get_haplo_freqs_pan_b(string_option, division)

                # If this part of the option is empty-the option is not good
                if div_dict == {}:
                    result = all(elem in missing for elem in division)
                    if (
                        result
                    ):  # and one:#len(division) == 1 and len(missing) == 1 and division[0] == missing[0]:
                        type = self.cypher_plan_b.findTypeFromIndexes(division)
                        div_dict = self.netGraph.haps_with_probs_by_label(type)
                        # one = False
                    else:
                        # There is no freqs
                        dict_all = {}

                        # No need to continue
                        break
                # Add the new data to old data and make it a new "dict_all"
                dict_all = self.open_option_(div_dict, dict_all)
        return dict_all

    def comp_hap_prob_plan_b(self, Hap, division, missing):
        if division[0] == list(set(self.index_dict.values())):  # [1, 2, 3, 4, 5]:
            return self.comp_hap_prob([Hap[0]], 0, 0, 0)
        haplo_probs = self.find_option_freq(division, Hap, missing)

        # Return {'Haps': haplos, 'Probs': probs} of the dict
        return self.open_dict_data(haplo_probs)

    def missing_from_data_to_string(self, hap, not_in_data):
        str_hap = ""
        str_not_in = []
        str_list = []

        for luci in hap:
            type_luci = luci.split("*")[0]
            if self.index_dict[type_luci] in not_in_data:
                str_not_in.append(luci)
            else:
                str_hap += "~" + str(luci)

        str_list.append(str_hap[1:])
        str_not_in = list(set(str_not_in))

        return [str_list, str_not_in]

    def find_option_freq_missing_data(self, option, haplos, missing, not_in_data):
        # Here we calculate the  freq of first the first division of the option
        # The first division

        all_the_data = set(self.index_dict.values())  # [1,2,3,4,5]

        # All data that we dont want to
        # all_missing = list(set(missing + not_in_data))
        all_missing = list(set(not_in_data))
        all_the_data = [x for x in all_the_data if x not in all_missing]

        dict_res = dict()
        for hap in haplos[0]:
            # The string of hap
            string_option = self.missing_from_data_to_string(hap, not_in_data)
            # Find the dict of hap and his freqs
            if len(string_option[0]) > 0 and string_option[0][0] != "":
                dict_all = self.get_haplo_freqs_pan_b(string_option[0], all_the_data)
                for key in dict_all.keys():
                    list_key = key.split("~")
                    list_key = (
                        list_key[: not_in_data[0] - 1]
                        + string_option[1]
                        + list_key[not_in_data[0] - 1 :]
                    )
                    dict_res["~".join(sorted(list_key))] = [
                        x * (self._factor_missing_data ** len(all_missing))
                        for x in dict_all[key]
                    ]

        return dict_res

    def comp_hap_prob_plan_b_missing_data(self, Hap, division, missing, not_in_data):

        haplo_probs = self.find_option_freq_missing_data(
            division, Hap, missing, not_in_data
        )
        # Return {'Haps': haplos, 'Probs': probs} of the dict
        return self.open_dict_data(haplo_probs)

    def read_matrix(self, index):
        # Init
        div_option = []

        # If the index is in the size of the matrix
        if len(self.matrix_planb) > index:
            div_option = self.matrix_planb[index]

        # Return the division
        return div_option

    def check_full_haplo(self, hap):
        # To get the array
        wanted_hap = hap[0][0]
        missing = []
        if len(wanted_hap[0][0]) < len(self.full_loci):
            # wanted_hap[1] is a list of locus
            missing = self.find_missing_indexes(wanted_hap[0][0])
        return missing

    def findLinkType(self, typing):
        loci = sorted([self.getLocus(a) for a in typing[0].split("~")])
        loci_lc = "".join([allele.lower() for allele in loci])
        return loci_lc

    def get_haplo_freqs_pan_b(self, haplos_string, division):
        haplo_probs = {}
        if len(haplos_string) > 0:
            # Check the type of the haplos
            type = self.cypher_plan_b.findTypeFromIndexes(division)
            label_haplo = self.cypher_plan_b.findLinkType(haplos_string)
            haplo_probs = self.netGraph.adjs_query_by_color(
                haplos_string, label_haplo, type
            )
        return haplo_probs

    def check_if_alleles_exist(self, allels):
        type = self.find_type_loci(allels[0])
        div = [self.index_dict[type]]
        probs = self.get_haplo_freqs_pan_b(allels, div)
        return probs

    def check_if_alleles_in_data(self, phase, index):
        size_haplo = phase[0][0]
        loci = []
        missing = []
        for t in range(len(size_haplo[0][0])):
            for i in range(len(phase)):
                for hapes in phase[i][index]:
                    for hap in hapes:
                        loci.append(hap[t])
            # Remove duplicates
            loci = list(set(loci))
            probs = self.check_if_alleles_exist(loci)
            if probs == {}:
                type = self.find_type_loci(loci[0])
                missing.append(self.index_dict[type])

            loci = []
        return missing

    def check_if_alleles_of_one_phase_in_data(self, phase):
        size_haplo = phase[0]
        loci = []
        missing = []
        for t in range(len(size_haplo[0])):
            for hapes in phase[0]:
                loci.append(hapes[t])
            # Remove duplicates
            loci = list(set(loci))
            probs = self.check_if_alleles_exist(loci)
            if probs == {}:
                type = self.find_type_loci(loci[0])
                missing.append(self.index_dict[type])

            loci = []
        return missing

    def allel_to_SR(self, allel_dict):
        for allel, probs in allel_dict.items():
            allel_dict[allel] = [sum(probs)]

    def comp_hap_prob_plan_c(self, phases, missing):
        dict_all = {}

        for phase in phases:
            dict_all_tmp = {}
            miss = []
            for i, loci in enumerate(phase):
                index = self.find_type_loci(loci)
                index_loci = self.index_dict[index]
                div_dict = self.get_haplo_freqs_pan_b([loci], [index_loci])
                self.allel_to_SR(div_dict)
                if div_dict == {}:
                    miss.append(loci)
                else:

                    if dict_all_tmp == {}:
                        dict_all_tmp = div_dict
                    else:
                        dict_all_tmp = self.open_option_(div_dict, dict_all_tmp, True)
                        if not dict_all_tmp:
                            break
            if len(miss) > 0:
                for key in dict_all_tmp:
                    list_key = key.split("~")
                    list_key = list_key + miss
                    dict_all["~".join(sorted(list_key))] = [
                        x * (self._factor_missing_data ** len(miss))
                        for x in dict_all_tmp[key]
                    ]
            else:
                for key in dict_all_tmp:
                    dict_all[key] = dict_all_tmp[key]

        type = self.cypher_plan_b.findTypeFromIndexes(missing)
        div_dict = self.netGraph.haps_with_probs_by_label(type)
        self.allel_to_SR(div_dict)
        if dict_all:
            if div_dict:
                dict_all = self.open_option_(div_dict, dict_all, True)
            else:
                for miss in missing:
                    type = self.cypher_plan_b.findTypeFromIndexes([miss])
                    div_dict = self.netGraph.haps_with_probs_by_label(type)
                    self.allel_to_SR(div_dict)
                    if div_dict:
                        dict_all = self.open_option_(div_dict, dict_all, True)

        return dict_all

    def comp_phase_prob_plan_c(self, phases, N_Loc, epsilon, MUUG_output):
        epsilon = 0
        hap_pairs_total = []
        pop_res_haplo = []
        geno_seen = set([])
        hap_total = {}
        p_total = []
        pop_res = {}
        Prob2 = []
        Haps2 = []
        maxProb = 0
        missing = self.check_full_haplo(phases)
        for i in range(len(phases)):
            # If we dont know about missing alleles

            P1 = self.open_dict_data(
                self.comp_hap_prob_plan_c(phases[i][0][0], missing)
            )

            # This will open locus ambiguities and comp probabilities for Hap1
            Haps1 = P1["Haps"]
            Prob1 = P1["Probs"]

            if len(Prob1) > 0:
                P2 = self.open_dict_data(
                    self.comp_hap_prob_plan_c(phases[i][1][0], missing)
                )

                Haps2 = P2["Haps"]
                Prob2 = P2["Probs"]

            Prob1WithIndexes = self.convert_list_to_one_dim(Prob1)

            Prob2WithIndexes = self.convert_list_to_one_dim(Prob2)

            if MUUG_output:
                maxProb = self.calc_haps_pairs(
                    Haps1,
                    Haps2,
                    Prob1WithIndexes,
                    Prob2WithIndexes,
                    epsilon,
                    hap_total,
                    pop_res,
                    maxProb,
                    geno_seen,
                )
            else:
                maxProb = self.calc_haps_pairs_haplotype(
                    Haps1,
                    Haps2,
                    Prob1WithIndexes,
                    Prob2WithIndexes,
                    epsilon,
                    hap_total,
                    pop_res,
                    maxProb,
                    hap_pairs_total,
                    geno_seen,
                    pop_res_haplo,
                    p_total,
                )
        if MUUG_output:
            sum_prob = sum(pop_res.values())
            pop_res_final = {"all_pops,all_pops": sum_prob}
            return {"MaxProb": maxProb, "Haps": hap_total, "Pops": pop_res_final}

        for i in range(len(pop_res_haplo)):
            pop_res_haplo[i][0] = "all_pops"
            pop_res_haplo[i][1] = "all_pops"

        return {
            "MaxProb": maxProb,
            "Haps": hap_pairs_total,
            "Probs": p_total,
            "Pops": pop_res_haplo,
        }

    # Calculate full plan b
    def comp_phase_prob_plan_b(self, phases, N_Loc, epsilon, MUUG_output):
        # receives a list of phases and computes haps and
        # probabilties and accumulate cartesian product
        hap_pairs_total = []
        pop_res_haplo = []
        geno_seen = set([])
        hap_total = {}
        p_total = []
        pop_res = {}
        Prob2 = []
        Haps2 = []
        matrix_index = 0
        maxProb = 0
        # Want to check if all in data base
        missing_data_1 = self.check_if_alleles_in_data(phases, 0)
        missing_data_2 = self.check_if_alleles_in_data(phases, 1)

        for i in range(len(phases)):
            phases[i][0].append(10)
            phases[i][1].append(10)

        # While we didnt get any haps
        while hap_total == {}:
            # Get the hapl division
            option = self.read_matrix(matrix_index)
            if option == []:
                break

            missing = self.check_full_haplo(phases)
            for i in range(len(phases)):
                # If we dont know about missing alleles
                if missing_data_1 == []:
                    index = min(matrix_index, phases[i][0][1])
                    option = self.read_matrix(index)
                    P1 = self.comp_hap_prob_plan_b(phases[i][0], option, missing)
                    if len(P1["Haps"]):
                        phases[i][0][1] = index
                else:
                    P1 = self.comp_hap_prob_plan_b_missing_data(
                        phases[i][0], option, missing, missing_data_1
                    )
                # This will open locus ambiguities and comp probabilities for Hap1
                Haps1 = P1["Haps"]
                Prob1 = P1["Probs"]
                # if len(Prob1) > 0:
                # If we dont know about missing alleles
                if missing_data_2 == []:
                    index = min(matrix_index, phases[i][1][1])
                    option = self.read_matrix(index)
                    P2 = self.comp_hap_prob_plan_b(phases[i][1], option, missing)
                    if len(P2["Haps"]):
                        phases[i][1][1] = index
                    Haps2 = P2["Haps"]
                    Prob2 = P2["Probs"]
                else:
                    if len(Prob1) > 0:
                        P2 = self.comp_hap_prob_plan_b_missing_data(
                            phases[i][1], option, missing, missing_data_2
                        )
                        Haps2 = P2["Haps"]
                        Prob2 = P2["Probs"]
                    # This will do the same for Hap 2;

                Prob1WithIndexes = self.convert_list_to_one_dim(Prob1)

                Prob2WithIndexes = self.convert_list_to_one_dim(Prob2)

                if MUUG_output:
                    maxProb = self.calc_haps_pairs(
                        Haps1,
                        Haps2,
                        Prob1WithIndexes,
                        Prob2WithIndexes,
                        epsilon,
                        hap_total,
                        pop_res,
                        maxProb,
                        geno_seen,
                    )
                else:
                    maxProb = self.calc_haps_pairs_haplotype(
                        Haps1,
                        Haps2,
                        Prob1WithIndexes,
                        Prob2WithIndexes,
                        epsilon,
                        hap_total,
                        pop_res,
                        maxProb,
                        hap_pairs_total,
                        geno_seen,
                        pop_res_haplo,
                        p_total,
                    )

            # Go to the next row in thr matrix
            matrix_index += 1

        matrix_index = 10
        matrix_index_curr = 0
        while hap_total == {} and matrix_index_curr < 6:
            for i in range(len(phases)):
                index_1 = min(matrix_index, phases[i][0][1])
                index_2 = min(matrix_index, phases[i][1][1])
                if not (index_1 == 10 and index_2 == 10):
                    if index_1 == 10 and len(phases[i][0][0]) > 0:
                        option = self.read_matrix(matrix_index_curr)
                        missing_data_1 = self.check_if_alleles_of_one_phase_in_data(
                            phases[i][0]
                        )
                        P1 = self.comp_hap_prob_plan_b_missing_data(
                            phases[i][0], option, missing, missing_data_1
                        )

                        option = self.read_matrix(index_2)
                        P2 = self.comp_hap_prob_plan_b(phases[i][1], option, missing)

                    if index_2 == 10 and len(phases[i][1][0]) > 0:
                        option = self.read_matrix(index_1)
                        P1 = self.comp_hap_prob_plan_b(phases[i][0], option, missing)

                        option = self.read_matrix(matrix_index_curr)
                        missing_data_2 = self.check_if_alleles_of_one_phase_in_data(
                            phases[i][1]
                        )
                        P2 = self.comp_hap_prob_plan_b_missing_data(
                            phases[i][1], option, missing, missing_data_2
                        )

                    Haps1 = P1["Haps"]
                    Prob1 = P1["Probs"]
                    Haps2 = P2["Haps"]
                    Prob2 = P2["Probs"]

                    Prob1WithIndexes = self.convert_list_to_one_dim(Prob1)

                    Prob2WithIndexes = self.convert_list_to_one_dim(Prob2)

                    if MUUG_output:
                        maxProb = self.calc_haps_pairs(
                            Haps1,
                            Haps2,
                            Prob1WithIndexes,
                            Prob2WithIndexes,
                            epsilon,
                            hap_total,
                            pop_res,
                            maxProb,
                            geno_seen,
                        )
                    else:
                        maxProb = self.calc_haps_pairs_haplotype(
                            Haps1,
                            Haps2,
                            Prob1WithIndexes,
                            Prob2WithIndexes,
                            epsilon,
                            hap_total,
                            pop_res,
                            maxProb,
                            hap_pairs_total,
                            geno_seen,
                            pop_res_haplo,
                            p_total,
                        )

            matrix_index_curr += 1

            # p_total returns an array of N*N (N is number of populations), hap_total - pairs of haplotypes.
            # pop_res are the names of the populations
        if MUUG_output:
            return {"MaxProb": maxProb, "Haps": hap_total, "Pops": pop_res}

        return {
            "MaxProb": maxProb,
            "Haps": hap_pairs_total,
            "Probs": p_total,
            "Pops": pop_res_haplo,
        }

    # haplotype - list of all loci options
    # return - type of tne haplotype
    def input_type(self, haplotype):
        type_list = []
        for locus_options in haplotype:
            locus = locus_options.split("*")[0]
            type_list.append(self.index_dict[locus])
        return type_list

    def open_ambiguities(self, hap, loc, split_loc):
        return open_ambiguities(hap, loc, split_loc)

    def comp_cand(
        self, gl_string, binary, epsilon, n, MUUG_output, haps_output, planb, em
    ):
        # receives a list of phases and computes haps and
        # probabilties and accumulate cartesian productEpsilon=0.0001
        chr = self.gl2haps(gl_string)
        if chr == []:
            return None, None
        # if we in 9-loci, check if the type input in valid format
        if self.nodes_for_plan_A:
            geno_type = self.input_type(chr["Genotype"][0])
            if not geno_type in self.nodes_for_plan_A:
                return None, None

        n_loci = chr["N_Loc"]

        # generate the 2^(N-1) phases
        pmags = self.gen_phases(chr["Genotype"], n_loci, binary)

        # return if the result is empty (why would that be?)
        if pmags == []:
            return None, None

        # res_muugs = {'Haps': 'NaN', 'Probs': 0}
        res_muugs = {"MaxProb": 0, "Haps": {}, "Pops": {}}
        res_haps = {"Haps": "Nan", "Probs": 0, "Pops": {}}
        # this step "opens" each phase by generating the cartesian product
        # NOTE: why do it here and not pass the ambiguity to the graph query?
        # [['A*01:01/A*01:02/A*01:03/A*01:16N', 'C*06:02']
        # to
        # [['A*01:01', 'C*06:02'],
        #  ['A*01:02', 'C*06:02'],
        #  ['A*01:03', 'C*06:02'],
        #  ['A*01:16', 'C*06:02']]
        #
        phases = self.open_phases(pmags, n_loci, gl_string)
        if not phases:
            print("in reduce_phase_to_valid_alleles")
            self.reduce_phase_to_valid_allels(pmags, n_loci)
            phases = self.open_phases(pmags, n_loci, gl_string)
        if not phases:
            print("in reduce_phase_to_commons_alleles")
            self.reduce_phase_to_commons_alleles(pmags, n_loci, commons_number=10)
            phases = self.open_phases(pmags, n_loci, gl_string)

        if phases:
            if MUUG_output:
                prior_matrix_orig = np.array(
                    self.priorMatrix, order="K", copy=True
                )  # copy.deepcopy(self.priorMatrix)
                res_muugs = self.call_comp_phase_prob(
                    epsilon, n, phases, chr, True, planb
                )
                if planb and len(res_muugs["Haps"]) == 0:
                    self.plan = "c"
                    self.reduce_phase_to_commons_alleles(pmags, n_loci, 1, True)
                    phases = self.open_phases(pmags, n_loci, gl_string)
                    res_muugs = self.comp_phase_prob_plan_c(
                        phases, n_loci, epsilon, True
                    )
                self.priorMatrix = prior_matrix_orig
            if haps_output:
                res_haps = self.call_comp_phase_prob(
                    epsilon, n, phases, chr, False, planb
                )
                if planb and len(res_haps["Haps"]) == 0 and not em:  # em
                    self.reduce_phase_to_commons_alleles(pmags, n_loci, 1, True)
                    phases = self.open_phases(pmags, n_loci, gl_string)
                    res_haps = self.comp_phase_prob_plan_c(
                        phases, n_loci, epsilon, False
                    )

        return res_muugs, res_haps

    def call_comp_phase_prob(self, epsilon, n, phases, chr, MUUG_output, planb):

        n_res = 0  # number of results
        min_res = 100  # minimum number of results (NOTE: why 10?)
        min_epsilon = 1.0e-9  # minimum threshold, then why the other?
        res = {"Haps": "NaN", "Probs": 0}
        lastRound = False
        while epsilon > 0:
            # each time through reduce epsilon by an order of magnitude
            # until you have seen "min_res" results
            epsilon /= 10
            #
            # if you have reduced epsilon below min_epsilon then just make it 0
            #
            if epsilon < min_epsilon:
                epsilon = 0.0

            # get results

            if MUUG_output:
                res = self.comp_phase_prob_genotype(phases, chr["N_Loc"], epsilon, n)
            else:
                res = self.comp_phase_prob_haplotype(phases, chr["N_Loc"], epsilon, n)
            haps = res["Haps"]
            n_res = len(haps)
            if n_res > 0 and epsilon > 0:
                # sorted_by_value = sorted(haps.items(), key=lambda kv: kv[1], reverse=True)
                epsilon = res["MaxProb"] / 100000
                lastRound = True
                break

        if lastRound:
            if MUUG_output:
                res = self.comp_phase_prob_genotype(phases, chr["N_Loc"], epsilon, n)
            else:
                res = self.comp_phase_prob_haplotype(phases, chr["N_Loc"], epsilon, n)

        # no plan b
        for level in range(2):
            if level == 1:
                self.priorMatrix = np.ones(
                    (len(self.populations), len(self.populations))
                )  ####
            if planb and len(res["Haps"]) == 0:
                self.plan = "b"
                epsilon = 1e-14
                n_res = 0
                min_res = 10
                min_epsilon = 1.0e-3
                # self.priorMatrix = np.ones((len(self.populations), len(self.populations)))
                while (epsilon > 0) & (n_res < min_res):
                    epsilon /= 10
                    if epsilon < min_epsilon:
                        epsilon = 0.0
                    phases_planb = deepcopy_list(phases)
                    # Find the option according to plan b
                    if MUUG_output:
                        res = self.comp_phase_prob_plan_b(
                            phases_planb, chr["N_Loc"], epsilon, True
                        )
                    else:
                        res = self.comp_phase_prob_plan_b(
                            phases_planb, chr["N_Loc"], epsilon, False
                        )
                    n_res = len(res["Haps"])

        return res

    # sum all probs for each first race and return the one with the highest prob
    def find_first_race(self, pops, probs):
        prob_by_race = dict()
        for i, pop in enumerate(pops):
            key = pop[0]
            if key in prob_by_race:
                prob_by_race[key] += probs[i][0] * probs[i][1]
            else:
                prob_by_race[key] = probs[i][0] * probs[i][1]

        sortedDict = sorted(prob_by_race.items(), key=lambda x: x[1], reverse=True)
        return sortedDict[0][0]

    # aum all probs of pairs- race1,second race
    # and return second race of the highest pair
    def find_second_race(self, race1, pops, probs):
        prob_by_race = dict()
        for i, pop in enumerate(pops):
            if pop[0] == race1:
                key = pop[1]
                if key in prob_by_race:
                    prob_by_race[key] += probs[i][0] * probs[i][1]
                else:
                    prob_by_race[key] = probs[i][0] * probs[i][1]

        sortedDict = sorted(prob_by_race.items(), key=lambda x: x[1], reverse=True)
        return sortedDict[0][0]

    # find the haps of race1-race2, and sum the probs
    def write_haplotype(
        self, race1, race2, haps, pops, probs, name_gl, foutHap, foutPop
    ):
        haps_dict = dict()
        for i, line in enumerate(pops):
            if race1 == line[0] and race2 == line[1]:
                key = haps[i][0] + "," + haps[i][1]
                if key in haps_dict:
                    haps_dict[key] += probs[i][0] * probs[i][1]
                else:
                    haps_dict[key] = probs[i][0] * probs[i][1]

        # write the 1000 highest to file
        sortedDict = sorted(haps_dict.items(), key=lambda x: x[1], reverse=True)
        minBestResult = min(1000, len(sortedDict))
        for k in range(minBestResult):
            foutHap.write(
                name_gl
                + ","
                + str(sortedDict[k][0])
                + ","
                + str(sortedDict[k][1])
                + ","
                + str(k)
                + ","
                + "\n"
            )
        foutPop.write(name_gl + "," + str(race1) + "," + str(race2) + "\n")

    # find pair of races with the highest prob
    def find_races(self, res, foutHap, foutPop, name_gl):
        if len(res["Pops"]) > 0:
            race1 = self.find_first_race(res["Pops"], res["Probs"])
            race2 = self.find_second_race(race1, res["Pops"], res["Probs"])
            self.write_haplotype(
                race1,
                race2,
                res["Haps"],
                res["Pops"],
                res["Probs"],
                name_gl,
                foutHap,
                foutPop,
            )

    """def calc_priority_matrix(self, list_race1, list_race2, priority ):
        list_race1 = list_race1.split(';')
        list_race2 = list_race2.split(';')
        popLen = len(self.populations)
        self.priorMatrix = np.zeros((popLen, popLen))
        # Identity matrix
        Imatrix = np.identity(popLen)

        for race_1 in list_race1:
            for race_2 in list_race2:
                    tmp_matrix = np.zeros((popLen, popLen))
                #if race1 in self.populations and race2 in self.populations:
                    race1 = self.populations.index(race_1)
                    race2 = self.populations.index(race_2)
                # calc priority
                    for i in range(popLen):
                        tmp_matrix[race1, i] = tmp_matrix[race1, i] + priority['gamma']
                        tmp_matrix[i, race2] = tmp_matrix[i, race2] + priority['gamma']

                    tmp_matrix[race1, race2] -= priority['gamma']

                    tmp_matrix[race1, race2] = tmp_matrix[race1, race2] + priority['alpha']

                    if race1 != race2:
                        tmp_matrix = tmp_matrix + tmp_matrix.transpose()
                        tmp_matrix[race1, race1] -= priority['gamma']
                        tmp_matrix[race2, race2] -= priority['gamma']

                    tmp_matrix[race1, race1] += priority['delta']
                    if race1 != race2:
                        tmp_matrix[race2, race2] += priority['delta']
                    tmp_matrix = priority['eta'] * np.ones((popLen, popLen)) + tmp_matrix + priority['beta'] * Imatrix

                    self.priorMatrix += tmp_matrix

        # multiply by pop ratio
        prior_sum = 0
        for i in range(popLen):
            for j in range(popLen):
                self.priorMatrix[i][j] = self.priorMatrix[i][j] * self.count_by_prob[i] * self.count_by_prob[j]
                prior_sum += self.priorMatrix[i][j]

        self.priorMatrix = self.priorMatrix / prior_sum"""

    def calc_priority_matrix(self, list_race1, list_race2, priority):
        popLen = len(self.populations)
        self.priorMatrix = np.zeros((popLen, popLen))
        # Identity matrix
        Imatrix = np.identity(popLen)

        for race_1 in list_race1:
            for race_2 in list_race2:
                if race_1 == "" and race_2 == "":
                    continue
                if race_1 == "" or race_2 == "":
                    race = (
                        self.populations.index(race_2)
                        if race_1 == ""
                        else self.populations.index(race_1)
                    )
                    """if race_1 == '':
                        race = self.populations.index(race_2)
                    else:
                        race = self.populations.index(race_1)"""

                    tmp_matrix = np.zeros((popLen, popLen))

                    # calc priority
                    for i in range(popLen):
                        tmp_matrix[race, i] = (
                            tmp_matrix[race, i] + priority["gamma"] * 2
                        )
                    tmp_matrix = tmp_matrix + tmp_matrix.transpose()
                    tmp_matrix[race, race] -= priority["gamma"] * 2

                    tmp_matrix = (
                        priority["eta"] * np.ones((popLen, popLen))
                        + tmp_matrix
                        + priority["beta"] * Imatrix
                    )

                    self.priorMatrix += tmp_matrix
                else:
                    tmp_matrix = np.zeros((popLen, popLen))
                    race1 = self.populations.index(race_1)
                    race2 = self.populations.index(race_2)
                    # calc priority
                    for i in range(popLen):
                        tmp_matrix[race1, i] = tmp_matrix[race1, i] + priority["gamma"]
                        tmp_matrix[i, race2] = tmp_matrix[i, race2] + priority["gamma"]

                    tmp_matrix[race1, race2] -= priority["gamma"]

                    tmp_matrix[race1, race2] = (
                        tmp_matrix[race1, race2] + priority["alpha"]
                    )

                    if race1 != race2:
                        tmp_matrix = tmp_matrix + tmp_matrix.transpose()
                        tmp_matrix[race1, race1] -= priority["gamma"]
                        tmp_matrix[race2, race2] -= priority["gamma"]

                    tmp_matrix[race1, race1] += priority["delta"]
                    if race1 != race2:
                        tmp_matrix[race2, race2] += priority["delta"]
                    tmp_matrix = (
                        priority["eta"] * np.ones((popLen, popLen))
                        + tmp_matrix
                        + priority["beta"] * Imatrix
                    )

                    self.priorMatrix += tmp_matrix

        # multiply by pop ratio
        prior_sum = 0
        for i in range(popLen):
            for j in range(popLen):
                self.priorMatrix[i][j] = (
                    self.priorMatrix[i][j]
                    * self.count_by_prob[i]
                    * self.count_by_prob[j]
                )
                prior_sum += self.priorMatrix[i][j]

        self.priorMatrix = self.priorMatrix / prior_sum

    def update_prob_by_priority(self, res, race1, race2, priority):
        # prob by priority
        pops = res["Pops"]
        for i, pop in enumerate(pops):
            race1 = self.populations.index(pop[0])
            race2 = self.populations.index(pop[1])
            res["Probs"][i][0] = float(res["Probs"][i][0]) * math.sqrt(
                self.priorMatrix[race1][race2]
            )
            res["Probs"][i][1] = float(res["Probs"][i][1]) * math.sqrt(
                self.priorMatrix[race1][race2]
            )
        return res

    def impute_one(
        self,
        subject_id,
        gl,
        binary,
        race1,
        race2,
        priority,
        epsilon,
        n,
        MUUG_output,
        haps_output,
        planb,
        em,
    ):  # em
        clean_gl = clean_up_gl(gl)
        if self.unk_priors == "MR":
            self.priorMatrix = np.ones((len(self.populations), len(self.populations)))
        else:
            self.priorMatrix = np.identity(len(self.populations))
        to_calc_prior_matrix = False
        if race1 or race2:
            race1 = race1.split(";")
            for i, race in enumerate(race1):
                if not race in self.populations:
                    race1[i] = ""
                else:
                    to_calc_prior_matrix = True
            race2 = race2.split(";")
            for i, race in enumerate(race2):
                if not race in self.populations:
                    race2[i] = ""
                else:
                    to_calc_prior_matrix = True
            if to_calc_prior_matrix:
                self.calc_priority_matrix(race1, race2, priority)
        # lookup based on genotype
        res_muugs = res_haps = None
        if gl:
            res_muugs, res_haps = self.comp_cand(
                clean_gl, binary, epsilon, n, MUUG_output, haps_output, planb, em
            )

        return subject_id, res_muugs, res_haps

    def impute_file(self, config, planb=None, em_mr=False, em=False):  ##em
        priority = config["priority"]
        MUUG_output = config["output_MUUG"]
        haps_output = config["output_haplotypes"]
        n = 1000
        epsilon = config["epsilon"]
        number_of_results = config["number_of_results"]
        number_of_pop_results = config["number_of_pop_results"]
        # planb = config["planb"]#em
        if planb is None:  # em
            planb = config["planb"]  # em

        # TODO: do the right thing if its a gzip
        if self.verbose:
            self.logger.info("Starting Imputation!")

        f_bin_exist = False
        if os.path.isfile(config["bin_imputation_input_file"]):
            with open(config["bin_imputation_input_file"]) as json_file:
                f_bin = json.load(json_file)
                f_bin_exist = True

        f = open(config["imputation_input_file"], "r")

        if MUUG_output:
            fout_hap_muug = open(config["imputation_out_umug_freq_file"], "w")
            fout_pop_muug = open(config["imputation_out_umug_pops_file"], "w")
        if haps_output:
            fout_hap_haplo = open(config["imputation_out_hap_freq_file"], "w")
            fout_pop_haplo = open(config["imputation_out_hap_pops_file"], "w")

        miss = open(config["imputation_out_miss_file"], "w")
        problem = open(config["imputation_out_problem_file"], "w")

        with f as lines:
            for i, name_gl in enumerate(lines):
                try:
                    name_gl = name_gl.rstrip()  # remove trailing whitespace
                    if "," in name_gl:
                        list_gl = name_gl.split(",")
                    else:
                        list_gl = name_gl.split("%")

                    subject_id = list_gl[0]
                    subject_gl = list_gl[1]
                    subject_bin = [1] * (len(self.full_loci) - 1)
                    if f_bin_exist:
                        subject_bin = f_bin[subject_id]
                    race1 = race2 = None
                    if len(list_gl) > 2:
                        race1 = list_gl[2]
                        race2 = list_gl[3]

                    start = timeit.default_timer()

                    #
                    # Impute One
                    # This method will impute one subject given the parameters
                    self.plan = "a"
                    self.option_1 = 0
                    self.option_2 = 0
                    subject_id, res_muugs, res_haps = self.impute_one(
                        subject_id,
                        subject_gl,
                        subject_bin,
                        race1,
                        race2,
                        priority,
                        epsilon,
                        n,
                        MUUG_output,
                        haps_output,
                        planb,
                        em,
                    )  # em

                    if res_muugs is None:
                        problem.write(str(i) + "," + str(subject_id) + "\n")
                        continue

                    if (
                        len(res_haps["Haps"]) == 0 or res_haps["Haps"] == "NaN"
                    ) and len(res_muugs["Haps"]) == 0:
                        miss.write(str(i) + "," + str(subject_id) + "\n")

                    if haps_output:
                        haps = res_haps["Haps"]
                        probs = res_haps["Probs"]
                        pops = res_haps["Pops"]
                        print(
                            "{index} Subject: {id} {hap_length} haplotypes".format(
                                index=i, id=subject_id, hap_length=len(haps)
                            )
                        )
                        if em_mr:
                            write_best_hap_race_pairs(
                                subject_id,
                                haps,
                                pops,
                                probs,
                                number_of_results,
                                fout_hap_haplo,
                            )
                            write_best_prob(subject_id, pops, probs, 1, fout_pop_haplo)
                        else:
                            write_best_prob(
                                subject_id,
                                haps,
                                probs,
                                number_of_results,
                                fout_hap_haplo,
                                "+",
                            )
                            write_best_prob(
                                subject_id,
                                pops,
                                probs,
                                number_of_pop_results,
                                fout_pop_haplo,
                            )
                    if MUUG_output:
                        haps = res_muugs["Haps"]
                        pops = res_muugs["Pops"]
                        print(
                            "{index} Subject: {id} {hap_length} haplotypes".format(
                                index=i, id=subject_id, hap_length=len(haps)
                            )
                        )
                        write_best_prob_genotype(
                            subject_id, haps, number_of_results, fout_hap_muug
                        )
                        write_best_prob_genotype(
                            subject_id, pops, number_of_pop_results, fout_pop_muug
                        )

                    if self.verbose:
                        self.logger.info(
                            "{index} Subject: {id} {hap_length} haplotypes".format(
                                index=i, id=subject_id, hap_length=len(haps)
                            )
                        )
                        self.logger.info(
                            "{index} Subject: {id} plan: {plan} oppen_phases - count of open regular option: {option1}, count of alternative opening: {option2}".format(
                                index=i,
                                id=subject_id,
                                plan=self.plan,
                                option1=self.option_1,
                                option2=self.option_2,
                            )
                        )

                    stop = timeit.default_timer()
                    time_taken = stop - start
                    print(time_taken)
                    if self.verbose:
                        self.logger.info("Time taken: " + str(time_taken))
                except:
                    print(f"{i} Subject: {subject_id} - Exception")
                    problem.write(str(name_gl) + "\n")
                    continue

            f.close()
            if MUUG_output:
                fout_hap_muug.close()
                fout_pop_muug.close()
            if haps_output:
                fout_hap_haplo.close()
                fout_pop_haplo.close()

            miss.close()
            problem.close()
