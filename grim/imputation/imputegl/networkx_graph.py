import networkx as nx
import csv


def missing(labelA, labelB):
    a = list(labelA)
    b = list(labelB)
    return [x for x in b if x not in a]


class Graph(object):
    __slots__ = 'graph', 'labelDict', 'whole_graph', 'full_loci', 'nodes_plan_a', 'nodes_plan_b'

    def __init__(self, config):
        self.graph = nx.DiGraph()
        self.labelDict = {}
        self.whole_graph = nx.DiGraph()
        self.full_loci = config["full_loci"]
        self.nodes_plan_a, self.nodes_plan_b = [], []
        if config["nodes_for_plan_A"]:
            path = '/'.join(config["node_file"].split('/')[:-1])

            # bug: dies if file doesn't exist
            # bug: list_f doesn't exist
            with open(path + "/nodes_for_plan_a.txt") as list_f:
                for item in list_f:
                    self.nodes_plan_a.append(item.strip())
            # bug: dies if file doesn't exist
            with open(path + "/nodes_for_plan_b.txt") as list_f:
                for item in list_f:
                    self.nodes_plan_b.append(item.strip())
            # self.nodes_plan_a = pickle.load(open( path + '/nodes_for_plan_a.pkl', "rb"))
            # self.nodes_plan_b = pickle.load(open( path + '/nodes_for_plan_b.pkl', "rb"))

    # build graph from files of nodes and edges between nodes with top relation
    def build_graph(self, nodesFile, edgesFile, allEdgesFile):
        nodesDict = dict()
        # add nodes from file
        with open(nodesFile) as nodesfile:
            readNodes = csv.reader(nodesfile, delimiter=",")
            firstLine = next(readNodes)
            for row in readNodes:
                if len(row) > 0:
                    if not self.nodes_plan_a or row[2] in self.nodes_plan_a:
                        self.graph.add_node(
                            row[1],
                            label=row[2],
                            freq=list(map(float, row[3].split(";"))),
                        )
                    if not self.nodes_plan_b or row[2] in self.nodes_plan_b:
                        self.whole_graph.add_node(
                            row[1],
                            label=row[2],
                            freq=list(map(float, row[3].split(";"))),
                        )
                    nodesDict[row[0]] = row[1]

        nodesfile.close()

        # add edges from file
        with open(edgesFile) as edgesfile:
            readEdges = csv.reader(edgesfile, delimiter=",")
            firstLine = next(readEdges)
            for row in readEdges:
                if len(row) > 0:
                    node1 = nodesDict[row[0]]
                    node2 = nodesDict[row[1]]
                    if node1 in self.graph and node2 in self.graph:
                        if self.graph.nodes[node1]["label"] == self.full_loci:
                            self.graph.add_edge(node2, node1)
                        else:
                            self.graph.add_edge(node1, node2)

        edgesfile.close()

        # add edges from file
        with open(allEdgesFile) as allEdgesfile:
            readEdges = csv.reader(allEdgesfile, delimiter=",")
            firstLine = next(readEdges)
            for row in readEdges:
                if len(row) > 0:
                    node1 = nodesDict[row[0]]
                    node2 = nodesDict[row[1]]
                    kind = "-".join(
                        sorted(
                            [
                                self.whole_graph.nodes[node1]["label"],
                                self.whole_graph.nodes[node2]["label"],
                            ],
                            key=len,
                        )
                    )
                    self.whole_graph.add_edge(node1, node2, color=kind)

        allEdgesfile.close()
        nodesDict.clear()

    # return all haplotype by specific label
    def haps_by_label(self, label):
        # cheak if already found
        if label in self.labelDict:
            return self.labelDict[label]
        # not found yet. serach and save in labelDict
        hapsList = []
        if not self.nodes_plan_a or label in self.nodes_plan_a:
            for key, key_data in self.graph.nodes(data=True):
                if key_data["label"] == label:
                    hapsList.append(key)
        elif label in self.nodes_plan_b:
            for key, key_data in self.whole_graph.nodes(data=True):
                if key_data["label"] == label:
                    hapsList.append(key)
        self.labelDict[label] = hapsList
        return hapsList

    def haps_with_probs_by_label(self, label):
        dictAlleles = {}

        listLabel = self.haps_by_label(label)
        if not self.nodes_plan_a or label in self.nodes_plan_a:
            for allele in listLabel:
                dictAlleles[allele] = self.graph.nodes[allele]["freq"]
        elif label in self.nodes_plan_b:
            for allele in listLabel:
                dictAlleles[allele] = self.whole_graph.nodes[allele]["freq"]

        return dictAlleles

    # find all adj of alleleList from label 'ABCQR'
    def adjs_query(self, alleleList):
        adjDict = dict()
        for allele in alleleList:
            if allele in self.graph.nodes():
                if self.graph.nodes[allele]["label"] == self.full_loci:  # 'ABCQR':
                    adjDict[allele] = self.graph.nodes[allele]["freq"]
                else:
                    adjs = self.graph.adj[allele]
                    for adj in adjs:
                        adjDict[adj] = self.graph.nodes[adj]["freq"]
        return adjDict

    # find all adj of alleleList by label
    def adjs_query_by_color(self, alleleList, labelA, labelB):
        # copyLabelA = labelA
        adjDict = dict()
        if labelA == labelB:
            return self.node_probs(alleleList, labelA)

        for allele in alleleList:
            if allele in self.whole_graph.nodes():
                copyLabelA = labelA
                newLabelA = labelA
                miss = missing(labelA, labelB)
                alleles = [allele]

                while len(miss) > 0:
                    tmpAllels = list()
                    newLabelA = copyLabelA + miss[0]
                    newLabelA = "".join(sorted(newLabelA))
                    del miss[0]
                    for oneAllel in alleles:
                        #    alleles.remove(oneAllel)
                        adjs = self.whole_graph.adj[oneAllel]
                        label = copyLabelA + "-" + newLabelA
                        for adj in adjs:
                            if adjs[adj]["color"] == label:
                                tmpAllels.append(adj)
                    alleles = tmpAllels
                    copyLabelA = newLabelA

                for adj in alleles:
                    adjDict[adj] = self.whole_graph.nodes[adj]["freq"]
        return adjDict

    # return dict of nodes and there proper freq
    def node_probs(self, nodes, label):
        nodesDict = dict()
        if not self.nodes_plan_b or label in self.nodes_plan_b:
            for node in nodes:
                if node in self.whole_graph.nodes():
                    nodesDict[node] = self.whole_graph.nodes[node]["freq"]
        elif label in self.nodes_plan_a:
            for node in nodes:
                if node in self.whole_graph.nodes():
                    nodesDict[node] = self.graph.nodes[node]["freq"]
        return nodesDict
