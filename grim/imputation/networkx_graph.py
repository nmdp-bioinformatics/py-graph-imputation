"""
Yuli Tshuva
Create a data structure of directed Graph suitable for our needs with numpy arrays and dictionaries.
"""

import csv
import numpy as np
import pandas as pd
import gc
from tqdm.auto import tqdm
import time


class Graph(object):
    def __init__(self, config):
        """
        Clean initiation.
        """
        self.Edges = []
        self.Vertices = []
        self.Vertices_attributes = {}
        self.Neighbors_start = []

        self.Whole_Edges = []
        self.Whole_Vertices = []
        self.Whole_Vertices_attributes = {}
        self.Whole_Neighbors_start = []

        self.labelDict = {}

        self.full_loci = config["full_loci"]
        self.nodes_plan_a, self.nodes_plan_b = [], []
        if config["nodes_for_plan_A"]:
            path = "/".join(config["node_file"].split("/")[:-1])
            with open(path + "/nodes_for_plan_a.txt") as list_f:
                for item in list_f:
                    self.nodes_plan_a.append(item.strip())
            with open(path + "/nodes_for_plan_b.txt") as list_f:
                for item in list_f:
                    self.nodes_plan_b.append(item.strip())

    def build_graph(self, nodesFile, edgesFile, allEdgesFile):
        """Build the graph by scanning the three files, line by line and filling up the class variables."""
        nodesDict = dict()
        with open(nodesFile) as nodesfile:
            readNodes = csv.reader(nodesfile, delimiter=",")
            next(readNodes)
            for row in tqdm(readNodes, desc="Vertices Creation:"):
                if len(row) > 0:
                    if not self.nodes_plan_a or row[2] in self.nodes_plan_a:
                        self.Vertices.append(row[1])
                        vertex_id = len(self.Vertices) - 1
                        self.Vertices_attributes[row[1]] = (
                            row[2],
                            list(map(float, row[3].split(";"))),
                            vertex_id,
                        )

                    if not self.nodes_plan_b or row[2] in self.nodes_plan_b:
                        self.Whole_Vertices.append(row[1])
                        vertex_id = len(self.Whole_Vertices) - 1
                        self.Whole_Vertices_attributes[row[1]] = (
                            row[2],
                            list(map(float, row[3].split(";"))),
                            vertex_id,
                        )

                    nodesDict[row[0]] = row[1]

        # Add edges from file
        with open(edgesFile) as edgesfile:
            readEdges = csv.reader(edgesfile, delimiter=",")
            next(readEdges)
            for row in tqdm(readEdges, desc="Edges Creation:"):
                if len(row) > 0:
                    node1_id = row[0]
                    node2_id = row[1]
                    node1 = nodesDict[node1_id]
                    node2 = nodesDict[node2_id]
                    if (
                        node1 in self.Vertices_attributes
                        and node2 in self.Vertices_attributes
                    ):
                        node1_label = self.Vertices_attributes[node1][0]
                        if node1_label == self.full_loci:
                            self.Edges.append([node2_id, node1_id])
                        else:
                            self.Edges.append([node1_id, node2_id])

        # add edges from file
        with open(allEdgesFile) as allEdgesfile:
            readEdges = csv.reader(allEdgesfile, delimiter=",")
            next(readEdges)
            for row in tqdm(readEdges, "All Edges Creation:"):
                if len(row) > 0:
                    node1_id = row[0]
                    node2_id = row[1]
                    node1 = nodesDict[node1_id]
                    node2 = nodesDict[node2_id]
                    node1_label = self.Whole_Vertices_attributes[node1][0]
                    node2_label = self.Whole_Vertices_attributes[node2][0]

                    if len(node1_label) < len(node2_label):
                        # Create a connector
                        connector = node2_label + node1

                        if connector not in self.Whole_Vertices_attributes:
                            self.Whole_Vertices.append(connector)
                            connector_id = len(self.Whole_Vertices) - 1
                            self.Whole_Vertices_attributes[connector] = connector_id

                            self.Whole_Edges.append([node1_id, connector_id])
                        else:
                            connector_id = self.Whole_Vertices_attributes[connector]

                        # Append the connector to the whole edges array
                        self.Whole_Edges.append([connector_id, node2_id])

                    else:
                        # Create a connector
                        connector = node1_label + node2

                        if connector not in self.Whole_Vertices_attributes:
                            self.Whole_Vertices.append(connector)
                            connector_id = len(self.Whole_Vertices) - 1
                            self.Whole_Vertices_attributes[connector] = connector_id

                        # Append the connector to the whole edges array
                        self.Whole_Edges.append([node2_id, connector_id])
                        self.Whole_Edges.append([connector_id, node1_id])

        nodesDict.clear()
        del nodesDict

        # Concat all the lists of the edges lists to a numpy array
        self.Edges = np.vstack(self.Edges).astype(np.uint32)
        self.Whole_Edges = np.vstack(self.Whole_Edges).astype(np.uint32)
        self.Vertices = np.array(self.Vertices, dtype=np.object_)
        self.Whole_Vertices = np.array(self.Whole_Vertices, dtype=np.object_)

        # Drop duplications in edges
        df_e = pd.DataFrame(self.Whole_Edges)
        df_e.drop_duplicates(inplace=True)
        del self.Whole_Edges
        self.Whole_Edges = df_e.to_numpy()
        del df_e

        # Sort the Edges arrays - Takes numpy to sort an array of size (10**8, 2) about 45 secs on Google Colab.
        sorted_indices = np.lexsort((self.Edges[:, 1], self.Edges[:, 0]))
        self.Edges = self.Edges[sorted_indices]
        sorted_indices = np.lexsort((self.Whole_Edges[:, 1], self.Whole_Edges[:, 0]))
        self.Whole_Edges = self.Whole_Edges[sorted_indices]

        # Save memory
        del sorted_indices

        # Create a list of the first appearance of a number in the 0 column in the matrix
        unique_values, first_occurrences_indices = np.unique(
            self.Edges[:, 0], return_index=True
        )

        j = 0
        for i in range(0, self.Vertices.shape[0]):
            if int(unique_values[j]) == i:
                self.Neighbors_start.append(int(first_occurrences_indices[j]))
                j += 1
            else:
                try:
                    self.Neighbors_start.append(self.Neighbors_start[-1])
                except:  # In case of the start of the list - empty list
                    self.Neighbors_start.append(0)

        # Free some memory
        del unique_values, first_occurrences_indices

        # Create a list of the first appearance of a number in the 0 column in the matrix
        unique_values, first_occurrences_indices = np.unique(
            self.Whole_Edges[:, 0], return_index=True
        )

        j = 0
        for i in range(0, self.Whole_Vertices.shape[0]):
            if int(unique_values[j]) == i:
                self.Whole_Neighbors_start.append(int(first_occurrences_indices[j]))
                j += 1
            else:
                try:
                    self.Whole_Neighbors_start.append(self.Whole_Neighbors_start[-1])
                except:  # In case of the start of the list - empty list
                    self.Whole_Neighbors_start.append(0)

        # Free some memory
        del unique_values, first_occurrences_indices

        self.Neighbors_start.append(int(len(self.Vertices)))
        self.Whole_Neighbors_start.append(int(len(self.Whole_Vertices)))

        self.Neighbors_start = np.array(self.Neighbors_start, dtype=np.uint32)
        self.Whole_Neighbors_start = np.array(
            self.Whole_Neighbors_start, dtype=np.uint32
        )

        # Take the first column out of the Edges arrays
        ### Do the following to massive save of memory
        Edges = self.Edges[:, 1].copy()
        del self.Edges
        self.Edges = Edges

        Whole_Edges = self.Whole_Edges[:, 1].copy()
        del self.Whole_Edges
        self.Whole_Edges = Whole_Edges

        gc.collect()

    def haps_by_label(self, label):
        """Find haplotypes by their labels.
        Does not use the graphical features of the haplotypes.
        Returns a list of haplotypes.
        Approved."""
        # Check if already found
        if label in self.labelDict:
            return self.labelDict[label]
        # If you get here, label hasn't been found yet. So, I should find it and save in labelDict.
        hapsList = []
        if not self.nodes_plan_a or label in self.nodes_plan_a:
            for haplotype, hap_label in self.Vertices_attributes.items():
                hap_label = hap_label[0]
                if hap_label == label:
                    hapsList.append(haplotype)
        elif label in self.nodes_plan_b:
            for haplotype, hap_label in self.Whole_Vertices_attributes.items():
                hap_label = hap_label[0]
                if hap_label == label:
                    hapsList.append(haplotype)
        self.labelDict[label] = hapsList
        return hapsList

    def haps_with_probs_by_label(self, label):
        """Find the haplotypes just like the above function but with the haplotypes' probabilities.
        Does not use the graphical features of the haplotypes.
        Returns a dictionary of haplotype to frequency.
        Approved."""
        dictAlleles = {}
        listLabel = self.haps_by_label(label)
        if not self.nodes_plan_a or label in self.nodes_plan_a:
            for allele in listLabel:
                dictAlleles[allele] = self.Vertices_attributes[allele][1]
        elif label in self.nodes_plan_b:
            for allele in listLabel:
                dictAlleles[allele] = self.Whole_Vertices_attributes[allele][1]
        return dictAlleles

    def adjs_query(self, alleleList):
        """A filtering query on the alleles in the graph.
        Does use the graph class.
        Returns a dictionary of haplotypes (can be partial) to frequencies.
        Approved."""
        adjDict = dict()
        for allele in alleleList:
            if allele in self.Vertices_attributes:
                allele_attributes = self.Vertices_attributes[allele][0:2]
                if allele_attributes[0] == self.full_loci:
                    adjDict[allele] = allele_attributes[1]
                else:
                    allele_id = self.Vertices_attributes[allele][2]
                    # Find the neighbors of the allele
                    allele_neighbors = self.Vertices[
                        self.Edges[
                            range(
                                self.Neighbors_start[allele_id],
                                self.Neighbors_start[allele_id + 1],
                            )
                        ]
                    ]
                    # The frequencies of the neighbors to the dictionary
                    for adj in allele_neighbors:
                        adjDict[adj] = self.Vertices_attributes[adj][1]
        return adjDict

    def adjs_query_by_color(self, alleleList, labelA, labelB):
        """A filtering query on the alleles in the graph.
        Does use the graph class.
        Returns a dictionary of haplotypes (can be partial) to frequencies.
        Approved."""
        adjDict = dict()
        if labelA == labelB:
            return self.node_probs(alleleList, labelA)

        for allele in alleleList:
            if allele in self.Whole_Vertices_attributes:
                alleles = []
                connector = labelB + allele

                if connector in self.Whole_Vertices_attributes:
                    connector_id = self.Whole_Vertices_attributes[connector]
                    alleles = self.Whole_Vertices[
                        self.Whole_Edges[
                            range(
                                self.Whole_Neighbors_start[connector_id],
                                self.Whole_Neighbors_start[connector_id + 1],
                            )
                        ]
                    ]

                for adj in alleles:
                    adjDict[adj] = self.Whole_Vertices_attributes[adj][1]
        return adjDict

    def node_probs(self, nodes, label):
        """Get a list of haplotypes and a label,
        Return a dictionary of nodes and their proper frequency."""
        nodesDict = {}
        if not self.nodes_plan_b or label in self.nodes_plan_b:
            for node in nodes:
                if node in self.Whole_Vertices_attributes:
                    nodesDict[node] = self.Whole_Vertices_attributes[node][1]
        elif label in self.nodes_plan_a:
            for node in nodes:
                if node in self.Whole_Vertices_attributes:
                    nodesDict[node] = self.Vertices_attributes[node][1]
        return nodesDict
