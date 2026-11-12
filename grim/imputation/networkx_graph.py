"""
Yuli Tshuva
Create a data structure of directed Graph suitable for our needs with numpy arrays and dictionaries.
"""

import csv
import gc
import sys
import time
from array import array
from contextlib import contextmanager

import numpy as np
from tqdm.auto import tqdm


@contextmanager
def _phase(desc):
    """Progress marker for the (otherwise silent) post-edge-load graph-build
    phases — the steps that run AFTER "All Edges Creation" and account for most
    of a large graph's load time. Prints '<desc> ...' then '<desc> done (Ns)'.

    Unlike the per-row tqdm bars (which are suppressed by TQDM_DISABLE because
    of their '\\r' spam), these emit ONE line per phase, so they're always shown
    — that's the whole point: to reveal where a slow load is spending time.
    Uses tqdm.write so it interleaves cleanly if any bars are active.
    """
    tqdm.write(f"[graph build] {desc} ...")
    t0 = time.time()
    try:
        yield
    finally:
        tqdm.write(f"[graph build] {desc} done ({time.time() - t0:.1f}s)")


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
        """Build the graph by scanning the three files, line by line and filling up the class variables.

        Memory notes: edges are accumulated as packed uint64 keys
        ``(src << 32) | dst`` in ``array('Q')`` buffers (8 bytes/edge) instead of
        Python lists of 2-item lists of strings (~150-200 bytes/edge), which
        also removes the string->uint32 ``np.vstack`` copy. Sorting the keys is
        identical to the previous ``np.lexsort((col1, col0))``, and ``np.unique``
        on the keys replaces pandas ``drop_duplicates`` + lexsort for the whole
        edges. The resulting arrays are identical to the previous implementation.
        """
        nodesDict = dict()  # csv node id (int) -> haplotype name
        with open(nodesFile) as nodesfile:
            readNodes = csv.reader(nodesfile, delimiter=",")
            next(readNodes)
            for row in tqdm(readNodes, desc="Vertices Creation:"):
                if len(row) > 0:
                    name = row[1]
                    # Labels repeat for every node of a loci combination; share one object.
                    label = sys.intern(row[2])
                    in_plan_a = not self.nodes_plan_a or label in self.nodes_plan_a
                    in_plan_b = not self.nodes_plan_b or label in self.nodes_plan_b
                    # Parse once; plan A and plan B attributes share the same (read-only) list.
                    freqs = (
                        list(map(float, row[3].split(";")))
                        if in_plan_a or in_plan_b
                        else None
                    )
                    if in_plan_a:
                        self.Vertices.append(name)
                        vertex_id = len(self.Vertices) - 1
                        self.Vertices_attributes[name] = (label, freqs, vertex_id)

                    if in_plan_b:
                        self.Whole_Vertices.append(name)
                        vertex_id = len(self.Whole_Vertices) - 1
                        self.Whole_Vertices_attributes[name] = (label, freqs, vertex_id)

                    nodesDict[int(row[0])] = name

        edge_keys = array("Q")
        # Add edges from file
        with open(edgesFile) as edgesfile:
            readEdges = csv.reader(edgesfile, delimiter=",")
            next(readEdges)
            for row in tqdm(readEdges, desc="Edges Creation:"):
                if len(row) > 0:
                    node1_id = int(row[0])
                    node2_id = int(row[1])
                    node1 = nodesDict[node1_id]
                    node2 = nodesDict[node2_id]
                    if (
                        node1 in self.Vertices_attributes
                        and node2 in self.Vertices_attributes
                    ):
                        node1_label = self.Vertices_attributes[node1][0]
                        if node1_label == self.full_loci:
                            edge_keys.append((node2_id << 32) | node1_id)
                        else:
                            edge_keys.append((node1_id << 32) | node2_id)

        whole_edge_keys = array("Q")
        # add edges from file
        with open(allEdgesFile) as allEdgesfile:
            readEdges = csv.reader(allEdgesfile, delimiter=",")
            next(readEdges)
            for row in tqdm(readEdges, "All Edges Creation:"):
                if len(row) > 0:
                    node1_id = int(row[0])
                    node2_id = int(row[1])
                    node1 = nodesDict[node1_id]
                    node2 = nodesDict[node2_id]
                    node1_attr = self.Whole_Vertices_attributes[node1]
                    node2_attr = self.Whole_Vertices_attributes[node2]
                    node1_label = node1_attr[0]
                    node2_label = node2_attr[0]
                    # Whole_Edges is consumed as positional indices into
                    # Whole_Vertices (self.Whole_Vertices[self.Whole_Edges[...]]),
                    # so edges must reference vertex POSITIONS, not the raw CSV
                    # node ids. The position of each real vertex is stored as the
                    # third element of its attribute tuple.
                    node1_pos = node1_attr[2]
                    node2_pos = node2_attr[2]

                    if len(node1_label) < len(node2_label):
                        # Create a connector
                        connector = node2_label + node1

                        if connector not in self.Whole_Vertices_attributes:
                            self.Whole_Vertices.append(connector)
                            connector_id = len(self.Whole_Vertices) - 1
                            self.Whole_Vertices_attributes[connector] = connector_id

                            whole_edge_keys.append((node1_pos << 32) | connector_id)
                        else:
                            connector_id = self.Whole_Vertices_attributes[connector]

                        # Append the connector to the whole edges array
                        whole_edge_keys.append((connector_id << 32) | node2_pos)

                    else:
                        # Create a connector
                        connector = node1_label + node2

                        if connector not in self.Whole_Vertices_attributes:
                            self.Whole_Vertices.append(connector)
                            connector_id = len(self.Whole_Vertices) - 1
                            self.Whole_Vertices_attributes[connector] = connector_id

                        # Append the connector to the whole edges array
                        whole_edge_keys.append((node2_pos << 32) | connector_id)
                        whole_edge_keys.append((connector_id << 32) | node1_pos)

        nodesDict.clear()
        del nodesDict

        # Sort the packed (src, dst) keys: same order as lexsort by (src, dst).
        with _phase("Sorting edges"):
            keys = np.frombuffer(edge_keys, dtype=np.uint64).copy()
            del edge_keys
            keys.sort()
            edges_src = (keys >> np.uint64(32)).astype(np.uint32)
            self.Edges = (keys & np.uint64(0xFFFFFFFF)).astype(np.uint32)
            del keys

        # Drop duplicate whole edges and sort them in one step.
        with _phase("Dropping duplicate whole edges + sorting"):
            keys = np.unique(np.frombuffer(whole_edge_keys, dtype=np.uint64))
            del whole_edge_keys
            whole_edges_src = (keys >> np.uint64(32)).astype(np.uint32)
            self.Whole_Edges = (keys & np.uint64(0xFFFFFFFF)).astype(np.uint32)
            del keys

        self.Vertices = np.array(self.Vertices, dtype=np.object_)
        self.Whole_Vertices = np.array(self.Whole_Vertices, dtype=np.object_)

        # Create a list of the first appearance of a number in the 0 column in the matrix
        with _phase("Indexing neighbors (edges)"):
            self.Neighbors_start = self._neighbors_start(
                edges_src, self.Vertices.shape[0], "Neighbors index:"
            )
            del edges_src

        with _phase("Indexing neighbors (whole edges)"):
            self.Whole_Neighbors_start = self._neighbors_start(
                whole_edges_src, self.Whole_Vertices.shape[0], "Whole neighbors index:"
            )
            del whole_edges_src

        gc.collect()

    @staticmethod
    def _neighbors_start(src, n_vertices, desc):
        """Start offset of each vertex's block in the (sorted) edge source column,
        plus a trailing ``n_vertices`` sentinel, as a uint32 array.

        Kept loop-for-loop identical to the original implementation (including
        how vertices without outgoing edges inherit the previous offset).
        """
        unique_values, first_occurrences_indices = np.unique(src, return_index=True)
        starts = []
        j = 0
        for i in tqdm(range(0, n_vertices), desc=desc):
            if int(unique_values[j]) == i:
                starts.append(int(first_occurrences_indices[j]))
                j += 1
            else:
                try:
                    starts.append(starts[-1])
                except IndexError:  # In case of the start of the list - empty list
                    starts.append(0)
        del unique_values, first_occurrences_indices
        starts.append(int(n_vertices))
        return np.array(starts, dtype=np.uint32)

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
                # Whole_Vertices_attributes holds two kinds of values: real
                # haplotype vertices map to a (label, probs, id) tuple, while
                # connector nodes added during graph build map to a bare int id.
                # Only the real haplotypes carry a label, so skip connectors.
                if not isinstance(hap_label, tuple):
                    continue
                if hap_label[0] == label:
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
