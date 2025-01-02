# -*- coding: utf-8 -*-

#
#    grim Graph Imputation
#    Copyright (c) 2021 Be The Match operated by National Marrow Donor Program. All Rights Reserved.
#
#    This library is free software; you can redistribute it and/or modify it
#    under the terms of the GNU Lesser General Public License as published
#    by the Free Software Foundation; either version 3 of the License, or (at
#    your option) any later version.
#
#    This library is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; with out even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
#    License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with this library;  if not, write to the Free Software Foundation,
#    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA.
#
#    > http://www.fsf.org/licensing/licenses/lgpl.html
#    > http://www.opensource.org/licenses/lgpl-license.php
#


from .imputation.impute import Imputation
from .imputation.networkx_graph import Graph

import sys
import os

# adding Folder_2 to the system path
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)).replace("/grim", ""))


from graph_generation import generate_neo4j_multi_hpf
from grim.run_impute_def import run_impute


def graph_freqs(conf_file="", for_em=False, em_pop=None):
    use_default_path = False
    if conf_file == "":
        use_default_path = True
        conf_file = (
            os.path.dirname(os.path.realpath(__file__)).replace("/grim", "")
            + "/conf/minimal-configuration.json"
        )

    generate_neo4j_multi_hpf.generate_graph(
        config_file=conf_file,
        em_pop=em_pop,
        em=for_em,
        use_default_path=use_default_path,
    )


def impute(conf_file="", hap_pop_pair = False, graph = None):

    project_dir_in_file, project_dir_graph = "", ""
    if conf_file == "":

        conf_file = (
            os.path.dirname(os.path.realpath(__file__)).replace("/grim", "")
            + "/conf/minimal-configuration.json"
        )
        project_dir_graph = (
            os.path.dirname(os.path.realpath(__file__)).replace("/grim", "")
            + "/graph_generation/"
        )
        project_dir_in_file = (
            os.path.dirname(os.path.realpath(__file__)).replace("/grim", "") + "/"
        )
    graph = run_impute(conf_file, project_dir_graph, project_dir_in_file, hap_pop_pair, graph)
    return graph


def impute_instance(config, graph, count_by_prob=None):
    imputation = Imputation(graph, config, count_by_prob)
    return imputation


def graph_instance(config):
    graph = Graph(config)
    graph.build_graph(
        config["node_file"], config["top_links_file"], config["edges_file"]
    )
    return graph
