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


from .validation import runfile
from .imputation.graph_generation import generate_neo4j_multi_hpf

from .imputation.imputegl import Imputation
from .imputation.imputegl.networkx_graph import Graph
import os

def graph_freqs(conf_file = "",  for_em = False,  em_pop=None ):
    if conf_file == "":
        conf_file = os.path.dirname(os.path.realpath(__file__))  + '/conf/minimal-configuration.json'


    generate_neo4j_multi_hpf.generate_graph(config_file=conf_file,  em_pop=em_pop,
                       em=for_em)

def impute(conf_file = ""):
    project_dir_in_file, project_dir_graph = "", ""
    if conf_file == "":
        conf_file = os.path.dirname(os.path.realpath(__file__)) + '/conf/minimal-configuration.json'
        project_dir_graph =  os.path.dirname(os.path.realpath(__file__)) + '/imputation/graph_generation/'
        project_dir_in_file =  os.path.dirname(os.path.realpath(__file__)) + '/'
    runfile.run_impute(conf_file, project_dir_graph, project_dir_in_file)

def impute_instance(config, graph, count_by_prob= None):
    imputation = Imputation(graph, config, count_by_prob)
    return imputation

def graph_instance(config):
    graph = Graph(config)
    graph.build_graph(config["node_file"], config["top_links_file"], config["edges_file"])
    return graph


