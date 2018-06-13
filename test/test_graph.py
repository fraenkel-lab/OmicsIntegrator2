#!/usr/bin/env python3

# Core python modules
import sys
import os
from pathlib import Path
import logging
# import pytest


# python external libraries
import numpy as np
import pandas as pd
import networkx as nx

sys.path.insert(0,'/Users/alex/Documents/OmicsIntegrator2/src')
import graph as oi
oi.logger.setLevel(logging.WARNING)

# try disconnected interactome
# try interactome with repeated edges
# try empty graph
# try interactome file with no connected components except pairs
# try interactome with more than 3 columns

# @pytest.mark.incremental
class Test_Oi2(object):
    """
    """

    tmp_interactome_filepath = Path.cwd() / 'tmp_test_graph.pickle'
    tmp_prize_filepath = Path.cwd() / 'tmp_test_prizes.pickle'
    tmp_files = [tmp_interactome_filepath, tmp_prize_filepath]

    number_of_nodes = 1000
    p_nodes_connected = 0.1
    number_of_prized_nodes = 100

    # create a random graph with networkx with random weights in 0-1
    # write it as an interactome pickle file
    def __init__(self):
        self.g = nx.gnp_random_graph(self.number_of_nodes, self.p_nodes_connected)
        self.df = nx.to_pandas_edgelist(self.g, 'protein1', 'protein2')
        self.number_of_edges = self.df.shape[0]
        self.df['cost'] = np.random.uniform(0, 1, self.number_of_edges)
        self.df.to_pickle(self.tmp_interactome_filepath)

        self.prizes = pd.Series(np.random.uniform(0, 3, self.number_of_nodes)).to_frame().sample(self.number_of_prized_nodes)
        self.terminals = self.prizes.index.values
        self.prizes.to_csv(self.tmp_prize_filepath, sep='\t')


    ###########################################################################
                #######          Initialization            #######
    ###########################################################################

    def test_init(self):
        self.graph = oi.Graph(self.tmp_interactome_filepath, {})

        assert hasattr(self.graph, "interactome_dataframe")  # test that this is deep equal to self.df"
        assert hasattr(self.graph, "interactome_graph")  # test that this is deep equal to self.g"
        assert len(self.graph.nodes) == self.number_of_nodes
        assert len(self.graph.edges) == self.number_of_edges
        assert len(self.graph.edge_costs) == self.number_of_edges
        assert len(self.graph.node_degrees) == self.number_of_nodes

        assert hasattr(self.graph, "params")
        assert hasattr(self.graph, "edge_penalties")
        assert hasattr(self.graph, "costs")

        print("pass test_init()")


    def test__reset_hyperparameters(self):
        params = {"w":5, "b":2, "g":2, "noise":0.1, "dummy_mode":"terminals", "seed":0}
        self.graph._reset_hyperparameters(params)

        assert self.graph.params.w == params['w']
        assert self.graph.params.b == params['b']
        assert self.graph.params.g == params['g']
        assert self.graph.params.noise == params['noise']
        assert self.graph.params.dummy_mode == params['dummy_mode']
        assert self.graph.params.seed == params['seed']

        assert hasattr(self.graph, "edge_penalties")
        assert hasattr(self.graph, "costs")

        print("pass test__reset_hyperparameters()")


    def test_prepare_prizes(self):
        self.graph.prepare_prizes(self.tmp_prize_filepath)

        assert hasattr(self.graph, "node_attributes")
        assert hasattr(self.graph, "bare_prizes")
        assert hasattr(self.graph, "prizes")
        assert hasattr(self.graph, "terminals")

        print("pass test_prepare_prizes()")


    ###########################################################################
                #######              PCSF               #######
    ###########################################################################

    def test__add_dummy_node(self):
        dummy_edges, dummy_costs, dummy_id, dummy_prize = self.graph._add_dummy_node(connected_to=self.terminals)

        assert dummy_id <= self.number_of_nodes
        assert np.array_equal(dummy_costs, np.array([self.graph.params.w] * self.number_of_prized_nodes))
        assert set(map(frozenset, dummy_edges.tolist())) == set([frozenset((dummy_id, node_id)) for node_id in self.terminals])

        print("pass test__add_dummy_node()")


    def test__check_validity_of_instance(self):

        # self.graph._check_validity_of_instance(edges, prizes, costs)
        assert True

        print("pass test__check_validity_of_instance()")


    def test_pcsf(self):
        self.vertex_indices, self.edge_indices = self.graph.pcsf()

        assert isinstance(self.vertex_indices, np.ndarray)
        assert isinstance(self.edge_indices, np.ndarray)

        # check that both are proper subsets of the originals
        # check that no dummy data is in the arrays
        # maybe check that the mean selected cost is less than the mean cost and the mean selected prize is higher than the mean prize?

        print("pass test_pcsf()")


    def test_output_forest_as_networkx(self):
        self.forest, self.augmented_forest = self.graph.output_forest_as_networkx(*self.graph.pcsf())

        assert isinstance(self.forest, nx.Graph)
        assert isinstance(self.augmented_forest, nx.Graph)

        # test that the forests contain all the nodes in the vertex_indices
        # test that the forsts have some basic set of attributes we care about

        print("pass test_output_forest_as_networkx()")


    def test_pcsf_objective_value(self):
        objective_value = self.graph.pcsf_objective_value(self.forest)

        assert objective_value >= 0

        print("pass test_pcsf_objective_value()")


    ###########################################################################
                #######          Randomziations         #######
    ###########################################################################

    def test_randomizations(self):
        forest, augmented_forest = self.graph.randomizations(noisy_edges_reps=3, random_terminals_reps=3)

        assert isinstance(forest, nx.Graph)
        assert isinstance(augmented_forest, nx.Graph)

        # test that robustness and specificity are set

        forest, augmented_forest = self.graph.randomizations(noisy_edges_reps=0, random_terminals_reps=3)

        assert isinstance(forest, nx.Graph)
        assert isinstance(augmented_forest, nx.Graph)

        # test that only one of the two is set

        forest, augmented_forest = self.graph.randomizations(noisy_edges_reps=3, random_terminals_reps=0)

        assert isinstance(forest, nx.Graph)
        assert isinstance(augmented_forest, nx.Graph)

        # test that only one of the two is set

        forest, augmented_forest = self.graph.randomizations(noisy_edges_reps=0, random_terminals_reps=0)

        assert isinstance(forest, nx.Graph)
        assert isinstance(augmented_forest, nx.Graph)

        # test that you get a regular pcsf solution from this.

        print("pass test_randomizations()")


    ###########################################################################
                #######          Grid Search          #######
    ###########################################################################

    def test_grid_randomization(self):
        Ws = [4,5]
        Bs = [1,2]
        Gs = [3,4]
        results = self.graph.grid_randomization(self.tmp_prize_filepath, Ws=Ws, Bs=Bs, Gs=Gs, noisy_edges_reps=2, random_terminals_reps=2)
        # unknown what should be tested here
        print("pass test_grid_randomization()")


    def test_grid_search(self):
        Ws = [4,5]
        Bs = [1,2]
        Gs = [3,4]
        results = self.graph.grid_search(self.tmp_prize_filepath, Ws=Ws, Bs=Bs, Gs=Gs)
        # unknown what should be tested here
        print("pass test_grid_search()")


    ###########################################################################
                #######        Subgraph Augmentation    #######
    ###########################################################################

    def test_betweenness(self):
        oi.betweenness(self.g)
        assert set(nx.get_node_attributes(self.g, "betweenness").keys()) == set(self.g.nodes())
        # assert nx.get_node_attributes(self.g, "betweenness").values() ==

        print("pass test_betweenness()")


    def test_louvain_clustering(self):
        oi.louvain_clustering(self.g)
        assert set(nx.get_node_attributes(self.g, "louvainClusters").keys()) == set(self.g.nodes())
        # assert nx.get_node_attributes(self.g, "louvainClusters").values() ==

        print("pass test_louvain_clustering()")


    def test_k_clique_clustering(self):
        oi.k_clique_clustering(self.g, 10)
        assert set(nx.get_node_attributes(self.g, "kCliqueClusters").keys()) == set(self.g.nodes())
        # assert nx.get_node_attributes(self.g, "kCliqueClusters").values() ==

        print("pass test_k_clique_clustering()")


    def test_spectral_clustering(self):
        oi.spectral_clustering(self.g, 10)
        assert set(nx.get_node_attributes(self.g, "spectralClusters").keys()) == set(self.g.nodes())
        # assert nx.get_node_attributes(self.g, "spectralClusters").values() ==

        print("pass test_spectral_clustering()")


    def test_augment_with_all_GO_terms(self):
        # oi.augment_with_all_GO_terms(nxgraph)

        oi.augment_with_all_GO_terms(self.g)

        # test that the graph has a location column columns
        assert True


        print("pass test_augment_with_all_GO_terms()")


    ###########################################################################
                #######            Results           #######
    ###########################################################################

    def test_summarize_grid_search(self):
        oi.summarize_grid_search(results, "membership")
        oi.summarize_grid_search(results, "robustness")
        oi.summarize_grid_search(results, "specificity")

        # unknown, but needs tests
        print("pass test_summarize_grid_search()")


    def test_get_robust_subgraph_from_randomizations(self):
        oi.get_robust_subgraph_from_randomizations(nxgraph, max_size=400, min_component_size=5)

        # unknown, but needs tests
        print("pass test_get_robust_subgraph_from_randomizations()")


    def test_filter_graph_by_component_size(self):
        oi.filter_graph_by_component_size(nxgraph, min_size=5)

        # unknown, but needs tests
        print("pass test_filter_graph_by_component_size()")


    ###########################################################################
                #######              Export             #######
    ###########################################################################

    def test_get_networkx_graph_as_dataframe_of_nodes(self):
        oi.get_networkx_graph_as_dataframe_of_nodes(nxgraph)

        # test that all features on the graph are also on the df
        # test that the length of the graph is the length of the df
        print("pass test_get_networkx_graph_as_dataframe_of_nodes()")


    def test_get_networkx_graph_as_dataframe_of_edges(self):
        oi.get_networkx_graph_as_dataframe_of_edges(nxgraph)

        # test that the number of edges in the graph is the length of the df

        print("pass test_get_networkx_graph_as_dataframe_of_edges()")


    def test_output_networkx_graph_as_pickle(self):
        path = oi.output_networkx_graph_as_pickle(nxgraph, ".", filename="pcsf_results.pickle")
        self.tmp_files.append(path)

        # stat the file, assert it exists
        print("pass test_output_networkx_graph_as_pickle()")


    def test_output_networkx_graph_as_graphml_for_cytoscape(self):
        path = oi.output_networkx_graph_as_graphml_for_cytoscape(nxgraph, ".", filename="pcsf_results.graphml.gz")
        self.tmp_files.append(path)

        # stat the file, assert it exists
        print("pass test_output_networkx_graph_as_graphml_for_cytoscape()")


    def test_output_networkx_graph_as_interactive_html(self):
        path = oi.output_networkx_graph_as_interactive_html(nxgraph, ".", filename="graph.html")
        self.tmp_files.append(path)

        # stat the file, assert it exists
        print("pass test_output_networkx_graph_as_interactive_html()")







    def tearDown(self):
        # remove the temporary files we created
        for file in self.tmp_files:
            os.remove(file)



if __name__ == '__main__':

    test = Test_Oi2()

    ###  Initialization
    test.test_init()
    test.test__reset_hyperparameters()      # note that this is a private method
    test.test_prepare_prizes()
    test.test__add_dummy_node()             # note that this is a private method

    ###  PCSF
    test.test__check_validity_of_instance() # note that this is a private method
    test.test_pcsf()
    test.test_output_forest_as_networkx()
    test.test_pcsf_objective_value()

    ###  Randomziations
    # test.test_randomizations()

    ###  Grid Search
    # test.test_grid_randomization()
    # test.test_grid_search()

    ###  Subgraph Augmentation
    test.test_betweenness()
    test.test_louvain_clustering()
    test.test_k_clique_clustering()
    test.test_spectral_clustering()
    # test.test_augment_with_all_GO_terms()

    ###  Results
    test.test_summarize_grid_search()
    test.test_get_robust_subgraph_from_randomizations()
    test.test_filter_graph_by_component_size()

    ###  Export
    test.test_get_networkx_graph_as_dataframe_of_nodes()
    test.test_get_networkx_graph_as_dataframe_of_edges()
    test.test_output_networkx_graph_as_pickle()
    test.test_output_networkx_graph_as_graphml_for_cytoscape()
    test.test_output_networkx_graph_as_interactive_html()


    test.tearDown()


