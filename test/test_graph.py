#!/usr/bin/env python3

# Core python modules
import sys
import os
from pathlib import Path
import logging
import pytest


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
        print("test_init", end='')
        self.graph = oi.Graph(self.tmp_interactome_filepath, {})

        assert hasattr(self.graph, "interactome_dataframe")
        assert hasattr(self.graph, "interactome_graph")
        assert len(self.graph.nodes) == self.number_of_nodes
        assert len(self.graph.edges) == self.number_of_edges
        assert len(self.graph.edge_costs) == self.number_of_edges
        assert len(self.graph.node_degrees) == self.number_of_nodes

        assert hasattr(self.graph, "params")
        assert hasattr(self.graph, "edge_penalties")
        assert hasattr(self.graph, "costs")

        print("...pass")


    def test__reset_hyperparameters(self):
        print("test__reset_hyperparameters", end='')
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

        print("...pass")


    def test_prepare_prizes(self):
        print("test_prepare_prizes", end='')
        self.graph.prepare_prizes(self.tmp_prize_filepath)

        assert hasattr(self.graph, "node_attributes")
        assert hasattr(self.graph, "bare_prizes")
        assert hasattr(self.graph, "prizes")
        assert hasattr(self.graph, "terminals")

        print("...pass")


    ###########################################################################
                #######              PCSF               #######
    ###########################################################################

    def test__add_dummy_node(self):
        print("test__add_dummy_node", end='')
        dummy_edges, dummy_costs, dummy_id, dummy_prize = self.graph._add_dummy_node(connected_to=self.terminals)

        assert dummy_id <= self.number_of_nodes
        assert np.array_equal(dummy_costs, np.array([self.graph.params.w] * self.number_of_prized_nodes))
        assert set(map(frozenset, dummy_edges.tolist())) == set([frozenset((dummy_id, node_id)) for node_id in self.terminals])

        print("...pass")


    def test__check_validity_of_instance(self):
        print("test__check_validity_of_instance", end='')

        edges = self.graph.edges
        prizes = self.graph.prizes
        costs = self.graph.costs
        root = 0
        num_clusters = 1
        pruning = "strong"
        verbosity_level = 0
        func_params = [edges, prizes, costs, root, num_clusters, pruning, verbosity_level]

        assert self.graph._check_validity_of_instance(*func_params)

        check = self.graph._check_validity_of_instance

        # Test malformed edges
        with pytest.raises(ValueError): check(                  edges.tolist(),         *func_params[1:])
        with pytest.raises(ValueError): check(                  np.expand_dims(edges,1),*func_params[1:])
        with pytest.raises(ValueError): check(                  edges[:,1],             *func_params[1:])

        # Test malformed prizes
        with pytest.raises(ValueError): check(edges,            prizes.tolist(),        *func_params[2:])
        with pytest.raises(ValueError): check(edges,            np.expand_dims(prizes,1),*func_params[2:])
        with pytest.raises(ValueError): check(edges,            prizes[:-1],            *func_params[2:])

        # Test malformed costs
        with pytest.raises(ValueError): check(*func_params[:2], costs.tolist(),         *func_params[3:])
        with pytest.raises(ValueError): check(*func_params[:2], np.expand_dims(costs,1),*func_params[3:])
        with pytest.raises(ValueError): check(*func_params[:2], costs[:-1],             *func_params[3:])

        # Test malformed root
        with pytest.raises(ValueError): check(*func_params[:3], -1,                     *func_params[4:])
        with pytest.raises(ValueError): check(*func_params[:3], self.number_of_nodes+2, *func_params[4:])
        with pytest.raises(ValueError): check(*func_params[:3], "0",                    *func_params[4:])

        # Test malformed num_clusters
        with pytest.raises(ValueError): check(*func_params[:4], 0,                      *func_params[5:])
        with pytest.raises(ValueError): check(*func_params[:4], -1,                     *func_params[5:])
        with pytest.raises(ValueError): check(*func_params[:4], prizes+1,               *func_params[5:])
        with pytest.raises(ValueError): check(*func_params[:4], "1",                    *func_params[5:])
        with pytest.raises(ValueError): check(*func_params[:4], None,                   *func_params[5:])

        # Test malformed pruning
        with pytest.raises(ValueError): check(*func_params[:5], "prune",                verbosity_level)
        with pytest.raises(ValueError): check(*func_params[:5], 7,                      verbosity_level)
        with pytest.raises(ValueError): check(*func_params[:5], None,                   verbosity_level)

        # Test malformed verbosity_level
        with pytest.raises(ValueError): check(*func_params[:6], 4)
        with pytest.raises(ValueError): check(*func_params[:6], -1)
        with pytest.raises(ValueError): check(*func_params[:6], "1")

        print("...pass")


    def test_pcsf(self):
        print("test_pcsf", end='')
        self.vertex_indices, self.edge_indices = self.graph.pcsf()


        assert isinstance(self.vertex_indices, np.ndarray)
        assert isinstance(self.edge_indices, np.ndarray)

        assert ((0 <= self.vertex_indices) & (self.vertex_indices < self.number_of_nodes)).all()
        assert set(self.vertex_indices) == set(np.unique(self.vertex_indices))

        assert ((0 <= self.edge_indices) & (self.edge_indices < self.number_of_edges)).all()
        assert set(self.edge_indices) == set(np.unique(self.edge_indices))

        print("...pass")


    def test_output_forest_as_networkx(self):
        print("test_output_forest_as_networkx", end='')
        self.forest, self.augmented_forest = self.graph.output_forest_as_networkx(*self.graph.pcsf())

        assert isinstance(self.forest, nx.Graph)
        assert isinstance(self.augmented_forest, nx.Graph)

        print(self.forest.nodes())

        assert set(self.forest.nodes()) == set(self.graph.nodes[self.vertex_indices])
        assert set(self.augmented_forest.nodes()) == set(self.graph.nodes[self.vertex_indices])

        print("...pass")


    def test_pcsf_objective_value(self):
        print("test_pcsf_objective_value", end='')
        objective_value = self.graph.pcsf_objective_value(self.forest)

        assert objective_value >= 0

        print("...pass")


    ###########################################################################
                #######          Randomziations         #######
    ###########################################################################

    def test_randomizations(self):
        print("test_randomizations", end='')
        forest, augmented_forest = self.graph.randomizations(noisy_edges_reps=3, random_terminals_reps=3)

        assert isinstance(forest, nx.Graph)
        assert isinstance(augmented_forest, nx.Graph)

        print(nx.get_node_attributes(forest,           "robustness"))
        print(set(self.graph.nodes))

        assert set(nx.get_node_attributes(forest,           "robustness").keys())  == set(self.graph.nodes)
        assert set(nx.get_node_attributes(augmented_forest, "robustness").keys())  == set(self.graph.nodes)

        assert set(nx.get_node_attributes(forest,           "specificity").keys()) == set(self.graph.nodes)
        assert set(nx.get_node_attributes(augmented_forest, "specificity").keys()) == set(self.graph.nodes)

        forest, augmented_forest = self.graph.randomizations(noisy_edges_reps=0, random_terminals_reps=3)

        assert isinstance(forest, nx.Graph)
        assert isinstance(augmented_forest, nx.Graph)

        assert set(nx.get_node_attributes(forest,           "robustness").keys())  == set(self.graph.nodes)
        assert set(nx.get_node_attributes(augmented_forest, "robustness").keys())  == set(self.graph.nodes)

        assert set(nx.get_node_attributes(forest,           "specificity").keys()) == set()
        assert set(nx.get_node_attributes(augmented_forest, "specificity").keys()) == set()

        forest, augmented_forest = self.graph.randomizations(noisy_edges_reps=3, random_terminals_reps=0)

        assert isinstance(forest, nx.Graph)
        assert isinstance(augmented_forest, nx.Graph)

        assert set(nx.get_node_attributes(forest,           "robustness").keys())  == set()
        assert set(nx.get_node_attributes(augmented_forest, "robustness").keys())  == set()

        assert set(nx.get_node_attributes(forest,           "specificity").keys()) == set(self.graph.nodes)
        assert set(nx.get_node_attributes(augmented_forest, "specificity").keys()) == set(self.graph.nodes)

        forest, augmented_forest = self.graph.randomizations(noisy_edges_reps=0, random_terminals_reps=0)

        assert isinstance(forest, nx.Graph)
        assert isinstance(augmented_forest, nx.Graph)

        assert set(nx.get_node_attributes(forest,           "robustness").keys())  == set()
        assert set(nx.get_node_attributes(augmented_forest, "robustness").keys())  == set()

        assert set(nx.get_node_attributes(forest,           "specificity").keys()) == set()
        assert set(nx.get_node_attributes(augmented_forest, "specificity").keys()) == set()

        print("...pass")


    ###########################################################################
                #######          Grid Search          #######
    ###########################################################################

    def test_grid_randomization(self):
        print("test_grid_randomization", end='')
        Ws = [4,5]
        Bs = [1,2]
        Gs = [3,4]
        self.results = self.graph.grid_randomization(self.tmp_prize_filepath, Ws=Ws, Bs=Bs, Gs=Gs, noisy_edges_reps=2, random_terminals_reps=2)
        # unknown what should be tested here
        print("...pass")


    def test_grid_search(self):
        print("test_grid_search", end='')
        Ws = [4,5]
        Bs = [1,2]
        Gs = [3,4]
        results = self.graph.grid_search(self.tmp_prize_filepath, Ws=Ws, Bs=Bs, Gs=Gs)
        # unknown what should be tested here
        print("...pass")


    ###########################################################################
                #######        Subgraph Augmentation    #######
    ###########################################################################

    def test_betweenness(self):
        print("test_betweenness", end='')
        oi.betweenness(self.g)
        assert set(nx.get_node_attributes(self.g, "betweenness").keys()) == set(self.g.nodes())
        assert all([isinstance(betweenness, float) for betweenness in nx.get_node_attributes(self.g, "betweenness").values()])

        print("...pass")


    def test_louvain_clustering(self):
        print("test_louvain_clustering", end='')
        oi.louvain_clustering(self.g)
        assert set(nx.get_node_attributes(self.g, "louvainClusters").keys()) == set(self.g.nodes())
        assert all([isinstance(louvainClusters, str) for louvainClusters in nx.get_node_attributes(self.g, "louvainClusters").values()])

        print("...pass")


    def test_k_clique_clustering(self):
        print("test_k_clique_clustering", end='')
        oi.k_clique_clustering(self.g, 3)
        assert set(nx.get_node_attributes(self.g, "kCliqueClusters").keys()) == set(self.g.nodes())
        assert all([isinstance(kCliqueClusters, str) for kCliqueClusters in nx.get_node_attributes(self.g, "kCliqueClusters").values()])

        print("...pass")


    def test_spectral_clustering(self):
        print("test_spectral_clustering", end='')
        oi.spectral_clustering(self.g, 10)
        assert set(nx.get_node_attributes(self.g, "spectralClusters").keys()) == set(self.g.nodes())
        assert all([isinstance(spectralClusters, str) for spectralClusters in nx.get_node_attributes(self.g, "spectralClusters").values()])

        print("...pass")


    def test_augment_with_all_GO_terms(self):
        print("test_augment_with_all_GO_terms", end='')
        oi.augment_with_all_GO_terms(self.g)
        # unknown what should be tested here
        print("...pass")


    ###########################################################################
                #######            Results           #######
    ###########################################################################

    def test_summarize_grid_search(self):
        print("test_summarize_grid_search", end='')
        node_summary_df = oi.summarize_grid_search(self.results, "membership")
        node_summary_df = oi.summarize_grid_search(self.results, "robustness")
        node_summary_df = oi.summarize_grid_search(self.results, "specificity")
        # unknown what should be tested here
        print("...pass")


    def test_get_robust_subgraph_from_randomizations(self):
        print("test_get_robust_subgraph_from_randomizations", end='')
        oi.get_robust_subgraph_from_randomizations(nxgraph, max_size=400, min_component_size=5)
        # unknown what should be tested here
        print("...pass")


    def test_filter_graph_by_component_size(self):
        print("test_filter_graph_by_component_size", end='')
        oi.filter_graph_by_component_size(nxgraph, min_size=5)
        # unknown what should be tested here
        print("...pass")


    ###########################################################################
                #######              Export             #######
    ###########################################################################

    def test_get_networkx_graph_as_dataframe_of_nodes(self):
        print("test_get_networkx_graph_as_dataframe_of_nodes", end='')
        df = oi.get_networkx_graph_as_dataframe_of_nodes(self.augmented_forest)

        assert len(df) == len(self.augmented_forest.nodes())

        df = oi.get_networkx_graph_as_dataframe_of_nodes(self.forest)

        # test that all features on the graph are also on the df
        # test that the length of the graph is the length of the df
        print("...pass")


    def test_get_networkx_graph_as_dataframe_of_edges(self):
        print("test_get_networkx_graph_as_dataframe_of_edges", end='')
        df = oi.get_networkx_graph_as_dataframe_of_edges(self.augmented_forest)



        df = oi.get_networkx_graph_as_dataframe_of_edges(self.forest)
        # test that the number of edges in the graph is the length of the df

        print("...pass")


    def test_output_networkx_graph_as_pickle(self):
        print("test_output_networkx_graph_as_pickle", end='')
        path = oi.output_networkx_graph_as_pickle(self.augmented_forest, ".", filename="pcsf_results.pickle")
        self.tmp_files.append(path)

        # stat the file, assert it exists
        print("...pass")


    def test_output_networkx_graph_as_graphml_for_cytoscape(self):
        print("test_output_networkx_graph_as_graphml_for_cytoscape", end='')
        path = oi.output_networkx_graph_as_graphml_for_cytoscape(self.augmented_forest, ".", filename="pcsf_results.graphml.gz")
        self.tmp_files.append(path)

        # stat the file, assert it exists
        print("...pass")


    def test_output_networkx_graph_as_interactive_html(self):
        print("test_output_networkx_graph_as_interactive_html", end='')
        path = oi.output_networkx_graph_as_interactive_html(self.augmented_forest, ".", filename="graph.html")
        self.tmp_files.append(path)

        # stat the file, assert it exists
        print("...pass")







    def tearDown(self):
        # remove the temporary files we created
        for file in self.tmp_files:
            os.remove(file)



if __name__ == '__main__':


    test = Test_Oi2()

    try:
        ###  Initialization
        test.test_init()
        test.test__reset_hyperparameters()
        test.test_prepare_prizes()
        test.test__add_dummy_node()

        ###  PCSF
        test.test__check_validity_of_instance()
        test.test_pcsf()
        test.test_output_forest_as_networkx()
        test.test_pcsf_objective_value()

        ###  Randomziations
        test.test_randomizations()

        ###  Grid Search
        test.test_grid_randomization()
        test.test_grid_search()

        ###  Subgraph Augmentation
        test.test_betweenness()
        test.test_louvain_clustering()
        test.test_k_clique_clustering()
        test.test_spectral_clustering()
        test.test_augment_with_all_GO_terms()

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

    finally:
        test.tearDown()


