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

# TODO: try disconnected interactome
# TODO: try interactome with repeated edges
# TODO: try empty graph
# TODO: try interactome file with no connected components except pairs
# TODO: try interactome with more than 3 columns

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
        self.df = nx.to_pandas_edgelist(self.g, 'protein1', 'protein2').astype(int)
        self.number_of_edges = self.df.shape[0]
        self.df['cost'] = np.random.uniform(0, 1, self.number_of_edges).astype(float)
        self.df.to_csv(self.tmp_interactome_filepath, sep='\t', index=False)

        self.prizes = pd.Series(np.random.uniform(0, 3, self.number_of_nodes)).to_frame().sample(self.number_of_prized_nodes).astype(float)
        self.terminals = self.prizes.index.values
        self.prizes.to_csv(self.tmp_prize_filepath, sep='\t')


    ###########################################################################
                #######          Initialization            #######
    ###########################################################################

    def test_init(self):
        print("test_init")
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
        print("test__reset_hyperparameters")
        reset = self.graph._reset_hyperparameters

        with pytest.raises(ValueError): reset({'w': -1})
        with pytest.raises(ValueError): reset({'w': "1"})

        with pytest.raises(ValueError): reset({'b': -1})
        with pytest.raises(ValueError): reset({'b': "1"})

        with pytest.raises(ValueError): reset({'g': -1})
        with pytest.raises(ValueError): reset({'g': "1"})

        with pytest.raises(ValueError): reset({'edge_noise': -1})
        with pytest.raises(ValueError): reset({'edge_noise': "1"})

        with pytest.raises(ValueError): reset({'dummy_mode': -1})
        with pytest.raises(ValueError): reset({'dummy_mode': "1"})
        with pytest.raises(ValueError): reset({'dummy_mode': " all"})

        with pytest.raises(ValueError): reset({'seed': [1]})
        with pytest.raises(ValueError): reset({'seed': "1"})

        params = {"w":5, "b":2, "g":2, "edge_noise":0.1, "dummy_mode":"terminals", "seed":0, "skip_checks":False}
        reset(params)

        assert self.graph.params.w == params['w']
        assert self.graph.params.b == params['b']
        assert self.graph.params.g == params['g']
        assert self.graph.params.edge_noise == params['edge_noise']
        assert self.graph.params.dummy_mode == params['dummy_mode']
        assert self.graph.params.seed == params['seed']

        assert hasattr(self.graph, "edge_penalties")
        assert hasattr(self.graph, "costs")

        print("...pass")


    def test_prepare_prizes(self):
        print("test_prepare_prizes")
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
        print("test__add_dummy_node")
        dummy_edges, dummy_costs, dummy_id, dummy_prize = self.graph._add_dummy_node(connected_to=self.terminals)

        assert dummy_id <= self.number_of_nodes
        assert np.array_equal(dummy_costs, np.array([self.graph.params.w] * self.number_of_prized_nodes))
        assert set(map(frozenset, dummy_edges.tolist())) == set([frozenset((dummy_id, node_id)) for node_id in self.terminals])

        print("...pass")


    def test__check_validity_of_instance(self):
        print("test__check_validity_of_instance")

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
        print("test_pcsf")
        self.vertex_indices, self.edge_indices = self.graph.pcsf()


        assert isinstance(self.vertex_indices, np.ndarray)
        assert isinstance(self.edge_indices, np.ndarray)

        assert ((0 <= self.vertex_indices) & (self.vertex_indices < self.number_of_nodes)).all()
        assert set(self.vertex_indices) == set(np.unique(self.vertex_indices))

        assert ((0 <= self.edge_indices) & (self.edge_indices < self.number_of_edges)).all()
        assert set(self.edge_indices) == set(np.unique(self.edge_indices))

        print("...pass")


    def test_output_forest_as_networkx(self):
        print("test_output_forest_as_networkx")
        self.forest, self.augmented_forest = self.graph.output_forest_as_networkx(*self.graph.pcsf())

        assert isinstance(self.forest, nx.Graph)
        assert isinstance(self.augmented_forest, nx.Graph)

        assert set(self.forest.nodes()) == set(self.graph.nodes[self.vertex_indices])
        assert set(self.augmented_forest.nodes()) == set(self.graph.nodes[self.vertex_indices])

        print("...pass")


    def test_pcsf_objective_value(self):
        print("test_pcsf_objective_value")
        objective_value = self.graph.pcsf_objective_value(self.forest)

        assert objective_value >= 0

        print("...pass")


    ###########################################################################
                #######          Randomziations         #######
    ###########################################################################

    def test_randomizations(self):
        print("test_randomizations")
        forest, augmented_forest = self.graph.randomizations(noisy_edges_reps=3, random_terminals_reps=3)

        assert isinstance(forest, nx.Graph)
        assert isinstance(augmented_forest, nx.Graph)

        assert set(nx.get_node_attributes(forest,           "robustness").keys())  != set()
        assert set(nx.get_node_attributes(augmented_forest, "robustness").keys())  != set()

        assert set(nx.get_node_attributes(forest,           "specificity").keys()) != set()
        assert set(nx.get_node_attributes(augmented_forest, "specificity").keys()) != set()

        forest, augmented_forest = self.graph.randomizations(noisy_edges_reps=0, random_terminals_reps=3)

        assert isinstance(forest, nx.Graph)
        assert isinstance(augmented_forest, nx.Graph)

        assert set(nx.get_node_attributes(forest,           "robustness").keys())  == set()
        assert set(nx.get_node_attributes(augmented_forest, "robustness").keys())  == set()

        assert set(nx.get_node_attributes(forest,           "specificity").keys()) != set()
        assert set(nx.get_node_attributes(augmented_forest, "specificity").keys()) != set()

        forest, augmented_forest = self.graph.randomizations(noisy_edges_reps=3, random_terminals_reps=0)

        assert isinstance(forest, nx.Graph)
        assert isinstance(augmented_forest, nx.Graph)

        assert set(nx.get_node_attributes(forest,           "robustness").keys())  != set()
        assert set(nx.get_node_attributes(augmented_forest, "robustness").keys())  != set()

        assert set(nx.get_node_attributes(forest,           "specificity").keys()) == set()
        assert set(nx.get_node_attributes(augmented_forest, "specificity").keys()) == set()

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
        print("test_grid_randomization")
        Ws = [1,2]
        Bs = [1,2]
        Gs = [3,4]
        self.results = self.graph.grid_randomization(self.tmp_prize_filepath, Ws=Ws, Bs=Bs, Gs=Gs, noisy_edges_reps=2, random_terminals_reps=2)
        # unknown what should be tested here
        print("...pass")


    def test_grid_search(self):
        print("test_grid_search")
        Ws = [1,2]
        Bs = [1,2]
        Gs = [3,4]
        results = self.graph.grid_search(self.tmp_prize_filepath, Ws=Ws, Bs=Bs, Gs=Gs)
        # unknown what should be tested here
        print("...pass")


    ###########################################################################
                #######        Subgraph Augmentation    #######
    ###########################################################################

    def test_betweenness(self):
        print("test_betweenness")
        oi.betweenness(self.g)
        assert set(nx.get_node_attributes(self.g, "betweenness").keys()) == set(self.g.nodes())
        assert all([isinstance(betweenness, float) for betweenness in nx.get_node_attributes(self.g, "betweenness").values()])

        print("...pass")


    def test_louvain_clustering(self):
        print("test_louvain_clustering")
        oi.louvain_clustering(self.g)
        assert set(nx.get_node_attributes(self.g, "louvain_clusters").keys()) == set(self.g.nodes())
        assert all([isinstance(louvainClusters, str) for louvainClusters in nx.get_node_attributes(self.g, "louvainClusters").values()])

        print("...pass")


    def test_k_clique_clustering(self):
        print("test_k_clique_clustering")
        oi.k_clique_clustering(self.g, 3)
        assert set(nx.get_node_attributes(self.g, "k_clique_clusters").keys()) == set(self.g.nodes())
        assert all([isinstance(k_clique_clusters, str) for k_clique_clusters in nx.get_node_attributes(self.g, "k_clique_clusters").values()])

        print("...pass")


    def test_spectral_clustering(self):
        print("test_spectral_clustering")
        oi.spectral_clustering(self.g, 10)
        assert set(nx.get_node_attributes(self.g, "spectral_clusters").keys()) == set(self.g.nodes())
        assert all([isinstance(spectral_clusters, str) for spectral_clusters in nx.get_node_attributes(self.g, "spectral_clusters").values()])

        print("...pass")


    def test_annotate_graph_nodes(self):
        print("test_annotate_graph_nodes")
        oi.annotate_graph_nodes(self.g)
        # unknown what should be tested here
        print("...pass")


    ###########################################################################
                #######            Results           #######
    ###########################################################################

    def test_summarize_grid_search(self):
        print("test_summarize_grid_search")
        node_summary_df = oi.summarize_grid_search(self.results, "membership")
        node_summary_df = oi.summarize_grid_search(self.results, "robustness")
        node_summary_df = oi.summarize_grid_search(self.results, "specificity")
        # unknown what should be tested here
        print("...pass")


    def test_get_robust_subgraph_from_randomizations(self):
        print("test_get_robust_subgraph_from_randomizations")
        oi.get_robust_subgraph_from_randomizations(nxgraph, max_size=400, min_component_size=5)
        # unknown what should be tested here
        print("...pass")


    def test_filter_graph_by_component_size(self):
        print("test_filter_graph_by_component_size")
        oi.filter_graph_by_component_size(nxgraph, min_size=5)
        # unknown what should be tested here
        print("...pass")


    ###########################################################################
                #######              Export             #######
    ###########################################################################

    def test_get_networkx_graph_as_dataframe_of_nodes(self):
        print("test_get_networkx_graph_as_dataframe_of_nodes")
        df = oi.get_networkx_graph_as_dataframe_of_nodes(self.augmented_forest)
        assert len(df) == len(self.augmented_forest.nodes())
        assert set(df.columns.tolist()) == set(list(self.augmented_forest.nodes(data=True))[0][1].keys())

        df = oi.get_networkx_graph_as_dataframe_of_nodes(self.forest)
        assert len(df) == len(self.augmented_forest.nodes())

        print("...pass")


    def test_get_networkx_graph_as_dataframe_of_edges(self):
        print("test_get_networkx_graph_as_dataframe_of_edges")
        df = oi.get_networkx_graph_as_dataframe_of_edges(self.augmented_forest)



        df = oi.get_networkx_graph_as_dataframe_of_edges(self.forest)
        # test that the number of edges in the graph is the length of the df

        print("...pass")


    def test_output_networkx_graph_as_pickle(self):
        print("test_output_networkx_graph_as_pickle")
        path = oi.output_networkx_graph_as_pickle(self.augmented_forest)
        self.tmp_files.append(path)

        # assert nx.read_gpickle(path) == self.augmented_forest  ## need to find graph deep equals function

        print("...pass")


    def test_output_networkx_graph_as_graphml_for_cytoscape(self):
        print("test_output_networkx_graph_as_graphml_for_cytoscape")
        path = oi.output_networkx_graph_as_graphml_for_cytoscape(self.augmented_forest)
        self.tmp_files.append(path)

        # assert nx.read_graphml(path) == self.augmented_forest  ## need to find graph deep equals function

        print("...pass")


    def test_output_networkx_graph_as_interactive_html(self):
        print("test_output_networkx_graph_as_interactive_html")
        path = oi.output_networkx_graph_as_interactive_html(self.augmented_forest)
        self.tmp_files.append(path)

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
        test.test_annotate_graph_nodes()

        ###  Results
        test.test_summarize_grid_search()  # TODO: @iamjli
        # test.test_get_robust_subgraph_from_randomizations()  # TODO: @iamjli
        # test.test_filter_graph_by_component_size()  # TODO: @iamjli

        ###  Export
        test.test_get_networkx_graph_as_dataframe_of_nodes()
        test.test_get_networkx_graph_as_dataframe_of_edges()
        test.test_output_networkx_graph_as_pickle()
        # test.test_output_networkx_graph_as_graphml_for_cytoscape()  # TODO: @zfrenchee requires networkx fix
        test.test_output_networkx_graph_as_interactive_html()


    finally:
        test.tearDown()


