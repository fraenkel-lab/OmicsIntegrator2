#!/usr/bin/env python3

# Core python modules
import sys
import os
import unittest

# python external libraries
import numpy as np
import pandas as pd
import networkx as nx

sys.path.insert(0,'/Users/alex/Documents/OmicsIntegrator2/src')
import graph as oi


# try disconnected interactome
# try interactome with repeated edges
# try empty graph
# try interactome file with no connected components except pairs
# try interactome with more than 3 columns

class TestOi2(unittest.TestCase):
	"""
	"""

	tmp_interactome_filename = 'tmp_test_graph.pickle'
	tmp_prize_filename = 'tmp_test_graph.pickle'

	number_of_nodes = 1000
	p_nodes_connected = 0.1
	number_of_prized_nodes = 100

	# create a random graph with networkx with random weights in 0-1
	# write it as an interactome pickle file
	def setUp(self):
		g = nx.gnp_random_graph(self.number_of_nodes, self.p_nodes_connected)
		df = nx.to_pandas_edgelist(g, 'protein1', 'protein2')
		number_of_edges = df.shape[0]
		df['cost'] = np.random.uniform(0, 1, number_of_edges)
		df.to_pickle(self.tmp_interactome_filename)

		prizes = pd.Series(np.random.uniform(0, 1, number_of_nodes)).to_frame().sample(number_of_prized_nodes)
		prizes.to_csv(tmp_prize_filename)


	def test_init(self):

		graph = oi.Graph(self.tmp_interactome_filename, {})

		# Test the following attributes:

			# graph.interactome_dataframe
			# graph.interactome_graph
			# graph.nodes
			# graph.edges
			# graph.edge_costs
			# graph.node_degrees


			# graph.params
			# graph.edge_penalties
			# graph.costs


		self.assertTrue(True)


	def test__reset_hyperparameters(self):
		graph._reset_hyperparameters(params={})



	def test_prepare_prizes(self):
		graph.prepare_prizes(prize_file)



	def test__prepare_prizes(self):
		graph._prepare_prizes(prizes_dataframe)



	def test__add_dummy_node(self):
		graph._add_dummy_node(connected_to=[])



	def test__check_validity_of_instance(self):
		graph._check_validity_of_instance(edges, prizes, costs)



	def test_pcsf(self):
		graph.pcsf(pruning="strong", verbosity_level=0)



	def test_output_forest_as_networkx(self):
		graph.output_forest_as_networkx(vertex_indices, edge_indices)



	def test_pcsf_objective_value(self):
		graph.pcsf_objective_value(forest)



	def test__noisy_edges(self):
		graph._noisy_edges()



	def test__random_terminals(self):
		graph._random_terminals()



	def test__aggregate_pcsf(self):
		graph._aggregate_pcsf(results, frequency_attribute_name="frequency")



	def test__noisy_edges_reps(self):
		graph._noisy_edges_reps(reps)



	def test__random_terminal_reps(self):
		graph._random_terminal_reps(reps)



	def test_randomizations(self):
		graph.randomizations(noisy_edges_reps=0, random_terminals_reps=0)



	def test__eval_PCSF_runs(self):
		graph._eval_PCSF_runs(params)



	def test_grid_randomization(self):
		graph.grid_randomization(prize_file, Ws, Bs, Gs, noisy_edges_reps, random_terminals_reps)



	def test_grid_search(self):
		graph.grid_search(prize_file, Ws, Bs, Gs)



	def test_betweenness(self):
		oi.betweenness(nxgraph)



	def test_louvain_clustering(self):
		oi.louvain_clustering(nxgraph)



	def test_edge_betweenness_clustering(self):
		oi.edge_betweenness_clustering(nxgraph)



	def test_k_clique_clustering(self):
		oi.k_clique_clustering(nxgraph, k)



	def test_spectral_clustering(self):
		oi.spectral_clustering(nxgraph, k)



	def test_augment_with_all_GO_terms(self):
		oi.augment_with_all_GO_terms(nxgraph)



	def test_augment_with_subcellular_localization(self):
		oi.augment_with_subcellular_localization(nxgraph)



	def test_augment_with_biological_process_terms(self):
		oi.augment_with_biological_process_terms(nxgraph)



	def test_augment_with_molecular_function_terms(self):
		oi.augment_with_molecular_function_terms(nxgraph)



	def test_perform_GO_enrichment_on_clusters(self):
		oi.perform_GO_enrichment_on_clusters(nxgraph, clustering)



	def test_summarize_grid_search(self):
		oi.summarize_grid_search(results, mode, top_n=np.Infinity)



	def test_get_robust_subgraph_from_randomizations(self):
		oi.get_robust_subgraph_from_randomizations(nxgraph, max_size=400, min_component_size=5)



	def test_filter_graph_by_component_size(self):
		oi.filter_graph_by_component_size(nxgraph, min_size=5)



	def test_get_networkx_graph_as_dataframe_of_nodes(self):
		oi.get_networkx_graph_as_dataframe_of_nodes(nxgraph)



	def test_get_networkx_graph_as_dataframe_of_edges(self):
		oi.get_networkx_graph_as_dataframe_of_edges(nxgraph)



	def test_output_networkx_graph_as_pickle(self):
		oi.output_networkx_graph_as_pickle(nxgraph, output_dir, filename="pcsf_results.pickle")



	def test_output_networkx_graph_as_graphml_for_cytoscape(self):
		oi.output_networkx_graph_as_graphml_for_cytoscape(nxgraph, output_dir, filename="pcsf_results.graphml.gz")



	def test_output_networkx_graph_as_interactive_html(self):
		oi.output_networkx_graph_as_interactive_html(nxgraph, output_dir, filename="graph.html")







	def tearDown(self):
		# remove the temporary files we created


		os.remove(self.tmp_interactome_filename)
		os.remove(self.tmp_prize_filename)








if __name__ == '__main__':
    unittest.main()
