#!/usr/bin/env python3

# Core python modules
import sys
import os

# Peripheral python modules
import argparse
import logging
import random
from collections import Counter
from itertools import product
from copy import copy
import json

# python external libraries
import numpy as np
import pandas as pd
import networkx as nx
import community    # pip install python-louvain
from py2cytoscape import util as cy

# Lab modules
from pcst_fast import pcst_fast

# list of classes and methods we'd like to export:
__all__ = [ "Graph",
			"output_networkx_graph_as_graphml_for_cytoscape",
			"output_networkx_graph_as_json_for_cytoscapejs",
			"get_networkx_graph_as_dataframe_of_nodes",
			"get_networkx_graph_as_dataframe_of_edges" ]


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
handler.setFormatter(logging.Formatter('%(asctime)s - Graph: %(levelname)s - %(message)s', "%I:%M:%S"))
logger.addHandler(handler)


# Helpers
def flatten(list_of_lists): return [item for sublist in list_of_lists for item in sublist]

def invert(list_of_lists): return {item: i for i, list in enumerate(list_of_lists) for item in list}

class Options(object):
	def __init__(self, options):
		self.__dict__.update(options)
	def __repr__(self):
		return dict(self.__dict__)


class Graph:
	"""
	A Graph object is a representation of a graph, with convenience methods for using the pcst_fast
	package, which computes an approximate minimization Prize-Collecting Steiner Forest objective.
	"""
	def __init__(self, interactome_file, params):
		"""
		Builds a representation of a graph from an interactome file.

		From the interactome_file, populates
		- `graph.interactome_dataframe` (pandas.DataFrame)
		- `graph.interactome_graph` (networkx.Graph)
		- `graph.nodes` (pandas.Index),
		- `graph.edges` (list of pairs),
		- `graph.costs` and `graph.edge_penalties` (lists, such that the ordering is the same as in graph.edges),
		- `graph.node_degrees` (list, such that the ordering is the same as in graph.nodes).

		Arguments:
			interactome_file (str or FILE): tab-delimited text file containing edges in interactome and their weights formatted like "ProteinA\tProteinB\tCost"
			params (dict): params with which to run the program
		"""

		interactome_fieldnames = ["source","target","cost"]
		self.interactome_dataframe = pd.read_csv(interactome_file, sep='\t', names=interactome_fieldnames)

		self.interactome_dataframe['temp'] = self.interactome_dataframe.apply(lambda row: ''.join(sorted([row['source'], row['target']])), axis=1)
		logger.info("Duplicated edges in the interactome file (we'll keep the max cost):")
		logger.info(self.interactome_dataframe[self.interactome_dataframe.set_index('temp').index.duplicated()][['source','target']].values.tolist())
		self.interactome_dataframe = self.interactome_dataframe.groupby('temp').max().reset_index()[["source","target","cost"]]

		self.interactome_graph = nx.from_pandas_dataframe(self.interactome_dataframe, 'source', 'target', edge_attr=['cost'])

		# We first take only the source and target columns from the interactome dataframe.
		# We then unstack them, which, unintuitively, stacks them into one column, allowing us to use factorize.
		# Factorize builds two datastructures, a unique pd.Index of each ID string to a numerical ID
		# and the datastructure we passed it with ID strings replaced with those numerical IDs.
		# We place those in self.nodes and self.edges respectively, but self.edges will need reshaping.
		(self.edges, self.nodes) = pd.factorize(self.interactome_dataframe[["source","target"]].unstack())

		# Here we do the inverse operation of "unstack" above, which gives us an interpretable edges datastructure
		self.edges = self.edges.reshape(self.interactome_dataframe[["source","target"]].shape, order='F')
		self.edge_costs = self.interactome_dataframe['cost'].astype(float).values

		# Numpy has a convenient counting function. However we're assuming here that each edge only appears once.
		# The indices into this datastructure are the same as those in self.nodes and self.edges.
		self.node_degrees = np.bincount(self.edges.flatten())

		self._reset_hyperparameters(params)


	def _reset_hyperparameters(self, params):
		"""
		Set the parameters on Graph and compute parameter-dependent features.

		Arguments:
			params (dict): params with which to run the program
		"""

		defaults = {"w": 6, "b": 1, "g": 20, "noise": 0.1, "exclude_terminals": False, "dummy_mode": "terminals", "knockout": [], "seed": None}

		self.params = Options({**defaults, **params})

		self._knockout(self.params.knockout)

		N = len(self.nodes)
		self.edge_penalties = self.params.g * np.array([self.node_degrees[a] * self.node_degrees[b] /
							((N - self.node_degrees[a] - 1) * (N - self.node_degrees[b] - 1) + self.node_degrees[a] * self.node_degrees[b]) for a, b in self.edges])

		self.costs = (self.edge_costs + self.edge_penalties)

		# In case specific nodes are penalized, we'd like to update the costs accordingly
		if hasattr(self, "additional_costs"): self.costs = (self.edge_costs + self.edge_penalties + self.additional_costs)

		# if this instance of graph has bare_prizes set, then presumably resetting the
		# hyperparameters should also reset the scaled prizes
		if hasattr(self, "bare_prizes"): self.prizes = self.bare_prizes * self.params.b


	def _penalize_nodes(self, nodes_to_penalize):
		"""
		Penalize a set of nodes by penalizing the edges connected to that node by some coefficient.

		Arguments:
			nodes_to_penalize (pandas.DataFrame): 2 columns: 'name' and 'penalty' with entries in [0, 1)
		"""

		# Here's we're indexing the penalized nodes by the indices we used for nodes during initialization
		nodes_to_penalize.index = graph.nodes.get_indexer(nodes_to_penalize['name'].values)

		# there will be some nodes in the penalty dataframe which we don't have in our interactome
		logger.info("Members of the penalty dataframe not present in the interactome:")
		logger.info(nodes_to_penalize[nodes_to_penalize.index == -1]['name'].tolist())
		nodes_to_penalize.drop(-1, inplace=True, errors='ignore')

		if not nodes_to_penalize['penalty_coefficient'].between(0, 1).all():
			logger.info("The node penalty coefficients must lie in [0, 1]. Passing..."); return

		nodes_to_knockout = nodes_to_penalize[nodes_to_penalize.penalty_coefficient == 1]
		logger.info("penalty coefficients of 1 are treated as knockouts. Proteins to remove from interactome:")
		logger.info(nodes_to_knockout.name.tolist())
		self._knockout(nodes_to_knockout.name.values)

		self.additional_costs = np.zeros(self.costs.shape)

		# Iterate through the rows of the nodes_to_penalize dataframe
		for index, name, penalty_coefficient in nodes_to_penalize.itertuples():
			# For each protein we'd like to penalize, get the indicies of the edges connected to that node
			edge_indices = self.interactome_dataframe[(self.interactome_dataframe.source == name) | (self.interactome_dataframe.target == name)].index
			# And compute an additional cost on those edges.
			self.additional_costs[edge_indices] += self.edge_costs[edge_indices] / (1 - penalty_coefficient)
		# Apply those additional costs by calling _reset_hyperparameters.
		self._reset_hyperparameters({})


	def _knockout(self, nodes_to_knockout):
		"""
		Knock out a set of nodes from the interactome, effectively removing them from results.

		Arguments:
			nodes_to_knockout (numpy.array): Array of string IDs of nodes to knock out.
		"""
		if len(nodes_to_knockout) > 0:
			logger.info("The knockout function has yet to be implemented, passing...");

		return


	def prepare_prizes(self, prize_file):
		"""
		Parses a prize file and returns an array of prizes, a list of terminal indices.

		This function logs duplicate assignments in the prize file and memebers of the prize file
		not found in the interactome.

		This file passed to this function must have at least two columns: node name and prize.
		Any additional columns will be assumed to be node attributes. However, in order to know
		the names of those attributes, this function now requires the input file contain headers,
		i.e. the first row of the tsv must be the names of the columns.

		Sets the graph attributes
		- `graph.prizes` (numpy.array): properly indexed
		- `graph.terminals` (numpy.array): their indices
		- `graph.node_attributes` (pandas.DataFrame)

		Arguments:
			prize_file (str or FILE): a filepath or file object containing a tsv **with headers**.
		"""

		prizes_dataframe = pd.read_csv(prize_file, sep='\t')
		prizes_dataframe.columns = ['name', 'prize'] + prizes_dataframe.columns[2:].tolist()  # TODO: error handling

		return self._prepare_prizes(prizes_dataframe)


	def _prepare_prizes(self, prizes_dataframe):

		# Strangely some files have duplicated genes, sometimes with different prizes. Keep the max prize.
		logger.info("Duplicated gene symbols in the prize file (we'll keep the max prize):")
		logger.info(prizes_dataframe[prizes_dataframe.set_index('name').index.duplicated()]['name'].tolist())
		prizes_dataframe = prizes_dataframe.groupby('name').max().reset_index()

		# Here's we're indexing the terminal nodes and associated prizes by the indices we used for nodes
		prizes_dataframe.set_index(self.nodes.get_indexer(prizes_dataframe['name']), inplace=True)

		# there will be some nodes in the prize file which we don't have in our interactome
		logger.info("Members of the prize file not present in the interactome:")
		logger.info(prizes_dataframe[prizes_dataframe.index == -1]['name'].tolist())
		prizes_dataframe.drop(-1, inplace=True, errors='ignore')

		self.node_attributes = prizes_dataframe.set_index('name').rename_axis(None)

		# Here we're making a dataframe with all the nodes as keys and the prizes from above or 0
		prizes_dataframe = pd.DataFrame(self.nodes, columns=["name"]).merge(prizes_dataframe, on="name", how="left").fillna(0)
		# Our return value is a 1D array, where each entry is a node's prize, indexed as above
		self.bare_prizes = prizes_dataframe['prize'].values
		self.prizes = self.bare_prizes * self.params.b

		self.terminals = pd.Series(self.prizes).nonzero()[0].tolist()


	def _add_dummy_node(self, connected_to=[]):

		dummy_id = len(self.nodes)
		dummy_prize = np.array([0])
		dummy_edges = np.array([(dummy_id, node_id) for node_id in connected_to])
		dummy_costs = np.array([self.params.w] * len(dummy_edges))

		return dummy_edges, dummy_costs, dummy_id, dummy_prize


	def _check_validity_of_instance(self, edges, prizes, costs):
		"""
		Assert that the parammeters and files passed to this program are valid, log useful error messages otherwise.
		"""
		pass


	def pcsf(self, pruning="strong", verbosity_level=0):
		"""
		Select the subgraph which approximately optimizes the Prize-Collecting Steiner Forest objective.

		This function mostly defers to pcst_fast, but does one important pre-processing step: it
		adds a dummy node which will serve as the PCSF root and connects that dummy node to either
		terminals, non-terminals, or all other nodes with edges weighted by self.params.w.

		In order to interpret the results of this function, use `output_forest_as_networkx` with
		the results of this function.

		Arguments:
			pruning (str): TODO
			verbosity_level (int): TODO

		Returns:
			numpy.array: indices of the selected vertices
			numpy.array: indices of the selected edges
		"""

		all = list(range(len(self.nodes)))
		others = list(set(all) - set(self.terminals))

		if self.params.dummy_mode == 'terminals': endpoints = self.terminals
		elif self.params.dummy_mode == 'other': endpoints = others
		elif self.params.dummy_mode == 'all': endpoints = all
		else: sys.exit("Invalid dummy mode")

		dummy_edges, dummy_costs, root, dummy_prize = self._add_dummy_node(connected_to=endpoints)

		# `edges`: a 2D int64 array. Each row (of length 2) specifies an undirected edge in the input graph. The nodes are labeled 0 to n-1, where n is the number of nodes.
		edges = np.concatenate((self.edges, dummy_edges))
		# `prizes`: the node prizes as a 1D float64 array.
		prizes = np.concatenate((self.prizes, dummy_prize))
		# `costs`: the edge costs as a 1D float64 array.
		costs = np.concatenate((self.costs, dummy_costs))
		# `root`: the root note for rooted PCST. For the unrooted variant, this parameter should be -1.
		# `num_clusters`: the number of connected components in the output.
		num_clusters = 1
		# `pruning`: a string value indicating the pruning method. Possible values are `'none'`, `'simple'`, `'gw'`, and `'strong'` (all literals are case-insensitive).
		# `verbosity_level`: an integer indicating how much debug output the function should produce.

		self._check_validity_of_instance(edges, prizes, costs)

		vertex_indices, edge_indices = pcst_fast(edges, prizes, costs, root, num_clusters, pruning, verbosity_level)
		# `vertex_indices`: indices of the vertices in the solution as a 1D int64 array.
		# `edge_indices`: indices of the edges in the output as a 1D int64 array. The list contains indices into the list of edges passed into the function.

		# Remove the dummy node and dummy edges for convenience
		vertex_indices = vertex_indices[vertex_indices != root]
		edge_indices = edge_indices[np.in1d(edge_indices, self.interactome_dataframe.index)]

		return vertex_indices, edge_indices


	def output_forest_as_networkx(self, vertex_indices, edge_indices):
		"""

		Arguments:
			vertex_indices (list): indices of the vertices selected in self.nodes
			edge_indices (list): indices of the edges selected in self.edges

		Returns:
			networkx.Graph: a networkx graph object
		"""

		# Replace the edge indices with the actual edges (source name, target name) by indexing into the interactome
		edges = self.interactome_dataframe.loc[edge_indices]
		forest = nx.from_pandas_dataframe(edges, 'source', 'target', edge_attr=True)
		# the above won't capture the singletons, so we'll add them here
		forest.add_nodes_from(list(set(self.nodes[vertex_indices]) - set(forest.nodes())))

		for attribute in self.node_attributes.columns.values:
			nx.set_node_attributes(forest, attribute, {node: attr for node, attr in self.node_attributes[attribute].to_dict().items() if node in forest.nodes()})

		# Add the degree as an attribute
		node_degree_dict = nx.degree(self.interactome_graph)
		nx.set_node_attributes(forest, 'degree', {node: degree for node, degree in node_degree_dict.items() if node in forest.nodes()})

		augmented_forest = nx.compose(self.interactome_graph.subgraph(forest.nodes()), forest)

		# Post-processing
		louvain_clustering(augmented_forest)

		return forest, augmented_forest


	def pcsf_objective_value(self, forest):
		"""
		Calculate PCSF objective function

		Arguments:
			forest (networkx.Graph): a forest like the one returned by output_forest_as_networkx -- Not an augmented forest!

		Returns:
			float: PCSF objective function score
		"""

		return (sum(self.prizes) - sum(nx.get_node_attributes(forest, 'prize').values())) + sum(nx.get_edge_attributes(forest, 'cost').values()) + (self.params.w * nx.number_connected_components(forest))


	def _noisy_edges(self):
		"""
		Adds gaussian noise to all edges in the graph

		Generate gaussian noise values, mean=0, stdev default=0.333 (edge values range between 0 and 1)

		Returns:
			numpy.array: edge weights with gaussian noise
		"""

		return np.clip(np.random.normal(self.costs, self.params.noise), 0.0001, None)  # None means don't clip above


	def _random_terminals(self):
		"""
		Switches the terminams with random nodes with a similar degree distribution.

		Returns:
			numpy.array: new prizes
			numpy.array: new terminals
		"""

		if len(self.edges) < 50: sys.exit("Cannot use random_terminals with such a small interactome.")

		nodes_sorted_by_degree = pd.Series(self.node_degrees).sort_values().index
		terminal_degree_rankings = np.array([nodes_sorted_by_degree.get_loc(terminal) for terminal in self.terminals])
		new_terminal_degree_rankings = np.clip(np.rint(np.random.normal(terminal_degree_rankings, 10)), 0, len(self.nodes)-1).astype(int)
		new_terminals = pd.Series(nodes_sorted_by_degree)[new_terminal_degree_rankings].values

		new_prizes = copy(self.prizes)

		for old_terminal, new_terminal in zip(self.terminals, new_terminals):
			new_prizes[old_terminal] = 0
			new_prizes[new_terminal] = self.prizes[old_terminal]

		return new_prizes, np.unique(new_terminals)


	def _aggregate_pcsf(self, results, frequency_attribute_name):
		"""
		Merge multiple PCSF results into one DataFrame

		Arguments:
			results (list): a list of [(vertex_indices, edge_indices),...] from multiple PCSF runs.
			frequency_attribute_name (str): Name of the attribute relating to the frequency of occurrence of components in the results.

		Returns:
			pandas.DataFrame: vertex indices and their fractional rate of occurrence in the PCSF results
			pandas.DataFrame: edge indices and their fractional rate of occurrence in the PCSF results
		"""

		# Transposes a list from [(vertex_indices, edge_indices),...] to ([vertex_indices,...], [edge_indices,...])
		vertex_indices, edge_indices = zip(*results)

		# These next steps are just data transformation/aggregation.
		# 1. Flatten the lists of lists of edge indices and vertex indices
		# 2. Count the occurrences of each edge and vertex index
		# 3. Transform from Counter object to DataFrame through list
		vertex_indices = pd.DataFrame(list(Counter(flatten(vertex_indices)).items()), columns=['node_index',frequency_attribute_name])
		edge_indices = pd.DataFrame(list(Counter(flatten(edge_indices)).items()), columns=['edge_index',frequency_attribute_name])
		# 4. Convert occurrences to fractions
		vertex_indices[frequency_attribute_name] /= len(results)
		edge_indices[frequency_attribute_name] /= len(results)

		return vertex_indices, edge_indices


	def randomizations(self, noisy_edges_reps, random_terminals_reps):
		"""
		Macro function which performs randomizations and merges the results

		Note that these are additive, not multiplicative:
		`noisy_edges_reps` = 5 and `random_terminals_reps` = 5 makes 10 PCSF runs, not 25.

		Arguments:
			noisy_edges_reps (int): Number of "Noisy Edges" type randomizations to perform
			random_terminals_reps (int): Number of "Random Terminals" type randomizations to perform

		Returns:
			networkx.Graph: forest
			networkx.Graph: augmented_forest
		"""

		if self.params.seed: random.seed(seed); numpy.random.seed(seed=seed)

		#### NOISY EDGES ####
		if noisy_edges_reps > 0:

			results = []

			true_edge_costs = copy(self.costs)

			for noisy_edge_costs in [self._noisy_edges() for rep in range(noisy_edges_reps)]:
				self.costs = noisy_edge_costs
				results.append(self.pcsf())

			self.costs = true_edge_costs

			robust_vertices, robust_edges = self._aggregate_pcsf(results, 'robustness')

		#### RANDOM TERMINALS ####
		if random_terminals_reps > 0:

			results = []

			true_prizes = self.prizes
			true_terminals = self.terminals

			for random_prizes, terminals in [self._random_terminals() for rep in range(random_terminals_reps)]:
				self.prizes = random_prizes
				self.terminals = terminals
				results.append(self.pcsf())

			self.prizes = true_prizes
			self.terminals = true_terminals

			specific_vertices, specific_edges = self._aggregate_pcsf(results, 'specificity')

		###########

		if random_terminals_reps == 0 and noisy_edges_reps > 0:
			vertex_indices = robust_vertices; edge_indices = robust_edges;

		elif noisy_edges_reps == 0 and random_terminals_reps > 0:
			vertex_indices = specific_vertices; edge_indices = specific_edges;

		elif noisy_edges_reps > 0 and random_terminals_reps > 0:
			vertex_indices = robust_vertices.merge(specific_vertices, how='outer', on='node_index')
			edge_indices = robust_edges.merge(specific_edges, how='outer', on='edge_index')

		else: sys.exit("Randomizations was called with invalid noisy_edges_reps and random_terminals_reps.")

		###########

		forest, augmented_forest = self.output_forest_as_networkx(vertex_indices.node_index.values, edge_indices.edge_index.values)

		vertex_indices.index = self.nodes[vertex_indices.node_index.values]

		if noisy_edges_reps > 0:
			nx.set_node_attributes(forest, 			 'robustness', vertex_indices['robustness'].to_dict())
			nx.set_node_attributes(augmented_forest, 'robustness', vertex_indices['robustness'].to_dict())

			edge_robustness_dic = edge_indices.set_index("edge_index")["robustness"].to_dict()
			nx.set_edge_attributes(forest, 'robustness', {tuple([self.nodes[x] for x in self.edges[edge]]): edge_robustness_dic[edge] for edge in edge_robustness_dic})

		if random_terminals_reps > 0:
			nx.set_node_attributes(forest, 			 'specificity', vertex_indices['specificity'].to_dict())
			nx.set_node_attributes(augmented_forest, 'specificity', vertex_indices['specificity'].to_dict())

		# TODO we aren't yet using edge_indices which contain robustness and specificity information.
		return forest, augmented_forest


	def _eval_pcsf(self, params):
		"""
		Convenience methods which sets parameters and performs PCSF
		"""

		self._reset_hyperparameters(params)
		paramstring = 'G_'+str(params['g'])+'_B_'+str(params['b'])+'_W_'+str(params['w'])
		logger.info(paramstring)
		return (paramstring, self.pcsf())


	def _grid_pcsf(self, prize_file, Gs, Bs, Ws):
		"""
		Internal function which executes pcsf at every point in a parameter grid.
		Subroutine of `grid_search`.

		Arguments:
			prize_file (str): filepath
			Gs (list): Values of gamma
			Bs (list): Values of beta
			Ws (list): Values of omega

		Returns:
			list: list of tuples of vertex indices and edge indices
		"""

		self.prepare_prizes(prize_file)

		parameter_permutations = [{'g':g,'b':b,'w':w} for (g, b, w) in product(Gs, Bs, Ws)]

		results = list(map(self._eval_pcsf, parameter_permutations))

		return results


	def grid_search(self, prize_file, Gs, Bs, Ws):
		"""
		Macro function which performs grid search and merges the results.

		This function is under construction and subject to change.

		Arguments:
			prize_file (str): filepath
			Gs (list): Values of gamma
			Bs (list): Values of beta
			Ws (list): Values of omega

		Returns:
			networkx.Graph: forest
			networkx.Graph: augmented_forest
			pd.DataFrame: parameters and node membership lists
		"""

		results = self._grid_pcsf(prize_file, Gs, Bs, Ws)

		### GET THE REGULAR OUTPUT ###
		vertex_indices, edge_indices = self._aggregate_pcsf(list(dict(results).values()), 'frequency')

		forest, augmented_forest = self.output_forest_as_networkx(vertex_indices.node_index.values, edge_indices.edge_index.values)

		vertex_indices.index = self.nodes[vertex_indices.node_index.values]

		nx.set_node_attributes(forest, 			 'frequency', vertex_indices['frequency'].to_dict())
		nx.set_node_attributes(augmented_forest, 'frequency', vertex_indices['frequency'].to_dict())

		### GET THE OUTPUT NEEDED BY TOBI'S VISUALIZATION ###
		params_by_nodes = pd.DataFrame({paramstring: dict(zip(self.nodes[vertex_indices], self.node_degrees[vertex_indices])) for paramstring, (vertex_indices, edge_indices) in results}).fillna(0)

		return forest, augmented_forest, params_by_nodes


def betweenness(nxgraph):
	nx.set_node_attributes(nxgraph, 'betweenness', nx.betweenness_centrality(nxgraph))


def louvain_clustering(nxgraph):
	"""
	"""
	nx.set_node_attributes(nxgraph, 'louvainClusters', community.best_partition(nxgraph))

# def edge_betweenness_clustering(nxgraph):  # is coming with NetworkX 2.0, to be released soon.
# 	"""
# 	"""
# 	nx.set_node_attributes(nxgraph, 'edgeBetweennessClusters', invert(nx.girvan_newman(nxgraph)))

def k_clique_clustering(nxgraph, k):
	"""
	"""
	nx.set_node_attributes(nxgraph, 'kCliqueClusters', invert(nx.k_clique_communities(nxgraph, k)))


def get_networkx_subgraph_from_randomizations(nxgraph, max_size=400, min_component_size=5): 
	"""
	Approach 1: from entire network, attempt to remove lowest robustness node. If removal results in a component 
	of size less than min_size, do not remove. 
	"""
	# TODO
	return 


def get_networkx_graph_as_dataframe_of_nodes(nxgraph):
	"""
	Arguments:
		nxgraph (networkx.Graph): any instance of networkx.Graph

	Returns:
		pd.DataFrame: nodes from the input graph and their attributes as a dataframe
	"""

	return pd.DataFrame.from_dict(dict(nxgraph.nodes(data=True))).transpose().fillna(0)


def get_networkx_graph_as_dataframe_of_edges(nxgraph):
	"""
	Arguments:
		nxgraph (networkx.Graph): any instance of networkx.Graph

	Returns:
		pd.DataFrame: edges from the input graph and their attributes as a dataframe
	"""

	intermediate = pd.DataFrame(nxgraph.edges(data=True))
	intermediate.columns = ['protein1', 'protein2'] + intermediate.columns[2:].tolist()
	# TODO: in the future, get the other attributes out into columns
	return intermediate[['protein1', 'protein2']]


def output_networkx_graph_as_graphml_for_cytoscape(nxgraph, output_dir, filename):
	"""
	Arguments:
		nxgraph (networkx.Graph): any instance of networkx.Graph
		output_dir (str): the directory in which to output the graph.
		filename (str): Filenames ending in .gz or .bz2 will be compressed.
	"""
	os.makedirs(os.path.abspath(output_dir), exist_ok=True)
	path = os.path.join(os.path.abspath(output_dir), filename)
	nx.write_graphml(nxgraph, path)


def output_networkx_graph_as_json_for_cytoscapejs(nxgraph, output_dir, filename="graph_json.json"):
	"""
	Arguments:
		nxgraph (networkx.Graph): any instance of networkx.Graph
		output_dir (str): the directory in which to output the file (named graph_json.json)
	"""
	path = os.path.join(os.path.abspath(output_dir), filename)
	njs = cy.from_networkx(nxgraph)
	with open(path,'w') as outf:
		outf.write(json.dumps(njs, indent=4))

