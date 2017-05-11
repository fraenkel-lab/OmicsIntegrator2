#!/usr/bin/env python3

# Core python modules
import sys
import os

# Peripheral python modules
import argparse
import pickle
import logging
import random
from collections import Counter
from copy import copy
import json

# Core python external libraries
import numpy as np
import pandas as pd

# Peripheral python external libraries
import networkx as nx
from py2cytoscape import util as cy

# Lab modules
from pcst_fast import pcst_fast

# list of classes and methods we'd like to export:
__all__ = [ "Graph",
			"output_networkx_graph_as_gml_for_cytoscape",
			"merge_two_prize_files",
			"get_networkx_graph_as_dataframe_of_nodes",
			"get_networkx_graph_as_dataframe_of_edges" ]


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
handler.setFormatter(logging.Formatter('%(asctime)s - Forest: %(levelname)s - %(message)s', "%I:%M:%S"))
logger.addHandler(handler)


# Some arbitrary helpers I believe should exist in the language anyhow
def flatten(list_of_lists): return [item for sublist in list_of_lists for item in sublist]

class Options(object):
	def __init__(self, options):
		self.__dict__.update(options)


class Graph:
	"""

	"""
	def __init__(self, interactome_file, params):
		"""
		Builds a representation of a graph from an interactome file.

		From the interactome_file, populates `self.nodes` (pd.Index), `self.edges` (list of pairs),
		`self.costs` (list, such that the ordering is the same as in self.edges), and
		`self.node_degrees` and self.negprizes (lists, such that the ordering is the same as in self.nodes).

		From the prize_file, populates self.terminals (list) and self.prizes (list which contains 0
		everywhere there isn't a terminal with an assigned prize).

		From the garnet_file, merge the TF terminals and prizes with `self.terminals` and `self.prizes`.

		Arguments:
			interactome_file (str or FILE): tab-delimited text file containing edges in interactome and their weights formatted like "ProteinA\tProteinB\tWeight"
			prize_file (str or FILE): tab-delimited text file containing all proteins with prizes formatted like "ProteinName\tPrizeValue"
			params (dict): params with which to run the program

		"""

		interactome_fieldnames = ["source","target","Weight"]
		self.interactome_dataframe = pd.read_csv(interactome_file, sep='\t', names=interactome_fieldnames)
		self.interactome_graph = nx.from_pandas_dataframe(self.interactome_dataframe, 'source', 'target', edge_attr=['Weight'])

		# We first take only the source and target columns from the interactome dataframe.
		# We then unstack them, which, unintuitively, stacks them into one column, allowing us to use factorize.
		# Factorize builds two datastructures, a unique pd.Index of each ID string to a numerical ID
		# and the datastructure we passed in with ID strings replaced with those numerical IDs.
		# We place those in self.nodes and self.edges respectively, but self.edges will need reshaping.
		(self.edges, self.nodes) = pd.factorize(self.interactome_dataframe[["source","target"]].unstack())

		# Here we do the inverse operation of "unstack" above, which gives us an interpretable edges datastructure
		self.edges = self.edges.reshape(self.interactome_dataframe[["source","target"]].shape, order='F')

		self.costs = self._determine_costs_from_interactome_file(self.interactome_dataframe['Weight'].astype(float).values)

		# Numpy has a convenient counting function. However we're assuming here that each edge only appears once.
		# The indices into this datastructure are the same as those in self.nodes and self.edges.
		self.node_degrees = np.bincount(self.edges.flatten())


		defaults = {"w": 6, "b": 1, "mu": 0, "a": 20, "noise": 0.1, "mu_squared": False, "exclude_terminals": False, "dummy_mode": "terminals", "seed": None}

		self.params = Options({**defaults, **params})

		self.negprizes = (self.node_degrees**2 if self.params.mu_squared else self.node_degrees) * self.params.mu # unless self.params.exclude_terminals TODO

		N = len(self.nodes)
		self.edge_penalties = self.params.a * np.array([self.node_degrees[a] * self.node_degrees[b] /
							((N - self.node_degrees[a] - 1) * (N - self.node_degrees[b] - 1) + self.node_degrees[a] * self.node_degrees[b]) for a, b in self.edges])

		self.costs = (self.costs + self.edge_penalties)


	def _determine_costs_from_interactome_file(self, native_costs_array): 
		flip_costs = 1 - native_costs_array
		pos_costs = flip_costs[flip_costs<0] = 0
		return pos_costs
	# def _determine_costs_from_interactome_file(self, native_costs_array): return native_costs_array


	def prepare_prizes(self, prize_file):
		"""
		Parses a prize file and returns an array of prizes, a list of terminal indices,
		and terminals missing from the interactome

		Arguments:
			prize_file (str or FILE): a filepath or file object containing a tsv of two columns: node name and prize

		Returns:
			list: prizes, properly indexed (ready for input to pcsf function)
			list: of indices of terminals
		"""

		prizes_dataframe = pd.read_csv(prize_file, sep='\t')
		prizes_dataframe.columns = ['name', 'prize'] + prizes_dataframe.columns[2:].tolist()

		# Strangely some files have duplicated genes, sometimes with different prizes. Keep the max prize.
		logger.info("Duplicated gene symbols in the prize file (we'll keep the max prize):")
		logger.info(prizes_dataframe[prizes_dataframe.set_index('name').index.duplicated()]['name'].tolist())
		prizes_dataframe = prizes_dataframe.groupby('name').max().reset_index()

		# Here's we're indexing the terminal nodes and associated prizes by the indices we used for nodes
		prizes_dataframe.set_index(self.nodes.get_indexer(prizes_dataframe['name']), inplace=True)

		# there will be some nodes in the prize file which we don't have in our interactome
		logger.info("Members of the prize file not present in the interactome:")
		logger.info(prizes_dataframe[prizes_dataframe.index == -1]['name'].tolist())
		prizes_dataframe.drop(-1, inplace=True)

		terminals = sorted(prizes_dataframe.index.values)
		terminal_attributes = prizes_dataframe.set_index('name').rename_axis(None)

		# Here we're making a dataframe with all the nodes as keys and the prizes from above or 0
		prizes_dataframe = pd.DataFrame(self.nodes, columns=["name"]).merge(prizes_dataframe, on="name", how="left").fillna(0)
		# Our return value is a 1D array, where each entry is a node's prize, indexed as above
		prizes = prizes_dataframe['prize'].values * self.params.b

		return prizes, terminals, terminal_attributes


	def _add_dummy_node(self, connected_to=[]):

		dummy_id = len(self.nodes)
		dummy_prize = np.array([0])
		dummy_edges = np.array([(dummy_id, connection) for connection in connected_to])
		dummy_costs = np.array([self.params.w] * len(dummy_edges))

		return dummy_edges, dummy_costs, dummy_id, dummy_prize


	def _check_validity_of_instance(self):
		"""
		Assert that the parammeters and files passed to this program are valid, log useful error messages otherwise.
		"""
		pass


	def pcsf(self, prizes, pruning="strong", verbosity_level=0):
		"""
		"""

		terminals = pd.Series(prizes).nonzero()[0].tolist()
		others = list(set(range(len(self.nodes))) - set(terminals))
		all = list(range(len(self.nodes)))

		if self.params.dummy_mode == 'terminals': endpoints = terminals
		elif self.params.dummy_mode == 'other': endpoints = others
		elif self.params.dummy_mode == 'all': endpoints = all
		else: sys.exit("Invalid dummy mode")

		dummy_edges, dummy_costs, root, dummy_prize = self._add_dummy_node(connected_to=endpoints)

		self._check_validity_of_instance()

		# `edges`: a 2D int64 array. Each row (of length 2) specifies an undirected edge in the input graph. The nodes are labeled 0 to n-1, where n is the number of nodes.
		edges = np.concatenate((self.edges, dummy_edges))
		# `prizes`: the node prizes as a 1D float64 array.
		prizes = np.concatenate((prizes, dummy_prize))
		# `costs`: the edge costs as a 1D float64 array.
		costs = np.concatenate((self.costs, dummy_costs))
		# `root`: the root note for rooted PCST. For the unrooted variant, this parameter should be -1.
		# `num_clusters`: the number of connected components in the output.
		num_clusters = 1
		# `pruning`: a string value indicating the pruning method. Possible values are `'none'`, `'simple'`, `'gw'`, and `'strong'` (all literals are case-insensitive).
		# `verbosity_level`: an integer indicating how much debug output the function should produce.
		vertex_indices, edge_indices = pcst_fast(edges, prizes, costs, root, num_clusters, pruning, verbosity_level)
		# `vertex_indices`: indices of the vertices in the solution as a 1D int64 array.
		# `edge_indices`: indices of the edges in the output as a 1D int64 array. The list contains indices into the list of edges passed into the function.

		return vertex_indices, edge_indices


	def _noisy_edges(self):
		"""
		Adds gaussian noise to all edges in the graph

		Generate gaussian noise values, mean=0, stdev default=0.333 (edge values range between 0 and 1)

		Returns:
			np.ndarray: edge weights with gaussian noise
		"""

		return np.clip(np.random.normal(self.costs, self.params.noise), 0.0001, None)  # None means don't clip above


	def _random_terminals(self, prizes, terminals):
		"""
		Succinct description of _random_terminals

		Selects nodes with a similar degree distribution to the original terminals, and assigns the prizes to them.

		Arguments:
			prizes ():
			terminals ():

		Returns:
			np.ndarray: new prizes
			np.ndarray: new terminals
		"""

		if len(self.edges) < 50: sys.exit("Cannot use random_terminals with such a small interactome.")

		nodes_sorted_by_degree = pd.Series(self.node_degrees).sort_values().index
		terminal_degree_rankings = np.array([nodes_sorted_by_degree.get_loc(terminal) for terminal in terminals])
		new_terminal_degree_rankings = np.clip(np.rint(np.random.normal(terminal_degree_rankings, 10)), 0, len(self.nodes)-1).astype(int)
		new_terminals = pd.Series(nodes_sorted_by_degree)[new_terminal_degree_rankings].values

		new_prizes = copy(prizes)

		for old_terminal, new_terminal in zip(terminals, new_terminals):
			new_prizes[old_terminal] = 0
			new_prizes[new_terminal] = prizes[old_terminal]

		return new_prizes, np.unique(new_terminals)


	def _aggregate_pcsf(self, results, frequency_attribute_name):
		"""
		Succinct description of _aggregate_pcsf

		Longer explanation of_aggregate_pcsf.

		Arguments:
			results (list): a list of [(vertex_indices, edge_indices),...] from multiple PCSF runs.
			frequency_attribute_name (str): Name of the attribute relating to the frequency of occurrence of components in the results.

		Returns:
			pd.DataFrame: new prizes
			pd.DataFrame: new terminals
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


	def randomizations(self, prizes, terminals, terminal_attributes, noisy_edges_reps, random_terminals_reps):
		"""

		Arguments:
			prizes ():
			terminals ():
			noisy_edges_reps (int):
			random_terminals_reps (int):

		Returns:
			networkx.Graph:
		"""

		if self.params.seed: random.seed(seed); numpy.random.seed(seed=seed)

		results = []

		# This is inelegant, but since the pcsf method relies on self.edges, we need to set self.edges
		# with randomized edges before each pcsf run. So we need to copy the true edges to reset later
		edge_costs = copy(self.costs)

		#### NOISY EDGES ####
		for noisy_edge_costs in [self._noisy_edges() for rep in range(noisy_edges_reps)]:
			self.costs = noisy_edge_costs
			results.append(self.pcsf(prizes))
		# Reset the true edges
		self.costs = edge_costs

		if len(results) > 0:
			robust_vertices, robust_edges = self._aggregate_pcsf(results, 'robustness')

		results = []

		#### RANDOM TERMINALS ####
		for random_prizes, terminals in [self._random_terminals(prizes, terminals) for rep in range(random_terminals_reps)]:
			results.append(self.pcsf(random_prizes))

		if len(results) > 0:
			specific_vertices, specific_edges = self._aggregate_pcsf(results, 'specificity')

		###########

		if random_terminals_reps == 0:  # but noisy_edges_reps != 0
			vertex_indices = robust_vertices; edge_indices = robust_edges;

		elif noisy_edges_reps == 0:  # but random_terminals_reps != 0
			vertex_indices = specific_vertices; edge_indices = specific_edges;

		else:  # noisy_edges_reps != 0 and random_terminals_reps != 0
			vertex_indices = robust_vertices.merge(specific_vertices, how='outer', on='node_index')
			edge_indices = robust_edges.merge(specific_edges, how='outer', on='edge_index')

		###########

		# Replace the edge indices with the actual edges (source name, target name) by merging with the interactome
		# By doing an inner join, we get rid of all the dummy node edges.
		edges = edge_indices.merge(self.interactome_dataframe, how='inner', left_on='edge_index', right_index=True)
		vertices = vertex_indices.merge(pd.DataFrame(self.nodes, columns=['name']), how='inner', left_on='node_index', right_index=True).set_index('name')

		forest = nx.from_pandas_dataframe(edges, 'source', 'target', edge_attr=True)
		forest_nodes = forest.nodes()

		if noisy_edges_reps > 0:
			nx.set_node_attributes(forest, 'robustness', {node: robustness for node, robustness in vertices['robustness'].to_dict().items() if node in forest_nodes})
		if random_terminals_reps > 0:
			nx.set_node_attributes(forest, 'specificity', {node: specificity for node, specificity in vertices['specificity'].to_dict().items() if node in forest_nodes})

		for attribute in terminal_attributes.columns.values:
			nx.set_node_attributes(forest, attribute, {node: attr for node, attr in terminal_attributes[attribute].to_dict().items() if node in forest_nodes})

		node_degree_dict = pd.DataFrame(list(zip(self.nodes, self.node_degrees)), columns=['name','degree']).set_index('name').to_dict()['degree']
		nx.set_node_attributes(forest, 'degree',  {node: degree for node, degree in node_degree_dict.items() if node in forest.nodes()})

		augmented_forest = nx.compose(self.interactome_graph.subgraph(vertices.index.tolist()), forest)

		return forest, augmented_forest


	def output_forest_as_networkx(self, vertex_indices, edge_indices, terminal_attributes):
		"""

		Arguments:
			vertex_indices (list): indices of the vertices selected in self.nodes
			edge_indices (list): indices of the edges selected in self.edges

		Returns:
			networkx.Graph: a networkx graph object
		"""

		vertex_indices = pd.DataFrame(vertex_indices, columns=['node_index'])
		edge_indices = pd.DataFrame(edge_indices, columns=['edge_index'])

		# Replace the edge indices with the actual edges (source name, target name) by merging with the interactome
		# By doing an inner join, we get rid of all the dummy node edges.
		edges = edge_indices.merge(self.interactome_dataframe, how='inner', left_on='edge_index', right_index=True)
		nodes = vertex_indices.merge(pd.DataFrame(self.nodes, columns=['name']), how='inner', left_on='node_index', right_index=True).set_index('name')

		forest = nx.from_pandas_dataframe(edges, 'source', 'target', edge_attr=['Weight'])

		for attribute in terminal_attributes.columns.values:
			nx.set_node_attributes(forest, attribute, {node: attr for node, attr in terminal_attributes[attribute].to_dict().items() if node in forest.nodes()})
			nx.set_node_attributes(forest, attribute, {node:0 for node in forest.nodes() if attribute not in forest.node[node]})

		#node_degree_dict = pd.DataFrame(list(zip(self.nodes, self.node_degrees)), columns=['name','degree']).set_index('name').to_dict()['degree']
		#nx.set_node_attributes(forest, 'degree', {node: degree for node, degree in node_degree_dict.items() if node in forest.nodes()})

		augmented_forest = nx.compose(self.interactome_graph.subgraph(nodes.index.tolist()), forest)

		return forest, augmented_forest


	def betweenness(self, nxgraph):
		"""
		Calculate betweenness centrality for all nodes in forest. The forest *should* be augmented

		Arguments:
			nxgraph (networkx.Graph): a networkx graph object

		Returns:
			networkx.Graph: a networkx graph object
		"""

		betweenness = nx.betweenness_centrality(nxgraph)
		nx.set_node_attributes(nxgraph, 'betweenness', betweenness)

		return nxgraph


	def pcsf_objective_value(self, prizes, forest):
		"""
		Calculate PCSF objective function

		Arguments:
			prizes (list): a list of prizes like the one returned by the prepare_prizes method.
			forest (networkx.Graph): a forest like the one returned by output_forest_as_networkx -- Not an augmented forest!

		Returns:
			float: PCSF objective function score
		"""

		return (sum(prizes) - sum(nx.get_node_attributes(forest, 'prize').values())) + sum(nx.get_edge_attributes(forest, 'Weight').values()) + (self.params.w * nx.number_connected_components(forest))


def output_networkx_graph_as_gml_for_cytoscape(nxgraph, output_dir, filename):
	"""
	Arguments:
		nxgraph (networkx.Graph): any instance of networkx.Graph
		output_dir (str): the directory in which to output the graph. Must already exist.
		filename (str): Filenames ending in .gz or .bz2 will be compressed.
	"""
	path = os.path.join(os.path.abspath(output_dir), filename)
	nx.write_gml(nxgraph, path)


def output_networkx_graph_as_json_for_cytoscapejs(nxgraph, output_dir):
	"""
	Arguments:
		nxgraph (networkx.Graph): any instance of networkx.Graph
		output_dir (str): the directory in which to output the file (named graph_json.json)
	"""
	path = os.path.join(os.path.abspath(output_dir), 'graph_json.json')
	njs = cy.from_networkx(nxgraph)
	with open(path,'w') as outf:
		outf.write(json.dumps(njs, indent=4))
	

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


def merge_two_prize_files(prize_file_1, prize_file_2, prize_file_1_node_type=None, prize_file_2_node_type=None):
	"""

	Arguments:
		prize_file_1 (str or FILE): a filepath or FILE object with a tsv of name(\t)prize(\t)more...
		prize_file_2 (str or FILE): a filepath or FILE object with a tsv of name(\t)prize(\t)more...
		prize_file_1_node_type
		prize_file_2_node_type

	Returns:
		pandas.DataFrame: a DataFrame of prizes with duplicates removed (first entry kept)
	"""

	prize_df1 = pd.read_csv(prize_file_1, sep='\t')
	prize_df1.columns = ['name', 'prize'] + prize_df1.columns[2:].tolist()
	if prize_file_1_node_type: prize_df1['type'] = prize_file_1_node_type
	prize_df2 = pd.read_csv(prize_file_2, sep='\t')
	prize_df2.columns = ['name', 'prize'] + prize_df2.columns[2:].tolist()
	if prize_file_2_node_type: prize_df2['type'] = prize_file_2_node_type

	return merge_two_prize_dataframes(prize_df1, prize_df2)


def merge_two_prize_dataframes(prize_df1, prize_df2):
	"""

	Arguments:
		prize_df1 (pandas.DataFrame): a dataframe with at least columns 'name' and 'prize'
		prize_df2 (pandas.DataFrame): a dataframe with at least columns 'name' and 'prize'

	Returns:
		pandas.DataFrame: a DataFrame of prizes with duplicates removed (first entry kept)
	"""

	prizes_dataframe = pd.concat((prize_df1, prize_df2))
	prizes_dataframe.drop_duplicates(subset=['name'], inplace=True) # Unclear if we should do this?

	return prizes_dataframe


def output_dataframe_to_tsv(dataframe, output_dir, filename):
	"""
	Output the dataframe to a csv
	"""
	path = os.path.join(os.path.abspath(output_dir), filename)
	dataframe.to_csv(path, sep='\t', header=True, index=False)


