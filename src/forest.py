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

# Core python external libraries
import numpy as np
import pandas as pd

# Peripheral python external libraries
import networkx as nx

# Lab modules
from pcst_fast.pcst_fast import pcst_fast

# list of classes and methods we'd like to export:
__all__ = ["Graph", "output_networkx_graph_as_gml_for_cytoscape", "merge_two_prize_files"]


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


parser = argparse.ArgumentParser(description="""
	Find multiple pathways within an interactome that are altered in a particular condition using
	the Prize Collecting Steiner Forest problem.""")

class FullPaths(argparse.Action):
	"""Expand user- and relative-paths"""
	def __call__(self,parser, namespace, values, option_string=None):
		setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def directory(dirname):
	if not os.path.isdir(dirname): raise argparse.ArgumentTypeError(dirname + " is not a directory")
	else: return dirname

# Input / Output parameters:
parser.add_argument("-e", "--edge", dest='edge_file', type=argparse.FileType('r'), required=True,
	help ='(Required) Path to the text file containing the edges. Should be a tab delimited file with 3 columns: "nodeA\tnodeB\tweight(between 0 and 1)"')
parser.add_argument("-p", "--prize", dest='prize_file', type=argparse.FileType('r'), required=True,
	help='(Required) Path to the text file containing the prizes. Should be a tab delimited file with lines: "nodeName\tprize"')
parser.add_argument('-o', '--output', dest='output_dir', action=FullPaths, type=directory, required=True,
	help='(Required) Output directory path')

# Command parameters (specify what the algorithm does):
parser.add_argument("--noisy_edges", dest='noisy_edges_repetitions', type=int, default=0,
	help='An integer specifying how many times you would like to add noise to the given edge values and re-run the algorithm. Results of these runs will be merged together and written in files with the word "_noisy_edges_" added to their names. The noise level can be controlled using the configuration file. [default: %default]')
parser.add_argument("--random_terminals", dest='random_terminals_repetitions', type=int, default=0,
	help='An integer specifying how many times you would like to apply your given prizes to random nodes in the interactome (with a similar degree distribution) and re-run the algorithm. Results of these runs will be merged together and written in files with the word "_random_terminals_" added to their names. [default: %default]')
# parser.add_argument("--knockout", dest='knockout', nargs='*', default=[],
# 	help='Protein(s) you would like to "knock out" of the interactome to simulate a knockout experiment. [default: %default]')

params = parser.add_argument_group('Parameters', 'Parameters description')

params.add_argument("-w", dest="w", type=float, required=False,
	help="   [default: 6]")
params.add_argument("-b", dest="b", type=float, required=False,
	help="   [default: 12]")
params.add_argument("-mu", dest="mu", type=float, required=False,
	help="   [default: 0.04]")
params.add_argument("-r", dest="r", type=float, required=False,
	help="   [default: None]")
params.add_argument("-noise", dest="noise", type=float, required=False,
	help="   [default: 0.33]")
params.add_argument("--dummyMode", dest='dummy_mode', choices=("terminals", "other", "all"), required=False,
	help='Tells the program which nodes in the interactome to connect the dummy node to. "terminals"= connect to all terminals, "others"= connect to all nodes except for terminals, "all"= connect to all nodes in the interactome. [default: terminals]')

params.add_argument("--muSquared", action='store_true', dest='mu_squared', required=False,
	help='Flag to add negative prizes to hub nodes proportional to their degree^2, rather than degree. Must specify a positive mu in conf file. [default: False]')

params.add_argument("--excludeTerminals", action='store_true', dest='exclude_terminals', required=False,
	help='Flag to exclude terminals when calculating negative prizes. Use if you want terminals to keep exact assigned prize regardless of degree. [default: False]')

params.add_argument("-s", "--seed", dest='seed', type=int, required=False,
	help='An integer seed for the pseudo-random number generators. If you want to reproduce exact results, supply the same seed. [default: None]')



if __name__ == '__main__':

	args = parser.parse_args()

	params = vars(args.params) # http://stackoverflow.com/questions/42400646/is-it-possible-to-denote-some-set-of-argparses-arguments-without-using-subparse

	graph = Graph(args.edge_file, params)

	prizes, terminals = graph.prepare_prizes(args.prize_file)

	if noisy_edges_repetitions + random_terminals_repetitions > 0:
		nxgraph = graph.randomizations(noisy_edges_repetitions, random_terminals_repetitions)

	else:
		vertices, edges = graph.pcsf(prizes)
		nxgraph = graph.output_forest_as_networkx(vertices, edges)

	output_networkx_graph_as_gml_for_cytoscape(nxgraph, output_dir+'/output.gml')



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

		interactome_fieldnames = ["source","target","cost"]
		self.interactome_dataframe = pd.read_csv(interactome_file, delimiter='\t', names=interactome_fieldnames)
		# self.interactome_graph = nx.from_pandas_dataframe(self.interactome_dataframe, 'source', 'target', edge_attr=['cost'])

		# We first take only the source and target columns from the interactome dataframe.
		# We then unstack them, which, unintuitively, stacks them into one column, which is a hack
		# 	allowing us to use factorize, which is a very convenient function.
		# Factorize builds two datastructures, a unique pd.Index of each ID string to a numerical ID
		# and the datastructure we passed in with ID strings replaced with those numerical IDs.
		# We place those in self.nodes and self.edges respectively, but self.edges will need reshaping.
		(self.edges, self.nodes) = pd.factorize(self.interactome_dataframe[["source","target"]].unstack())


		# Here we do the inverse operation of "unstack" above, which gives us an interpretable edges datastructure
		self.edges = self.edges.reshape(self.interactome_dataframe[["source","target"]].shape, order='F')

		self.costs = self.interactome_dataframe['cost'].values

		# Numpy has a convenient counting function. However we're assuming here that each edge only appears once.
		# The indices into this datastructure are the same as those in self.nodes and self.edges.
		self.node_degrees = np.bincount(self.edges.flatten())


		defaults = {"w": 6, "b": 12, "D": 6, "mu": 0.04, "r": None, "noise": 0.1, "mu_squared": False, "exclude_terminals": False, "dummy_mode": "terminals", "seed": None}

		self.params = Options({**defaults, **params})

		self.negprizes = (self.node_degrees**2 if self.params.mu_squared else self.node_degrees) * self.params.mu


	def prepare_prizes(self, prize_file):
		"""
		Parses a prize file and returns an array of prizes, a list of terminal indices,
		and terminals missing from the interactome

		Arguments:
			prize_file (str or FILE): a filepath or file object containing a tsv of two columns: node name and prize

		Returns:
			list: prizes, properly indexed (ready for input to pcsf function)
			list: of indices of terminals
			dataframe: terminals missing from interactome
		"""

		prizes_dataframe = pd.read_csv(prize_file, delimiter='\t', names=["name","prize"])

		# Here's we're indexing the terminal nodes and associated prizes by the indices we used for nodes
		prizes_dataframe.set_index(self.nodes.get_indexer(prizes_dataframe['name']), inplace=True)

		# there will be some nodes in the prize file which we don't have in our interactome
		terminals_missing_from_interactome = prizes_dataframe[prizes_dataframe.index == -1]
		logger.info("Members of the prize file not present in the interactome:")
		logger.info(terminals_missing_from_interactome)
		prizes_dataframe.drop(-1, inplace=True)

		terminals = sorted(prizes_dataframe.index.values)

		# Here we're making a dataframe with all the nodes as keys and the prizes from above or 0
		prizes_dataframe = prizes_dataframe.merge(pd.DataFrame(self.nodes, columns=["name"]), on="name", how="outer").fillna(0)
		# We re-index again, making sure we have a consistent index scheme
		prizes_dataframe.set_index(self.nodes.get_indexer(prizes_dataframe['name']), inplace=True)
		prizes_dataframe.sort_index(inplace=True)

		# Our return value is a list, where each entry is a node's prize, indexed as above
		prizes = prizes_dataframe['prize'].values

		return prizes, terminals


	def _add_dummy_node(self, connected_to=[]):

		dummy_id = len(self.nodes)
		dummy_edges = np.array([(dummy_id, connection) for connection in connected_to])
		dummy_costs = np.array([self.params.w] * len(dummy_edges))

		return dummy_edges, dummy_costs, dummy_id


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

		dummy_edges, dummy_costs, root = self._add_dummy_node(connected_to=endpoints)

		self._check_validity_of_instance()

		# `edges`: a 2D int64 array. Each row (of length 2) specifies an undirected edge in the input graph. The nodes are labeled 0 to n-1, where n is the number of nodes.
		edges = np.concatenate((self.edges, dummy_edges))
		# `prizes`: the node prizes as a 1D float64 array.
		# `costs`: the edge costs as a 1D float64 array.
		costs = np.concatenate((self.costs, dummy_costs))
		# `root`: the root note for rooted PCST. For the unrooted variant, this parameter should be -1.
		# `num_clusters`: the number of connected components in the output.
		num_clusters = 0
		# `pruning`: a string value indicating the pruning method. Possible values are `'none'`, `'simple'`, `'gw'`, and `'strong'` (all literals are case-insensitive).
		# `verbosity_level`: an integer indicating how much debug output the function should produce.
		vertex_indices, edge_indices = pcst_fast(edges, prizes, costs, root, num_clusters, pruning, verbosity_level)
		# `vertices`: the vertices in the solution as a 1D int64 array.
		# `edges`: the edges in the output as a 1D int64 array. The list contains indices into the list of edges passed into the function.

		return vertex_indices, edge_indices


	def _noisy_edges(self):
		"""
		Adds gaussian noise to all edges in the graph

		Generate gaussian noise values, mean=0, stdev default=0.333 (edge values range between 0 and 1)

		Returns:
			np.ndarray: edge weights with gaussian noise
		"""

		return np.clip(np.random.normal(self.costs, self.params.noise), 0, 1)


	def _random_terminals(self, prizes, terminals):
		"""
		Succinct description of random_terminals

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
		new_terminal_degree_rankings = np.clip(terminal_degree_rankings + numpy.random.normal(0,100.0,len(terminals)), 0, len(self.nodes)).tolist()
		new_terminals = pd.Series(nodes_sorted_by_degree)[new_terminal_degree_rankings][0]

		new_prizes = copy(prizes)

		for old_terminal, new_terminal in zip(terminals, new_terminals):
			new_prizes[old_terminal] = 0
			new_prizes[new_terminal] = prizes[old_terminal]

		return new_prizes, np.unique(new_terminals)


	def randomizations(self, prizes, noisy_edges_reps, random_terminals_reps):
		"""

		Arguments:
			prizes ():
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

		#### RANDOM TERMINALS ####
		for random_prizes, terminals in [self._random_terminals(prizes) for rep in range(random_terminals_reps)]:

			results.append(self.pcsf(random_prizes))

		denominator = len(results)
		# Transposes a list from [(vertex_indices, edge_indices),...] to ([vertex_indices,...], [edge_indices,...])
		vertex_indices, edge_indices = zip(*results)

		# These next steps are just data transformation/aggregation.
		# 1. Flatten the lists of lists of edge indices and vertex indices
		# 2. Count the occurrences of each edge and vertex index
		# 3. Transform from Counter object to DataFrame through list
		vertex_indices = pd.DataFrame(list(Counter(flatten(vertex_indices)).items()), columns=['node_index','occurrence'])
		edge_indices = pd.DataFrame(list(Counter(flatten(edge_indices)).items()), columns=['edge_index','occurrence'])
		# 4. Convert occurrences to fractions
		vertex_indices['occurrence'] /= denominator
		edge_indices['occurrence'] /= denominator

		# Replace the edge indices with the actual edges (source name, target name) by merging with the interactome
		# By doing an inner join, we get rid of all the dummy node edges.
		edges = edge_indices.merge(self.interactome_dataframe, how='inner', left_on='edge_index', right_index=True)
		vertices = vertex_indices.merge(pd.DataFrame(self.nodes, columns=['name']), how='inner', left_on='node_index', right_index=True).set_index('name')

		forest = nx.from_pandas_dataframe(edges, 'source', 'target', edge_attr=['cost','occurrence'])
		nx.set_node_attributes(forest, 'occurrence', vertices['occurrence'].to_dict())

		augmented_forest = nx.compose(forest, self.interactome_graph.subgraph(vertices.index.tolist()))

		return forest, augmented_forest


	def output_forest_as_networkx(self, vertex_indices, edge_indices):
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

		forest = nx.from_pandas_dataframe(edges, 'source', 'target', edge_attr=['cost'])

		augmented_forest = nx.compose(forest, self.interactome_graph.subgraph(nodes.index.tolist()))

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


def output_networkx_graph_as_gml_for_cytoscape(nxgraph, filename):
	"""
	Arguments:
		nxgraph (networkx.Graph): any instance of networkx.Graph
		filename (str): A full filepath. Filenames ending in .gz or .bz2 will be compressed.
	"""
	nx.write_gml(nxgraph, filename)


def merge_two_prize_files(prize_file_1, prize_file_2):
	"""

	Arguments:
		prize_file_1 ():
		prize_file_2 ():
	"""
	pass




