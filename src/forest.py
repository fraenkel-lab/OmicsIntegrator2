#!/usr/bin/env python3

# Core python modules
import sys
import os

# Peripheral python modules
import argparse
import pickle
import logging
import random

# Core python external libraries
import numpy as np
import pandas as pd

# Peripheral python external libraries
import networkx as nx

# Lab modules
import garnet
import pcst_fast.pcst_fast as pcst_fast

# list of public methods:
__all__ = ["Graph"]


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
handler.setFormatter(logging.Formatter('%(asctime)s - Forest: %(levelname)s - %(message)s', "%I:%M:%S"))
logger.addHandler(handler)


parser = argparse.ArgumenParser(description="""
	Find multiple pathways within an interactome that are altered in a particular condition using
	the Prize Collecting Steiner Forest problem.""")

class FullPaths(argparse.Action):
	"""Expand user- and relative-paths"""
	def __call__(self,parser, namespace, values, option_string=None):
		setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def directory(dirname):
	if not os.path.isdir(dirname): raise argparse.ArgumentTypeError(dirname + " is not a directory")
	else: return dirname

# File arguments
parser.add_argument("-i", "--edge", dest='edge_file', type=argparse.FileType('r'), required=True,
	help ='(Required) Path to the text file containing the edges. Should be a tab delimited file with 3 columns: "ProteinA\tProteinB\tWeight(between 0 and 1)"')
parser.add_argument("-p", "--prize", dest='prize_file', type=argparse.FileType('r'), required=True,
	help='(Required) Path to the text file containing the prizes. Should be a tab delimited file with lines: "ProteinName\tPrizeValue"')
parser.add_argument("-tf", "--tfprizes" dest='garnet_file', type=argparse.FileType('r'), required=False,
	help='tsv file containing the output of the GARNET module regression.')
parser.add_argument('-o', '--output', dest='output_dir', action=FullPaths, type=directory, required=True,
	help='(Required) Output directory path')

# Optional arguments

# confFile (): text file containing values for all parameters. Should include the lines "w=<value>", "D=<value>", and "b=<value>".
parser.add_argument("-c", "--conf", dest='config_file', type=argparse.FileType('r'), default='conf.txt',
	help='Path to the text file containing the parameters. Should be several lines that looks like: "ParameterName = ParameterValue". Must contain values for w, b, D.  May contain values for optional parameters mu, garnetBeta, noise, r, g. [default: %default]')

# dummyMode (): a string that indicates which nodes in the interactome to connect the dummy node to. 'terminals'=connect to all terminals (default), 'others'=connect to all nodes except for terminals, 'all'=connect to all nodes in the interactome, or a the path to a text file containing a list of proteins to connect to.
parser.add_argument("-d","--dummyMode", dest='dummy_mode', choices=("terminals", "other", "all"), default='terminals',
	help='Tells the program which nodes in the interactome to connect the dummy node to. "terminals"= connect to all terminals, "others"= connect to all nodes except for terminals, "all"= connect to all nodes in the interactome. If you wish you supply your own list of proteins, dummyMode could also be the path to a text file containing a list of proteins (one per line). [default: %default]')

parser.add_argument("--muSquared", action='store_true', dest='mu_squared', default=False,
	help='Flag to add negative prizes to hub nodes proportional to their degree^2, rather than degree. Must specify a positive mu in conf file. [default: %default]')
parser.add_argument("--excludeTerminals", action='store_true', dest='exclude_terminals', default=False,
	help='Flag to exclude terminals when calculating negative prizes. Use if you want terminals to keep exact assigned prize regardless of degree. [default: %default]')

parser.add_argument("--noisy_edges", dest='noisy_edges_repetitions', default=0, type=int,
	help='An integer specifying how many times you would like to add noise to the given edge values and re-run the algorithm. Results of these runs will be merged together and written in files with the word "_noisy_edges_" added to their names. The noise level can be controlled using the configuration file. [default: %default]')
parser.add_argument("--random_terminals", dest='random_terminals_repetitions', default=0, type=int,
	help='An integer specifying how many times you would like to apply your given prizes to random nodes in the interactome (with a similar degree distribution) and re-run the algorithm. Results of these runs will be merged together and written in files with the word "_random_terminals_" added to their names. [default: %default]')
parser.add_argument("--knockout", dest='knockout', nargs='*', default=[],
	help='Protein(s) you would like to "knock out" of the interactome to simulate a knockout experiment. [default: %default]')

parser.add_argument("-s", "--seed", dest='seed', type=int, default=None,
	help='An integer seed for the pseudo-random number generators. If you want to reproduce exact results, supply the same seed. [default: %default]')


if __name__ == '__main__':

	args = parser.parse_args()

	options = Options({"w": 6 , "b": 12, :"gb": 0.1, "D": 6, "mu": 0.04})

	graph = Graph(edge_file, options)

	prizes, terminals, terminals_missing_from_interactome = graph.prepare_prizes(prize_file)

	vertices, edges = graph.pcsf(prizes)

	nxgraph = graph.output_forest_as_networkx(vertices, edges)

	output_networkx_graph_as_gml_for_cytoscape(nxgraph)


class Options:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


class Graph:
	"""

	"""
	def __init__(self, interactome_file, options):
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
			options (dict): options with which to run the program

		"""

		interactome_fieldnames = ["source","target","weight"]
		interactome_dataframe = pd.read_csv(interactome_file, delimiter='\t', names=interactome_fieldnames)

		# We first take only the source and target columns from the interactome dataframe.
		# We then unstack them, which, unintuitively, stacks them into one column, which is a hack
		# 	allowing us to use factorize, which is a very convenient function.
		# Factorize builds two datastructures, a unique pd.Index of each ID string to a numerical ID
		# and the datastructure we passed in with ID strings replaced with those numerical IDs.
		# We place those in self.nodes and self.edges respectively, but self.edges will need reshaping.
		(self.edges, self.nodes) = pd.factorize(interactome_dataframe[["source","target"]].unstack())

		# Here we do the inverse operation of "unstack" above, which gives us an interpretable edges datastructure
		self.edges = self.edges.reshape(interactome_dataframe[["source","target"]].shape, order='F') #.tolist() # maybe needs to be list of tuples. in that case map(tuple, this)

		self.costs = interactome_dataframe['weight'].tolist()

		# Numpy has a convenient counting function. However we're assuming here that each edge only appears once.
		# The indices into this datastructure are the same as those in self.nodes and self.edges.
		self.node_degrees = np.bincount(self.edges.flatten())

		self.edges = self.edges.tolist()

		# self._add_dummy_node(connected_to=options.dummy_mode)

		# self._check_validity_of_instance()

		defaults = {"w": 6 , "b": 12, :"gb": 0.1, "D": 6, "mu": 0.04, "r": None, "noise": 0.33, "mu_squared": False, "exclude_terminals": False}

		self.params = Options(defaults.update(options))

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
		terminals_missing_from_interactome = prizes_dataframe[prizes_dataframe.index == -1]  # if this only returns a view, the next line will cause problems
		prizes_dataframe.drop(-1, inplace=True)

		terminals = sorted(prizes_dataframe.index.tolist())

		# Here we're making a dataframe with all the nodes as keys and the prizes from above or 0
		prizes_dataframe = prizes_dataframe.merge(pd.DataFrame(self.nodes, columns=["name"]), on="name", how="outer").fillna(0)
		# We re-index again, making sure we have a consistent index scheme
		prizes_dataframe.set_index(self.nodes.get_indexer(prizes_dataframe['name']), inplace=True)
		prizes_dataframe.sort_index(inplace=True)

		# Our return value is a list, where each entry is a node's prize, indexed as above
		prizes = prizes_dataframe['prize'].tolist()

		return prizes, terminals, terminals_missing_from_interactome


	def _add_dummy_node(self, connected_to="terminals"):

		self.dummy_id = len(self.nodes)

		# if connected_to == 'terminals' or connected_to == 'all':

		# if connected_to == 'others' or connected_to == 'all':


	def _check_validity_of_instance(self):
		"""
		Assert that the parammeters and files passed to this program are valid, log useful error messages otherwise.
		"""
		pass

		# does there exist an assert module in python? Is that what we would want here? does it play nice with logger?


	def pcsf(prizes, pruning="strong", verbosity_level=0):
		"""
		"""

		# `edges`: a list of pairs of integers. Each pair specifies an undirected edge in the input graph. The nodes are labeled 0 to n-1, where n is the number of nodes.
		# `prizes`: the list of node prizes.
		# `costs`: the list of edge costs.
		# `root`: the root note for rooted PCST. For the unrooted variant, this parameter should be -1.
		# `num_clusters`: the number of connected components in the output.
		num_clusters = 1
		# `pruning`: a string value indicating the pruning method. Possible values are `'none'`, `'simple'`, `'gw'`, and `'strong'` (all literals are case-insensitive).
		# `verbosity_level`: an integer indicating how much debug output the function should produce.
		vertices, edges = pcst_fast(self.edges, prizes, self.costs, self.root, num_clusters, pruning, verbosity_level)

		# `vertices`: a list of vertices in the output.
		# `edges`: a list of edges in the output. The list contains indices into the list of edges passed into the function.
		return vertices, edges


	def noisy_edges(self, seed=None):
		"""
		Adds gaussian noise to all edges in the graph

		Generate gaussian noise values, mean=0, stdev default=0.333 (edge values range between 0 and 1)

		Arguments:
			seed (): a random seed

		Returns:
			list: edge weights with gaussian noise
		"""

		if seed: random.seed(seed)

		std_dev = self.params.noise
		costs = [cost + random.gauss(0,std_dev) for cost in self.costs]

		return costs


	def random_terminals(self, terminals, prizes, seed=None):
		"""
		Succinct description of random_terminals

		Selects nodes with a similar degree distribution to the original terminals, and assigns the prizes to them.

		Arguments:
			seed (): a random seed

		Returns:
			list: new terminal nodes
		"""

		if len(self.edges) < 50: sys.exit("Cannot use --random_terminals with such a small interactome.")

		if seed: random.seed(seed)

		nodes_sorted_by_degree = pd.Series(self.node_degrees).sort_values()

		new_prizes = prizes.copy()
		new_terminals = []

		for terminal in terminals:
			prize = prizes.loc(terminal)
			new_prizes[terminal] = 0
			new_terminal = clip(terminal_index_in nodes_sorted_by_degree + random.gauss(0,100.0), lower=0, upper=len(self.nodes))
			new_prizes[new_terminal] = prize
			new_terminals.append(new_terminal)

		return new_prizes, new_terminals


	def randomizations(self, noisy_edges_repetitions, random_terminals_repetitions):
		"""
		"""

		results = []

		edge_costs = self.costs.copy()

		for noisy_edge_costs in [self.noisy_edges() for rep in range(noisy_edges_repetitions)]:
			self.costs = noisy_edge_costs
			results.append(self.pcsf(prizes))

		self.costs = edge_costs

		for random_prizes, terminals in [self.random_terminals() for rep in range(random_terminals_repetitions)]:

			results.append(self.pcsf(random_prizes))

		for vertices, edges in results:
			True
			# maybe a counter  to count occurrences?


	def output_forest_as_networkx(self, vertices, edges):

		# OUTPUT: self.optForest - a networkx digraph storing the forest returned by msgsteiner
		#         self.augForest - a networkx digraph storing the forest returned by msgsteiner, plus
		#                          all of the edges in the interactome between nodes in that forest
		#         self.dumForest - a networkx digraph storing the dummy node edges in the optimal
		#                          forest
		pass
			   # betweenness - a boolean flag indicating whether we should do the costly betweenness
			   #               calculation


def mergeOutputs(PCSFOutputObj1, PCSFOutputObj2, betweenness, n1=1, n2=1):
    """
    Merges two PCSFOutput objects together. Creates a new PCSFOutput object whose graphs contain
    all edges found in either original object, with updated fracOptContaining values and
    betweenness values.

    INPUT: Two PCSFOutput objects, either individual objects or themselves results of merges.
                Ideally, these output objects were created using the same interactome (though the
                prizes or algorithm parameters may have been different). This is not enforced.
           betweenness - a T/F flag indicating whether to do the costly betweenness calculation
           n1,n2 - integers, the number of msgsteiner runs each PCSFOutput object is a result of
                   (if one of PCSFOutputObj is the result of a merge, this should be >1).

    RETURNS: A new PCSFOutput object, with all edges found in either original object, with updated
                fracOptContaining and betweenness values. If a node or edge is found in both
                original objects, the prize or weight in this object is copied from
                PCSFOutputObj1. The inputObj reference in this object is the same as
                PCSFOutputObj1.
    """
    print 'Merging outputs to give summary over %i algorithm runs...'%(n1+n2)
    mergedObj = copy.deepcopy(PCSFOutputObj1)
    #Update fracOptContaining for all edges in outputObj1
    for (node1,node2,data) in PCSFOutputObj1.optForest.edges(data=True):
        numRuns1 = data['fracOptContaining']*n1
        try:
            #if the edge is not in outputObj2 this will return a KeyError
            numRuns2 = PCSFOutputObj2.optForest[node1][node2]['fracOptContaining'] * n2
        except KeyError:
            numRuns2 = 0.0
        mergedObj.optForest[node1][node2]['fracOptContaining'] = (numRuns1 + numRuns2)/(n1+n2)
    #Update fracOptContaining for all nodes in outputObj1
    for (node, data) in PCSFOutputObj1.optForest.nodes(data=True):
        numRuns1 = data['fracOptContaining']*n1
        try:
            #if the node is not in outputObj2 this will return a KeyError
            numRuns2 = PCSFOutputObj2.optForest.node[node]['fracOptContaining'] * n2
        except KeyError:
            numRuns2 = 0.0
        mergedObj.optForest.node[node]['fracOptContaining'] = (numRuns1 + numRuns2)/(n1+n2)
    #Add optForest edges to mergedObj that appear in outputObj2 but not in outputObj1
    for (node1, node2, data) in PCSFOutputObj2.optForest.edges(data=True):
        try:
            dataM = mergedObj.optForest[node1][node2]
        except KeyError:
            numRuns2 = data['fracOptContaining'] * n2
            #If there are nodes in outputObj2 not included in 1, they will be added to
            #mergedObj without error
            if node1 not in mergedObj.optForest.nodes():
                mergedObj.optForest.add_node(node1,
                                             prize=PCSFOutputObj2.optForest.node[node1]['prize'],
                                             TerminalType=PCSFOutputObj2.optForest.node[node1]['TerminalType'],
                                             fracOptContaining=numRuns2/(n1+n2))
            if node2 not in mergedObj.optForest.nodes():
                mergedObj.optForest.add_node(node2,
                                            prize=PCSFOutputObj2.optForest.node[node2]['prize'],
                                            TerminalType=PCSFOutputObj2.optForest.node[node2]['TerminalType'],
                                            fracOptContaining=numRuns2/(n1+n2))
            mergedObj.optForest.add_edge(node1, node2, weight=data['weight'],
                                         fracOptContaining=numRuns2/(n1+n2))
    #Add dumForest edges to mergedObj that appear in outputObj2 but not in outputObj1
    for (node1, node2, data) in PCSFOutputObj2.dumForest.edges(data=True):
        try:
            dataM = mergedObj.dumForest[node1][node2]
        except KeyError:
            mergedObj.dumForest.add_edge(node1, node2)

    #Create augForest based on new optForest
    #Need to first copy optForest in case an edge previously included in augForest
    #has a new fracOptContaining
    mergedObj.augForest = copy.deepcopy(mergedObj.optForest)
    for node in mergedObj.augForest.nodes():
        edges = {}
        try:
            edges.update(mergedObj.inputObj.undirEdges[node])
        except KeyError:
            pass
        try:
           edges.update(mergedObj.inputObj.dirEdges[node])
        except KeyError:
            #If a node found in mergedObj.optForest is not found in PCSFInputObj1's interactome,
            #it is quietly ignored in making augForest
            pass
        for node2 in edges:
            if node2 in mergedObj.augForest.nodes():
                if (node, node2) not in mergedObj.optForest.edges():
                    mergedObj.augForest.add_edge(node, node2, weight=edges[node2],
                                                 fracOptContaining=0.0)

    #Calculate betweenness centrality for all nodes in augmented forest
    if betweenness:
        betweenness = nx.betweenness_centrality(mergedObj.augForest)
        nx.set_node_attributes(mergedObj.augForest, 'betweenness', betweenness)

    print 'Outputs were successfully merged.\n'
    return mergedObj








def output_networkx_graph_as_gml_for_cytoscape(nxgraph):
	pass

