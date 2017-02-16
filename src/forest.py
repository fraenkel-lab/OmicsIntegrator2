#!/usr/bin/env python3

# Core python modules
import sys
import os

# Peripheral python modules
import argparse
import pickle
import logging

# Core python external libraries
import numpy as np
import pandas as pd

# Peripheral python external libraries

# Lab modules
import garnet
import pcst_fast

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

# Required arguments
parser.add_argument("-p", "--prize", dest='prize_file', type=argparse.FileType('r'), required=True,
	help='(Required) Path to the text file containing the prizes. Should be a tab delimited file with lines: "ProteinName\tPrizeValue"')
parser.add_argument("-e", "--edge", dest='edge_file', type=argparse.FileType('r'), required=True,
	help ='(Required) Path to the text file containing the interactome edges. Should be a tab delimited file with 3 or 4 columns: "ProteinA\tProteinB\tWeight(between 0 and 1)\tDirectionality(U or D, optional)"')
parser.add_argument('-o', '--output', dest='output_dir', action=FullPaths, type=directory, required=True,
	help='(Required) Output directory path')

# Optional arguments
parser.add_argument("-c", "--conf", dest='config_file', type=argparse.FileType('r'), default='conf.txt',
	help='Path to the text file containing the parameters. Should be several lines that looks like: "ParameterName = ParameterValue". Must contain values for w, b, D.  May contain values for optional parameters mu, garnetBeta, noise, r, g. [default: %default]')
parser.add_argument("-d","--dummyMode", dest='dummy_mode', choices=("terminals", "other", "all") default='terminals',
	help='Tells the program which nodes in the interactome to connect the dummy node to. "terminals"= connect to all terminals, "others"= connect to all nodes except for terminals, "all"= connect to all nodes in the interactome. If you wish you supply your own list of proteins, dummyMode could also be the path to a text file containing a list of proteins (one per line). [default: %default]')

# do we really want this?
parser.add_argument("--garnet", dest='garnet_directory_path_or_None', action=FullPaths, type=directory, default=None,
	help='Path to the text file containing the output of the GARNET module regression. Should be a tab delimited file with 2 columns: "TranscriptionFactorName\tScore". [default: %default]')

parser.add_argument("--muSquared", action='store_true', dest='mu_squared', default=False,
	help='Flag to add negative prizes to hub nodes proportional to their degree^2, rather than degree. Must specify a positive mu in conf file. [default: %default]')
parser.add_argument("--excludeTerms", action='store_true', dest='exclude_terms', default=False,
	help='Flag to exclude terminals when calculating negative prizes. Use if you want terminals to keep exact assigned prize regardless of degree. [default: %default]')
parser.add_argument("--outlabel", dest = 'output_label', default='result',
	help='A string to put at the beginning of the names of files output by the program. [default: %default]')
parser.add_argument("--noisyEdges", dest='noiseNum', default=0, type=int,
	help='An integer specifying how many times you would like to add noise to the given edge values and re-run the algorithm. Results of these runs will be merged together and written in files with the word "_noisyEdges_" added to their names. The noise level can be controlled using the configuration file. [default: %default]')
parser.add_argument("--randomTerminals", dest='termNum', default=0, type=int,
	help='An integer specifying how many times you would like to apply your given prizes to random nodes in the interactome (with a similar degree distribution) and re-run the algorithm. Results of these runs will be merged together and written in files with the word "_randomTerminals_" added to their names. [default: %default]')
parser.add_argument("--knockout", dest='knockout', nargs='*', default=[],
	help='Protein(s) you would like to "knock out" of the interactome to simulate a knockout experiment. [default: %default]')
parser.add_argument("-s", "--seed", dest='seed', type=int, default=None,
	help='An integer seed for the pseudo-random number generators. If you want to reproduce exact results, supply the same seed. [default: %default]')


if __name__ == '__main__':

	args = parser.parse_args()
	options = {}

	graph = Graph(edge_file)

	main(graph, prize_file, output_dir, options)


class Graph():
	"""docstring for Graph"""
	def __init__(self, interactome_file):

		interactome_fieldnames = ["source","target","weight"]
		interactome_dataframe = pd.read_csv(interactome_file, delimiter='\t', names=interactome_fieldnames)

		self.nodes = pd.Index(list(set(interactome_dataframe['source'].unique()).union(
								   set(interactome_dataframe['target'].unique()) ) ))
		# node_at_index = lambda k: self.nodes.get_value(k,0)

		self.edges =

		self.costs = interactome_dataframe['weight']

		# self.root  = root
		# self.pruning = pruning
		# self.verbosity_level = verbosity_level



	def pcsf(prizes):
		"""
		"""

		# `edges`: a list of pairs of integers. Each pair specifies an undirected edge in the input graph. The nodes are labeled 0 to n-1, where n is the number of nodes.
		# `prizes`: the list of node prizes.
		# `costs`: the list of edge costs.
		# `root`: the root note for rooted PCST. For the unrooted variant, this parameter should be -1.
		# `num_clusters`: the number of connected components in the output.
		# `pruning`: a string value indicating the pruning method. Possible values are `'none'`, `'simple'`, `'gw'`, and `'strong'` (all literals are case-insensitive).
		# `verbosity_level`: an integer indicating how much debug output the function should produce.
		vertices, edges = pcst_fast(self.edges, prizes, self.costs, self.root, num_clusters, self.pruning, self.verbosity_level)
		# `vertices`: a list of vertices in the output.
		# `edges`: a list of edges in the output. The list contains indices into the list of edges passed into the function.

		return vertices




def main(graph, prize_file, output_dir, options):
	pass




