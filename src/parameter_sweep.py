#!/usr/bin/env python3

# Core python modules
import sys, os

# Peripheral python modules
import argparse

# import this module
# from . import Graph, output_networkx_graph_as_graphml_for_cytoscape, output_networkx_graph_as_json_for_cytoscapejs
from graph import Graph, output_networkx_graph_as_graphml_for_cytoscape, output_networkx_graph_as_json_for_cytoscapejs, output_networkx_graph_as_pickle, get_networkx_graph_as_node_edge_dataframes, get_networkx_subgraph_from_randomizations

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
parser.add_argument("-e", "--edge", dest='edge_file', type=str, required=True,
	help ='(Required) Path to the text file containing the edges. Should be a tab delimited file with 3 columns: "nodeA\tnodeB\tcost"')
parser.add_argument("-p", "--prize", dest='prize_file', type=str, required=True,
	help='(Required) Path to the text file containing the prizes. Should be a tab delimited file with lines: "nodeName(tab)prize"')
parser.add_argument('-o', '--output', dest='output_dir', action=FullPaths, type=directory, required=True,
	help='(Required) Output directory path')

# Command parameters (specify what the algorithm does):
parser.add_argument("--noisy_edges", dest='noisy_edges_repetitions', type=int, default=10,
	help='An integer specifying how many times you would like to add noise to the given edge values and re-run the algorithm. Results of these runs will be merged together and written in files with the word "_noisy_edges_" added to their names. The noise level can be controlled using the configuration file. [default: %default]')
parser.add_argument("--random_terminals", dest='random_terminals_repetitions', type=int, default=10,
	help='An integer specifying how many times you would like to apply your given prizes to random nodes in the interactome (with a similar degree distribution) and re-run the algorithm. Results of these runs will be merged together and written in files with the word "_random_terminals_" added to their names. [default: %default]')
parser.add_argument("--knockout", dest='knockout', nargs='*', default=[],
	help='Protein(s) you would like to "knock out" of the interactome to simulate a knockout experiment. [default: []]')

params = parser.add_argument_group('Parameters', 'Parameters description')

params.add_argument("-w", dest="w", nargs="*", type=float, required=False, default=[6],
	help="Omega: the weight of the edges connecting the dummy node to the nodes selected by dummyMode [default: 6]")
params.add_argument("-b", dest="b", nargs="*", type=float, required=False, default=[1],
	help="Beta: scaling factor of prizes [default: 1]")
params.add_argument("-g", dest="g", nargs="*", type=float, required=False, default=[20],
	help="Gamma: multiplicative edge penalty from degree of endpoints [default: 20]")
params.add_argument("-noise", dest="noise", type=float, required=False, default=0.1,
	help="Standard Deviation of the gaussian noise added to edges in Noisy Edges Randomizations [default: 0.1]")
params.add_argument("--dummyMode", dest='dummy_mode', choices=("terminals", "other", "all"), required=False,
	help='Tells the program which nodes in the interactome to connect the dummy node to. "terminals"= connect to all terminals, "others"= connect to all nodes except for terminals, "all"= connect to all nodes in the interactome. [default: terminals]')
params.add_argument("--excludeTerminals", action='store_true', dest='exclude_terminals', required=False,
	help='Flag to exclude terminals when calculating negative prizes. Use if you want terminals to keep exact assigned prize regardless of degree. [default: False]')
params.add_argument("-s", "--seed", dest='seed', type=int, required=False,
	help='An integer seed for the pseudo-random number generators. If you want to reproduce exact results, supply the same seed. [default: None]')


def output_dataframe_to_tsv(dataframe, output_dir, filename):
	"""
	Output the dataframe to a csv
	"""

	path = os.path.join(os.path.abspath(output_dir), filename)
	dataframe.to_csv(path, sep='\t', header=True, index=False)


def main():

	args = parser.parse_args()

	params = {"w":args.w, "b":args.b, "g":args.g, "noise":args.noise, "dummy_mode":args.dummy_mode, "exclude_terminals":args.exclude_terminals, "seed":args.seed,
			  "noisy_edges_repetitions": args.noisy_edges_repetitions, "random_terminals_repetitions": args.random_terminals_repetitions}
	params = {param: value for param, value in params.items() if value}

	graph = Graph(args.edge_file, {})

	print(params)

	# Parameter search
	results = graph.grid_search_randomizations(args.prize_file, params)

	for tag, forest, augmented_forest in results: 

		augmented_nodes_df, augmented_edges_df = get_networkx_graph_as_node_edge_dataframes(augmented_forest)

		# Get top 400 nodes as subnetwork of augmented forest
		robust_net = get_networkx_subgraph_from_randomizations(augmented_forest, max_size=400)
		robust_net_nodes_df, robust_net_edges_df = get_networkx_graph_as_node_edge_dataframes(robust_net)

		# Save augmented forest as pickled networkx object
		output_networkx_graph_as_pickle(augmented_forest, args.output_dir, tag+".augmented_forest.gpickle")

		# Write node and edge attributes for augmented network to files
		output_dataframe_to_tsv(augmented_nodes_df, args.output_dir, tag+".augmented_forest.nodes.tsv")
		output_dataframe_to_tsv(augmented_edges_df, args.output_dir, tag+".augmented_forest.edges.tsv")

		# Write json for robust network
		output_networkx_graph_as_json_for_cytoscapejs(robust_net, args.output_dir, tag+".robust_network.json")


main()