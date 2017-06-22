#!/usr/bin/env python3

# Core python modules
import sys, os

# Peripheral python modules
import argparse

# import this module
# from . import Graph, output_networkx_graph_as_gml_for_cytoscape, merge_two_prize_files
from . import Graph, output_networkx_graph_as_gml_for_cytoscape, merge_two_prize_files, get_networkx_graph_as_dataframe_of_nodes, get_networkx_graph_as_dataframe_of_edges, output_networkx_graph_as_json_for_cytoscapejs

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
	help='(Required) Path to the text file containing the prizes. Should be a tab delimited file with lines: "nodeName(tab)prize"')
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
	help="Omega: the weight of the edges connecting the dummy node to the nodes selected by dummyMode [default: 6]")
params.add_argument("-b", dest="b", type=float, required=False,
	help="Beta: scaling factor of prizes [default: 1]")
params.add_argument("-mu", dest="mu", type=float, required=False,
	help="Mu: multiplicative node penalty from degree of node [default: 0]")
params.add_argument("-a", dest="a", type=float, required=False,
	help="Alpha: multiplicative edge penalty from degree of endpoints [default: 20]")
params.add_argument("-noise", dest="noise", type=float, required=False,
	help="Standard Deviation of the gaussian noise added to edges in Noisy Edges Randomizations [default: 0.1]")
params.add_argument("--dummyMode", dest='dummy_mode', choices=("terminals", "other", "all"), required=False,
	help='Tells the program which nodes in the interactome to connect the dummy node to. "terminals"= connect to all terminals, "others"= connect to all nodes except for terminals, "all"= connect to all nodes in the interactome. [default: terminals]')
params.add_argument("--muSquared", action='store_true', dest='mu_squared', required=False,
	help='Flag to add negative prizes to hub nodes proportional to their degree^2, rather than degree. Must specify a positive mu in conf file. [default: False]')
params.add_argument("--excludeTerminals", action='store_true', dest='exclude_terminals', required=False,
	help='Flag to exclude terminals when calculating negative prizes. Use if you want terminals to keep exact assigned prize regardless of degree. [default: False]')
params.add_argument("-s", "--seed", dest='seed', type=int, required=False,
	help='An integer seed for the pseudo-random number generators. If you want to reproduce exact results, supply the same seed. [default: None]')



def main():

	args = parser.parse_args()

	params = {"w":args.w, "b":args.b, "mu":args.mu, "a":args.a, "noise":args.noise, "dummy_mode":args.dummy_mode, "mu_squared":args.mu_squared, "exclude_terminals":args.exclude_terminals, "seed":args.seed}
	params = {param: value for param, value in params.items() if value}
	# gross code. http://stackoverflow.com/questions/42400646/is-it-possible-to-denote-some-set-of-argparses-arguments-without-using-subparse

	graph = Graph(args.edge_file, params)

	prizes, terminals, terminal_attributes = graph.prepare_prizes(args.prize_file)

	if args.noisy_edges_repetitions + args.random_terminals_repetitions > 0:
		forest, augmented_forest = graph.randomizations(prizes, terminals, terminal_attributes, args.noisy_edges_repetitions, args.random_terminals_repetitions)

	else:
		vertices, edges = graph.pcsf(prizes)
		forest, augmented_forest = graph.output_forest_as_networkx(vertices, edges, terminal_attributes)

	#output_networkx_graph_as_gml_for_cytoscape(augmented_forest, args.output_dir, 'output.gml')
	output_networkx_graph_as_json_for_cytoscapejs(augmented_forest, args.output_dir)

if __name__ == '__main__':
	main()
