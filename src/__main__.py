#!/usr/bin/env python3

# Core python modules
import sys, os

# Peripheral python modules
import argparse

# import this module
import OmicsIntegrator.graph as oi


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
io_params = parser.add_argument_group("Input / Output Files")

io_params.add_argument("-e", "--edge", dest='edge_file', type=str, required=True,
    help ='(Required) Path to the text file containing the edges. Should be a tab delimited file (no header) with 3 columns: "nodeA\tnodeB\tcost"')
io_params.add_argument("-p", "--prize", dest='prize_file', type=str, required=True,
    help='(Required) Path to the text file containing the prizes. Should be a tab delimited file (*with* header) with lines: "nodeName(tab)prize"')
io_params.add_argument('-o', '--output', dest='output_dir', action=FullPaths, type=directory, required=True,
    help='(Required) Output directory path')

# Command parameters (specify what the algorithm does):
pcsf_params = parser.add_argument_group("PCSF Parameters")

pcsf_params.add_argument("-w", dest="w", nargs="*", type=float, required=False, default=6,
    help="Omega: the weight of the edges connecting the dummy node to the nodes selected by dummyMode [default: 6]")
pcsf_params.add_argument("-b", dest="b", nargs="*", type=float, required=False, default=1,
    help="Beta: scaling factor of prizes [default: 1]")
pcsf_params.add_argument("-g", dest="g", nargs="*", type=float, required=False, default=20,
    help="Gamma: multiplicative edge penalty from degree of endpoints [default: 20]")
pcsf_params.add_argument("-noise", dest="edge_noise", type=float, required=False, default=0.1,
    help="Standard Deviation of the gaussian noise added to edges in Noisy Edges Randomizations [default: 0.1]")
pcsf_params.add_argument("--noisy_edges", dest='noisy_edges_repetitions', type=int, default=0,
    help='An integer specifying how many times you would like to add noise to the given edge values and re-run the algorithm. Results of these runs will be merged together and written in files with the word "_noisy_edges_" added to their names. The noise level can be controlled using the configuration file. [default: 0]')
pcsf_params.add_argument("--random_terminals", dest='random_terminals_repetitions', type=int, default=0,
    help='An integer specifying how many times you would like to apply your given prizes to random nodes in the interactome (with a similar degree distribution) and re-run the algorithm. Results of these runs will be merged together and written in files with the word "_random_terminals_" added to their names. [default: 0]')
pcsf_params.add_argument("--dummyMode", dest='dummy_mode', choices=("terminals", "other", "all"), required=False,
    help='Tells the program which nodes in the interactome to connect the dummy node to. "terminals"= connect to all terminals, "others"= connect to all nodes except for terminals, "all"= connect to all nodes in the interactome. [default: terminals]')
pcsf_params.add_argument("-s", "--seed", dest='seed', type=int, required=False,
    help='An integer seed for the pseudo-random number generators. If you want to reproduce exact results, supply the same seed. [default: None]')


def output_dataframe_to_tsv(dataframe, output_dir, filename):
    path = os.path.join(os.path.abspath(output_dir), filename)
    dataframe.to_csv(path, sep='\t', header=True, index=False)


def main():

    args = parser.parse_args()

    params = {"w":args.w, "b":args.b, "g":args.g, "edge_noise":args.edge_noise, "dummy_mode":args.dummy_mode, "seed":args.seed,
              "noisy_edges_repetitions": args.noisy_edges_repetitions, "random_terminals_repetitions": args.random_terminals_repetitions}

    params = {param: value for param, value in params.items() if value is not None}

    graph = oi.Graph(args.edge_file, params)
    graph.prepare_prizes(args.prize_file)

    if args.noisy_edges_repetitions + args.random_terminals_repetitions > 0:
        forest, augmented_forest = graph.randomizations(args.noisy_edges_repetitions, args.random_terminals_repetitions)

    else:
        vertex_indices, edge_indices = graph.pcsf()
        forest, augmented_forest = graph.output_forest_as_networkx(vertex_indices, edge_indices)

    #oi.output_networkx_graph_as_graphml_for_cytoscape(augmented_forest, args.output_dir)
    oi.output_networkx_graph_as_interactive_html(augmented_forest, args.output_dir)

if __name__ == '__main__':
    main()
