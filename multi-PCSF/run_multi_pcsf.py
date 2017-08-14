#!/usr/bin/env python3

#Run dendrogram-based multi-sample PCSF
#constraining similar samples to result in similar networks

#Inputs: 
#   -prize file per sample
#   -One interactome
#   -Dendrogram representing heirarchical clustering of samples

import argparse
from graph import Graph, output_networkx_graph_as_json_for_cytoscapejs

def run_single_PCSF(edgeFile, paramDict, outdir):
    #One straightforward application of PCSF
    graph = Graph(edgeFile, paramDict)
    vertex_indices, edge_indices = graph.pcsf()
    forest, augmented_forest = graph.output_forest_as_networkx(vertex_indices, edge_indices)
    output_networkx_graph_as_json_for_cytoscapejs(augmented_forest, outdir)
    return forest

def run_param_screen(prizeFile, edgeFile, w_list, b_list, a_list, outdir):
    #For a single clade & sample, run a parameter screen and return the best resulting forest
    #Did Gitter et al vary parameters per iteration? - they tested multiple params but mostly kept things unchanging.
    #TODO Decide how we want to do this


def run_multi_PCSF(dendrogram, prizefiles, edgeFile):
    #Iterate through dendrogram, and at each clade, run param screen for each sample, adding artificial prizes



def main():
    parser = argparse.ArgumentParser(description="""
	Run the PCSF algorithm several times on hierarchically clustered samples, refining the tree for each sample so that similar
    samples result in similar trees.""")

    parser.add_argument("-e", "--edge", dest='edge_file', type=argparse.FileType('r'), required=True,
	    help ='(Required) Path to the text file containing the edges. Should be a tab delimited file with 3 columns: "nodeA\tnodeB\tweight(between 0 and 1)"')
    parser.add_argument("-p", "--prizefiles", dest='prize_files', type=argparse.FileType('r'), required=True, nargs="+"
	    help='(Required, one or more) Paths to the text files containing the prizes. Should be tab delimited files with lines: "nodeName(tab)prize"')
    parset.add_argument("-d", "--dendrogram", dest="dendrogram", required=True,
        help='(Required) object denoting hierarchical clustering of samples, of the type returned by scipy.cluster.heirarchy\'s linkage().. Should be an array of length n-1, where dendrogram[i] indicates which clusters are merged at the i-th iteration.')
    parser.add_argument('-o', '--output', dest='output_dir', action=FullPaths, type=directory, required=True,
	    help='(Required) Output directory path')
    parser. add_argument("-ws", "--w_list", dest="w_list", nargs="+", default=[5,10], 
        help="A list of integers for the parameter w (number of trees). default='5 10'")
    parser. add_argument("-bs", "--b_list", dest="b_list", nargs="+", default=[5,10], 
        help="A list of integers for the parameter b (size of network). default='5 10'")
    parser. add_argument("-as", "--a_list", dest="a_list", nargs="+", default=[0,10000,100000], 
        help="A list of integers for the parameter a (negative prize on hubs). default='0 10000 100000'")
    params.add_argument("-s", "--seed", dest='seed', type=int, required=False,
	    help='An integer seed for the pseudo-random number generators. If you want to reproduce exact results, supply the same seed. [default: None]')

    args = parser.parse_args()




if __name__ == '__main__':
	main()
