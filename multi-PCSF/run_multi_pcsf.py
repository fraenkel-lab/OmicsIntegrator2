#!/usr/bin/env python3

#Run dendrogram-based multi-sample PCSF
#constraining similar samples to result in similar networks

#Inputs: 
#   -prize file per sample
#   -One interactome
#   -Dendrogram representing heirarchical clustering of samples

import argparse
import sys, os
from graph import Graph, output_networkx_graph_as_json_for_cytoscapejs

def run_single_PCSF(prizeFile, edgeFile, paramDict, outdir):
    #One straightforward application of PCSF
    graph = Graph(edgeFile, paramDict)
    graph.prepare_prizes(prizeFile)
    vertex_indices, edge_indices = graph.pcsf()
    forest, augmented_forest = graph.output_forest_as_networkx(vertex_indices, edge_indices)
    output_networkx_graph_as_json_for_cytoscapejs(augmented_forest, outdir)
    return forest

def run_param_screen(prizeFile, edgeFile, w_list, b_list, a_list, outdir):
    #For a single clade & sample, run a parameter screen and return the best resulting forest
    #Did Gitter et al vary parameters per iteration? - they tested multiple params but mostly kept things unchanging.
    #TODO Decide how we want to do this


def run_multi_PCSF(dendrogram, prizefiles, edgeFile, paramDict, outdir):
    #Iterate through dendrogram, and at each clade, run forest for each sample, adding artificial prizes

    N = len(prizefiles)
    lastF = []*N #list of N lists of nodes in latest version of each sample
    artificial_prize_dicts = {}*N #list of N dicts of latest artificial prizes in each sample

    #Run the first iteration with unaltered prizes
    #Store list of potential Steiner nodes for each sample
    os.makedirs(outdir + '/initial')
    for i,p in enumerate(prizefiles):
        unadjusted_forest = run_single_PCSF(p, edgeFile, paramDict, outdir + '/initial')  #change this to param screen later?
        lastF[i] = unadjusted_forest.nodes() #NOTE this currently doesn't include a way to distinguish between what was originally a Steiner node or terminal

    #now interate over dendrogram, and at each merge, re-run PCSF for samples in that merge
    for i,c in enumerate(dendrogram):
        os.makedirs(outdir = '/iter%i'%i)
        num_samples_in_clade = int(c[3])
        s_c = float(c[2])
        d_c = 1-s_c
        adj1 = int(c[0])
        adj2 = int(c[1])
        if adj1 < N:
            #then this is an initial sample
        else:
            #need to recursively find all the initial samples in this clade
        #repeat for adj2

        #after this calculation, we'll have gotten a list of samples_in_clade
        if len(samples_in_clade) != num_samples_in_clade: sys.exit()
        #get frequency of nodes in these networks
        forestFreq = nodeFrequency([lastF[k] for k in samples_in_clade])
        for s in samples_in_clade:
            nodes = lastF[s]
            artificial_prizes = artificial_prize_dicts[s]
            for node in forestFreq:
                if node in artificial_prizes:
                    artificial_prizes[node] = artificial_prizes[node] + (s_c*forestFreq[node])^-1*d_c #TODO add lambda and alpha
                else:
                    artificial_prizes[node] = (s_c*forestFreq[node])^-1*d_c #TODO add lambda and alpha

            #submit new artificial prizes + orig prize list to run_single_pcsf
            #update lastF and artificial_prize_dicts
        
def nodeFrequency(list_of_node_lists):
    #return dict with node:freq of all nodes in these lists
    n = float(len(list_of_node_lists)) #want floating point division
    countDict = {}
    for forest in list_of_node_lists:
        for node in forest:
            if node in countDict:
                countDict[node] += 1
            else:
                countDict[node] = 1
    freqDict = {}
    for node in countDict:
        freqDict[node] = countDict[node] / n
    return freqDict
        

def main():
    parser = argparse.ArgumentParser(description="""
	Run the PCSF algorithm several times on hierarchically clustered samples, refining the tree for each sample so that similar
    samples result in similar trees.""")

    parser.add_argument("-e", "--edge", dest='edge_file', type=argparse.FileType('r'), required=True,
	    help ='(Required) Path to the text file containing the edges. Should be a tab delimited file with 3 columns: "nodeA\tnodeB\tweight(between 0 and 1)"')
    parser.add_argument("-p", "--prizefiles", dest='prize_files', type=argparse.FileType('r'), required=True, nargs="+"
	    help='(Required, one or more) Paths to the text files containing the prizes. The list should be in the same order as was provided to create the dendrogram. Should be tab delimited files with lines: "nodeName(tab)prize"')
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
