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
parser.add_argument("-i", "--interactome", dest='interactome_file', type=argparse.FileType('r'), required=True,
	help ='(Required) Path to the text file containing the interactome edges. Should be a tab delimited file with 3 or 4 columns: "ProteinA\tProteinB\tWeight(between 0 and 1)\tDirectionality(U or D, optional)"')
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
parser.add_argument("--excludeTerms", action='store_true', dest='exclude_terms', default=False,
	help='Flag to exclude terminals when calculating negative prizes. Use if you want terminals to keep exact assigned prize regardless of degree. [default: %default]')

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

	options = Options({"w": 6 , "b": 12, :"gb": 0.1, "D": 6, "mu": 0.04})

	graph = Graph(edge_file, options)

	main(graph, prize_file, output_dir, options)


class Options:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


class Graph:
	"""docstring for Graph"""
	def __init__(self, interactome_file, options):
		"""
		Builds a representation of a graph from an interactome file.

		From the interactome_file, populates `self.nodes` (pd.Index), `self.edges` (list of pairs),
		`self.costs` (list, such that the ordering is the same as in self.edges), and
		`self.node_degrees` (list, such that the ordering is the same as in self.nodes).

		From the prize_file, populates self.terminals (list) and self.prizes (list which contains 0
		everywhere there isn't a terminal with an assigned prize).

		From the garnet_file, merge the TF terminals and prizes with `self.terminals` and `self.prizes`.

		Arguments:
			interactome_file (str or FILE): tab-delimited text file containing edges in interactome and their weights formatted like "ProteinA\tProteinB\tWeight"
			prize_file (str or FILE): tab-delimited text file containing all proteins with prizes formatted like "ProteinName\tPrizeValue"
			options (dict): options with which to run the program

		self.w, self.b, self.D, self.gb, self.mu, self.g, self.r, self.noise
		"""

		interactome_fieldnames = ["source","target","weight"]
		interactome_dataframe = pd.read_csv(interactome_file, delimiter='\t', names=interactome_fieldnames)

		(self.edges, self.nodes) = pd.factorize(interactome_dataframe[["source","target"]].unstack())

		self.edges = self.edges.reshape(interactome_dataframe[["source","target"]].shape, order='F').tolist()  # maybe needs to be list of tuples. in that case map(tuple, this)

		self.costs = interactome_dataframe['weight'].tolist()

		self.node_degrees = [ ]
		# possible to flatten edges, count occurrences / 2 of each number, and then sort by ID number (v easy)
		# might also be nice to sort by index, so that you can find the distribution by gene symbol or whatever


		prizes_dataframe = pd.read_csv(prizes_file, delimiter='\t', names=["name","prize"])

		self.terminals = prizes_dataframe["name"]

		self.prizes = [ for node in self.nodes if ]  # I still don't know how to merge pandas indices


		self._add_dummy_node(connected_to=options.dummy_mode)

		self._check_validity_of_instance()

		defaults = {"w": 6 , "b": 12, :"gb": 0.1, "D": 6, "mu": 0.04}
		self.params = Options(defaults.update(options))




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
		vertices, edges = pcst_fast(self.edges, prizes, self.costs, self.root, num_clusters, pruning, verbosity_level)
		# `vertices`: a list of vertices in the output.
		# `edges`: a list of edges in the output. The list contains indices into the list of edges passed into the function.

		return vertices, edges

			   # betweenness - a boolean flag indicating whether we should do the costly betweenness
			   #               calculation
		# OUTPUT: self.optForest - a networkx digraph storing the forest returned by msgsteiner
		#         self.augForest - a networkx digraph storing the forest returned by msgsteiner, plus
		#                          all of the edges in the interactome between nodes in that forest
		#         self.dumForest - a networkx digraph storing the dummy node edges in the optimal
		#                          forest


	def noiseEdges(self, seed=None):
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
		edges = [edge + random.gauss(0,std_dev) for edge in self.edges]

		return edges


	def randomTerminals(self, seed, excludeT):
		"""
		Succinct description of randomTerminals

		Selects nodes with a similar degree distribution to the original terminals, and assigns the prizes to them.

		Arguments:
			seed (): a random seed

		Returns:
			list: new terminal nodes
		"""

		if len(PCSFInputObj.undirEdges) + len(PCSFInputObj.dirEdges) < 50: sys.exit("Cannot use --randomTerminals with such a small interactome.")



		#Find index of current terminal in degrees list
		for k,terminal in enumerate(PCSFInputObj.origPrizes):
			for i,value in enumerate(degrees):
				if terminal == value[0]:
					index = i
					break
			#Choose an index offset to select new terminal (distance from orig terminal in degrees list)
			#Make sure newly chosen terminal is not already chosen on a previous round
			newTerm = ''
			i = -1
			while newTerm in newPCSFInputObj.origPrizes and i<=10000:
				i+=1
				if seed != None: random.seed(seed+k+i)
				offset = int(random.gauss(0.0,100.0))
				newIndex = index + offset
				#if offset causes the index to wraparound to the other side of the list, try again
				if newIndex<0: continue
				try:
					newNode = degrees[newIndex]
				except IndexError:
					#if offset points outside list, try loop again
					continue
				#To make truly random, need to choose randomly between all nodes with the same degree
				#Otherwise, ordering of dict iteration matters
				nodesWithSameDegree = []
				for node in degrees[newIndex:]:
					if node[1] == newNode[1]:
						nodesWithSameDegree.append(node)
					else:
						break
				for node in degrees[newIndex-1::-1]:
					if node[1] == newNode[1]:
						nodesWithSameDegree.append(node)
					else:
						break
				newTerm = random.choice(nodesWithSameDegree)[0]
			#if we've tried 10000 times, throw error to avoid infinite loop
			if newTerm in newPCSFInputObj.origPrizes:
				sys.exit('There was a problem with --randomTerminals. Aborting.')
			#Assign prize to newly chosen terminal
			newPCSFInputObj.origPrizes[newTerm] = PCSFInputObj.origPrizes[terminal]
		del newPCSFInputObj.origPrizes['']
		newPCSFInputObj.assignNegPrizes(newPCSFInputObj.musquared,excludeT)
		print 'New degree-matched terminals have been chosen.\n'
		return newPCSFInputObj


def changeValuesAndMergeResults(func, seed, inputObj, numRuns, msgpath, outputpath, outputlabel, excludeT):
	"""
	Changes the prizes/edges in the PCSFInput object according to func and runs the msgsteiner
	algorithm, then merges the results together with the given PCSFOutput object. Writes
	cytoscape files for final merged results.

	INPUT: func - the function which takes inputObj and changes the prize/edge values
				  (i.e. shuffles or adds noise)
		   numRums - the number of times to change the values and re-run msgsteiner
		   msgpath - path to the directory where msgsteiner is kept
		   outputpath - path to the directory where output files should be stored
		   outputlabel - a label with which to name all of the output files for this run

	OUTPUT: <outputlabel>_changed_#_info.txt - a text file FOR EACH RUN containing the
					  contents of stderr for all msgsteiner runs
	RETURNS: merged - the PCSFOutput object that is a result of all the merges

	"""
	print 'Preparing to change values %i times and get merged results of running the '\
		  'algorithm on new values.\n' %numRuns
	#Create multiprocessing Pool
	if inputObj.processes == None:
		pool = mp.Pool()
	else:
		pool = mp.Pool(inputObj.processes)
	if seed != None:
		#For each run, create process, change prize/edge values and run msgsteiner
		#Note that each run will create a info file
		results = [pool.apply_async(PCSF, args=(func(inputObj, seed+i, excludeT),
				  msgpath, seed+i,)) for i in xrange(numRuns)]
	else:
		results = [pool.apply_async(PCSF, args=(func(inputObj, seed, excludeT),
				  msgpath,seed,)) for i in xrange(numRuns)]
	output = [p.get() for p in results]
	i = 0
	#Merge output of new msgsteiner runs together
	while i <= numRuns-1:
		(newEdgeList, newInfo) = output[i]
		#By creating the output object with inputObj instead of changedInputObj,
		#the prizes stored in the networkx graphs will be the ORIGINAL CORRECT prizes,
		#not the changed prizes.
		if str(func)[10:23]  == 'shufflePrizes':
			changedOutputObj = PCSFOutput(inputObj, newEdgeList, newInfo, outputpath, outputlabel+'_shuffledPrizes_%i'%i, 0)
		elif str(func)[10:20] == 'noiseEdges':
			changedOutputObj = PCSFOutput(inputObj, newEdgeList, newInfo, outputpath, outputlabel+'_noisyEdges_%i'%i, 0)
		elif str(func)[10:25] == 'randomTerminals':
			changedOutputObj = PCSFOutput(inputObj, newEdgeList, newInfo, outputpath, outputlabel+'_randomTerminals_%i'%i, 0)
		if i == 0:
			#first run
			merged = changedOutputObj
		elif i == numRuns-1:
			#last run, merge results and calculate betweenness
			merged = mergeOutputs(merged, changedOutputObj, 1, i, 1)
		else:
			#Merge results of runs with the merged object containing all results so far
			merged = mergeOutputs(merged, changedOutputObj, 0, i, 1)
		i += 1
	#return merged outputobj
	return merged


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




def main(graph, prize_file, output_dir, options):
	pass




