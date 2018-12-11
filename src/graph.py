#!/usr/bin/env python3

# Core python modules
import sys
import os
import multiprocessing
import logging
import random
import numbers
import math

# Peripheral python modules
import argparse
from collections import Counter
from itertools import product
from pathlib import Path
from copy import copy
import json
from pkg_resources import resource_filename as get_path

# python external libraries
import numpy as np
import pandas as pd
import networkx as nx
from networkx.readwrite import json_graph as nx_json
import community    # pip install python-louvain
from sklearn.cluster import SpectralClustering

# Lab modules
from pcst_fast import pcst_fast
from axial import axial

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
handler.setFormatter(logging.Formatter('%(asctime)s - OI2: %(levelname)s - %(message)s', "%I:%M:%S"))
logger.addHandler(handler)

# Count the number of available CPUs for potential use in multiprocessing code
try: n_cpus = int(os.environ["SLURM_JOB_CPUS_PER_NODE"])
except KeyError: n_cpus = multiprocessing.cpu_count()

# Helpers

def flatten(list_of_lists): return [item for sublist in list_of_lists for item in sublist]

def invert(list_of_lists): return {item: i for i, list in enumerate(list_of_lists) for item in list}

def safe_string(unsafe_string): return ''.join(e for e in unsafe_string if e.isalnum())

class Options(object):
    def __init__(self, options):
        self.__dict__.update(options)
    def __repr__(self):
        return dict(self.__dict__)


class Graph:
    """
    A Graph object is a representation of a graph, with convenience methods for using the pcst_fast
    package, which approximately minimizes the Prize-Collecting Steiner Forest objective.
    """

    ###########################################################################
                #######          Initialization            #######
    ###########################################################################

    def __init__(self, interactome_file, params=dict()):
        """
        Builds a representation of a graph from an interactome file.

        From the interactome_file, populates
        - `graph.interactome_dataframe` (pandas.DataFrame)
        - `graph.interactome_graph` (networkx.Graph)
        - `graph.nodes` (pandas.Index),
        - `graph.edges` (list of pairs),
        - `graph.costs` and `graph.edge_penalties` (lists, such that the ordering is the same as in graph.edges),
        - `graph.node_degrees` (list, such that the ordering is the same as in graph.nodes).

        Arguments:
            interactome_file (str or FILE): tab-delimited text file containing edges in interactome and their weights formatted like "ProteinA\tProteinB\tCost"
            params (dict): params with which to run the program
        """

        self.interactome_dataframe = pd.read_csv(interactome_file, sep='\t')
        self.interactome_graph = nx.from_pandas_edgelist(self.interactome_dataframe, 'protein1', 'protein2', edge_attr=self.interactome_dataframe.columns[2:].tolist())

        # Convert the interactome dataframe from string interactor IDs to integer interactor IDs.
        # Do so by selecting the protein1 and protein2 columns from the interactome dataframe,
        # then unstacking them, which (unintuitively) stacks them into one column, allowing us to use factorize.
        # Factorize builds two datastructures, a unique pd.Index which maps each ID string to an integer ID,
        # and the datastructure we passed in with string IDs replaced with those integer IDs.
        # We place those in self.nodes and self.edges respectively, but self.edges will need reshaping.
        (self.edges, self.nodes) = pd.factorize(self.interactome_dataframe[["protein1","protein2"]].unstack())
        # Here we do the inverse operation of "unstack" above, which gives us an interpretable edges datastructure
        self.edges = self.edges.reshape(self.interactome_dataframe[["protein1","protein2"]].shape, order='F')

        self.edge_costs = self.interactome_dataframe['cost'].astype(float).values

        # Count the number of incident edges into each node.
        # The indices into this datastructure are the same as those in self.nodes which are the IDs in self.edges.
        self.node_degrees = np.bincount(self.edges.flatten())

        # The rest of the setup work is occasionally repeated, so use another method to complete setup.
        self._reset_hyperparameters(params=params)


    def _reset_hyperparameters(self, params=dict()):
        """
        Set the parameters on Graph and compute parameter-dependent features.

        Arguments:
            params (dict): params with which to run the program
        """

        defaults = {"w": 5, "b": 1, "g": 3, "edge_noise": 0.1, "dummy_mode": "terminals", "seed": 0, "skip_checks": False}

        # Overwrite the defaults with any user-specified parameters.
        self.params = Options({**defaults, **params})

        if not self.params.skip_checks: self._check_validity_of_hyperparameters()

        # Add costs to each edge, proportional to the degrees of the nodes it connects, modulated by parameter g.
        N = len(self.nodes)
        self.edge_penalties = (10**self.params.g) * np.array([self.node_degrees[a] * self.node_degrees[b] /
                            ((N - self.node_degrees[a] - 1) * (N - self.node_degrees[b] - 1) + self.node_degrees[a] * self.node_degrees[b]) for a, b in self.edges])

        self.costs = (self.edge_costs + self.edge_penalties)

        # If this instance of graph has bare_prizes set, then presumably resetting the
        # hyperparameters should also reset the scaled prizes
        if hasattr(self, "bare_prizes"): self.prizes = self.bare_prizes * self.params.b


    def _check_validity_of_hyperparameters(self):
        """
        Assert that the hyperparameters passed to this program are valid, otherwise raise helpful error messages.
        """

        if not (isinstance(self.params.w, numbers.Number) and (self.params.w >= 0)): raise ValueError("parameter w must be a positive number. Was "+str(self.params.w))
        if not (isinstance(self.params.b, numbers.Number) and (self.params.b >= 0)): raise ValueError("parameter b must be a positive number. Was "+str(self.params.b))
        if not (isinstance(self.params.g, numbers.Number) and (self.params.g >= 0)): raise ValueError("parameter g must be a positive number. Was "+str(self.params.g))
        if not (isinstance(self.params.edge_noise, numbers.Number) and (self.params.edge_noise >= 0)): raise ValueError("parameter edge_noise must be a positive number. Was "+str(self.params.edge_noise))
        if not (self.params.dummy_mode in ['terminals', 'other', 'all']): raise ValueError("parameter dummy_mode must be one of 'terminals', 'other', or 'all'. Was "+str(self.params.dummy_mode))
        if not (isinstance(self.params.seed, int)): raise ValueError("parameter seed must be a int. Was type "+str(type(self.params.seed)))


    def prepare_prizes(self, prize_file):
        """
        Parses a prize file and adds prizes and other attributes to the graph object.

        The file passed to this function must have at least two columns: node name and prize.
        Any additional columns will be assumed to be node attributes. However, in order to know
        the names of those attributes, this function requires the input file to contain headers,
        i.e. the first row of the tsv must be the names of the columns.

        Sets the graph attributes
        - `graph.bare_prizes` (numpy.array): properly indexed (same as `graph.nodes`) prizes from the file.
        - `graph.prizes` (numpy.array): properly indexed prizes, scaled by beta (`graph.params.b`)
        - `graph.terminals` (numpy.array): their indices
        - `graph.node_attributes` (pandas.DataFrame) Any node attributes passed in with the prize file (columns 3, ...)

        Arguments:
            prize_file (str or FILE): a filepath or file object containing a tsv **with column headers**.
        """

        prizes_dataframe = pd.read_csv(prize_file, sep='\t')
        prizes_dataframe.columns = ['name', 'prize'] + prizes_dataframe.columns[2:].tolist()
        prizes_dataframe['prize'] = pd.to_numeric(prizes_dataframe['prize'])

        return self._prepare_prizes(prizes_dataframe)


    def _prepare_prizes(self, prizes_dataframe):

        # Some files have duplicated genes with different prizes. Keep the max prize.
        logger.info("Duplicated gene symbols in the prize file (we'll keep the max prize):")
        logger.info(prizes_dataframe[prizes_dataframe.set_index('name').index.duplicated()]['name'].tolist())
        prizes_dataframe = prizes_dataframe.groupby('name').max().reset_index()

        # Find the indices of the nodes being assigned prizes
        prizes_dataframe.set_index(self.nodes.get_indexer(prizes_dataframe['name']), inplace=True)

        # there will be some nodes in the prize file which we don't have in our interactome
        logger.info("Members of the prize file not present in the interactome:")
        logger.info(prizes_dataframe[prizes_dataframe.index == -1]['name'].tolist())
        prizes_dataframe.drop(-1, inplace=True, errors='ignore')

        # All nodes in the prizes dataframe are defined to be terminals.
        prizes_dataframe["terminal"] = True
        prizes_dataframe["type"] = prizes_dataframe.get("type", default="protein")

        # Node attributes dataframe for all proteins in self.nodes.
        self.node_attributes = prizes_dataframe.set_index('name').rename_axis(None).reindex(self.nodes)
        self.node_attributes["degree"] = self.node_degrees
        self.node_attributes["prize"].fillna(0, inplace=True)
        self.node_attributes["type"].fillna("protein", inplace=True)
        self.node_attributes["terminal"].fillna(False, inplace=True)

        self.bare_prizes = self.node_attributes["prize"].values
        self.prizes = self.bare_prizes * self.params.b

        self.terminals = np.where(self.node_attributes["terminal"] == True)[0]


    ###########################################################################
                #######              PCSF               #######
    ###########################################################################

    def _add_dummy_node(self, connected_to=[]):

        dummy_id = len(self.nodes)
        dummy_prize = np.array([0])
        dummy_edges = np.array([(dummy_id, node_id) for node_id in connected_to])
        dummy_costs = np.array([self.params.w] * len(dummy_edges))

        return dummy_edges, dummy_costs, dummy_id, dummy_prize


    def _check_validity_of_instance(self, edges, prizes, costs, root, num_clusters, pruning, verbosity_level):
        """
        Assert that the data passed to this program are valid, otherwise raise helpful error messages.
        """

        if not (isinstance(edges, np.ndarray)): raise ValueError("edges must be a numpy array, type was: "+str(type(edges)))
        if not (len(edges.shape) == 2): raise ValueError("edges must be an array of dimension 2, dimension was: "+str(len(edges.shape)))
        if not (edges.shape[1] == 2): raise ValueError("edges array must have two columns, number of columns was: "+str(edges.shape[1]))

        if not (isinstance(prizes, np.ndarray)): raise ValueError("prizes must be a numpy array, type was: "+str(type(prizes)))
        if not (len(prizes.shape) == 1): raise ValueError("prizes must be an array of dimension 1, dimension was: "+str(len(prizes.shape)))
        if not (len(prizes) == len(np.unique(edges.flatten()))): raise ValueError("there must be as many prizes as nodes. # prizes: "+str(len(prizes))+", # nodes: "+str(len(np.unique(edges.flatten()))))

        if not (isinstance(costs, np.ndarray)): raise ValueError("costs must be a numpy array, type was: "+str(type(costs)))
        if not (len(costs.shape) == 1): raise ValueError("costs must be an array of dimension 1, dimension was: "+str(len(costs.shape)))
        if not (len(costs) == len(edges)): raise ValueError("there must be as many costs as edges. # costs: "+str(len(costs))+", # edges: "+str(len(edges)))

        if not (isinstance(root, int)): raise ValueError("root must be an int, type was: "+str(type(root)))
        if not (0 <= root < len(prizes)): raise ValueError("root must be one of the nodes in the graph. root: "+str(root)+", nodes: [0, "+str(len(prizes-1))+"]")

        if not (isinstance(num_clusters, int)): raise ValueError("num_clusters must be an int, type was: "+str(type(num_clusters)))
        if not (0 < num_clusters < len(prizes)): raise ValueError("num_clusters must be greater than 0, and less than the number of nodes. num_clusters was: "+str(num_clusters)+"# nodes was: "+str(len(prizes)))

        if not (pruning in ['none', 'simple', 'gw', 'strong']): raise ValueError("pruning must be one of ['none', 'simple', 'gw', 'strong']. pruning was: "+str(pruning))

        if not (verbosity_level in [0, 1, 2, 3]): raise ValueError("verbosity_level must be an integer, any of [0, 1, 2, 3]. verbosity_level was: "+str(verbosity_level))

        return True


    def pcsf(self, pruning="strong", verbosity_level=0):
        """
        Select the subgraph which approximately optimizes the Prize-Collecting Steiner Forest objective.

        This function mostly defers to pcst_fast, but does one important pre-processing step: it
        adds a dummy node which will serve as the PCSF root and connects that dummy node to either
        terminals, non-terminals, or all other nodes with edges weighted by self.params.w.

        In order to interpret the results of this function, use `output_forest_as_networkx` with
        the results of this function.

        Arguments:
            pruning (str): a string value indicating the pruning method. Possible values are `'none'`, `'simple'`, `'gw'`, and `'strong'` (all literals are case-insensitive).
            verbosity_level (int): an integer indicating how much debug output the function should produce.

        Returns:
            numpy.array: indices of the selected vertices
            numpy.array: indices of the selected edges
        """

        all = list(range(len(self.nodes)))
        others = list(set(all) - set(self.terminals))

        if self.params.dummy_mode == 'terminals': endpoints = self.terminals
        elif self.params.dummy_mode == 'other': endpoints = others
        elif self.params.dummy_mode == 'all': endpoints = all
        else: raise ValueError("Improper input to PCSF: dummy_mode must be one of 'terminals', 'other', or 'all'")

        dummy_edges, dummy_costs, root, dummy_prize = self._add_dummy_node(connected_to=endpoints)

        # `edges`: a 2D int64 array. Each row (of length 2) specifies an undirected edge in the input graph. The nodes are labeled 0 to n-1, where n is the number of nodes.
        edges = np.concatenate((self.edges, dummy_edges))
        # `prizes`: the node prizes as a 1D float64 array.
        prizes = np.concatenate((self.prizes, dummy_prize))
        # `costs`: the edge costs as a 1D float64 array.
        costs = np.concatenate((self.costs, dummy_costs))
        # `root`: the root note for rooted PCST. For the unrooted variant, this parameter should be -1.
        # `num_clusters`: the number of connected components in the output.
        num_clusters = 1
        # `pruning`: a string value indicating the pruning method. Possible values are `'none'`, `'simple'`, `'gw'`, and `'strong'` (all literals are case-insensitive).
        # `verbosity_level`: an integer indicating how much debug output the function should produce.

        if not self.params.skip_checks: self._check_validity_of_instance(edges, prizes, costs, root, num_clusters, pruning, verbosity_level)

        vertex_indices, edge_indices = pcst_fast(edges, prizes, costs, root, num_clusters, pruning, verbosity_level)
        # `vertex_indices`: indices of the vertices in the solution as a 1D int64 array.
        # `edge_indices`: indices of the edges in the output as a 1D int64 array. The list contains indices into the list of edges passed into the function.

        # Remove the dummy node and dummy edges for convenience
        vertex_indices = vertex_indices[vertex_indices != root]
        edge_indices = edge_indices[np.in1d(edge_indices, self.interactome_dataframe.index)]

        return vertex_indices, edge_indices


    def output_forest_as_networkx(self, vertex_indices, edge_indices):
        """
        Construct a networkx graph from a set of vertex and edge indices (i.e. a pcsf output)

        Arguments:
            vertex_indices (list): indices of the vertices selected in self.nodes. Note, this list must be of type int or boolean. Errors may arise from empty numpy arrays of dtype='object'
            edge_indices (list): indices of the edges selected in self.edges

        Returns:
            networkx.Graph: a networkx graph object
        """

        if len(vertex_indices) == 0:
            logger.warning("The resulting Forest is empty. Try different parameters.")
            return nx.empty_graph(0), nx.empty_graph(0)

        # Replace the edge indices with the actual edges (protein1 name, protein2 name) by indexing into the interactome
        edges = self.interactome_dataframe.loc[edge_indices]
        forest = nx.from_pandas_edgelist(edges, 'protein1', 'protein2', edge_attr=True)
        # the above won't capture the singletons, so we'll add them here
        forest.add_nodes_from(list(set(self.nodes[vertex_indices]) - set(forest.nodes())))

        # Set all the attributes on graph
        nx.set_node_attributes(forest, self.node_attributes.reindex(list(forest.nodes())).dropna(how='all').to_dict(orient='index'))
        # Set a flag on all the edges which were selected by PCSF (before augmenting the forest)
        nx.set_edge_attributes(forest, True, name='in_solution')
        # Create a new graph including all edges between all selected nodes, not just those edges selected by PCSF.
        augmented_forest = nx.compose(self.interactome_graph.subgraph(forest.nodes()), forest)

        # Post-processing
        betweenness(augmented_forest)
        louvain_clustering(augmented_forest)
        annotate_graph_nodes(augmented_forest)

        return forest, augmented_forest


    def pcsf_objective_value(self, forest):
        """
        Calculate PCSF objective function

        Arguments:
            forest (networkx.Graph): a forest like the one returned by output_forest_as_networkx -- Not an augmented forest!

        Returns:
            float: PCSF objective function score
        """

        return ((sum(self.prizes) - sum(nx.get_node_attributes(forest, 'prize').values())) +
                 sum(nx.get_edge_attributes(forest, 'cost').values()) +
                 (self.params.w * nx.number_connected_components(forest)))


    ###########################################################################
                #######          Randomziations         #######
    ###########################################################################

    def _noisy_edges(self):
        """
        Adds gaussian edge_noise to all edge costs in the graph, modulated by parameter `edge_noise`

        Returns:
            numpy.array: edge weights with added gaussian edge_noise
        """

        return np.clip(np.random.normal(self.costs, self.params.edge_noise), 0.0001, None)  # None means don't clip above


    def _random_terminals(self):
        """
        Switches the terminals with random nodes with a similar degree.

        Returns:
            numpy.array: new prizes
            numpy.array: new terminals
        """

        nodes_sorted_by_degree = pd.Series(self.node_degrees).sort_values().index
        terminal_degree_rankings = np.array([nodes_sorted_by_degree.get_loc(terminal) for terminal in self.terminals])
        new_terminal_degree_rankings = np.clip(np.rint(np.random.normal(terminal_degree_rankings, 10)), 0, len(self.nodes)-1).astype(int)
        new_terminals = pd.Series(nodes_sorted_by_degree)[new_terminal_degree_rankings].values

        new_prizes = copy(self.prizes)

        for old_terminal, new_terminal in zip(self.terminals, new_terminals):
            new_prizes[old_terminal] = 0
            new_prizes[new_terminal] = self.prizes[old_terminal]

        return new_prizes, np.unique(new_terminals)


    def _noisy_edges_reps(self, reps):
        """
        Perform PCSF and collect results for some number of noisy edges randomizations
        """

        results = []

        true_edge_costs = copy(self.costs)

        for noisy_edge_costs in [self._noisy_edges() for rep in range(reps)]:
            self.costs = noisy_edge_costs
            results.append(self.pcsf())

        self.costs = true_edge_costs

        return self._aggregate_pcsf(results, 'robustness')


    def _random_terminal_reps(self, reps):
        """
        Perform PCSF and collect results for some number of random_terminal randomizations
        """
        results = []

        true_prizes = self.prizes
        true_terminals = self.terminals

        for random_prizes, terminals in [self._random_terminals() for rep in range(reps)]:
            self.prizes = random_prizes
            self.terminals = terminals
            results.append(self.pcsf())

        self.prizes = true_prizes
        self.terminals = true_terminals

        return self._aggregate_pcsf(results, 'specificity')


    def _aggregate_pcsf(self, results, frequency_attribute_name="frequency"):
        """
        Merge multiple PCSF results into one DataFrame

        Arguments:
            results (list): a list of [(vertex_indices, edge_indices),...] from multiple PCSF runs.
            frequency_attribute_name (str): Name of the attribute relating to the frequency of occurrence of components in the results.

        Returns:
            pandas.DataFrame: vertex indices and their fractional rate of occurrence in the PCSF results
            pandas.DataFrame: edge indices and their fractional rate of occurrence in the PCSF results
        """

        if len(results) == 0: return pd.DataFrame(), pd.DataFrame()

        # Transposes a list from [(vertex_indices, edge_indices),...] to ([vertex_indices,...], [edge_indices,...])
        vertex_indices, edge_indices = zip(*results)

        # These next steps are just data transformation/aggregation.
        # 1. Flatten the lists of lists of edge indices and vertex indices
        # 2. Count the occurrences of each edge and vertex index
        # 3. Transform from Counter object to DataFrame through list
        vertex_indices_df = pd.DataFrame(list(Counter(flatten(vertex_indices)).items()), columns=['node_index',frequency_attribute_name]).set_index('node_index')
        edge_indices_df = pd.DataFrame(list(Counter(flatten(edge_indices)).items()), columns=['edge_index',frequency_attribute_name]).set_index('edge_index')
        # 4. Convert occurrences to fractions
        vertex_indices_df[frequency_attribute_name] /= len(results)
        edge_indices_df[frequency_attribute_name] /= len(results)

        return vertex_indices_df, edge_indices_df


    def randomizations(self, noisy_edges_reps=0, random_terminals_reps=0):
        """
        Macro function which performs randomizations and merges the results

        Note that thee parameters are additive, not multiplicative:
        `noisy_edges_reps` = 5 and `random_terminals_reps` = 5 makes 10 PCSF runs, not 25.

        Arguments:
            noisy_edges_reps (int): Number of "Noisy Edges" type randomizations to perform
            random_terminals_reps (int): Number of "Random Terminals" type randomizations to perform

        Returns:
            networkx.Graph: forest
            networkx.Graph: augmented_forest
        """

        if self.params.seed: random.seed(self.params.seed); np.random.seed(seed=self.params.seed)

        # For single PCSF run
        if noisy_edges_reps == random_terminals_reps == 0:
            return self.output_forest_as_networkx(*self.pcsf())

        robust_vertices,   robust_edges   = self._noisy_edges_reps(noisy_edges_reps)
        specific_vertices, specific_edges = self._random_terminal_reps(random_terminals_reps)

        vertex_indices = pd.concat([robust_vertices, specific_vertices], axis=1).fillna(0)
        edge_indices = pd.concat([robust_edges, specific_edges], axis=1).fillna(0)

        forest, augmented_forest = self.output_forest_as_networkx(vertex_indices.index.values, edge_indices.index.values)

        # Skip attribute setting if solution is empty.
        if forest.number_of_nodes() == 0: return forest, augmented_forest

        # reindex `vertex_indices` by name: basically we "dereference" the vertex indices to vertex names
        vertex_indices.index = self.nodes[vertex_indices.index.values]

        nx.set_node_attributes(forest,           vertex_indices.reindex(list(forest.nodes())).dropna(how='all').to_dict(orient='index'))
        nx.set_node_attributes(augmented_forest, vertex_indices.reindex(list(augmented_forest.nodes())).dropna(how='all').to_dict(orient='index'))

        return forest, augmented_forest


    ###########################################################################
                #######          Grid Search          #######
    ###########################################################################

    def _eval_PCSF_runs(self, params):
        """
        Convenience method which sets parameters and performs PCSF randomizations.

        Arguments:
            params (dict): dictionary with regular OI2 parameters _AND_ noisy_edge_reps and random_terminals_reps

        Returns:
            str: Parameter values in string format
            networkx.Graph: forest
            networkx.Graph: augmented_forest
        """

        self._reset_hyperparameters(params=params)
        paramstring = "W_{:04.2f}_B_{:04.2f}_G_{:04.2f}".format(self.params.w, self.params.b, self.params.g)

        if params["noisy_edge_reps"] == params["random_terminals_reps"] == 0: logger.info("Single PCSF run for " + paramstring)
        else: logger.info("Randomizations for " + paramstring)

        forest, augmented_forest = self.randomizations(params["noisy_edge_reps"], params["random_terminals_reps"])

        return paramstring, forest, augmented_forest


    def grid_randomization(self, prize_file, Ws=[5], Bs=[1], Gs=[3], noisy_edges_reps=0, random_terminals_reps=0):
        """
        Macro function which performs grid search or randomizations or both.

        Arguments:
            prize_file (str): filepath
            Gs (list): Values of gamma
            Bs (list): Values of beta
            Ws (list): Values of omega
            noisy_edges_reps (int): Number of robustness experiments
            random_terminals_reps (int): Number of specificity experiments

        Returns:
            dict: Forest and augmented forest networkx graphs, keyed by parameter string
        """

        pool = multiprocessing.Pool(n_cpus)

        self.prepare_prizes(prize_file)

        param_sets = [{'w': w, 'b': b, 'g': g, 'noisy_edge_reps': noisy_edges_reps, 'random_terminals_reps': random_terminals_reps} for (w, b, g) in product(Ws, Bs, Gs)]

        results = pool.map(self._eval_PCSF_runs, param_sets)
        # Convert to dictionary format
        results = {paramstring: {"forest": forest, "augmented_forest": augmented_forest} for paramstring, forest, augmented_forest in results }

        return results


    def grid_search(self, prize_file, Ws, Bs, Gs):
        """
        Macro function which performs grid search.

        Arguments:
            prize_file (str): filepath
            Gs (list): Values of gamma
            Bs (list): Values of beta
            Ws (list): Values of omega

        Returns:
            networkx.Graph: forest
            networkx.Graph: augmented_forest
            pd.DataFrame: parameters and node membership lists
        """

        return self.grid_randomization(prize_file, Ws, Bs, Gs, 0, 0)



###############################################################################
            #######          Subgraph Augmentation        #######
###############################################################################

def betweenness(nxgraph):
    """
    Compute and add as an attribute the betweenness of each node.

    Betweenness centrality of a node v is the sum of the fraction of all-pairs shortest paths that pass through v.

    Arguments:
        nxgraph (networkx.Graph): a networkx graph, usually the augmented_forest.
    """
    nx.set_node_attributes(nxgraph, {node: {'betweenness':betweenness} for node,betweenness in nx.betweenness_centrality(nxgraph).items()})


def louvain_clustering(nxgraph):
    """
    Compute "Louvain"/"Community" clustering on a networkx graph, and add the cluster labels as attributes on the nodes.


    Arguments:
        nxgraph (networkx.Graph): a networkx graph, usually the augmented_forest.
    """
    nx.set_node_attributes(nxgraph, {node: {'louvain_clusters':str(cluster)} for node,cluster in community.best_partition(nxgraph).items()})


def k_clique_clustering(nxgraph, k):
    """
    Compute "k-Clique" clustering on a networkx graph, and add the cluster labels as attributes on the nodes.

    See the [networkx docs](https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.community.kclique.k_clique_communities.html#networkx.algorithms.community.kclique.k_clique_communities)

    Arguments:
        nxgraph (networkx.Graph): a networkx graph, usually the augmented_forest.
    """

    if k < 2: logger.critical("K-Clique Clustering requires that k be an integer larger than 1."); raise ValueError("Improper input to k_clique_clustering")

    clustering = pd.Series(invert(nx.algorithms.community.kclique.k_clique_communities(nxgraph, k)), name='k_clique_clusters').astype(str).reindex(nxgraph.nodes())
    nx.set_node_attributes(nxgraph, clustering.to_frame().to_dict(orient='index'))


def spectral_clustering(nxgraph, k):
    """
    Compute "spectral" clustering on a networkx graph, and add the cluster labels as attributes on the nodes.


    Arguments:
        nxgraph (networkx.Graph): a networkx graph, usually the augmented_forest.
    """
    adj_matrix = nx.to_pandas_adjacency(nxgraph)
    clustering =  SpectralClustering(k, affinity='precomputed', n_init=100, assign_labels='discretize').fit_predict(adj_matrix.values)
    nx.set_node_attributes(nxgraph, {node: {'spectral_clusters':str(cluster)} for node,cluster in zip(adj_matrix.index, clustering)})


def annotate_graph_nodes(nxgraph):
    """
    Arguments:
        nxgraph (networkx.Graph): a networkx graph, usually the augmented_forest.
    """

    try:
        annotation = pd.read_pickle(get_path('OmicsIntegrator', 'annotation/final_annotation.pickle'))
    except:
        annotation = pd.read_pickle(Path.cwd().parent / 'annotation' / 'final_annotation.pickle')

    nx.set_node_attributes(nxgraph, annotation.reindex(list(nxgraph.nodes())).dropna(how='all').to_dict(orient='index'))


###############################################################################
            #######            Results           #######
###############################################################################

def summarize_grid_search(results, mode, top_n=np.Infinity):
    """
    Summarizes results of `grid_randomization` or `grid_search` into a matrix where each row is a gene
    and each column is a parameter run. If summarizing "membership", entries will be 0 or 1
    indicating whether or not a node appeared in each experiment. If summarizing "robustness"
    or "specificity", entries indicate robustness or specificity values for each experiment.

    Arguments:
        results (list of tuples): Results of `grid_randomization` or `grid_search` of form `{'paramstring': { 'forest': object, 'augmented forest': object}}`
        mode (str): Reported values "membership", "robustness", "specificity"
        top_n (int): Takes the top_n values of the summary dataframe. top_n=-1 sets no threshold

    Returns:
        pd.DataFrame: Columns correspond to each parameter experiment, indexed by nodes
    """

    # Exclude any degenerate results
    results = {paramstring: graphs for paramstring, graphs in results.items() if graphs["augmented_forest"].number_of_nodes() > 0}

    if mode == "membership": # Summarize single-run parameter search
        series = [pd.Series(1, index=graphs["augmented_forest"].nodes(), name=paramstring) for paramstring, graphs in results.items()]
    elif mode == "robustness": # Summarize randomized robustness
        series = [get_networkx_graph_as_dataframe_of_nodes(graphs["augmented_forest"])["robustness"].rename(paramstring) for paramstring, graphs in results.items()]
    elif mode == "specificity": # Summarize randomized specificity
        series = [get_networkx_graph_as_dataframe_of_nodes(graphs["augmented_forest"])["specificity"].rename(paramstring) for paramstring, graphs in results.items()]
    else:
        logger.critical("`mode` must be one of the following: 'membership', 'robustness', or 'specificity'."); raise ValueError("Improper input to summarize_grid_search")

    node_summary_df = pd.concat(series, axis=1).fillna(0)

    # df can get quite large with many sparse entries, so let's filter for the top_n entries
    if len(node_summary_df) > top_n:
        cutoff = sorted(node_summary_df.sum(axis=1).tolist(), reverse=True)[top_n]
        node_summary_df = node_summary_df[node_summary_df.sum(axis=1) > cutoff]

    return node_summary_df


def get_robust_subgraph_from_randomizations(nxgraph, max_size=400, min_component_size=5):
    """
    Given a graph with robustness attributes, take the top `max_size` robust nodes and
    prune any "small" components.

    Arguments:
        nxgraph (networkx.Graph): Network from randomization experiment
        max_size (int): Max size of robust network

    Returns:
        networkx.Graph: Robust network
    """

    # TODO: Potential alternative approach - from entire network, attempt to remove lowest robustness node.
    # If removal results in a component of size less than min_size, do not remove.

    if nxgraph.number_of_nodes() == 0:
        logger.warning("Augmented forest is empty.")
        return nxgraph

    # Get indices of top nodes sorted by high robustness, then low specificity. Don't allow nodes with robustness = 0.
    node_attributes_df = get_networkx_graph_as_dataframe_of_nodes(nxgraph)
    node_attributes_df = node_attributes_df[node_attributes_df["robustness"] > 0]
    node_attributes_df.sort_values(["robustness", "specificity"], ascending=[False, True], inplace=True)
    top_hits = node_attributes_df.index[:min(max_size,len(node_attributes_df))]
    # Get robust subnetwork and remove small components.
    robust_network = nxgraph.subgraph(top_hits)
    robust_network = filter_graph_by_component_size(robust_network, min_component_size)

    return robust_network


def filter_graph_by_component_size(nxgraph, min_size=5):
    """
    Removes any components that are less than `min_size`.

    Arguments:
        nxgraph (networkx.Graph): Network from randomization experiment
        min_size (int): Min size of components in `nxgraph`. Set to 2 to remove singletons only.

    Returns:
        networkx.Graph: Network with components less than specified size removed.
    """

    filtered_subgraph = nxgraph.copy()

    small_components = [g.nodes() for g in nx.connected_component_subgraphs(nxgraph, copy=False) if g.number_of_nodes() < min_size]
    filtered_subgraph.remove_nodes_from(flatten(small_components))

    return filtered_subgraph


###############################################################################
            #######              Export             #######
###############################################################################

def get_networkx_graph_as_dataframe_of_nodes(nxgraph):
    """
    Arguments:
        nxgraph (networkx.Graph): any instance of networkx.Graph
    Returns:
        pd.DataFrame: nodes from the input graph and their attributes as a dataframe
    """

    return pd.DataFrame.from_dict(dict(nxgraph.nodes(data=True)), orient='index')


def get_networkx_graph_as_dataframe_of_edges(nxgraph):
    """
    Arguments:
        nxgraph (networkx.Graph): any instance of networkx.Graph
    Returns:
        pd.DataFrame: edges from the input graph and their attributes as a dataframe
    """

    return nx.to_pandas_edgelist(nxgraph, 'protein1', 'protein2')


def output_networkx_graph_as_pickle(nxgraph, output_dir=".", filename="pcsf_results.pickle"):
    """
    Arguments:
        nxgraph (networkx.Graph): any instance of networkx.Graph
        output_dir (str): the directory in which to output the graph.
        filename (str): Filenames ending in .gz or .bz2 will be compressed.
    Returns:
        Path: the filepath which was outputted to
    """

    path = Path(output_dir)
    path.mkdir(exist_ok=True, parents=True)
    path = path / filename
    nx.write_gpickle(nxgraph, open(path, 'wb'))

    return path.resolve()


def output_networkx_graph_as_graphml_for_cytoscape(nxgraph, output_dir=".", filename="pcsf_results.graphml.gz"):
    """
    Arguments:
        nxgraph (networkx.Graph): any instance of networkx.Graph
        output_dir (str): the directory in which to output the graph.
        filename (str): Filenames ending in .gz or .bz2 will be compressed.
    Returns:
        Path: the filepath which was outputted to
    """
    path = Path(output_dir)
    path.mkdir(exist_ok=True, parents=True)
    path = path / filename
    nx.write_graphml(nxgraph, path)

    return path.resolve()


def output_networkx_graph_as_interactive_html(nxgraph, attribute_metadata=dict(), output_dir=".", filename="graph.html"):
    """
    Arguments:
        nxgraph (networkx.Graph): any instance of networkx.Graph
        output_dir (str): the directory in which to output the file
        filename (str): the filename of the output file
    Returns:
        Path: the filepath which was outputted to
    """
    return axial.graph(nxgraph,
        title='OmicsIntegrator2 Results',
        scripts_mode="inline",
        data_mode="inline",
        output_dir=output_dir,
        filename=filename)

