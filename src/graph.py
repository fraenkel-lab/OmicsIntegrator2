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
from copy import copy
import json
from pkg_resources import resource_filename as get_path

# python external libraries
import numpy as np
import pandas as pd
import networkx as nx
from networkx.readwrite import json_graph
import community    # pip install python-louvain
from sklearn.cluster import SpectralClustering
import jinja2

# Lab modules
from pcst_fast import pcst_fast

# list of classes and methods we'd like to export:
__all__ = [ "Graph",
            "output_networkx_graph_as_graphml_for_cytoscape",
            "output_networkx_graph_as_json_for_cytoscapejs",
            "output_networkx_graph_as_interactive_html",
            "get_networkx_graph_as_dataframe_of_nodes",
            "get_networkx_graph_as_dataframe_of_edges",
            "summarize_grid_search", 
            "get_robust_subgraph_from_randomizations" ]

templateLoader = jinja2.FileSystemLoader(os.path.dirname(os.path.abspath(__file__)))
templateEnv = jinja2.Environment(loader=templateLoader)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
handler.setFormatter(logging.Formatter('%(asctime)s - Graph: %(levelname)s - %(message)s', "%I:%M:%S"))
logger.addHandler(handler)

# Count the number of available CPUs for potential use in multiprocessing code
try: n_cpus = int(os.environ["SLURM_JOB_CPUS_PER_NODE"])
except KeyError: n_cpus = multiprocessing.cpu_count()

# Helpers
def flatten(list_of_lists): return [item for sublist in list_of_lists for item in sublist]

def invert(list_of_lists): return {item: i for i, list in enumerate(list_of_lists) for item in list}

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

    def __init__(self, interactome_file, params={}, skip_checks=False):
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
        self.interactome_dataframe.columns = ["source","target","cost"] + self.interactome_dataframe.columns[3:].tolist()  # TODO: error handling

        if not skip_checks:
            # Handle the case of possible duplicate edges.
            # Do so by creating a string column of both interactors, e.g. "ABCD1EFGR" and remove duplicates
            # This operation is time consuming, especially if there do exist duplicates. TODO: optimize this code.
            self.interactome_dataframe['temp'] = self.interactome_dataframe.apply(lambda row: ''.join(sorted([row['source'], row['target']])), axis=1)
            duplicated_edges = self.interactome_dataframe[self.interactome_dataframe.set_index('temp').index.duplicated()][['source','target']].values.tolist()
            logger.info("Duplicated edges in the interactome file (we'll keep the max cost):")
            logger.info(duplicated_edges)
            if len(duplicated_edges) > 0:
                self.interactome_dataframe = self.interactome_dataframe.groupby('temp').max().reset_index()[["source","target","cost"]]
            else: del self.interactome_dataframe['temp']

        self.interactome_graph = nx.from_pandas_edgelist(self.interactome_dataframe, 'source', 'target', edge_attr=self.interactome_dataframe.columns[2:].tolist())

        # Convert the interactome dataframe from string interactor IDs to integer interactor IDs.
        # Do so by selecting the source and target columns from the interactome dataframe,
        # then unstacking them, which (unintuitively) stacks them into one column, allowing us to use factorize.
        # Factorize builds two datastructures, a unique pd.Index which maps each ID string to an integer ID,
        # and the datastructure we passed in with string IDs replaced with those integer IDs.
        # We place those in self.nodes and self.edges respectively, but self.edges will need reshaping.
        (self.edges, self.nodes) = pd.factorize(self.interactome_dataframe[["source","target"]].unstack())
        # Here we do the inverse operation of "unstack" above, which gives us an interpretable edges datastructure
        self.edges = self.edges.reshape(self.interactome_dataframe[["source","target"]].shape, order='F')

        self.edge_costs = self.interactome_dataframe['cost'].astype(float).values

        # Count the number of incident edges into each node.
        # The indices into this datastructure are the same as those in self.nodes which are the IDs in self.edges.
        self.node_degrees = np.bincount(self.edges.flatten())

        # The rest of the setup work is occasionally repeated, so use another method to complete setup.
        self._reset_hyperparameters(params=params)


    def _reset_hyperparameters(self, params={}):
        """
        Set the parameters on Graph and compute parameter-dependent features.

        Arguments:
            params (dict): params with which to run the program
        """

        defaults = {"w": 6, "b": 1, "g": 1000, "noise": 0.1, "exclude_terminals": False, "dummy_mode": "terminals", "knockout": [], "seed": None}

        # Overwrite the defaults with any user-specified parameters.
        self.params = Options({**defaults, **params})
        # Knockout any proteins from the interactome
        self._knockout(self.params.knockout)
        # Add costs to each edge, proportional to the degrees of the nodes it connects, modulated by parameter g.
        N = len(self.nodes)
        self.edge_penalties = self.params.g * np.array([self.node_degrees[a] * self.node_degrees[b] /
                            ((N - self.node_degrees[a] - 1) * (N - self.node_degrees[b] - 1) + self.node_degrees[a] * self.node_degrees[b]) for a, b in self.edges])

        self.costs = (self.edge_costs + self.edge_penalties)

        # In case specific nodes are penalized, we'd like to update the costs accordingly
        if hasattr(self, "additional_costs"): self.costs = (self.edge_costs + self.edge_penalties + self.additional_costs)

        # if this instance of graph has bare_prizes set, then presumably resetting the
        # hyperparameters should also reset the scaled prizes
        if hasattr(self, "bare_prizes"): self.prizes = self.bare_prizes * self.params.b


    def _penalize_nodes(self, nodes_to_penalize):
        """
        Penalize a set of nodes by penalizing the edges connected to that node by some coefficient.

        Arguments:
            nodes_to_penalize (pandas.DataFrame): 2 columns: 'name' and 'penalty' with entries in [0, 1)
        """

        # Here's we're indexing the penalized nodes by the indices we used for nodes during initialization
        nodes_to_penalize.index = graph.nodes.get_indexer(nodes_to_penalize['name'].values)

        # there will be some nodes in the penalty dataframe which we don't have in our interactome
        logger.info("Members of the penalty dataframe not present in the interactome (we'll need to drop these):")
        logger.info(nodes_to_penalize[nodes_to_penalize.index == -1]['name'].tolist())
        nodes_to_penalize.drop(-1, inplace=True, errors='ignore')

        if not nodes_to_penalize['penalty_coefficient'].between(0, 1).all():
            logger.info("The node penalty coefficients must lie in [0, 1]. Skipping penalization..."); return

        nodes_to_knockout = nodes_to_penalize[nodes_to_penalize.penalty_coefficient == 1]
        logger.info("penalty coefficients of 1 are treated as knockouts. Proteins to knock out from interactome:")
        logger.info(nodes_to_knockout.name.tolist())
        self._knockout(nodes_to_knockout.name.values)
        nodes_to_penalize= nodes_to_penalize[nodes_to_penalize.penalty_coefficient < 1]

        self.additional_costs = np.zeros(self.costs.shape)

        # Iterate through the rows of the nodes_to_penalize dataframe
        for index, name, penalty_coefficient in nodes_to_penalize.itertuples():
            # For each protein we'd like to penalize, get the indicies of the edges connected to that node
            edge_indices = self.interactome_dataframe[(self.interactome_dataframe.source == name) | (self.interactome_dataframe.target == name)].index
            # And compute an additional cost on those edges.
            self.additional_costs[edge_indices] += self.edge_costs[edge_indices] / (1 - penalty_coefficient)
        # Apply those additional costs by calling _reset_hyperparameters.
        self._reset_hyperparameters()


    def _knockout(self, nodes_to_knockout):
        """
        Knock out a set of nodes from the interactome, effectively removing them from results.

        Arguments:
            nodes_to_knockout (numpy.array): Array of string IDs of nodes to knock out.
        """
        if len(nodes_to_knockout) > 0:
            logger.info("The knockout function has yet to be implemented, passing...");

        return


    def prepare_prizes(self, prize_file):
        """
        Parses a prize file and adds prize-related attributes to the graph object.

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
        prizes_dataframe.columns = ['name', 'prize'] + prizes_dataframe.columns[2:].tolist()  # TODO: error handling

        if "type" not in prizes_dataframe.columns: prizes_dataframe["type"] = "terminal"

        return self._prepare_prizes(prizes_dataframe)


    def _prepare_prizes(self, prizes_dataframe):

        # Some files have duplicated genes, sometimes with different prizes. Keep the max prize.
        logger.info("Duplicated gene symbols in the prize file (we'll keep the max prize):")
        logger.info(prizes_dataframe[prizes_dataframe.set_index('name').index.duplicated()]['name'].tolist())
        prizes_dataframe = prizes_dataframe.groupby('name').max().reset_index()

        # Here's we're indexing the terminal nodes and associated prizes by the indices we used for nodes
        prizes_dataframe.set_index(self.nodes.get_indexer(prizes_dataframe['name']), inplace=True)

        # there will be some nodes in the prize file which we don't have in our interactome
        logger.info("Members of the prize file not present in the interactome:")
        logger.info(prizes_dataframe[prizes_dataframe.index == -1]['name'].tolist())
        prizes_dataframe.drop(-1, inplace=True, errors='ignore')

        # Node attributes dataframe for all proteins in self.nodes. Default for non-prized nodes is Steiner.
        self.node_attributes = prizes_dataframe.set_index('name').rename_axis(None).reindex(self.nodes)
        self.node_attributes["degree"] = self.node_degrees
        self.node_attributes["prize"].fillna(0, inplace=True)
        self.node_attributes["type"].fillna("steiner", inplace=True)

        # Here we're making a dataframe with all the nodes as keys and the prizes from above or 0
        prizes_dataframe = pd.DataFrame(self.nodes, columns=["name"]).merge(prizes_dataframe, on="name", how="left").fillna(0)
        # Our return value is a 1D array, where each entry is a node's prize, indexed as above
        self.bare_prizes = prizes_dataframe['prize'].values
        self.prizes = self.bare_prizes * self.params.b

        self.terminals = pd.Series(self.prizes).nonzero()[0].tolist()


    ###########################################################################
                #######              PCSF               #######
    ###########################################################################

    def _add_dummy_node(self, connected_to=[]):

        dummy_id = len(self.nodes)
        dummy_prize = np.array([0])
        dummy_edges = np.array([(dummy_id, node_id) for node_id in connected_to])
        dummy_costs = np.array([self.params.w] * len(dummy_edges))

        return dummy_edges, dummy_costs, dummy_id, dummy_prize


    def _check_validity_of_instance(self, edges, prizes, costs):
        """
        Assert that the parammeters and files passed to this program are valid, log useful error messages otherwise.
        """
        pass


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
        else: sys.exit("Invalid dummy mode")

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

        self._check_validity_of_instance(edges, prizes, costs)

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

        if np.array(vertex_indices).dtype == "O": 
            logger.warning("`vertex_indices` is of dtype Object. This may be caused by an empty dataframe with default Object types.")

        # Replace the edge indices with the actual edges (source name, target name) by indexing into the interactome
        edges = self.interactome_dataframe.loc[edge_indices]
        forest = nx.from_pandas_edgelist(edges, 'source', 'target', edge_attr=True)
        # the above won't capture the singletons, so we'll add them here 
        forest.add_nodes_from(list(set(self.nodes[vertex_indices]) - set(forest.nodes())))

        # Set all the attributes on graph
        nx.set_node_attributes(forest, self.node_attributes.loc[list(forest.nodes())].dropna(how='all').to_dict(orient='index'))
        # Set a flag on all the edges which were selected by PCSF (before augmenting the forest)
        nx.set_edge_attributes(forest, True, name='in_solution')
        # Create a new graph including all edges between all selected nodes, not just those edges selected by PCSF.
        augmented_forest = nx.compose(self.interactome_graph.subgraph(forest.nodes()), forest)

        # Post-processing
        betweenness(augmented_forest)
        louvain_clustering(augmented_forest)
        augment_with_subcellular_localization(augmented_forest)

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
        Adds gaussian noise to all edge costs in the graph, modulated by parameter `noise`

        Returns:
            numpy.array: edge weights with added gaussian noise
        """

        return np.clip(np.random.normal(self.costs, self.params.noise), 0.0001, None)  # None means don't clip above


    def _random_terminals(self):
        """
        Switches the terminams with random nodes with a similar degree.

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


    def _aggregate_pcsf(self, results, frequency_attribute_name):
        """
        Merge multiple PCSF results into one DataFrame

        Arguments:
            results (list): a list of [(vertex_indices, edge_indices),...] from multiple PCSF runs.
            frequency_attribute_name (str): Name of the attribute relating to the frequency of occurrence of components in the results.

        Returns:
            pandas.DataFrame: vertex indices and their fractional rate of occurrence in the PCSF results
            pandas.DataFrame: edge indices and their fractional rate of occurrence in the PCSF results
        """

        # Transposes a list from [(vertex_indices, edge_indices),...] to ([vertex_indices,...], [edge_indices,...])
        vertex_indices, edge_indices = zip(*results)

        # These next steps are just data transformation/aggregation.
        # 1. Flatten the lists of lists of edge indices and vertex indices
        # 2. Count the occurrences of each edge and vertex index
        # 3. Transform from Counter object to DataFrame through list
        vertex_indices_df = pd.DataFrame(list(Counter(flatten(vertex_indices)).items()), columns=['node_index',frequency_attribute_name])
        edge_indices_df = pd.DataFrame(list(Counter(flatten(edge_indices)).items()), columns=['edge_index',frequency_attribute_name])
        # 4. Convert occurrences to fractions
        vertex_indices_df[frequency_attribute_name] /= len(results)
        edge_indices_df[frequency_attribute_name] /= len(results)

        return vertex_indices_df, edge_indices_df


    def randomizations(self, noisy_edges_reps, random_terminals_reps):
        """
        Macro function which performs randomizations and merges the results

        Note that these are additive, not multiplicative:
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

        #### NOISY EDGES ####
        if noisy_edges_reps > 0:

            results = []

            true_edge_costs = copy(self.costs)

            for noisy_edge_costs in [self._noisy_edges() for rep in range(noisy_edges_reps)]:
                self.costs = noisy_edge_costs
                results.append(self.pcsf())

            self.costs = true_edge_costs

            robust_vertices, robust_edges = self._aggregate_pcsf(results, 'robustness')

        #### RANDOM TERMINALS ####
        if random_terminals_reps > 0:

            results = []

            true_prizes = self.prizes
            true_terminals = self.terminals

            for random_prizes, terminals in [self._random_terminals() for rep in range(random_terminals_reps)]:
                self.prizes = random_prizes
                self.terminals = terminals
                results.append(self.pcsf())

            self.prizes = true_prizes
            self.terminals = true_terminals

            specific_vertices, specific_edges = self._aggregate_pcsf(results, 'specificity')

        ###########

        if random_terminals_reps == 0 and noisy_edges_reps > 0:
            vertex_indices = robust_vertices; edge_indices = robust_edges;

        elif noisy_edges_reps == 0 and random_terminals_reps > 0:
            vertex_indices = specific_vertices; edge_indices = specific_edges;

        elif noisy_edges_reps > 0 and random_terminals_reps > 0:
            vertex_indices = robust_vertices.merge(specific_vertices, how='outer', on='node_index').fillna(0)
            edge_indices = robust_edges.merge(specific_edges, how='outer', on='edge_index').fillna(0)

        else: sys.exit("Randomizations was called with invalid noisy_edges_reps and random_terminals_reps.")

        ###########

        forest, augmented_forest = self.output_forest_as_networkx(vertex_indices.node_index.values.astype(int), edge_indices.edge_index.values.astype(int))

        # Skip attribute setting if solution is empty. 
        if forest.number_of_nodes() == 0: return forest, augmented_forest

        # reindex `vertex_indices_df` by name: basically we "dereference" the vertex indices to vertex names
        vertex_indices.index = self.nodes[vertex_indices.node_index.values]

        nx.set_node_attributes(forest,           vertex_indices.loc[list(forest.nodes())].dropna(how='all').to_dict(orient='index'))
        nx.set_node_attributes(augmented_forest, vertex_indices.loc[list(augmented_forest.nodes())].dropna(how='all').to_dict(orient='index'))

        return forest, augmented_forest


    ###########################################################################
                #######          Grid Search          #######
    ###########################################################################

    def _eval_PCSF_runs(self, params):
        """
        Convenience method which sets parameters and performs PCSF randomizations.

        Arguments: 
            params (dict): params with which to run the program

        Returns: 
            str: Parameter values in string format 
            networkx.Graph: forest
            networkx.Graph: augmented_forest
        """

        self._reset_hyperparameters(params=params)
        paramstring = "W_{:04.2f}_B_{:04.2f}_G_{:d}".format(self.params.w, self.params.b, self.params.g)

        if params["noisy_edge_reps"] + params["random_terminals_reps"] == 0:
            logger.info("Single PCSF run for " + paramstring)
        else: 
            logger.info("Randomizations for " + paramstring)

        forest, augmented_forest = self.randomizations(params["noisy_edge_reps"], params["random_terminals_reps"])
        
        return paramstring, forest, augmented_forest


    def grid_randomization(self, prize_file, Ws, Bs, Gs, noisy_edges_reps, random_terminals_reps):
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

        # get number of cpus available to job
        try:
            n_cpus = int(os.environ["SLURM_JOB_CPUS_PER_NODE"])
        except KeyError:
            n_cpus = multiprocessing.cpu_count()

        pool = multiprocessing.Pool(n_cpus)

        self.prepare_prizes(prize_file)

        model_params = [{'w': w, 'b': b, 'g': g, 'noisy_edge_reps': noisy_edges_reps, 'random_terminals_reps': random_terminals_reps} for (w, b, g) in product(Ws, Bs, Gs)]

        results = pool.map(self._eval_PCSF_runs, model_params)
        # Convert to dictionary format
        results = { paramstring: {"forest": forest, "augmented_forest": augmented_forest} for paramstring, forest, augmented_forest in results }

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
    nx.set_node_attributes(nxgraph, {node: {'louvainClusters':str(cluster)} for node,cluster in community.best_partition(nxgraph).items()})

def edge_betweenness_clustering(nxgraph):
    """
    Compute "Edge-betweenness"/"Girvan-Newman" clustering on a networkx graph, and add the cluster labels as attributes on the nodes.

    The Girvan–Newman algorithm detects communities by progressively removing edges from the original graph.
    The algorithm removes the “most valuable” edge, traditionally the edge with the highest betweenness centrality, at each step.
    As the graph breaks down into pieces, the tightly knit community structure is exposed and the result can be depicted as a dendrogram.

    TODO: currently, we're removing a single edge, which has no effect, so this isn't a real clustering method yet.

    Arguments:
        nxgraph (networkx.Graph): a networkx graph, usually the augmented_forest.
    """
    nx.set_node_attributes(nxgraph, {node: {'edgeBetweennessClusters':str(cluster)} for node,cluster in invert(next(nx.algorithms.community.centrality.girvan_newman(nxgraph)))})

def k_clique_clustering(nxgraph, k):
    """
    Compute "k-Clique" clustering on a networkx graph, and add the cluster labels as attributes on the nodes.


    Arguments:
        nxgraph (networkx.Graph): a networkx graph, usually the augmented_forest.
    """
    nx.set_node_attributes(nxgraph, {node: {'kCliqueClusters':str(cluster)} for node,cluster in invert(nx.algorithms.community.kclique.k_clique_communities(nxgraph, k)).items()})

def spectral_clustering(nxgraph, k):
    """
    Compute "spectral" clustering on a networkx graph, and add the cluster labels as attributes on the nodes.


    Arguments:
        nxgraph (networkx.Graph): a networkx graph, usually the augmented_forest.
    """
    adj_matrix = nx.to_pandas_adjacency(nxgraph)
    clustering =  SpectralClustering(k, affinity='precomputed', n_init=100, assign_labels='discretize').fit_predict(adj_matrix.values)
    nx.set_node_attributes(nxgraph, {node: {'spectral_clusters':str(cluster)} for node,cluster in zip(adj_matrix.index, clustering)})


###############################################################################
            #######            GO Enrichment          #######
###############################################################################

def augment_with_all_GO_terms(nxgraph):
    """
    Arguments:
        nxgraph (networkx.Graph): a networkx graph, usually the augmented_forest.
    """
    augment_with_subcellular_localization(nxgraph)
    augment_with_biological_process_terms(nxgraph)
    augment_with_molecular_function_terms(nxgraph)


def augment_with_subcellular_localization(nxgraph):
    """
    Arguments:
        nxgraph (networkx.Graph): a networkx graph, usually the augmented_forest.
    """

    try:
        subcellular = pd.read_pickle(get_path('OmicsIntegrator', 'subcellular_compartments/subcellular.pickle'))
    except:
        # maybe need os.path.realpath(__file__)
        subcellular = pd.read_pickle(os.path.dirname(os.path.realpath(__file__))+'/../subcellular/subcellular.pickle')

    try:
        nx.set_node_attributes(nxgraph, subcellular.loc[list(nxgraph.nodes())].dropna(how='all').to_dict(orient='index'))
    except KeyError:
        logger.info('Not assigning subcellular locations. For that function use Gene Symbols.')

def augment_with_biological_process_terms(nxgraph):
    """
    Arguments:
        nxgraph (networkx.Graph): a networkx graph, usually the augmented_forest.
    """
    pass

def augment_with_molecular_function_terms(nxgraph):
    """
    Arguments:
        nxgraph (networkx.Graph): a networkx graph, usually the augmented_forest.
    """
    pass

def perform_GO_enrichment_on_clusters(nxgraph, clustering):
    """
    Arguments:
        nxgraph (networkx.Graph): a networkx graph, usually the augmented_forest.
        clustering (str): the column name of the clustering to perform GO enrichment with.
    """
    pass


###############################################################################
            #######            Results           #######
###############################################################################

def summarize_grid_search(results, mode, top_n=False): 
    """
    Summarizes results of `_grid_randomization` into a matrix where each row is a gene 
    and each column is a parameter run. If summarizing "membership", entries will be 0 or 1
    indicating whether or not a node appeared in each experiment. If summarizing "robustness"
    or "specificity", entries indicate robustness or specificity values for each experiment. 
    Also, node attributes columns are added for plotting purposes.

    Arguments:
        results (list of tuples): Results of `_grid_randomization` (paramstring, forest, augmented forest)
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
        logger.info("`mode` must be one of the following: 'membership', 'robustness', or 'specificity'.")
        return
    
    df = pd.concat(series, axis=1).fillna(0)
    
    # df can get quite large with many sparse entries, so let's filter for the top_n entries
    if not top_n: return df
    
    if len(df) > top_n: 
        cutoff = sorted(df.sum(axis=1).tolist(), reverse=True)[top_n]
        df = df[df.sum(axis=1) > cutoff]
    
    return df


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
        logger.info("Augmented forest is empty.")
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

    return pd.DataFrame.from_dict(dict(nxgraph.nodes(data=True))).transpose()


def get_networkx_graph_as_dataframe_of_edges(nxgraph):
    """
    Arguments:
        nxgraph (networkx.Graph): any instance of networkx.Graph
    Returns:
        pd.DataFrame: edges from the input graph and their attributes as a dataframe
    """

    intermediate = pd.DataFrame(nxgraph.edges(data=True))
    intermediate.columns = ['protein1', 'protein2'] + intermediate.columns[2:].tolist()
    # TODO: in the future, get the other attributes out into columns
    return intermediate[['protein1', 'protein2']]


def output_networkx_graph_as_pickle(nxgraph, output_dir, filename="pcsf_results.pickle"):
    """
    Arguments:
        nxgraph (networkx.Graph): any instance of networkx.Graph
        output_dir (str): the directory in which to output the graph.
        filename (str): Filenames ending in .gz or .bz2 will be compressed.
    Returns:
        str: filepath to output
    """
    os.makedirs(os.path.abspath(output_dir), exist_ok=True)
    path = os.path.join(os.path.abspath(output_dir), filename)
    nx.write_gpickle(nxgraph, path)

    return path


def output_networkx_graph_as_graphml_for_cytoscape(nxgraph, output_dir, filename="pscf_results.graphml.gz"):
    """
    Arguments:
        nxgraph (networkx.Graph): any instance of networkx.Graph
        output_dir (str): the directory in which to output the graph.
        filename (str): Filenames ending in .gz or .bz2 will be compressed.
    Returns:
        str: filepath to output
    """
    os.makedirs(os.path.abspath(output_dir), exist_ok=True)
    path = os.path.join(os.path.abspath(output_dir), filename)
    nx.write_graphml(nxgraph, path)

    return path

class Encoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.int64):
            return str(obj)
        return json.JSONEncoder.default(self, obj)

def output_networkx_graph_as_json_for_cytoscapejs(nxgraph, output_dir, filename="graph_json.json"):
    """
    Arguments:
        nxgraph (networkx.Graph): any instance of networkx.Graph
        output_dir (str): the directory in which to output the file (named graph_json.json)
    Returns:
        str: filepath to output
    """
    from py2cytoscape import util as cy

    njs = cy.from_networkx(nxgraph)
    njs["data"]["name"] = filename.replace(".json", "")

    os.makedirs(os.path.abspath(output_dir), exist_ok=True)
    path = os.path.join(os.path.abspath(output_dir), filename)
    with open(path,'w') as output_file:
        output_file.write(json.dumps(njs, cls=Encoder))

    return path


def output_networkx_graph_as_interactive_html(nxgraph, output_dir, filename="graph.html"):
    """
    Arguments:
        nxgraph (networkx.Graph): any instance of networkx.Graph
        output_dir (str): the directory in which to output the file
    Returns:
        str: filepath to output
    """

    try:
        vizjinja = get_path('OmicsIntegrator', 'viz.jinja')
    except:
        # maybe need os.path.realpath(__file__)
        print('except')
        vizjinja = os.path.dirname(os.path.realpath(__file__)) + '/viz.jinja'

    graph_json = json_graph.node_link_data(nxgraph, attrs=dict(source='source_name', target='target_name', name='id', key='key', link='links'))

    def indexOf(node_id): return [i for (i,node) in enumerate(graph_json['nodes']) if node['id'] == node_id][0]
    graph_json["links"] = [{**link, **{"source":indexOf(link['source_name']), "target":indexOf(link['target_name'])}} for link in graph_json["links"]]
    graph_json = json.dumps(graph_json, cls=Encoder)

    nodes = nxgraph.nodes()

    os.makedirs(os.path.abspath(output_dir), exist_ok=True)
    path = os.path.join(os.path.abspath(output_dir), filename)

    if len(nodes) > 0:
        numerical_node_attributes = list(set().union(*[[attribute_key for attribute_key,attribute_value in node[1].items() if isinstance(attribute_value, numbers.Number)] for node in nxgraph.node(data=True)]))
        #print(numerical_node_attributes)
        non_numerical_node_attributes = list(set().union(*[[attribute_key for attribute_key,attribute_value in node[1].items() if attribute_key not in numerical_node_attributes] for node in nxgraph.node(data=True)]))
        min_max = lambda l: (min(l),max(l))
        numerical_node_attributes = {attribute: min_max(nx.get_node_attributes(nxgraph, attribute).values()) for attribute in numerical_node_attributes}
        #print(numerical_node_attributes)
        html_output = templateEnv.get_template('viz.jinja').render(graph_json=graph_json, nodes=nodes, numerical_node_attributes=numerical_node_attributes, non_numerical_node_attributes=non_numerical_node_attributes)
        with open(path,'w') as output_file:
            output_file.write(html_output)

    return path

def output_networkx_graph_as_edgelist(nxgraph, output_dir, filename="graph_json.json"):
    """
    Arguments:
        nxgraph (networkx.Graph): any instance of networkx.Graph
        output_dir (str): the directory in which to output the file (named graph_edgelist.txt)
    Returns:
        str: filepath to output
    """
    os.makedirs(os.path.abspath(output_dir), exist_ok=True)
    path = os.path.join(os.path.abspath(output_dir), filename)
    nx.write_edgelist(nxgraph, path, data=False)

    return path
