import numpy as np
import pandas as pd
import networkx as nx
from sklearn.preprocessing import LabelEncoder

from gslr import gslr

from graph import Graph
from graph import betweenness, louvain_clustering, augment_with_subcellular_localization

def gslr_macro(graph, dataset, sparsity_low=100, sparsity_high=200, num_steps=25, pruning="strong", verbosity_level=1):
    """
    Graph-Sparse Logistic Regression

    Arguments:
        graph (OmicsIntegrator.Graph): an OmicsIntegrator Graph object
        dataset (pd.DataFrame): the n x d data matrix (n examples in dimension d). The column headers are geneSymbols and the index is the class label.
        sparsity_low (int): the (approximate) lower bound for the output sparsity
        sparsity_high (int): the (approximate) upper bound for the output sparsity
        num_steps (int): the number of iterations to take
        pruning (str): a string value indicating the pruning method. Possible values are `'none'`, `'simple'`, `'gw'`, and `'strong'` (all literals are case-insensitive).
        verbosity_level (int): an integer indicating how much debug output the function should produce.
    """

    # Create root and dummy edges
    all_nodes = list(range(len(graph.nodes)))
    dummy_edges, dummy_costs, root, dummy_prize = graph._add_dummy_node(connected_to=all_nodes)

    # X (np.array): the n x d data matrix (n examples in dimension d)
    dataset = dataset.reindex(columns=graph.nodes).fillna(0)
    X = np.concatenate([dataset.values, np.zeros((len(dataset),1))], axis=1)
    # y (np.array): the n-dimensional label vector. Each entry is an integer between 0 and c-1 (inclusive), where c is the number of classes.
    labeler = LabelEncoder()
    y = labeler.fit_transform(dataset.index.tolist())
    classes = labeler.classes_
    num_classes = len(classes)
    # W0 (np.array): the classes x features weight matrix
    num_nodes = len(graph.nodes)+1
    W0 = np.zeros((num_classes, num_nodes))
    # sparsity_low (int): the (approximate) lower bound for the output sparsity
    # sparsity_high (int): the (approximate) upper bound for the output sparsity
    # graph_opts (GraphOptions): passed directly to `pcst_fast`
    edges = np.concatenate((graph.edges, dummy_edges))
    num_clusters = 1
    graph_opts = gslr.GraphOptions(edges=edges, root=root, num_clusters=num_clusters, pruning=pruning)
    # steps (np.array): the step size schedule, represented by a matrix of size num_steps x num_choices. In each iteration, the algorithm tries all current choices for the step size and chooses the one that makes largest progress.
    possible_steps = np.array([0.03, 0.1, 0.3])
    steps = np.tile(possible_steps, (num_steps, 1))
    # verbosity_level (int): indicates whether intermediate output should be printed verbosity_level - 1 is being passed to pcst_fast
    # graph_proj_max_num_iter (int): the maximum number of iterations in the graph-sparsity projection.
    # edge_costs (np.array): a real vector with non-negative edge costs
    costs = np.concatenate((graph.costs, dummy_costs))
    # edge_costs_multiplier (np.array): a factor weighing edge costs vs prizes

    W_hat, losses = gslr.gslr(X, y, W0, sparsity_low, sparsity_high, graph_opts, steps, verbosity_level, edge_costs=costs, edge_costs_multiplier=(1/graph.params.b))

    # remove last column of What
    # W_hat = W_hat[:-1,:]
    W_hat = np.delete(W_hat, -1, 1)

    # Output preparation
    if num_classes == 2: W_hat = [W_hat[0]]

    class_networks = []
    for class_parameter_vector in W_hat:

        network = graph.interactome_graph.subgraph(graph.nodes[np.nonzero(class_parameter_vector)].tolist())

        nx.set_node_attributes(network, pd.DataFrame(graph.node_degrees, index=graph.nodes, columns=['degree']).loc[list(network.nodes())].astype(int).to_dict(orient='index'))

        # Post-processing
        betweenness(network)
        louvain_clustering(network)
        annotate_graph_nodes(network)

        class_networks.append(network)

    return class_networks, W_hat, losses
