pcsf_fast
==============

A library for clustering graph-structured data based on the prize-collecting Steiner forest (PCSF) problem. See [this paper](http://www.jmlr.org/proceedings/papers/v37/hegde15.html) for details about the algorithm. The code for the sparse recovery experiments in the paper can be found [here](https://github.com/ludwigschmidt/graph_sparsity_experiments).

Installation
------------

The core library has no dependencies besides a basic build system for C++11.
Both g++ and clang are currently supported.
The Python wrapper requires a functioning Python build system.
The Matlab wrapper requires the mex compiler.

### Python

Compile the python wrapper with the supplied makefile:

    make pcst_fast_py

You can then import the package via `import pcst_fast`.

### Matlab

Compile the mex wrapper with the supplied makefile:

    make mexfiles

The makefile assumes that Matlab's mex compiler is available via `mex`.
The resulting mex function has the name `cluster_grid_pcst`.


Usage
-----

The Matlab and Python interfaces currently support somewhat different use cases.
The Matlab wrapper was written specifically for sparse recovery experiments.
The Python wrapper gives direct access to the PCSF solver.

### Python

The `pcst_fast` package contains the following function:

    vertices, edges = pcst_fast(edges, prizes, costs, root, num_clusters, pruning, verbosity_level)
    
The parameters are:
* `edges`: a 2D int64 array. Each row (of length 2) specifies an undirected edge in the input graph. The nodes are labeled 0 to n-1, where n is the number of nodes.
* `prizes`: the node prizes as a 1D float64 array.
* `costs`: the edge costs as a 1D float64 array.
* `root`: the root note for rooted PCST. For the unrooted variant, this parameter should be -1.
* `num_clusters`: the number of connected components in the output.
* `pruning`: a string value indicating the pruning method. Possible values are `'none'`, `'simple'`, `'gw'`, and `'strong'` (all literals are case-insensitive). `'none'` and `'simple'` return intermediate stages of the algorithm and do not have approximation guarantees. They are only intended for development. The standard GW pruning method is `'gw'`, which is also the default. `'strong'` uses "strong pruning", which was introduced in [\[JMP00\]](http://dl.acm.org/citation.cfm?id=338637). It has the same theoretical guarantees as GW pruning but better empirical performance in some cases. For the PCSF problem, the output of strong pruning is at least as good as the output of GW pruning.
* `verbosity_level`: an integer indicating how much debug output the function should produce.

The output variables are:
* `vertices`: the vertices in the solution as a 1D int64 array.
* `edges`: the edges in the output as a 1D int64 array. The list contains indices into the list of edges passed into the function.


### Matlab

The basic function signature is

    support, sparsity = cluster_grid_pcst(values, num_clusters, lambda)

The parameters have the following meaning:
* `values`: a 2D double matrix containing the value at each grid node (or equivalently, pixel). Each element must be non-negative.
* `num_clusters`: the target number of clusters in the output. This is a hard constraint which the output will always satisfy.
* `lambda`: the cost of the edges in the grid. This can be used as a Lagrangian relaxation parameter for controling the sparsity of the output: larger lambda means smaller support. Lambda must always be non-negative.

The output variables are:
* `support`: a 2D double matrix indicating the resulting support. Supported elements are marked with a 1.0, other elements with a 0.0.
* `sparsity`: the sparsity of the returned support. This output variable is optional.

Moreover, `cluster_grid_pcst` also accepts a fourth, optional paramter `opts` with the following fields, all of which are optional:
* `gamma`: If this double value is present, the underlying PCST problem is solved on a *rooted* graph. This means that the input graph contains an extra root node connected to all grid nodes. The cost of the edges from the root to the grid nodes is `lambda (1 + gamma)`. If the clustering problem should correspond to the relaxed cluster model, `num_clusters` should be set to 0. This has the effect that the number of clusters in the output is not explicitly controlled: after removing the root node, the support on the grid nodes can contain many clusters. `gamma` gives indirect control over the number of clusters in this case: larger gamma means a smaller number of clusters. Note that the overall sparsity and number of clusters depends on both `gamma` and `lambda`, the two variables are not orthogonal.
* `pruning`: a string value indicating the pruning method. Possible values are `'none'`, `'simple'`, `'gw'`, and `'strong'` (all literals are case-insensitive). `'none'` and `'simple'` return intermediate stages of the algorithm and do not have approximation guarantees. They are only intended for development. The standard GW pruning method is `'gw'`, which is also the default. `'strong'` uses "strong pruning", which was introduced in [\[JMP00\]](http://dl.acm.org/citation.cfm?id=338637). It has the same theoretical guarantees as GW pruning but better empirical performance in some cases. For the PCSF problem, the output of strong pruning is at least as good as the output of GW pruning.
* `verbose`: a double value indicating the level of log messages written to the Matlab shell. Currently, the following values are supported:
    - `0` means no output (unless errors occur).
    - `2` prints a large amount of debug output from the inner PCST algorithm.
