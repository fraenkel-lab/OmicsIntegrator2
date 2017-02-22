#ifndef __CLUSTER_GRID_H__
#define __CLUSTER_GRID_H__

#include <utility>
#include <vector>

#include "pcst_fast.h"

namespace cluster_approx {

typedef std::vector<std::vector<double> > Matrix2d;
typedef std::vector<std::vector<bool> > Matrix2b;
typedef std::pair<int, int> EdgePair;


void build_grid_graph(const Matrix2d& values,
                      bool include_root,
                      double gamma,
                      int* n,
                      std::vector<EdgePair>* edges,
                      std::vector<double>* prizes,
                      std::vector<double>* costs,
                      int* root);


void convert_forest_to_support(const std::vector<int>& forest_node_indices,
                               int root,
                               int width,
                               int height,
                               Matrix2b* result,
                               int* num_clusters);


bool cluster_grid_pcst(const Matrix2d& values,
                       int target_num_clusters,
                       double lambda,
                       bool include_root,
                       double gamma,
                       PCSTFast::PruningMethod pruning,
                       int verbosity_level,
                       void (*output_function_)(const char*),
                       Matrix2b* result,
                       int* result_sparsity,
                       int* result_num_clusters);


bool cluster_grid_pcst_binsearch(const Matrix2d& values,
                                 int target_num_clusters,
                                 int sparsity_low,
                                 int sparsity_high,
                                 int max_num_iter,
                                 PCSTFast::PruningMethod pruning,
                                 int verbosity_level,
                                 void (*output_function_)(const char*),
                                 Matrix2b* result,
                                 int* result_sparsity);


}  // namespace cluster_approx


#endif
