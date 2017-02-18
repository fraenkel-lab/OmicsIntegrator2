#ifndef __PCST_FAST_SWIG_H__
#define __PCST_FAST_SWIG_H__

#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "pcst_fast.h"

namespace cluster_approx {

void output_function(const char* output) {
  printf(output);
  fflush(stdout);
}

std::pair<std::vector<int>, std::vector<int> > pcst_fast(
    const std::vector<std::pair<int, int> >& edges,
    const std::vector<double>& prizes,
    const std::vector<double>& costs,
    int root,
    int num_clusters,
    const std::string& pruning,
    int verbosity_level) {
  PCSTFast::PruningMethod pruning_method =
      PCSTFast::parse_pruning_method(pruning);
  PCSTFast algo(edges, prizes, costs, root, num_clusters,
                pruning_method, verbosity_level, output_function);
  std::vector<int> result_nodes;
  std::vector<int> result_edges;
  algo.run(&result_nodes, &result_edges);
  return std::make_pair(result_nodes, result_edges);
}

}  // namespace cluster_approx

#endif
