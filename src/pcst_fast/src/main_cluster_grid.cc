#include <cstdio>
#include <utility>
#include <vector>

#include "cluster_grid.h"
#include "pcst_fast.h"

using cluster_approx::cluster_grid_pcst;
using cluster_approx::PCSTFast;
using std::vector;

void output_function(const char* output) {
  printf(output);
}


int main() {
  vector<vector<double> > values;
  values.resize(5);
  values[0].push_back(1.0);
  values[0].push_back(1.0);
  values[0].push_back(0.0);
  values[0].push_back(0.0);
  values[0].push_back(0.0);

  values[1].push_back(1.0);
  values[1].push_back(1.0);
  values[1].push_back(0.0);
  values[1].push_back(0.0);
  values[1].push_back(0.0);

  values[2].push_back(0.0);
  values[2].push_back(0.0);
  values[2].push_back(0.0);
  values[2].push_back(0.0);
  values[2].push_back(0.0);

  values[3].push_back(0.0);
  values[3].push_back(0.0);
  values[3].push_back(0.0);
  values[3].push_back(0.0);
  values[3].push_back(0.0);
  
  values[4].push_back(0.0);
  values[4].push_back(0.0);
  values[4].push_back(1.0);
  values[4].push_back(1.0);
  values[4].push_back(1.0);

  int target_num_clusters = 2;
  double lambda = 0.5;
  bool include_root = false;
  double gamma = 1.0;
  PCSTFast::PruningMethod pruning = PCSTFast::kGWPruning;
  int verbosity_level = 2;


  /*vector<vector<double> > values;
  values.resize(1);
  values[0].push_back(2.1);
  values[0].push_back(0.0);
  values[0].push_back(0.0);
  values[0].push_back(0.0);
  values[0].push_back(2.2);

  int target_num_clusters = 0;
  double lambda = 1.5;
  bool include_root = true;
  double gamma = 0.1;
  PCSTFast::PruningMethod pruning = PCSTFast::kGWPruning;
  int verbosity_level = 2;*/

  vector<vector<bool> > result;
  int result_sparsity;
  int result_num_clusters;

  bool res = cluster_grid_pcst(values, target_num_clusters, lambda,
                               include_root, gamma, pruning, verbosity_level,
                               output_function, &result, &result_sparsity,
                               &result_num_clusters);

  /*int k_low = 6;
  int k_high = 8;
  int max_num_iter = 10;
  vector<vector<bool> > result;
  int result_sparsity;
  bool res = cluster_grid_pcst_binsearch(values, target_num_clusters, k_low,
      k_high, max_num_iter, pruning, verbosity_level, output_function,
      &result, &result_sparsity);*/
  
  if (!res) {
    printf("Error returned by cluster_grid_pcst\n");
    return 0;
  }

  printf("Result (sparsity %d, num_clusters %d):\n",
         result_sparsity, result_num_clusters);
  for (int yy = 0; yy < static_cast<int>(result.size()); ++yy) {
    for (int xx = 0; xx < static_cast<int>(result[0].size()); ++xx) {
      printf("%d", static_cast<int>(result[yy][xx]));
    }
    printf("\n");
  }
  
  return 0;
}
