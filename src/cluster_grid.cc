#include "cluster_grid.h"

#include <algorithm>
#include <cstdio>
#include <utility>
#include <vector>

#include "pcst_fast.h"

using std::make_pair;
using std::max;
using std::vector;

namespace cluster_approx {

const int kOutputBufferSize = 10000;
char output_buffer[kOutputBufferSize];


void build_grid_graph(const Matrix2d& values,
                      bool include_root,
                      double gamma,
                      vector<EdgePair>* edges,
                      vector<double>* prizes,
                      vector<double>* costs,
                      int* root) {
  edges->clear();
  prizes->clear();
  costs->clear();
  *root = PCSTFast::kNoRoot;
  int height = static_cast<int>(values.size());
  int width;
  if (height == 0) {
    width = 0;
    return;
  } else {
    width = static_cast<int>(values[0].size());
  }
  int n = width * height;
  if (include_root) {
    n += 1;
  }

  prizes->resize(n);

  for (int yy = 0; yy < height; ++yy) {
    for (int xx = 0; xx < width; ++xx) {
      int cur_index = yy * width + xx;
      (*prizes)[cur_index] = values[yy][xx];

      if (xx != width - 1) {
        int next_right = cur_index + 1;
        edges->push_back(make_pair(cur_index, next_right));
        costs->push_back(1.0);
      }

      if (yy != height - 1) {
        int next_down = cur_index + width;
        edges->push_back(make_pair(cur_index, next_down));
        costs->push_back(1.0);
      }
    }
  }

  if (include_root) {
    *root = n - 1;
    (*prizes)[*root] = 0.0;
    double root_edge_cost = 1.0 + gamma;

    for (int ii = 0; ii < *root; ++ii) {
      edges->push_back(make_pair(*root, ii));
      costs->push_back(root_edge_cost);
    }
  }
}


void convert_forest_to_support(const vector<int>& forest_node_indices,
                               int root,
                               int width,
                               int height,
                               Matrix2b* result,
                               int* num_clusters) {
  result->clear();
  result->resize(height);
  for (int yy = 0; yy < height; ++yy) {
    (*result)[yy].resize(width, false);
  }

  for (size_t ii = 0; ii < forest_node_indices.size(); ++ii) {
    int node = forest_node_indices[ii];
    if (node != root) {
      int node_xx = node % width;
      int node_yy = node / width;
      (*result)[node_yy][node_xx] = true;
    }
  }

  if (num_clusters != NULL) {
    *num_clusters = 0;
    vector<int> label(width * height, 0);
    vector<int> q;

    for (size_t ii = 0; ii < forest_node_indices.size(); ++ii) {
      int node = forest_node_indices[ii];
      if (node != root && label[node] == 0) {
        *num_clusters += 1;
        label[node] = *num_clusters;

        q.resize(0);
        q.push_back(node);
        size_t q_pointer = 0;
        while (q_pointer < q.size()) {
          int cur_node = q[q_pointer];
          q_pointer += 1;
          int xx = cur_node % width;
          int yy = cur_node / width;

          int next_node = 0;
          if (xx > 0) {
            next_node = yy * width + xx - 1;
            if (label[next_node] == 0 && (*result)[yy][xx]) {
              label[next_node] = *num_clusters;
              q.push_back(next_node);
            }
          }
          if (xx < width - 1) {
            next_node = yy * width + xx + 1;
            if (label[next_node] == 0 && (*result)[yy][xx]) {
              label[next_node] = *num_clusters;
              q.push_back(next_node);
            }
          }
          if (yy > 0) {
            next_node = (yy - 1) * width + xx;
            if (label[next_node] == 0 && (*result)[yy][xx]) {
              label[next_node] = *num_clusters;
              q.push_back(next_node);
            }
          }
          if (yy < height - 1) {
            next_node = (yy + 1) * width + xx;
            if (label[next_node] == 0 && (*result)[yy][xx]) {
              label[next_node] = *num_clusters;
              q.push_back(next_node);
            }
          }
        }
      }
    }
  }
}


bool cluster_grid_pcst(const Matrix2d& values,
                       int target_num_clusters,
                       double lambda,
                       bool include_root,
                       double gamma,
                       PCSTFast::PruningMethod pruning,
                       int verbosity_level,
                       void (*output_function)(const char*),
                       Matrix2b* result,
                       int* result_sparsity,
                       int* result_num_clusters) {
  vector<EdgePair> edges;
  vector<double> prizes;
  vector<double> costs;
  int root;

  build_grid_graph(values, include_root, gamma, &edges, &prizes, &costs, &root);

  for (size_t ii = 0; ii < costs.size(); ++ii) {
    costs[ii] *= lambda;
  }

  PCSTFast algo(edges, prizes, costs, root, target_num_clusters, pruning,
                verbosity_level, output_function);

  vector<int> forest_node_indices;
  vector<int> forest_edge_indices;
  bool res = algo.run(&forest_node_indices, &forest_edge_indices);
  if (!res) {
    return false;
  }
  
  convert_forest_to_support(forest_node_indices, root, values[0].size(),
                            values.size(), result, result_num_clusters);
  
  *result_sparsity = static_cast<int>(forest_node_indices.size());
  if (include_root) {
    *result_sparsity -= 1;
  }

  return true;
}


bool cluster_grid_pcst_binsearch(const Matrix2d& values,
                                 int target_num_clusters,
                                 int sparsity_low,
                                 int sparsity_high,
                                 int max_num_iter,
                                 PCSTFast::PruningMethod pruning,
                                 int verbosity_level,
                                 void (*output_function)(const char*),
                                 Matrix2b* result,
                                 int* result_sparsity) {
  vector<EdgePair> edges;
  vector<double> prizes;
  vector<double> costs;
  int root;
  
  build_grid_graph(values, false, 0.0, &edges, &prizes, &costs, &root);

  vector<double> cur_costs(costs);

  double lambda_low = 0.0;
  vector<double> sorted_prizes(prizes);
  int guess_pos = sorted_prizes.size() - sparsity_high;
  std::nth_element(sorted_prizes.begin(),
                   sorted_prizes.begin() + guess_pos,
                   sorted_prizes.end());
  double lambda_high = 2.0 * sorted_prizes[guess_pos];

  bool using_sparsity_low = false;
  bool using_max_value = false;
  if (lambda_high == 0.0) {
    guess_pos = sorted_prizes.size() - sparsity_low;
    std::nth_element(sorted_prizes.begin(),
                     sorted_prizes.begin() + guess_pos,
                     sorted_prizes.end());
    lambda_high = 2.0 * sorted_prizes[guess_pos]; 
    if (lambda_high != 0.0) {
      using_sparsity_low = true;
    } else {
      using_max_value = true;
      lambda_high = prizes[0];
      for (size_t ii = 1; ii < prizes.size(); ++ii) {
        lambda_high = max(lambda_high, prizes[ii]);
      }
      lambda_high *= 2.0;
    }
  }

  if (verbosity_level >= 1) {
    const char* sparsity_low_text = "k_low";
    const char* sparsity_high_text = "k_high";
    const char* max_value_text = "max value";
    const char* guess_text = sparsity_high_text;
    if (using_sparsity_low) {
      guess_text = sparsity_low_text;
    } else if (using_max_value) {
      guess_text = max_value_text;
    }
    snprintf(output_buffer, kOutputBufferSize, "n = %lu  c: %d  k_low: %d  "
        "k_high: %d  l_low: %e  l_high: %e  max_num_iter: %d  (using %s for "
        "initial guess).\n",
        prizes.size(), target_num_clusters, sparsity_low, sparsity_high,
        lambda_low, lambda_high, max_num_iter, guess_text);
    output_function(output_buffer);
  }
  
  int num_iter = 0;
  int cur_k = sparsity_high + 1;
  lambda_high /= 2.0;
  vector<int> forest_node_indices;
  vector<int> forest_edge_indices;

  do {
    num_iter += 1;
    lambda_high *= 2.0;

    // run algo
    for (size_t ii = 0; ii < costs.size(); ++ii) {
      cur_costs[ii] = lambda_high * costs[ii];
    }

    PCSTFast algo(edges, prizes, cur_costs, root, target_num_clusters,
                  pruning, verbosity_level, output_function);

    bool res = algo.run(&forest_node_indices, &forest_edge_indices);
    if (!res) {
      snprintf(output_buffer, kOutputBufferSize, "Algorithm returned false.\n");
      output_function(output_buffer);
      return false;
    }

    cur_k = forest_node_indices.size();

    if (verbosity_level >= 1) {
      snprintf(output_buffer, kOutputBufferSize, "increase:   l_high: %e  "
          "k: %d\n", lambda_high, cur_k);
      output_function(output_buffer);
    }
  } while (cur_k > sparsity_high && num_iter < max_num_iter);

  if (num_iter < max_num_iter && cur_k >= sparsity_low) {
    *result_sparsity = cur_k;
    convert_forest_to_support(forest_node_indices, root, values[0].size(),
                              values.size(), result, NULL);
    if (verbosity_level >= 1) {
      snprintf(output_buffer, kOutputBufferSize, "Found good lambda "
          "in exponential increase phase, returning.\n");
      output_function(output_buffer);
    }
    return true;
  }

  double lambda_mid = 0.0;

  while (num_iter < max_num_iter) {
    num_iter += 1;
    lambda_mid = (lambda_low + lambda_high) / 2.0;

    // run algo
    for (size_t ii = 0; ii < costs.size(); ++ii) {
      cur_costs[ii] = lambda_mid * costs[ii];
    }

    PCSTFast algo(edges, prizes, cur_costs, root, target_num_clusters,
                  pruning, verbosity_level, output_function);

    bool res = algo.run(&forest_node_indices, &forest_edge_indices);
    if (!res) {
      snprintf(output_buffer, kOutputBufferSize, "Algorithm returned false.\n");
      output_function(output_buffer);
      return false;
    }

    cur_k = forest_node_indices.size();
    
    if (verbosity_level >= 1) {
      snprintf(output_buffer, kOutputBufferSize, "bin_search: l_mid:  %e  "
          "k: %d  (lambda_low: %e  lambda_high: %e)\n", lambda_mid, cur_k,
          lambda_low, lambda_high);
      output_function(output_buffer);
    }

    if (cur_k <= sparsity_high && cur_k >= sparsity_low) {
      *result_sparsity = cur_k;
      convert_forest_to_support(forest_node_indices, root, values[0].size(),
                                values.size(), result, NULL);
      if (verbosity_level >= 1) {
        snprintf(output_buffer, kOutputBufferSize, "Found good lambda "
            "in binary search phase, returning.\n");
        output_function(output_buffer);
      }
      return true;
    }

    if (cur_k > sparsity_high) {
      lambda_low = lambda_mid;
    } else {
      lambda_high = lambda_mid;
    }
  }

  // run algo
  for (size_t ii = 0; ii < costs.size(); ++ii) {
    cur_costs[ii] = lambda_high * costs[ii];
  }

  PCSTFast algo(edges, prizes, cur_costs, root, target_num_clusters,
                pruning, verbosity_level, output_function);

  bool res = algo.run(&forest_node_indices, &forest_edge_indices);
  if (!res) {
    snprintf(output_buffer, kOutputBufferSize, "Algorithm returned false.\n");
    output_function(output_buffer);
    return false;
  }
  *result_sparsity = forest_node_indices.size();
  convert_forest_to_support(forest_node_indices, root, values[0].size(),
                            values.size(), result, NULL);

  if (verbosity_level >= 1) {
    snprintf(output_buffer, kOutputBufferSize, "Reached the maximum number of "
        "iterations, using the last l_high: %e  k: %d\n",
        lambda_high, *result_sparsity);
    output_function(output_buffer);
  }

  return true;
}


}  // namespace cluster_approx
