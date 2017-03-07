#include <algorithm>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <math.h>
#include <matrix.h>
#include <mex.h>

#include "mex_helper.h"
#include "pcst_fast.h"

using cluster_approx::PCSTFast;
using std::make_pair;
using std::pair;
using std::set;
using std::sort;
using std::string;
using std::vector;

void output_function(const char* s) {
  mexPrintf(s);
  mexEvalString("drawnow;");
}

bool edges_unique(const vector<pair<int, int> >& edges) {
  vector<pair<int, int> > tmp_edges(edges);
  for (size_t ii = 0; ii < tmp_edges.size(); ++ii) {
    if (tmp_edges[ii].first > tmp_edges[ii].second) {
      int tmp = tmp_edges[ii].first;
      tmp_edges[ii].first = tmp_edges[ii].second;
      tmp_edges[ii].second = tmp;
    }
  }
  sort(tmp_edges.begin(), tmp_edges.end());
  for (int ii = 1; ii < static_cast<int>(tmp_edges.size()); ++ii) {
    if (tmp_edges[ii].first == tmp_edges[ii - 1].first
        && tmp_edges[ii].second == tmp_edges[ii - 1].second) {
      return false;
    }
  }
  return true;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs < 3) {
    mexErrMsgTxt("At least three input arguments required (edges, node prizes,"
        " edge_costs).");
  }
  if (nrhs > 4) {
    mexErrMsgTxt("Too many input arguments, at most four: edges, node prizes,"
        " edge_costs, and the options struct.");
  }
  if (nlhs > 2) {
    mexErrMsgTxt("Too many output arguments.");
  }

  // Get node costs
  bool prizes_are_column_vector;
  vector<double> prizes;
  if (!get_double_vector(prhs[1], &prizes, &prizes_are_column_vector)) {
    mexErrMsgTxt("Could not read the node prizes as a double vector.");
  }
  for (size_t ii = 0; ii < prizes.size(); ++ii) {
    if (prizes[ii] < 0.0) {
      mexErrMsgTxt("All node prizes must be non-negative.");
    }
  }
  int n = prizes.size();

  // Get edges
  vector<vector<int> > edge_matrix;
  if (!get_double_matrix_as_int(prhs[0], &edge_matrix)) {
    mexErrMsgTxt("Could not read the edges as a double matrix.");
  }
  int m = edge_matrix.size();
  if (m > 0 && edge_matrix[0].size() != 2) {
    mexErrMsgTxt("The edge matrix must have exactly two columns.");
  }
  vector<pair<int, int> > edges(m);
  for (int ii = 0; ii < m; ++ii) {
    int node1 = edge_matrix[ii][0];
    int node2 = edge_matrix[ii][1];
    if (node1 < 1 || node1 > n || node2 < 1 || node2 > n) {
      mexErrMsgTxt("Edge endpoints must be between 1 and n (both inclusive).");
    }
    if (node1 == node2) {
      mexErrMsgTxt("Loops are not allowed.");
    }
    edges[ii] = make_pair(node1 - 1, node2 - 1);
  }
  if (!edges_unique(edges)) {
    mexErrMsgTxt("Duplicate edges are not allowed.");
  }

  // Get edge costs
  bool costs_are_column_vector;
  vector<double> costs;
  if (!get_double_vector(prhs[2], &costs, &costs_are_column_vector)) {
    mexErrMsgTxt("Could not read the edge costs as a double vector.");
  }
  if (static_cast<int>(costs.size()) != m) {
    mexErrMsgTxt("Incorrect number of edge costs.");
  }
  for (size_t ii = 0; ii < costs.size(); ++ii) {
    if (costs[ii] < 0.0) {
      mexErrMsgTxt("All edge costs must be non-negative.");
    }
  }

  int verbosity_level = 0;
  int target_num_active_clusters = 1;
  int root = PCSTFast::kNoRoot;
  PCSTFast::PruningMethod pruning = PCSTFast::kGWPruning;

  if (nrhs == 4) {
    const mxArray* opts = prhs[3];
    // set of accepted options
    set<string> known_options;
    known_options.insert("target_num_active_clusters");
    known_options.insert("pruning");
    known_options.insert("root");
    known_options.insert("verbosity_level");

    // retrieve option struct fields
    vector<string> options;
    if (!get_fields(opts, &options)) {
      mexErrMsgTxt("Cannot get fields from options argument.");
    }

    for (size_t ii = 0; ii < options.size(); ++ii) {
      if (known_options.find(options[ii]) == known_options.end()) {
        const size_t tmp_size = 1000;
        char tmp[tmp_size];
        snprintf(tmp, tmp_size, "Unknown option \"%s\"\n", options[ii].c_str());
        mexErrMsgTxt(tmp);
      }
    }

    if (has_field(opts, "verbosity_level")
        && !get_double_field_as_int(opts, "verbosity_level",
            &verbosity_level)) {
      mexErrMsgTxt("The verbosity_level field must be a double value.");
    }

    if (has_field(opts, "root")) {
      if (!get_double_field_as_int(opts, "root", &root)) {
        mexErrMsgTxt("The root field must be a double value.");
      }
      if (root < 1 || root > n) {
        mexErrMsgTxt("The root must be between 1 and n (both inclusive).");
      }
      root -= 1;
    }

    if (has_field(opts, "target_num_active_clusters")) {
      if (!get_double_field_as_int(opts, "target_num_active_clusters",
          &target_num_active_clusters)) {
        mexErrMsgTxt("The target_num_active_clusters field must be a double "
            "value.");
      }
      if (target_num_active_clusters < 0 || target_num_active_clusters > n) {
        mexErrMsgTxt("target_num_active_clusters must be between 0 and n "
            "(both inclusive).");
      }
    }

    if (has_field(opts, "pruning")) {
      string pruning_string;
      if (!get_string_field(opts, "pruning", &pruning_string)) {
        mexErrMsgTxt("The pruning field must be a string.");
      }

      pruning = PCSTFast::parse_pruning_method(pruning_string);
      if (pruning == PCSTFast::kUnknownPruning) {
        mexErrMsgTxt("Unknown pruning method.");
      }
    }
  }

  PCSTFast algo(edges, prizes, costs, root, target_num_active_clusters, pruning,
      verbosity_level, output_function);
  
  vector<int> result_nodes;
  vector<int> result_edges;
  bool res = algo.run(&result_nodes, &result_edges);

  if (!res) {
    mexErrMsgTxt("PCSTFast returned false.");
  }

  if (nlhs >= 1) {
    for (size_t ii = 0; ii < result_nodes.size(); ++ii) {
      result_nodes[ii] += 1;
    }
    set_double_row_vector(&(plhs[0]), result_nodes);
  }

  if (nlhs >= 2) {
    for (size_t ii = 0; ii < result_edges.size(); ++ii) {
      result_edges[ii] += 1;
    }
    set_double_row_vector(&(plhs[1]), result_edges);
  }
  
  return;
}
