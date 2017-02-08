#include <set>
#include <string>
#include <utility>
#include <vector>

#include <math.h>
#include <matrix.h>
#include <mex.h>

#include "cluster_grid.h"
#include "mex_helper.h"
#include "pcst_fast.h"

using cluster_approx::cluster_grid_pcst;
using cluster_approx::PCSTFast;
using std::make_pair;
using std::set;
using std::string;
using std::vector;

void output_function(const char* s) {
  mexPrintf(s);
  mexEvalString("drawnow;");
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs < 4) {
    mexErrMsgTxt("At least four input argument required (values, num_clusters,"
        " k_low, k_high).");
  }
  if (nrhs > 5) {
    mexErrMsgTxt("Too many input arguments, at most four: values, num_clusters,"
        " k_low, k_high, and the options struct.");
  }
  if (nlhs > 2) {
    mexErrMsgTxt("Too many output arguments.");
  }

  // required parameters
  vector<vector<double> > values;
  int target_num_clusters;
  int k_low;
  int k_high;

  if (!get_double_matrix(prhs[0], &values)) {
    mexErrMsgTxt("Could not read the grid values as a double matrix.");
  }

  if (!get_double_as_int(prhs[1], &target_num_clusters)) {
    mexErrMsgTxt("Could not read num_clusters as a double value.");
  }

  if (!get_double_as_int(prhs[2], &k_low)) {
    mexErrMsgTxt("Could not read k_low as a double value.");
  }

  if (!get_double_as_int(prhs[3], &k_high)) {
    mexErrMsgTxt("Could not read k_high as a double value.");
  }

  int max_num_iter = 10;
  int verbosity_level = 0;
  PCSTFast::PruningMethod pruning = PCSTFast::kGWPruning;

  if (nrhs == 5) {
    const mxArray* opts = prhs[4];

    // set of accepted options
    set<string> known_options;
    known_options.insert("max_num_iter");
    known_options.insert("pruning");
    known_options.insert("verbose");

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

    if (has_field(opts, "max_num_iter")
        && !get_double_field_as_int(opts, "max_num_iter", &max_num_iter)) {
      mexErrMsgTxt("The max_num_iter field must be a double value.");
    }
    
    if (has_field(opts, "verbose")
        && !get_double_field_as_int(opts, "verbose", &verbosity_level)) {
      mexErrMsgTxt("The verbose field must be a double value.");
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
  
  vector<vector<bool> > support;
  int result_sparsity;

  bool res = cluster_grid_pcst_binsearch(values, target_num_clusters, k_low,
                                         k_high, max_num_iter, pruning,
                                         verbosity_level, output_function,
                                         &support, &result_sparsity);
  
  if (!res) {
    mexErrMsgTxt("cluster_grid_pcst returned false.");
  }

  if (nlhs >= 1) {
    set_double_matrix(&(plhs[0]), support);
  }

  if (nlhs >= 2) {
    set_double(&(plhs[1]), result_sparsity);
  }

  return;
}

