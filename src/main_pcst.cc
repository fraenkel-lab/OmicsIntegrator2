#include <cstdio>
#include <iostream>
#include <utility>
#include <vector>

#include "pcst_fast.h"

using cluster_approx::PCSTFast;
using std::cout;
using std::endl;
using std::make_pair;
using std::pair;
using std::vector;

void output_function(const char* s) {
  fprintf(stderr, "%s", s);
  fflush(stderr);
}

int main() {
  // Now included in tests
  /*int n = 3;
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  vector<double> prizes;
  prizes.push_back(0);
  prizes.push_back(5);
  prizes.push_back(6);
  vector<double> costs;
  costs.push_back(3);
  costs.push_back(4);
  int root = 0;
  int target_num_active_clusters = 0;
  PCSTFast::PruningMethod pruning = PCSTFast::kNoPruning;
  int verbosity_level = 2;*/

  // Now included in tests
  /*int n = 4;
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  edges.push_back(make_pair(2, 3));
  vector<double> prizes;
  prizes.push_back(10);
  prizes.push_back(0);
  prizes.push_back(1);
  prizes.push_back(10);
  vector<double> costs;
  costs.push_back(10);
  costs.push_back(4);
  costs.push_back(3);
  int root = 0;
  int target_num_active_clusters = 0;
  PCSTFast::PruningMethod pruning = PCSTFast::kNoPruning;
  int verbosity_level = 2;*/

  // Now included in tests
  /*int n = 4;
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  edges.push_back(make_pair(2, 3));
  vector<double> prizes;
  prizes.push_back(10);
  prizes.push_back(10);
  prizes.push_back(1);
  prizes.push_back(10);
  vector<double> costs;
  costs.push_back(10);
  costs.push_back(6);
  costs.push_back(5);
  int root = 0;
  int target_num_active_clusters = 0;
  PCSTFast::PruningMethod pruning = PCSTFast::kNoPruning;
  int verbosity_level = 2;*/

  // Now included in tests
  /*int n = 3;
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  vector<double> prizes;
  prizes.push_back(10);
  prizes.push_back(3);
  prizes.push_back(3);
  vector<double> costs;
  costs.push_back(100);
  costs.push_back(2);
  int root = 0;
  int target_num_active_clusters = 0;
  PCSTFast::PruningMethod pruning = PCSTFast::kNoPruning;
  int verbosity_level = 2;*/

  // Now included in tests
  /*int n = 3;
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  vector<double> prizes;
  prizes.push_back(10);
  prizes.push_back(3);
  prizes.push_back(3);
  vector<double> costs;
  costs.push_back(100);
  costs.push_back(2);
  int root = -1;
  int target_num_active_clusters = 2;
  PCSTFast::PruningMethod pruning = PCSTFast::kGWPruning;
  int verbosity_level = 2;*/

  // Now included in tests
  /*int n = 3;
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  vector<double> prizes;
  prizes.push_back(10);
  prizes.push_back(0);
  prizes.push_back(3);
  vector<double> costs;
  costs.push_back(100);
  costs.push_back(2);
  int root = -1;
  int target_num_active_clusters = 2;
  PCSTFast::PruningMethod pruning = PCSTFast::kGWPruning;
  int verbosity_level = 2;*/

  // Now included in tests
  /*int n = 4;
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  edges.push_back(make_pair(2, 3));
  vector<double> prizes;
  prizes.push_back(10);
  prizes.push_back(0);
  prizes.push_back(6);
  prizes.push_back(6);
  vector<double> costs;
  costs.push_back(100);
  costs.push_back(2);
  costs.push_back(5);
  int root = -1;
  int target_num_active_clusters = 2;
  PCSTFast::PruningMethod pruning = PCSTFast::kGWPruning;
  int verbosity_level = 2;*/

  // Now included in tests
  /*int n = 10;
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  edges.push_back(make_pair(2, 3));
  edges.push_back(make_pair(0, 9));
  edges.push_back(make_pair(0, 2));
  edges.push_back(make_pair(0, 3));
  edges.push_back(make_pair(0, 5));
  edges.push_back(make_pair(1, 9));
  edges.push_back(make_pair(1, 3));
  edges.push_back(make_pair(1, 5));
  edges.push_back(make_pair(1, 7));
  edges.push_back(make_pair(2, 8));
  edges.push_back(make_pair(2, 3));
  edges.push_back(make_pair(3, 4));
  edges.push_back(make_pair(3, 5));
  edges.push_back(make_pair(3, 6));
  edges.push_back(make_pair(3, 7));
  edges.push_back(make_pair(3, 8));
  edges.push_back(make_pair(3, 9));
  edges.push_back(make_pair(4, 5));
  edges.push_back(make_pair(4, 6));
  edges.push_back(make_pair(4, 7));
  edges.push_back(make_pair(5, 8));
  edges.push_back(make_pair(6, 8));
  vector<double> prizes;
  prizes.push_back(0.032052554364677466);
  prizes.push_back(0.32473378289799926);
  prizes.push_back(0.069699345546302638);
  prizes.push_back(0);
  prizes.push_back(0.74867253235151754);
  prizes.push_back(0.19804330340026255);
  prizes.push_back(0.85430521133171622);
  prizes.push_back(0.83819939651391351);
  prizes.push_back(0.71744625276884877);
  prizes.push_back(0.016798567754083948);
  vector<double> costs;
  costs.push_back(100);
  costs.push_back(2);
  costs.push_back(5);
  costs.push_back(0.8);
  costs.push_back(0.8);
  costs.push_back(0.8800000000000001);
  costs.push_back(0.8);
  costs.push_back(0.8);
  costs.push_back(0.8800000000000001);
  costs.push_back(0.8);
  costs.push_back(0.8);
  costs.push_back(0.8);
  costs.push_back(0.8800000000000001);
  costs.push_back(0.8800000000000001);
  costs.push_back(0.8800000000000001);
  costs.push_back(0.8800000000000001);
  costs.push_back(0.8800000000000001);
  costs.push_back(0.8800000000000001);
  costs.push_back(0.8800000000000001);
  costs.push_back(0.8);
  costs.push_back(0.8);
  costs.push_back(0.8);
  costs.push_back(0.8);
  costs.push_back(0.8);
  int root = 3;
  int target_num_active_clusters = 0;
  PCSTFast::PruningMethod pruning = PCSTFast::kGWPruning;
  int verbosity_level = 2;*/

  // Now included in tests
  /*int n = 8;
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(1, 2));
  edges.push_back(make_pair(2, 3));
  edges.push_back(make_pair(3, 4));
  edges.push_back(make_pair(4, 5));
  edges.push_back(make_pair(5, 6));
  edges.push_back(make_pair(6, 7));
  vector<double> prizes;
  prizes.push_back(100.0);
  prizes.push_back(0.0);
  prizes.push_back(0.0);
  prizes.push_back(1.0);
  prizes.push_back(0.0);
  prizes.push_back(0.0);
  prizes.push_back(0.0);
  prizes.push_back(100.0);
  vector<double> costs;
  costs.push_back(0.9);
  costs.push_back(0.9);
  costs.push_back(0.9);
  costs.push_back(0.9);
  costs.push_back(0.9);
  costs.push_back(0.9);
  costs.push_back(0.9);
  int root = -1;
  int target_num_active_clusters = 1;
  PCSTFast::PruningMethod pruning = PCSTFast::kGWPruning;
  int verbosity_level = 2;*/

  // Now included in tests
  int n = 5;
  vector<pair<int, int> > edges;
  edges.push_back(make_pair(0, 1));
  edges.push_back(make_pair(0, 2));
  edges.push_back(make_pair(2, 3));
  edges.push_back(make_pair(3, 4));
  vector<double> prizes;
  prizes.push_back(0.0);
  prizes.push_back(2.2);
  prizes.push_back(0.0);
  prizes.push_back(0.0);
  prizes.push_back(2.1);
  vector<double> costs;
  costs.push_back(1.0);
  costs.push_back(1.0);
  costs.push_back(1.0);
  costs.push_back(1.0);
  int root = -1;
  int target_num_active_clusters = 1;
  PCSTFast::PruningMethod pruning = PCSTFast::kStrongPruning;
  int verbosity_level = 2;

  PCSTFast algo(n, edges, prizes, costs, root, target_num_active_clusters,
                pruning, verbosity_level, output_function);
  
  vector<int> node_result;
  vector<int> edge_result;
  algo.run(&node_result, &edge_result);

  cout << "Result: " << endl;
  cout << "Nodes: " << endl;
  for (size_t ii = 0; ii < node_result.size(); ++ii) {
    cout << node_result[ii] << " ";
  }
  cout << endl;
  cout << "Edges: " << endl;
  for (size_t ii = 0; ii < edge_result.size(); ++ii) {
    cout << edge_result[ii] << " ";
  }
  cout << endl;

  return 0;
}
