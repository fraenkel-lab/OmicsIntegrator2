%module pcst_fast
%{
#include "pcst_fast_swig.h"
%}

%include "std_pair.i"
%include "std_string.i"
%include "std_vector.i"

namespace std {
  %template(EdgePair) pair<int, int>;
  %template(DoubleVector) vector<double>;
  %template(IntVector) vector<int>;
  %template(EdgeVector) vector<pair<int, int> >;
  %template(ResultPair) pair<vector<int>, vector<int> >;
}

%ignore output_function;

%include "pcst_fast_swig.h"
