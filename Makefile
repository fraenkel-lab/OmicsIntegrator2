# Ludwig Schmidt (ludwigschmidt2@gmail.com) 2014
#
# This makefile is based on http://make.paulandlesley.org/autodep.html .

CXX = clang++
MEX = mex
CXXFLAGS = -std=c++11 -Wall -Wextra -O3 -fPIC
MEXCXXFLAGS = -Wall -Wextra -O3
GTESTDIR = external/googletest/googletest

SRCDIR = src
DEPDIR = .deps
OBJDIR = obj

SRCS = cluster_grid.cc pcst_fast.cc cluster_grid_pcst_mex_wrapper.cc pcst_fast_test.cc

.PHONY: clean archive

clean:
	rm -rf $(OBJDIR)
	rm -rf $(DEPDIR)
	rm -f cluster_grid_pcst.mex*
	rm -f cluster_grid_pcst_binsearch.mex*
	rm -f pcst_fast_wrap.cxx
	rm -f _pcst_fast.so
	rm -f pcst_fast.py
	rm -f pcst_fast.pyc
	rm -f pcst_fast_test

run_tests: run_pcst_fast_test

mexfiles: cluster_grid_pcst_mexfile cluster_grid_pcst_binsearch_mexfile

# gtest
$(OBJDIR)/gtest-all.o: $(GTESTDIR)/src/gtest-all.cc
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -I $(GTESTDIR)/include -I $(GTESTDIR) -c -o $@ $<

$(OBJDIR)/gtest_main.o: $(GTESTDIR)/src/gtest_main.cc
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -I $(GTESTDIR)/include -c -o $@ $<


PCST_FAST_OBJS = pcst_fast.o

PCST_FAST_TEST_OBJS = pcst_fast.o gtest-all.o gtest_main.o
pcst_fast_test: $(PCST_FAST_TEST_OBJS:%=$(OBJDIR)/%) $(SRCDIR)/pcst_fast_test.cc
	$(CXX) $(CXXFLAGS) -I $(GTESTDIR)/include -c -o $(OBJDIR)/pcst_fast_test.o $(SRCDIR)/pcst_fast_test.cc
	$(CXX) $(CXXFLAGS) -o $@ $(PCST_FAST_TEST_OBJS:%=$(OBJDIR)/%) $(OBJDIR)/pcst_fast_test.o -pthread

run_pcst_fast_test: pcst_fast_test
	./pcst_fast_test


CLUSTER_GRID_PCST_MEXFILE_SRC = cluster_grid_pcst_mex_wrapper.cc cluster_grid.cc pcst_fast.cc
CLUSTER_GRID_PCST_MEXFILE_SRC_DEPS = $(CLUSTER_GRID_PCST_MEXFILE_SRC) mex_helper.h cluster_grid.h pcst_fast.h

cluster_grid_pcst_mexfile: $(CLUSTER_GRID_PCST_MEXFILE_SRC_DEPS:%=$(SRCDIR)/%)
	$(MEX) -v CXXFLAGS="\$$CXXFLAGS $(MEXCXXFLAGS)" -output cluster_grid_pcst $(CLUSTER_GRID_PCST_MEXFILE_SRC:%=$(SRCDIR)/%)


CLUSTER_GRID_PCST_BINSEARCH_MEXFILE_SRC = cluster_grid_pcst_binsearch_mex_wrapper.cc cluster_grid.cc pcst_fast.cc
CLUSTER_GRID_PCST_BINSEARCH_MEXFILE_SRC_DEPS = $(CLUSTER_GRID_PCST_BINSEARCH_MEXFILE_SRC) mex_helper.h cluster_grid.h pcst_fast.h

cluster_grid_pcst_binsearch_mexfile: $(CLUSTER_GRID_PCST_BINSEARCH_MEXFILE_SRC_DEPS:%=$(SRCDIR)/%)
	$(MEX) -v CXXFLAGS="\$$CXXFLAGS $(MEXCXXFLAGS)" -output cluster_grid_pcst_binsearch $(CLUSTER_GRID_PCST_BINSEARCH_MEXFILE_SRC:%=$(SRCDIR)/%)


PCST_FAST_MEXFILE_SRC = pcst_fast_mex_wrapper.cc pcst_fast.cc
PCST_FAST_MEXFILE_SRC_DEPS = $(PCST_FAST_MEXFILE_SRC) mex_helper.h pcst_fast.h

pcst_fast_mexfile: $(PCST_FAST_MEXFILE_SRC_DEPS:%=$(SRCDIR)/%)
	$(MEX) -v CXXFLAGS="\$$CXXFLAGS $(MEXCXXFLAGS)" -output pcst_fast $(PCST_FAST_MEXFILE_SRC:%=$(SRCDIR)/%)


PCST_FAST_PY_SRC = pcst_fast_pybind.cc
PCST_FAST_PY_SRC_DEPS = $(PCST_FAST_PY_SRC) pcst_fast.h pcst_fast.cc
pcst_fast_py: $(PCST_FAST_PY_SRC_DEPS:%=$(SRCDIR)/%)
	$(CXX) $(CXXFLAGS) -shared -I $(SRCDIR) -I external/pybind11/include `python-config --cflags --ldflags` $(SRCDIR)/pcst_fast_pybind.cc $(SRCDIR)/pcst_fast.cc -o pcst_fast.so
