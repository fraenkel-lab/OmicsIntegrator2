# Ludwig Schmidt (ludwigschmidt2@gmail.com) 2014
#
# This makefile is based on http://make.paulandlesley.org/autodep.html .

CXX = g++
MEX = mex
CXXFLAGS = -Wall -Wextra -O3 -fPIC
MEXCXXFLAGS = -Wall -Wextra -O3
GTESTDIR = /usr/src/gtest

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
	$(CXX) $(CXXFLAGS) -I $(GTESTDIR) -c -o $@ $<


PCST_FAST_OBJS = pcst_fast.o

PCST_FAST_TEST_OBJS = $(PCST_FAST_OBJS) pcst_fast_test.o gtest-all.o
pcst_fast_test: $(PCST_FAST_TEST_OBJS:%=$(OBJDIR)/%)
	$(CXX) $(CXXFLAGS) -o $@ $^ -pthread

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


PCST_FAST_SWIG_SRC = pcst_fast.cc
PCST_FAST_SWIG_SRC_DEPS = $(PCST_FAST_SWIG_SRC) pcst_fast_swig.h pcst_fast.h pcst_fast.i

pcst_fast_swig: $(PCST_FAST_SWIG_SRC_DEPS:%=$(SRCDIR)/%)
	swig -c++ -python -builtin -outcurrentdir $(SRCDIR)/pcst_fast.i
	mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $(SRCDIR)/pcst_fast.cc -o $(OBJDIR)/pcst_fast.o
	$(CXX) $(CXXFLAGS) `python-config --includes` -c pcst_fast_wrap.cxx -I $(SRCDIR) -o $(OBJDIR)/pcst_fast_wrap.o
	$(CXX) -shared $(OBJDIR)/pcst_fast.o $(OBJDIR)/pcst_fast_wrap.o -o _pcst_fast.so `python-config --ldflags` 
	rm -f pcst_fast_wrap.cxx



$(OBJDIR)/%.o: $(SRCDIR)/%.cc
	# Create the directory the current target lives in.
	@mkdir -p $(@D)
	# Compile and generate a dependency file.
	# See http://gcc.gnu.org/onlinedocs/gcc/Preprocessor-Options.html .
	$(CXX) $(CXXFLAGS) -MMD -MP -c -o $@ $<
	# Move dependency file to dependency file directory.
	# Create the dependency file directory if necessary.
	@mkdir -p $(DEPDIR)
	@mv $(OBJDIR)/$*.d $(DEPDIR)/$*.d

# Include the generated dependency files.
# The command replaces each file name in SRCS with its dependency file.
# See http://www.gnu.org/software/make/manual/html_node/Substitution-Refs.html#Substitution-Refs for the GNU make details.
-include $(SRCS:%.cc=$(DEPDIR)/%.d)
