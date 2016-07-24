### Makefile --- for unittesting

## Author: wkretzsch@gmail.com
## Version: $Id: Makefile,v 0.0 2013/06/11 20:46:11 winni Exp $
## Keywords: 
## X-URL: 

# where to find src code to be tested
SRC_DIR := ..

# Points to the root of Google Test, relative to where this file is.
# Remember to tweak this if you move this file.
GTEST_DIR = $(SRC_DIR)/gtest-1.7.0

# google test throwing warnings, workaround is to mark the dir as a system header dir
export CXXFLAGS += -isystem $(GTEST_DIR)

export CPATH += :$(GTEST_DIR)/include:$(SRC_DIR)/u-wkretzsch-snptools/samtools-0.1.16:$(SRC_DIR)/u-wkretzsch-snptools/tabix-0.2.5:$(GTEST_DIR)
export CPLUS_INCLUDE_PATH += :$(SRC_DIR):$(SRC_DIR)/olivier:$(GTEST_DIR)/src:$(SRC_DIR)/tabixpp:$(SRC_DIR)/vcfPhasedBreakMulti/src
export LIBRARY_PATH += :$(SRC_DIR)/u-wkretzsch-snptools/samtools-0.1.16:$(SRC_DIR)/tabixpp
AR= ar
LIBS += -lpthread

# All tests produced by this Makefile.  Remember to add new tests you
# created to the list.
TESTS := insti_unittests other_unittests glPack_unittests hmm_unittests
ifndef NCUDA
TESTS += cuda_tests
endif
RUNTESTS:=$(addprefix run_,$(TESTS))
O_FILES :=  $(addprefix $(SRC_DIR)/, $(OBJECT_FILES))

# hmmLike.o contains some gpu code
INSTI_O_FILES :=  $(filter-out %main.o, $(O_FILES))
CUDA_O_FILES := $(filter-out %allCuda.link.o, $(INSTI_O_FILES))

all: $(RUNTESTS)


debug:
	$(CXX) -E -x c++ - -v < /dev/null

clean:
	rm -f $(TESTS) gtest.a gtest_main.a *.o *.hpp.gch *.h.gch insti



# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h


# Builds gtest.a and gtest_main.a.

# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CXXFLAGS) -c $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CXXFLAGS) -c $(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

# Builds a sample test.  A test should link with either gtest.a or
# gtest_main.a, depending on whether it defines its own main()
# function.

insti_unittests.o : insti_unittests.cpp $(SRC_DIR)/*.hpp.gch $(GTEST_HEADERS)
	$(CXX) $(CXXFLAGS) -c $<

insti_unittests : insti_unittests.o gtest_main.a $(INSTI_O_FILES) 
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

run_insti_unittests: insti_unittests
	./$^

other_unittests.o : other_unittests.cpp $(SRC_DIR)/*.hpp.gch $(GTEST_HEADERS)
	$(CXX) $(CXXFLAGS) -c $<

other_unittests : other_unittests.o gtest_main.a $(INSTI_O_FILES)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

run_other_unittests: other_unittests
	./$^

glPack_unittests.o : glPack_unittests.cpp $(SRC_DIR)/glPack.hpp.gch $(SRC_DIR)/glVQ.hpp.gch $(GTEST_HEADERS)
	$(CXX) $(CXXFLAGS) -c $<

glPack_unittests : glPack_unittests.o gtest_main.a $(SRC_DIR)/glPack.o $(SRC_DIR)/glVQ.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

run_glPack_unittests: glPack_unittests
	./$^

hmm_unittests.o : hmm_unittests.cpp $(SRC_DIR)/*.hpp.gch $(GTEST_HEADERS)
	$(CXX) $(CXXFLAGS) -c $<

hmm_unittests : hmm_unittests.o gtest_main.a $(INSTI_O_FILES)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

run_hmm_unittests: hmm_unittests
	./$^


ifndef NCUDA
# compile cuda object file
# ok, need to turn on sm_20 option for separate .cu file compilation
# see: http://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#using-separate-compilation-in-cuda
# order of nvcc options is important! options should always go last!

NVCCFLAGS += -L $(GTEST_DIR)

cuda_tests.cu.o: cuda_tests.cu
	$(NVCC) $< $(NVCCFLAGS) -dc -o $@ 

# compile g++ object file
cuda_tests.o: cuda_tests.cpp $(GTEST_HEADERS) $(SRC_DIR)/*.hpp.gch
	$(CXX) $(CXXFLAGS) -c $<

#$(SRC_DIR)/hmmLike.cu.o
cuda_tests.link.o: cuda_tests.cu.o $(SRC_DIR)/hmmLike.cu.o 
	$(NVCC) $^ $(NVCCFLAGS) -dlink -o $@

cuda_tests: cuda_tests.link.o gtest_main.a cuda_tests.o cuda_tests.cu.o $(SRC_DIR)/hmmLike.cu.o $(SRC_DIR)/hmmLike.o $(CUDA_O_FILES)
	$(CXX) $^ $(CXXFLAGS) -o $@ $(LIBS)

run_cuda_tests:cuda_tests
	./$^
endif

.PHONY: clean all debug $(RUNTESTS)
.SUFFIXES:
### Makefile ends here
