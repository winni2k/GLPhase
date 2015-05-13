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
export CXXFLAGS += -I$(SRC_DIR) -I$(SRC_DIR)/olivier -I$(GTEST_DIR)/src -I$(SRC_DIR)/vcfPhasedBreakMulti/src -I$(SRC_DIR)/htslibpp -I$(GTEST_DIR)/include -I$(GTEST_DIR)
export CPATH += :
AR= ar
LIBS += -lpthread -L$(SRC_DIR)/htslibpp

# All tests produced by this Makefile.  Remember to add new tests you
# created to the list.
TESTS := insti_unittests.bin other_unittests.bin glPack_unittests.bin glReader_tests.bin
ALLTESTS := $(TESTS)
ifndef NCUDA
TESTSCU:= cuda_tests.bin
ALLTESTS += $(TESTSCU)
endif
RUNTESTS=$(ALLTESTS:.bin=.run)
INSTILIB := $(SRC_DIR)/$(INSTILIB)
O_FILES :=  $(addprefix $(SRC_DIR)/, $(OBJECT_FILES))

# hmmLike.o contains some gpu code
INSTI_O_FILES :=  $(filter-out %main.o, $(O_FILES))
CUDA_O_FILES := $(filter-out %allCuda.link.o, $(INSTI_O_FILES))


all: $(RUNTESTS)

debug:
	$(CXX) -E -x c++ - -v < /dev/null

clean:
	rm -f $(ALLTESTS) gtest.a gtest_main.a *.o *.hpp.gch *.h.gch insti



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

# instilib is not needed here, but I'm using it as a proxy
# for updates to the headers that go into creating the library
%.o: %.cpp $(INSTILIB)
	$(CXX) $(CXXFLAGS) -c $<

$(TESTS): %.bin: %.o gtest_main.a $(INSTILIB)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

$(RUNTESTS): %.run: %.bin
	./$^

ifndef NCUDA
# compile cuda object file
# ok, need to turn on sm_20 option for separate .cu file compilation
# see: http://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#using-separate-compilation-in-cuda
# order of nvcc options is important! options should always go last!

NVCCFLAGS += -L $(GTEST_DIR) -I$(SRC_DIR)

cuda_tests.cu.o: cuda_tests.cu
	$(NVCC) $< $(NVCCFLAGS) -dc -o $@ 

# compile g++ object file
cuda_tests.o: cuda_tests.cpp $(GTEST_HEADERS) $(SRC_DIR)/*.hpp.gch
	$(CXX) $(CXXFLAGS) -c $<

#$(SRC_DIR)/hmmLike.cu.o
cuda_tests.link.o: cuda_tests.cu.o $(SRC_DIR)/hmmLike.cu.o 
	$(NVCC) $^ $(NVCCFLAGS) -dlink -o $@

cuda_tests.bin: cuda_tests.link.o gtest_main.a cuda_tests.o cuda_tests.cu.o $(SRC_DIR)/hmmLike.cu.o $(SRC_DIR)/hmmLike.o $(CUDA_O_FILES)
	$(CXX) $^ $(CXXFLAGS) -o $@ $(LIBS)

run_cuda_tests:cuda_tests
	./$^
endif

.PHONY: clean all debug $(RUNTESTS)
.SUFFIXES: 
### Makefile ends here
