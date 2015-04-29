
#targets
.PHONY: clean test all debug static release disttest


BDIR := ../bin
VERSION := $(shell grep VERSION_MAJOR version.hpp | cut -f3 -d' ').$(shell grep VERSION_MINOR version.hpp | cut -f3 -d' ').$(shell grep VERSION_REVISION version.hpp | cut -f3 -d' ')

# object and header files are defined here
SOURCES=impute.cpp insti.cpp emcchain.cpp olivier/utils.cpp haplotype.cpp kMedoids.cpp kNN.cpp hapPanel.cpp sampler.cpp glPack.cpp tabix.cpp glVQ.cpp geneticMap.cpp bio.cpp glReader.cpp
OBJECT_FILES=$(SOURCES:.cpp=.o)
HEADER_FILES=$(SOURCES:.cpp=.hpp.gch)
ifndef NCUDA
OBJECT_FILES += hmmLike.cu.o allCuda.link.o hmmLike.o
endif

ifndef NCUDA
HEADER_FILES += hmmLike.hpp.gch
endif

BINBASE := insti.$(VERSION)
ifdef NCUDA
BINBASE := $(BINBASE).ncuda
CXXFLAGS += -DNCUDA
endif

#system specific choices
NOMP = 
SYSLINKFLAGS_START :=
SYSLINKFLAGS_END :=
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    SYSLINKFLAGS_START += -Wl,-whole-archive
    SYSLINKFLAGS_END += -Wl,-no-whole-archive
endif
# Darwin does not support OMP, so disable
ifeq ($(UNAME_S),Darwin)
    NOMP = 1
endif


# disables OMP
ifndef NOMP
CXXFLAGS += -fopenmp
else
BINBASE := $(BINBASE).nomp
CXXFLAGS += -DNOMP
endif

BIN := $(BDIR)/$(BINBASE)
INSTILIB=$(BINBASE).a
HTSLIBPP_DIR:=htslibpp
HTSLIBPP:=$(HTSLIBPP_DIR)/libhtspp.a


release: BIN += .release


# targets
all: $(BIN)
release: $(BIN)
debug: $(BIN).debug
test: debug
test: $(INSTILIB) $(BIN).debug $(HEADER_FILES)
	$(MAKE) -C testing

disttest: test debug
	../t/regression_tests.pl $(BIN).debug

# define CXXFLAGS
all: CXXFLAGS += -O3
release: CXXFLAGS += -O3 -DNDEBUG
static: CXXFLAGS += -O3 -DNDEBUG
debug: CXXFLAGS += -O0 -ggdb
test: CXXFLAGS += -O0 -ggdb

# gcc
CXX := g++
CXXFLAGS += -Wall -Wextra -std=c++11 -mpopcnt -Ihtslibpp
BOOSTLIB := -lboost_iostreams #alternatively: /opt/boost/boost-1.54.0-intel/lib/libboost_iostreams.a
LIBS += -lgsl -lgslcblas -lhts -lbz2 $(BOOSTLIB) -lm -lstdc++ -lz -lhtspp

# nvcc
ifndef NCUDA
all:NVCCFLAGS += -O3
release: NVCCFLAGS += -O3

NVCCFLAGS += -arch=sm_35 -m64 #--ptxas-options=-v
NVCC := nvcc
CUDALIBS := -lcudart
LIBS += $(CUDALIBS)
test: NVCCFLAGS += -g -G
endif

# path variables
CPLUS_INCLUDE_PATH += :olivier:vcfPhasedBreakMulti/src
LIBRARY_PATH := $(LIBRARY_PATH)
LD_LIBRARY_PATH := $(LIBRARY_PATH):$(LD_LIBRARY_PATH)


# perform all exporting here
export CPLUS_INCLUDE_PATH
export LIBRARY_PATH
export LD_LIBRARY_PATH
export CXXFLAGS
export CXX
export LIBS
export OBJECT_FILES
export INSTILIB

ifndef NCUDA
export NVCCFLAGS
export NVCC
endif

static: KERNEL := $(shell uname -r)
static: SLIBS := u-wkretzsch-snptools/samtools-0.1.16/libbam.a /usr/lib64/libm.a ~/lib/libz.a

# I'm deprecating this
#testDebug: $(OBJECT_FILES) $(HEADER_FILES) debug
#	$(MAKE) -C testing debug

clean:
	rm -f *.a *.o *.hpp.gch olivier/utils.hpp.gch
	$(MAKE) -C htslibpp clean

.DELETE_ON_ERROR:

### normal compilation

## precompile headers
%.hpp.gch: %.hpp
	$(CXX) $(CXXFLAGS) -x c++-header $<

ifndef NCUDA
hmmLike.hpp.gch: hmmLike.hpp sampler.hpp.gch olivier/utils.hpp.gch glPack.hpp.gch hmmLike.h
	$(CXX) $(CXXFLAGS) -x c++-header $<
endif

glReader.hpp.gch: glReader.hpp.gch olivier/utils.hpp.gch bio.hpp.gch

glPack.hpp.gch: glPack.hpp olivier/utils.hpp.gch glVQ.hpp.gch
	$(CXX) $(CXXFLAGS) -x c++-header $<

hapPanel.hpp.gch: hapPanel.hpp bio.hpp.gch
	$(CXX) $(CXXFLAGS) -x c++-header $<

olivier/utils.hpp.gch: olivier/utils.hpp
	$(CXX) $(CXXFLAGS) -x c++-header $<

kNN.hpp.gch: kNN.hpp haplotype.hpp.gch sampler.hpp.gch
	$(CXX) $(CXXFLAGS) -x c++-header $<

kMedoids.hpp.gch: kMedoids.hpp haplotype.hpp.gch sampler.hpp.gch
	$(CXX) $(CXXFLAGS) -x c++-header $<

MHSampler.hpp.gch: MHSampler.hpp olivier/utils.hpp.gch
	$(CXX) $(CXXFLAGS) -x c++-header $<

geneticMap.hpp.gch: geneticMap.hpp olivier/utils.hpp.gch
	$(CXX) $(CXXFLAGS) -x c++-header $<

INSTIHEADERPRE :=insti.hpp impute.hpp.gch emcchain.hpp.gch olivier/utils.hpp.gch sampler.hpp.gch version.hpp.gch bio.hpp.gch hapPanel.hpp.gch MHSampler.hpp.gch kNN.hpp.gch kMedoids.hpp.gch glPack.hpp.gch tabix.hpp.gch geneticMap.hpp.gch

ifndef NCUDA
INSTIHEADERPRE += hmmLike.hpp.gch
endif

insti.hpp.gch: $(INSTIHEADERPRE)
	$(CXX) $(CXXFLAGS) -x c++-header $<

## object files
%.o: %.cpp %.hpp.gch
	$(CXX) $(CXXFLAGS) -c $<

glReader.o: glReader.cpp glReader.hpp.gch bio.hpp.gch olivier/utils.hpp.gch
	$(CXX) $(CXXFLAGS) -c $< -o $@

olivier/utils.o: olivier/utils.cpp olivier/utils.hpp.gch
	$(CXX) $(CXXFLAGS) -Wno-sign-compare -c $< -o $@

hmmLike.o: hmmLike.cpp hmmLike.hpp.gch hmmLike.h
	$(CXX) $(CXXFLAGS) -c $<

ifndef NCUDA
hmmLike.cu.o: hmmLike.cu hmmLike.hcu hmmLike.h
	$(NVCC) $< $(NVCCFLAGS) -dc -o $@

allCuda.link.o: hmmLike.cu.o
	$(NVCC) $^ $(NVCCFLAGS) -dlink -o $@
endif

main.o: main.cpp version.hpp.gch impute.hpp.gch insti.hpp.gch
	$(CXX) $(CXXFLAGS) -c $<

## library
$(INSTILIB): $(OBJECT_FILES) $(HTSLIBPP)
	$(AR) $(ARFLAGS) $@ $(OBJECT_FILES)

$(HTSLIBPP):
	$(MAKE) -C htslibpp

## binaries
## need to include --whole-archive option for static member variables
$(BIN): $(INSTILIB) main.o
	$(CXX) $(CXXFLAGS) $(SYSLINKFLAGS_START) $^ $(SYSLINKFLAGS_END) -o $@ $(LIBS) -Lhtslibpp
	cd $(BDIR) && rm -f insti && ln -s $@ insti

static: $(INSTILIB) main.o
	$(CXX) $(CXXFLAGS) $^ -o $(BDIR)/insti.$(VERSION).static-$(KERNEL) $(SLIBS)

$(BIN).debug: $(INSTILIB) main.o
	$(CXX) $(CXXFLAGS) $(SYSLINKFLAGS_START) $^ $(SYSLINKFLAGS_END) -o $@ $(LIBS)  -Lhtslibpp 


## installation
oxford: all
	cd $(BDIR) && cp insti.$(VERSION) ~/bin/
	cd ~/bin && rm -f insti && ln -s insti.$(VERSION) insti

.SUFFIXES:
