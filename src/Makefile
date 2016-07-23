
#targets
.PHONY: clean test all debug static release


BDIR := ../bin
VERSION := $(shell grep VERSION_MAJOR version.hpp | cut -f3 -d' ').$(shell grep VERSION_MINOR version.hpp | cut -f3 -d' ').$(shell grep VERSION_REVISION version.hpp | cut -f3 -d' ')

# object and header files are defined here
OBJECT_FILES := main.o impute.o insti.o emcchain.o olivier/utils.o haplotype.o kMedoids.o kNN.o hapPanel.o sampler.o glPack.o tabixpp/tabix.o glVQ.o geneticMap.o hmm.o
ifndef NCUDA
OBJECT_FILES += hmmLike.cu.o allCuda.link.o hmmLike.o
endif

HEADER_FILES := insti.hpp.gch impute.hpp.gch emcchain.hpp.gch olivier/utils.hpp.gch version.hpp.gch bio.hpp.gch hapPanel.hpp.gch MHSampler.hpp.gch sampler.hpp.gch tabixpp/tabix.hpp.gch globals.h glVQ.hpp geneticMap.hpp.gch hmm.hpp.gch
ifndef NCUDA
HEADER_FILES += hmmLike.hpp.gch
endif

ifndef NCUDA
BIN := $(BDIR)/insti.$(VERSION)
else
BIN :=  $(BDIR)/insti.$(VERSION).ncuda
CXXFLAGS += -DNCUDA
endif

release: BIN += .release

# targets
all: $(BIN)
release: $(BIN)
debug: $(BIN).debug
test: debug
test: $(OBJECT_FILES) $(HEADER_FILES)
	$(MAKE) -C testing

# define CXXFLAGS
all: CXXFLAGS += -O3
release: CXXFLAGS += -O3 -DNDEBUG
static: CXXFLAGS += -O3 -DNDEBUG
debug: CXXFLAGS += -O0 -ggdb
test: CXXFLAGS += -O0 -ggdb

# gcc
CXX := g++-4.7+
CXXFLAGS += -Wall -Wextra -std=c++11 -mpopcnt -fopenmp
BOOSTLIB := -lboost_iostreams-mt #alternatively: /opt/boost/boost-1.54.0-intel/lib/libboost_iostreams.a
LIBS += -lgsl -lgslcblas -ltabix -lbz2 $(BOOSTLIB) -lm -lstdc++ -lz

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
CPLUS_INCLUDE_PATH += :olivier:tabixpp:vcfPhasedBreakMulti/src
LIBRARY_PATH := tabixpp:$(LIBRARY_PATH)
LD_LIBRARY_PATH := $(LIBRARY_PATH):$(LD_LIBRARY_PATH)


# perform all exporting here
export CPLUS_INCLUDE_PATH
export LIBRARY_PATH
export LD_LIBRARY_PATH
export CXXFLAGS
export CXX
export LIBS
export OBJECT_FILES

ifndef NCUDA
export NVCCFLAGS
export NVCC
endif

static: KERNEL := $(shell uname -r)
static: SLIBS := u-wkretzsch-snptools/samtools-0.1.16/libbam.a u-wkretzsch-snptools/tabix-0.2.5/libtabix.a /usr/lib64/libm.a ~/lib/libz.a ~/li

# I'm deprecating this
#testDebug: $(OBJECT_FILES) $(HEADER_FILES) debug
#	$(MAKE) -C testing debug

clean:
	rm -f *.o *.hpp.gch olivier/utils.hpp.gch


.DELETE_ON_ERROR:

### normal compilation

## precompile headers
%.hpp.gch: %.hpp
	$(CXX) $(CXXFLAGS) -x c++-header $<

ifndef NCUDA
hmmLike.hpp.gch: hmmLike.hpp sampler.hpp.gch olivier/utils.hpp.gch glPack.hpp.gch hmmLike.h
	$(CXX) $(CXXFLAGS) -x c++-header $<
endif

glPack.hpp.gch: glPack.hpp olivier/utils.hpp.gch glVQ.hpp.gch
	$(CXX) $(CXXFLAGS) -x c++-header $<

hapPanel.hpp.gch: hapPanel.hpp bio.hpp.gch
	$(CXX) $(CXXFLAGS) -x c++-header $<

hmm.hpp.gch: hmm.hpp impute.hpp.gch
	$(CXX) $(CXXFLAGS) -x c++-header $<

olivier/utils.hpp.gch: olivier/utils.hpp
	$(CXX) $(CXXFLAGS) -x c++-header $<

tabixpp/tabix.hpp.gch: tabixpp/tabix.hpp
	$(CXX) $(CXXFLAGS) -x c++-header $<

kNN.hpp.gch: kNN.hpp haplotype.hpp.gch sampler.hpp.gch
	$(CXX) $(CXXFLAGS) -x c++-header $<

kMedoids.hpp.gch: kMedoids.hpp haplotype.hpp.gch sampler.hpp.gch
	$(CXX) $(CXXFLAGS) -x c++-header $<

MHSampler.hpp.gch: MHSampler.hpp olivier/utils.hpp.gch
	$(CXX) $(CXXFLAGS) -x c++-header $<

geneticMap.hpp.gch: geneticMap.hpp olivier/utils.hpp.gch
	$(CXX) $(CXXFLAGS) -x c++-header $<

INSTIHEADERPRE :=insti.hpp impute.hpp.gch emcchain.hpp.gch olivier/utils.hpp.gch sampler.hpp.gch version.hpp.gch bio.hpp.gch hapPanel.hpp.gch MHSampler.hpp.gch kNN.hpp.gch kMedoids.hpp.gch glPack.hpp.gch tabixpp/tabix.hpp.gch geneticMap.hpp.gch hmm.hpp.gch

ifndef NCUDA
INSTIHEADERPRE += hmmLike.hpp.gch
endif

insti.hpp.gch: $(INSTIHEADERPRE)
	$(CXX) $(CXXFLAGS) -x c++-header $<

## object files
%.o: %.cpp %.hpp.gch
	$(CXX) $(CXXFLAGS) -c $< $(LIBS)

olivier/utils.o: olivier/utils.cpp olivier/utils.hpp.gch
	$(CXX) $(CXXFLAGS) -Wno-sign-compare -c $< -o $@ $(LIBS)

tabixpp/tabix.o: 
	$(MAKE) -C tabixpp

hmmLike.o: hmmLike.cpp hmmLike.hpp.gch hmmLike.h
	$(CXX) $(CXXFLAGS) -c $< $(LIBS)

ifndef NCUDA
hmmLike.cu.o: hmmLike.cu hmmLike.hcu hmmLike.h
	$(NVCC) $< $(NVCCFLAGS) -dc -o $@

allCuda.link.o: hmmLike.cu.o
	$(NVCC) $^ $(NVCCFLAGS) -dlink -o $@
endif

main.o: main.cpp version.hpp.gch impute.hpp.gch insti.hpp.gch
	$(CXX) $(CXXFLAGS) -c $< $(LIBS)

## binaries
$(BIN): $(OBJECT_FILES)
	$(CXX) $(CXXFLAGS) $(OBJECT_FILES) -o $@ $(LIBS)
	cd $(BDIR) && rm -f insti && ln -s $@ insti

static: $(OBJECT_FILES)
	$(CXX) $(CXXFLAGS) $(OBJECT_FILES) -o $(BDIR)/insti.$(VERSION).static-$(KERNEL) $(SLIBS)

$(BIN).debug: $(OBJECT_FILES)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)


## installation
oxford: all
	cd $(BDIR) && cp insti.$(VERSION) ~/bin/
	cd ~/bin && rm -f insti && ln -s insti.$(VERSION) insti

.SUFFIXES:
