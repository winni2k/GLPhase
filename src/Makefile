export CXX= g++-4.7+
export DEBUGFLAGS= -ggdb
export CXXFLAGS= -Wall -Wextra -std=c++11 -std=gnu++11 -mpopcnt
BOOSTLIB = -lboost_iostreams-mt #alternatively: /opt/boost/boost-1.54.0-intel/lib/libboost_iostreams.a
export LIBS += -lgsl -lgslcblas -lbam -ltabix -lbz2 $(BOOSTLIB) -lm -lstdc++ -lz

# path variables
export CPATH += :u-wkretzsch-snptools/samtools-0.1.16:u-wkretzsch-snptools/tabix-0.2.5
export CPLUS_INCLUDE_PATH += :olivier
export LIBRARY_PATH += :u-wkretzsch-snptools/samtools-0.1.16:u-wkretzsch-snptools/tabix-0.2.5
export LD_LIBRARY_PATH += $(LIBRARY_PATH)

BDIR := ../bin
VERSION := $(shell grep VERSION_MAJOR version.hpp | cut -f3 -d' ').$(shell grep VERSION_MINOR version.hpp | cut -f3 -d' ').$(shell grep VERSION_REVISION version.hpp | cut -f3 -d' ')
OBJECT_FILES := main.o impute.o insti.o emcchain.o relationshipGraph.o olivier/utils.o haplotype.o kMedoids.o kNN.o hapPanel.o
HEADER_FILES := insti.hpp.gch impute.hpp.gch emcchain.hpp.gch olivier/utils.hpp.gch relationship.hpp.gch version.hpp.gch enums.hpp.gch snp.hpp.gch hapPanel.hpp.gch MHSampler.hpp.gch

all: $(BDIR)/insti.$(VERSION)

$(BDIR)/insti.$(VERSION): CXXFLAGS += -O3

debug: CXXFLAGS += -O1
debug: CXXFLAGS += $(DEBUGFLAGS)
#	$(CXX) -E -x c++ - -v < /dev/null

debug: $(BDIR)/insti.$(VERSION).debug

test: debug

test : $(OBJECT_FILES) $(HEADER_FILES)
	$(MAKE) -C testing

testDebug: $(OBJECT_FILES) $(HEADER_FILES) debug
	$(MAKE) -C testing debug

clean:
	rm -f *.o *.hpp.gch olivier/utils.hpp.gch

.PHONY: clean test all testDebug debug

.DELETE_ON_ERROR:

### normal compilation

## precompile headers
%.hpp.gch: %.hpp
	$(CXX) $(CXXFLAGS) -x c++-header $<

hapPanel.hpp.gch: hapPanel.hpp snp.hpp.gch
	$(CXX) $(CXXFLAGS) -x c++-header $<

olivier/utils.hpp.gch: olivier/utils.hpp
	$(CXX) $(CXXFLAGS) -x c++-header $<

kNN.hpp.gch: kNN.hpp haplotype.hpp.gch
	$(CXX) $(CXXFLAGS) -x c++-header $<

kMedoids.hpp.gch: kMedoids.hpp haplotype.hpp.gch
	$(CXX) $(CXXFLAGS) -x c++-header $<

relationshipGraph.hpp.gch: relationshipGraph.hpp olivier/utils.hpp.gch enums.hpp.gch
	$(CXX) $(CXXFLAGS) -x c++-header $<

relationship.hpp.gch: relationship.hpp relationshipGraph.hpp.gch kMedoids.hpp.gch kNN.hpp.gch hapPanel.hpp.gch enums.hpp.gch
	$(CXX) $(CXXFLAGS) -x c++-header $<

MHSampler.hpp.gch: MHSampler.hpp olivier/utils.hpp.gch
	$(CXX) $(CXXFLAGS) -x c++-header $<

insti.hpp.gch: insti.hpp impute.hpp.gch emcchain.hpp.gch olivier/utils.hpp.gch relationship.hpp.gch version.hpp.gch enums.hpp.gch snp.hpp.gch hapPanel.hpp.gch MHSampler.hpp.gch
	$(CXX) $(CXXFLAGS) -x c++-header $<

## object files
%.o: %.cpp %.hpp.gch
	$(CXX) $(CXXFLAGS) -c $< $(LIBS)

olivier/utils.o: olivier/utils.cpp olivier/utils.hpp.gch
	$(CXX) $(CXXFLAGS) -Wno-sign-compare -c $< -o $@ $(LIBS)

main.o: main.cpp version.hpp.gch impute.hpp.gch insti.hpp.gch
	$(CXX) $(CXXFLAGS) -c $< $(LIBS)


## binaries
$(BDIR)/insti.$(VERSION): $(OBJECT_FILES)
	$(CXX) $(CXXFLAGS) $(OBJECT_FILES) -o $@ $(LIBS)
	cd $(BDIR) && rm -f insti && ln -s insti.$(VERSION) insti

$(BDIR)/insti.$(VERSION).debug: $(OBJECT_FILES)
	$(CXX) $(CXXFLAGS) $(OBJECT_FILES) -o $@ $(LIBS)


## installation
oxford: all
	cd $(BDIR) && cp insti.$(VERSION) ~/bin/
	cd ~/bin && rm -f insti && ln -s insti.$(VERSION) insti

