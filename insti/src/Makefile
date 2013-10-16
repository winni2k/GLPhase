export CXX= g++-4.7+
export DEBUGFLAGS= -ggdb
export CXXFLAGS= -Wall -Wextra -std=c++11 -std=gnu++11 -mpopcnt $(DEBUGFLAGS)
BOOSTLIB=-lboost_iostreams-mt #alternatively: /opt/boost/boost-1.54.0-intel/lib/libboost_iostreams.a
export LIBS += -lpthread -lgsl -lgslcblas -lbam -ltabix -lbz2 $(BOOSTLIB) -lm -lstdc++ -lz
export CPATH := u-wkretzsch-snptools/samtools-0.1.16:u-wkretzsch-snptools/tabix-0.2.5:olivier
export LIBRARY_PATH := u-wkretzsch-snptools/samtools-0.1.16:u-wkretzsch-snptools/tabix-0.2.5
export LD_LIBRARY_PATH := $(LIBRARY_PATH)

all: insti

insti: CXXFLAGS += -O1

instiO0: CXXFLAGS += -O0


debug:
	echo $(CPATH)

test : insti.o impute.o emcchain.o olivier/utils.o
	cd testing && make

clean:
	rm -f *.o *.hpp.gch insti instiO0 olivier/utils.hpp.gch

.PHONY: clean

.DELETE_ON_ERROR:

### normal compilation

## precompile headers
%.hpp.gch: %.hpp
	$(CXX) $(CXXFLAGS) -x c++-header $<

olivier/utils.hpp.gch: olivier/utils.hpp
	$(CXX) $(CXXFLAGS) -x c++-header $<

kMedoids.hpp.gch: kMedoids.hpp haplotype.hpp.gch
	$(CXX) $(CXXFLAGS) -x c++-header $<

relationshipGraph.hpp.gch: relationshipGraph.hpp olivier/utils.hpp.gch
	$(CXX) $(CXXFLAGS) -x c++-header $<

relationship.hpp.gch: relationship.hpp relationshipGraph.hpp.gch kMedoids.hpp.gch
	$(CXX) $(CXXFLAGS) -x c++-header $<

insti.hpp.gch: insti.hpp impute.hpp.gch emcchain.hpp.gch olivier/utils.hpp.gch relationship.hpp.gch version.hpp.gch
	$(CXX) $(CXXFLAGS) -x c++-header $<

## object files
%.o: %.cpp %.hpp.gch
	$(CXX) $(CXXFLAGS) -c $< $(LIBS)

olivier/utils.o: olivier/utils.cpp olivier/utils.hpp.gch
	$(CXX) $(CXXFLAGS) -Wno-sign-compare -c $< -o $@ $(LIBS)

main.o: main.cpp version.hpp.gch impute.hpp.gch insti.hpp.gch
	$(CXX) $(CXXFLAGS) -c $< $(LIBS)

OBJECT_FILES = main.o impute.o insti.o emcchain.o relationshipGraph.o olivier/utils.o haplotype.o kMedoids.o

insti: $(OBJECT_FILES)
	$(CXX) $(CXXFLAGS) $(OBJECT_FILES)  -o $@ $(LIBS)

instiO0: $(OBJECT_FILES)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJECT_FILES) $(LIBS)


