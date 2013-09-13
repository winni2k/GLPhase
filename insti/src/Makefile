export CXX= g++-4.7+
export DEBUGFLAGS= -ggdb
export CXXFLAGS= -Wall -Wextra -std=c++11 -std=gnu++11 -mpopcnt $(DEBUGFLAGS)
export LIBS += -lpthread -lgsl -lgslcblas -lz -lbam -ltabix -lbz2 -lboost_iostreams-mt -lm -lstdc++
export CPATH := snptools.cpp/samtools-0.1.16:snptools.cpp/tabix-0.2.5:olivier
export LIBRARY_PATH := snptools.cpp/samtools-0.1.16:snptools.cpp/tabix-0.2.5
export LD_LIBRARY_PATH := snptools.cpp/samtools-0.1.16:snptools.cpp/tabix-0.2.5

all: insti

insti: CXXFLAGS += -O3

instiO0: CXXFLAGS += -O0


debug:
	echo $(CPATH)

test : insti.o impute.o emcchain.o olivier/utils.o
	cd testing && make

clean:
	rm -f *.o *.h.gch insti

.PHONY: clean

.DELETE_ON_ERROR:

### normal compilation

## precompile headers
%.h.gch: %.h 
	$(CXX) $(CXXFLAGS) -x c++-header $^

olivier/utils.h.gch: olivier/utils.h
	$(CXX) $(CXXFLAGS) -x c++-header $^

HEADER_GCH = impute.h.gch insti.h.gch emcchain.h.gch relationshipGraph.h.gch olivier/utils.h.gch haplotype.h.gch version.h.gch

## object files
%.o: %.cpp $(HEADER_GCH)
	$(CXX) $(CXXFLAGS) -c $< $(LIBS)

olivier/utils.o: olivier/utils.cpp $(HEADER_GCH)
	$(CXX) -Wno-sign-compare $(CXXFLAGS) -c $< $(LIBS)

OBJECT_FILES = main.o impute.o insti.o emcchain.o relationshipGraph.o olivier/utils.o haplotype.o 

insti: $(OBJECT_FILES) $(HEADER_GCH)
	$(CXX) $(CXXFLAGS) -o insti $(OBJECT_FILES) $(LIBS)

instiO0: $(OBJECT_FILES) $(HEADER_GCH)
	$(CXX) $(CXXFLAGS) -o instiO0 $(OBJECT_FILES) $(LIBS)


