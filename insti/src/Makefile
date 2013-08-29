export CXX= g++-4.7+
export DEBUGFLAGS= -ggdb
export CXXFLAGS= -Wall -Wextra -std=c++11 -std=gnu++11 $(DEBUGFLAGS)
export LIBS += -lgomp -lpthread -lgsl -lgslcblas -lz -lbam -ltabix -lbz2 -lboost_iostreams-mt -lm
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

## headers
version.h.gch: version.h
	$(CXX) $(CXXFLAGS) -x c++-header version.h


## object files
impute.o: impute.cpp impute.h
	$(CXX) $(CXXFLAGS) -c impute.cpp impute.h $(LIBS)

insti.o: insti.cpp insti.h impute.h emcchain.h relationshipGraph.h olivier/utils.h
	$(CXX) $(CXXFLAGS) -c $^ $(LIBS)

relationshipGraph.o: relationshipGraph.cpp relationshipGraph.h olivier/utils.h
	$(CXX) $(CXXFLAGS) -c $^ $(LIBS)

emcchain.o: emcchain.cpp emcchain.h
	$(CXX) $(CXXFLAGS) -c $^ $(LIBS)

main.o: main.cpp impute.h insti.h version.h.gch
	$(CXX) $(CXXFLAGS) -c $^ $(LIBS)

olivier/utils.o: olivier/utils.cpp olivier/utils.h 
	$(CXX) -Wno-sign-compare $(CXXFLAGS) -c $^ $(LIBS)

insti: main.o impute.o insti.o emcchain.o relationshipGraph.o olivier/utils.o
	$(CXX) $(CXXFLAGS) -o insti $^ $(LIBS)

instiO0: main.o impute.o insti.o emcchain.o relationshipGraph.o olivier/utils.o
	$(CXX) $(CXXFLAGS) -o instiO0 $^ $(LIBS)


