export CXX= g++-4.7+
export DEBUGFLAGS= -ggdb
export CXXFLAGS= -Wall -Wextra -std=c++11 -std=gnu++11 $(DEBUGFLAGS)
export LIBS += -lgomp -lpthread -lgsl -lgslcblas -lz -lbam -ltabix -lbz2 -lboost_iostreams-mt -lm -lstdc++
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
	$(CXX) $(CXXFLAGS) -x c++-header $^
insti.h.gch: insti.h
	$(CXX) $(CXXFLAGS) -x c++-header $^
impute.h.gch: impute.h
	$(CXX) $(CXXFLAGS) -x c++-header $^
olivier/utils.h.gch: olivier/utils.h
	$(CXX) $(CXXFLAGS) -x c++-header $^
relationshipGraph.h.gch: relationshipGraph.h
	$(CXX) $(CXXFLAGS) -x c++-header $^
emcchain.h.gch: emcchain.h
	$(CXX) $(CXXFLAGS) -x c++-header $^

## object files
impute.o: impute.cpp
	$(CXX) $(CXXFLAGS) -c $^ $(LIBS)

insti.o: insti.cpp
	$(CXX) $(CXXFLAGS) -c $^ $(LIBS)

relationshipGraph.o: relationshipGraph.cpp 
	$(CXX) $(CXXFLAGS) -c $^ $(LIBS)

emcchain.o: emcchain.cpp
	$(CXX) $(CXXFLAGS) -c $^ $(LIBS)

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c $^ $(LIBS)

olivier/utils.o: olivier/utils.cpp
	$(CXX) -Wno-sign-compare $(CXXFLAGS) -c $^ $(LIBS)

insti: main.o impute.o insti.o emcchain.o relationshipGraph.o olivier/utils.o
	$(CXX) $(CXXFLAGS) -o insti $^ $(LIBS)

instiO0: main.o impute.o insti.o emcchain.o relationshipGraph.o olivier/utils.o
	$(CXX) $(CXXFLAGS) -o instiO0 $^ $(LIBS)


