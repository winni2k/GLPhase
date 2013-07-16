export CC= g++-4.7+
export CXXFLAGS= -g -Wall -Wextra -O0 -fopenmp -std=c++11 -std=gnu++11 -g
export LIBS += -lgomp -lpthread -lgsl -lgslcblas -lz -lbam -ltabix -lbz2 -lboost_iostreams-mt
export CPATH := snptools.cpp/samtools-0.1.16:snptools.cpp/tabix-0.2.5:olivier
export LIBRARY_PATH := snptools.cpp/samtools-0.1.16:snptools.cpp/tabix-0.2.5
export LD_LIBRARY_PATH := snptools.cpp/samtools-0.1.16:snptools.cpp/tabix-0.2.5

all: insti

debug:
	echo $(CPATH)

test : insti.o impute.o emcchain.o utils.o
	cd testing && make

clean:
	rm -f *.o *.h.gch insti

.PHONY: clean

.DELETE_ON_ERROR:

### normal compilation
impute.o: impute.cpp impute.h
	$(CC) $(CXXFLAGS) -c impute.cpp impute.h $(LIBS)

insti.o: insti.cpp insti.h impute.h emcchain.h relationshipGraph.h olivier/utils.h
	$(CC) $(CXXFLAGS) -c $^ $(LIBS)

relationshipGraph.o: relationshipGraph.cpp relationshipGraph.h
	$(CC) $(CXXFLAGS) -c $^ $(LIBS)

emcchain.o: emcchain.cpp emcchain.h
	$(CC) $(CXXFLAGS) -c $^ $(LIBS)

main.o: main.cpp impute.h insti.h
	$(CC) $(CXXFLAGS) -c $^ $(LIBS)

utils.o: olivier/utils.cpp olivier/utils.h 
	$(CC) -Wno-sign-compare $(CXXFLAGS) -c $^ $(LIBS)

insti: main.o impute.o insti.o emcchain.o utils.o relationshipGraph.o
	$(CC) $(CXXFLAGS) -o insti $^ $(LIBS)


