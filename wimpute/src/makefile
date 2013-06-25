export CC= g++-4.7+
export CFLAGS= -g -Wall -Wextra -O0 -fopenmp -std=c++11 -std=gnu++11 -g
export LIBS += -lgomp -lpthread -lgsl -lgslcblas -lz -lbam -ltabix -lbz2 -lboost_iostreams-mt
CPATH := $(GTEST_DIR)/include:snptools.cpp/samtools-0.1.16:snptools.cpp/tabix-0.2.5:olivier
LIBRARY_PATH := snptools.cpp/samtools-0.1.16:snptools.cpp/tabix-0.2.5
LD_LIBRARY_PATH := snptools.cpp/samtools-0.1.16:snptools.cpp/tabix-0.2.5

all: wimpute

test : wimpute.o impute.o emcchain.o utils.o
	cd testing && make

clean:
	rm -f *.o *.h.gch wimpute

.PHONY: clean

### normal compilation
impute.o: impute.cpp impute.h
	$(CC) $(CFLAGS) -c impute.cpp impute.h $(LIBS)

wimpute.o: wimpute.cpp wimpute.h impute.h emcchain.h relationshipGraph.h olivier/utils.h
	$(CC) $(CFLAGS) -c $^ $(LIBS)

relationshipGraph.o: relationshipGraph.cpp relationshipGraph.h
	$(CC) $(CFLAGS) -c $^ $(LIBS)

emcchain.o: emcchain.cpp emcchain.h
	$(CC) $(CFLAGS) -c $^ $(LIBS)

main.o: main.cpp impute.h wimpute.h
	$(CC) $(CFLAGS) -c $^ $(LIBS)

utils.o: olivier/utils.cpp olivier/utils.h 
	$(CC) -Wno-sign-compare $(CFLAGS) -c $^ $(LIBS)

wimpute: main.o impute.o wimpute.o emcchain.o utils.o relationshipGraph.o
	$(CC) $(CFLAGS) -o wimpute $^ $(LIBS)


