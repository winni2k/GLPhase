export CC= g++-4.7+
export CFLAGS= -g -Wall -Wextra -O3 -fopenmp -std=c++11 -std=gnu++11 -g
export LIBS += -lgomp -lpthread -lgsl -lgslcblas -lz -lbam -ltabix -lbz2 -lboost_iostreams-mt
CPLUS_INCLUDE_PATH := $(GTEST_DIR)/include:snptools.cpp/samtools-0.1.16:snptools.cpp/tabix-0.2.5:olivier
LIBRARY_PATH := snptools.cpp/samtools-0.1.16:snptools.cpp/tabix-0.2.5

all: wimpute

test : wimpute.o impute.o emcchain.o utils.o
	cd testing && make

clean:
	rm -f *.o *.h.gch wimpute

.PHONY: clean

### normal compilation
impute.o: impute.cpp impute.h
	$(CC) $(CFLAGS) -c impute.cpp impute.h $(LIBS)

wimpute.o: wimpute.cpp wimpute.h impute.h emcchain.h
	$(CC) $(CFLAGS) -c wimpute.cpp wimpute.h olivier/utils.h $(LIBS)

emcchain.o: emcchain.cpp emcchain.h
	$(CC) $(CFLAGS) -c emcchain.cpp emcchain.h $(LIBS)

main.o: main.cpp impute.h wimpute.h
	$(CC) $(CFLAGS) -c main.cpp impute.h wimpute.h $(LIBS)

utils.o: olivier/utils.h olivier/utils.cpp
	$(CC) -Wno-sign-compare $(CFLAGS) -c olivier/utils.cpp olivier/utils.h $(LIBS)


wimpute: main.o impute.o wimpute.o emcchain.o utils.o
	$(CC) $(CFLAGS) -o wimpute main.o impute.o wimpute.o emcchain.o utils.o $(LIBS)


