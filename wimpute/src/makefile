CC=g++-4.7+
CFLAGS=-Wall -Wextra -O3 -fopenmp -std=c++11 -std=gnu++11
INCLUDE=-Isnptools.cpp/samtools-0.1.16 -Isnptools.cpp/tabix-0.2.5
LIBS=-lgomp -lpthread -lgsl -lgslcblas -lz -lbam -ltabix -Lsnptools.cpp/samtools-0.1.16 -Lsnptools.cpp/tabix-0.2.5 -lbz2


all: wimpute

impute.o: impute.cpp impute.h
	$(CC) $(CFLAGS) $(INCLUDE) -c impute.cpp impute.h $(LIBS)

wimpute.o: wimpute.cpp wimpute.h
	$(CC) $(CFLAGS) $(INCLUDE) -c wimpute.cpp wimpute.h $(LIBS)

main.o: main.cpp impute.h wimpute.h
	$(CC) $(CFLAGS) $(INCLUDE) -c main.cpp impute.h wimpute.h $(LIBS)

wimpute: main.o impute.o wimpute.o
	$(CC) $(CFLAGS) $(INCLUDE) -o wimpute main.o impute.o wimpute.o $(LIBS)