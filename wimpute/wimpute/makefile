CC=g++
CFLAGS=-Wall -Wextra -O3 -fopenmp
INCLUDE=-Isnptools.cpp/samtools-0.1.16 -Isnptools.cpp/tabix-0.2.5
LIBS=-lgomp -lpthread -lgsl -lgslcblas -lz -lbam -ltabix -Lsnptools.cpp/samtools-0.1.16 -Lsnptools.cpp/tabix-0.2.5 -lbz2


all: wimpute

impute.o: impute.cpp impute.h
	$(CC) $(CFLAGS) $(INCLUDE) -c impute.cpp impute.h $(LIBS)

main.o: main.cpp impute.h
	$(CC) $(CFLAGS) $(INCLUDE) -c main.cpp impute.h $(LIBS)

wimpute: main.o impute.o
	$(CC) $(CFLAGS) $(INCLUDE) -o wimpute main.o impute.o $(LIBS)