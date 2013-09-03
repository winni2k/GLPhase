CC=g++
CFLAGS=-Wall -Wextra -O3 -fopenmp
INCLUDE=-Isamtools-0.1.16 -Itabix-0.2.5
LIBS=-lgomp -lpthread -lgsl -lgslcblas -lz -lbam -ltabix -Lsamtools-0.1.16 -Ltabix-0.2.5 -lbz2

all:	document pileup varisite bamodel poprob probin impute hapfuse prob2vcf

document:	document.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o document document.cpp $(LIBS)
	./document >readme.txt

pileup:	pileup.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o pileup pileup.cpp $(LIBS)

varisite:	varisite.cpp ebd.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o varisite varisite.cpp $(LIBS)

bamodel:	bamodel.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o bamodel bamodel.cpp $(LIBS)

poprob:	poprob.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o poprob poprob.cpp $(LIBS)

probin:	probin.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o probin probin.cpp $(LIBS)

impute:	impute.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o impute impute.cpp $(LIBS)

hapfuse:	hapfuse.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o hapfuse hapfuse.cpp $(LIBS)

prob2vcf:	prob2vcf.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o prob2vcf prob2vcf.cpp $(LIBS)

