all: wimpute

wimpute.o: wimpute.cpp
	g++ -c wimpute.cpp
add.o: add.cpp
	g++ -c add.cpp
wimpute: wimpute.o add.o add.h
	g++ -o wimpute wimpute.o add.o add.h