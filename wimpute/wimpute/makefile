all: wimpute

main.o: main.cpp
	g++ -c main.cpp
add.o: add.cpp
	g++ -c add.cpp
wimpute: main.o add.o add.h
	g++ -o wimpute main.o add.o add.h