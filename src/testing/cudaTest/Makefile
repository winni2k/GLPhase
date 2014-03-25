### Makefile --- 

## Author: wkretzsch@gmail.com
## Version: $Id: Makefile,v 0.0 2014/03/24 17:40:53 winni Exp $
## Keywords: 
## X-URL: 

all:test

%.o: %.cu
	nvcc -dc -arch=sm_35 -m64 $<

%.o: %.cpp
	g++ $< -c -o $@

link.o: a.o b.o
	nvcc -arch=sm_35 -m64 -dlink $^ -o link.o

test: link.o a.o b.o c.o
	g++ $^ -o $@ -lcudart
	./$@



### Makefile ends here
