all: single cluster
CC=g++
CPPFLAGS=-O3
.SUFFIXES: .c .cc .h. .o

single: main.cc
	$(CC) $(CPPFLAGS) $< -o $@

cluster: main.cc
	$(CC) $(CPPFLAGS) -DCLUSTER $< -o $@

single.dat: single
	./single >single.dat

cluster.dat: cluster
	./cluster >cluster.dat

graph: single.dat cluster.dat
	gnuplot plot.plt

clean:
	rm -f single cluster single.dat cluster.dat
