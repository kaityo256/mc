all: site bond
CC=g++
CPPFLAGS=-O3 -march=native -std=c++11 -Wall -Wextra


bond: bond.cpp
	$(CC) $(CPPFLAGS) $< -o $@

site: site.cpp
	$(CC) $(CPPFLAGS) $< -o $@

bond.dat: bond
	./bond > bond.dat

site.dat: site
	./site > site.dat

graph: bond.dat site.dat
	gnuplot crossing.plt
	gnuplot percolation.plt

clean:
	rm -f site bond
