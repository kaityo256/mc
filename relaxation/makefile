all: a.out

CC=g++
CPPFLAGS=-O3 -std=c++11 -march=native

a.out: main.cpp
	$(CC) $(CPPFLAGS) $< -o $@

relaxation.dat: a.out
	./a.out > relaxation.dat


graph: relaxation.dat
	gnuplot relaxation.plt

clean:
	rm -f a.out
