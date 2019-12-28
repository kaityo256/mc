all: a.out
CC=g++
CPPFLAGS=-O3 -march=native -std=c++11 -Wall -Wextra
DATA=Q2_L16.dat Q2_L32.dat Q3_L16.dat Q3_L32.dat


a.out: main.cpp
	$(CC) $(CPPFLAGS) $< -o $@


run: a.out
	./a.out 2 16 0.5 1.2 > Q2_L16.dat
	./a.out 2 32 0.5 1.2 > Q2_L32.dat
	./a.out 3 16 0.8 1.2 > Q3_L16.dat
	./a.out 3 32 0.8 1.2 > Q3_L32.dat

graph: $(DATA)
	gnuplot q2.plt
	gnuplot q3.plt

clean:
	rm -f a.out
