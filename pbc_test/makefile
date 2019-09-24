all: speed.out mc.out

CC=g++
CPPFLAGS=-O3 -std=c++11

speed.out: main.cpp
	$(CC) $(CPPFLAGS) $< -o $@

mc.out: main.cpp
	$(CC) $(CPPFLAGS) -DMC $< -o $@

graph: mc.out
	./mc.out > data.dat
	gnuplot graph.plt

clean:
	rm -f speed.out mc.out
