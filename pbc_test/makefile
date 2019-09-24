all: a.out

CC=g++
CPPFLAGS=-O3 -std=c++11

a.out: main.cpp
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	rm -f a.out
