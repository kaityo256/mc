EXE: single sw wolff
all: $(EXE)

CC=g++
CPPFLAGS=-O3
.SUFFIXES: .c .cc .h. .o .dat .acf

single: main.cc
	$(CC) $(CPPFLAGS) $< -o $@

sw: main.cc
	$(CC) $(CPPFLAGS) -DSW $< -o $@

wolff: main.cc
	$(CC) $(CPPFLAGS) -DWOLFF $< -o $@

single.dat: single
	./single > $@

sw.dat: sw
	./sw > $@

wolff.dat: wolff
	./wolff > $@

.dat.acf:
	ruby auto.rb $< > $@

acf.png: single.acf sw.acf wolff.acf
	gnuplot acf.plt

graph: acf.png

clean:
	rm -f $(EXE) *.dat *.png *.acf
