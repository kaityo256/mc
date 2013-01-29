all: nucleation.png

nucleation.png: nucleation.dat nucleation.plt
	gnuplot nucleation.plt

nucleation.dat:
	ruby nucleation.rb > $@

clean:
	rm -f nucleation.dat nucleation.png
