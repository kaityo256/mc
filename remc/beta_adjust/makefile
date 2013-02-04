all: beta.png

beta.png: beta.dat beta.plt
	gnuplot beta.plt

beta.dat: exmc.rb
	ruby exmc.rb > $@

clean:
	rm -f beta.dat beta.png
