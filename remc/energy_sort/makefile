all: index.png

index.png: index.dat index.plt
	gnuplot index.plt

index.dat: exmc.rb
	ruby exmc.rb > $@

clean:
	rm -f index.dat index.png
