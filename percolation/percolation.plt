set term pngcairo
set out "percolation.png"

set xlabel "p"
set ylabel "Percolation Probability"

set key at 0.9,0.2
set ytics 0.2
p "bond.dat" u 1:3 w linespoints  pt 6 t "Bond"\
, "site.dat" u 1:3 w linespoints  pt 6 t "Site"\
