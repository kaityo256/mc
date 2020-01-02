set term pngcairo
set out "crossing.png"

set xlabel "p"
set ylabel "Crossing Probability"

set ytics 0.2
p "bond.dat" w linespoints  pt 6 t "Bond"\
, "site.dat" w linespoints  pt 6 t "Site"\
