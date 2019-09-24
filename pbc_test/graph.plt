set term pngcairo
set out "graph.png"
set xla "beta"
set yla "m^2"
set ytics 0.2
unset key
p "data.dat" pt 7 
