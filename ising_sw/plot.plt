set term png
set out "plot.png"

set xla "Beta"
set yla "Magnetization^2"
set key at 0.40,0.8
p "single.dat" pt 6 t "Single-Flip", "cluster.dat" pt 7 t "Cluster Update"
