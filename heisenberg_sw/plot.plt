set term png
set out "magnetization.png"
set xla "Temperature"
set yla "Magnetization^2"
set key at 2,0.8
set style data linespoints
p "single.dat" pt 6 t "Single-Flip", "cluster.dat" pt 7 t "Cluster Update"

set out "heatcapacity.png"
set xla "Temperature"
set yla "Heat Capacity"
set key at 2.5,1.2
set style data linespoints
p "single.dat" u 1:3 pt 6 t "Single-Flip", "cluster.dat" u 1:3 pt 7 t "Cluster Update"


