set term png 
set out "beta.png"

set xla "Index"
set yla "Beta"
set ytics (1,5,10)
set log y
set pointsize 2.0
p [][1:10] "beta.dat" pt 6 t"Data", 1.23616**x lt 1 lc 0

