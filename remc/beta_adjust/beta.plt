set term png 
set out "beta.png"

set xla "Index"
set yla "Beta"
set ytics (1,5,10,15)
set log y
set pointsize 2.0
p [][1:15] "beta.dat" pt 6 t"Data", 1.34361**x lt 1 lc 0

