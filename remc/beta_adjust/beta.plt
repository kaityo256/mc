set term png 
set out "beta.png"

set xla "Index"
set yla "Beta"
set log y
set pointsize 2.0
P = 0.6
c = -log(P)
a = ((2+c)+(c**2+4*c)**0.5)*0.5

p [][1:] "beta.dat" pt 6 t"Data", a**x lt 1 lc 0

