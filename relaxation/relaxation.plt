set term pngcairo
set out "relaxation.png"

f(x) = (a-b)*exp(-x/c)+b
a = 1.0
b = 0.05
c = 40
fit [50:500] f(x) "relaxation.dat" via a,b,c
unset key
set xlabel "MCs"
set ylabel "m^2"
p [0:500] "relaxation.dat" pt 6, f(x) lw 2.0 lc "black"
