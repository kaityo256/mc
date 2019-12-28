set term pngcairo
set out "q2.png"

set xla "Beta"
set yla "U"
set yra [1:3]
set style data line
bc = log(1.0+sqrt(2.0))
set arrow 1 from bc,1.5 to bc, 1
p "Q2_L16.dat" u 1:3 t "L=16"\
, "Q2_L32.dat" u 1:3 t "L=32"\
